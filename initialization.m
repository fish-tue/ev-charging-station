%% INITIALIZATION - Generate framework for EV charging station
% Description:
%   Runs all parameters for the model and computes the discomfort for
%   priority and non-priority users according to Section III in [2].
% Outputs:
%   - initialization.mat ( 'r', 'N', 'M', 'par', 'T', 'u_min', 'u_bar', ...
%    'u_max', 'd_xp', 'd_xnp')
% Assumptions and limitations:
%   - Urgency distribution is uniform
% Other m-files required:
%   - none
% MAT-files required:
%   - none
% Toolboxes required:
%    - none
% Authors: Lotte Hollander, Leonardo Pedroso, Mauro Salazar, 
%    Ricardo de Castro
% Revision history:
%   22/08/2024 - Lotte Hollander
%       * Implementation for MSc thesis [1]
% References: 
%   [1] L. Hollander, L. Pedroso, R. de Castro and M. Salazar, "Fair 
%   Incentive Mechanisms for Differentiated Services in a Public Electric 
%   Vehicle Charging Station", Eindhoven University of Technology, 2024.

%% Parameters
clear;

r = 1 + 4;                      % timeslots INCLUDING r0. order = [r0, r_day, r_night]
N = 100;                        % users
M = 10;                         % maximum number of charging spots

par.P_veh_max = 11;             % max charging power of Peugeot e-208 [kW]
par.P_station_max = 19.2*5;     % max charging power L2 dual charger [kW]. 
                                % Note: dual chargers split power equally

% Probabilities
par.P_p = 0.2;                  % proportion of belonging to priority group
par.P_go = 0.5;                 % proportion of charging (someplace)
par.P_home = 1-par.P_go;        % proportion of not charging that day

tau = linspace(0, 24, 100);     % continuous time of day

% Urgency is an uniform distribution
u_min = 0;
u_bar = 1;
u_max = 2;

T = 4;                          % individual decision window 

% Energy mix
par.E0 = 0.25;                  % estimated guess 
par.Egrid = 0.23;               % based on CAISO data

%% Discomfort function - for different sized timeslots
alpha = [0; 0.03 * ones(r-1,1)];      % scaling variable, r0 NOT dependend on x

% Discomfort of charging elsewhere
d0_xp = 1.2*ones(1, length(tau));     % discomfort r0 (charging elsewhere)
d0_xnp = 0.6*ones(1, length(tau));    % discomfort r0 (charging elsewhere)

% Generate plugin discomfort function
d_plugin = zeros(1,length(tau));      % initalize plugin discomfort vector

for i = 1:length(tau)
    if tau(i) <= 6 || tau(i) > 21
        d_plugin(1,i) = 1;
    else
        d_plugin(1,i) = exp((-tau(i).^2)/11) + 0.4 * exp((-(tau(i)-15).^2)/6) + exp((-(tau(i)-23).^2)/13);
    end
end

% find values between 6AM-9PM
indx_day = find(tau > 6 & tau <= 21.1);          
t_day = buffer(tau(indx_day), length(indx_day)/(r-2));
d_day = buffer(d_plugin(indx_day), length(indx_day)/(r-2));

% find values between 9PM-6AM
indx_night = find(tau <= 6 | tau > 21.1);
t_night1 = tau(tau<=6)';                     % 0AM-6AM
t_night2 = tau(tau >= 21.1)';                % 9PM-0AM
d_night = d_plugin(indx_night);

% Mean plug-in discomfort: % r_0 , r_day , r_night
dplug_xp = [mean(d0_xp), mean(d_day), mean(d_night)];     
dplug_xnp = [mean(d0_xnp), mean(d_day), mean(d_night)];

% Combine with function for charging power
d_xp = @(x) dplug_xp.' + alpha .* max(par.P_veh_max - (par.P_station_max./(x*N)),0);
d_xnp = @(x) dplug_xnp.' + alpha .* max(par.P_veh_max - (par.P_station_max./(x*N)),0);

% Compute duration of each timeslot
par.delta_r = zeros(1,(r-1));
for i=1:(r-2)
par.delta_r(1,i) = round(t_day(end,i)-t_day(1,i));
end
par.delta_r(1,(r-1)) = round(t_night1(end)-t_night1(1)) + round(t_night2(end) - t_night2(1));

%% Function curve of renewable energy
% Load renewable dataset
addpath('data\');               % add data folder to path

% Table with time, solar, wind, geothermal, biomass, biogass, small hydro
data.Allrenew = readmatrix('CAISO-renew-spring.csv','Range','B2:G289');      % Select season
data.Prenew = data.Allrenew(:,1);                                            % Solar energy [MW]

% Generate maximum renewable power curve 
P_renew_max = data.Prenew / (28*4);                                          % Convert range [kW]
P_renew_max = interp1(linspace(0,24,length(P_renew_max)), P_renew_max, tau); % Interpolate with tau considering 0-24h

P_renew_day = buffer(P_renew_max(indx_day), length(indx_day)/(r-2));         % Values between 6AM-9PM
P_renew_night = P_renew_max(indx_night);                                     % Values between 9PM-6AM

% Average P_renew
par.P_renew_mean = [mean(P_renew_day).'; mean(P_renew_night).'];             % [day slots, night slot]

%% Plots
%  magenta, pink, blue, orange, yellow
newcolors = ["#785EF0";"#DC267F";"#4EC3F9";"#FE6100";"#FFB000"];

% Plot timeslots of d_plugin
figure()
hold on
box on
set(gca,'FontSize',14);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
set(gca,'xminorgrid','on','yminorgrid','on')
% plot timeslot 0-6AM
plot(tau(find(tau <= 6.1)), d_plugin(find(tau <= 6.1)), Color="#FFB000", LineWidth=2.0);
% plot timeslots during the day
plot(t_day, d_day, LineWidth=2.0);  
% plot timeslot 21-0PM
plot(tau(find(tau > 21)), d_plugin(find(tau > 21)), Color="#FFB000", LineWidth=2.0);     
plot(tau, d0_xp, 'Color', 'black', 'LineStyle', ':', linewidth=2.0)
plot(tau, d0_xnp, 'Color', 'black', 'LineStyle','--', LineWidth=2.0)
grid on
xlim([0 24]);
ylabel('$\mathrm{d_{plugin}}(\tau)$','Interpreter','latex', 'FontSize', 14);
xlabel('time $\tau$ [h]','Interpreter','latex', 'FontSize', 14);
lgd = legend({'Night time slot','Morning time slot','Afternoon time slot','Evening time slot', '','$d_0$ for $\mathbf{x^p}$', '$d_0$ for $\mathbf{x^{np}}$'},'Location','northeast','Interpreter','latex');
fontsize(lgd,12,'points');
xlim([0 24]);
ylim([0 1.8]);
hold off;
colororder(newcolors)
% Save figure to .fig and .svg formats
savefig('./fig/03_discomfort.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/03_discomfort.svg');

% Plot renewable energy power curve
figure()
box on;
hold on;
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
set(gca,'xminorgrid','on','yminorgrid','on')
plot(tau,P_renew_max, 'LineWidth', 2.0)
grid on;
lgd = legend({'April 30, 2023'},'Location','northwest','Interpreter','latex');
fontsize(lgd,11.5,'points');
ylabel('renewable power [kW]','Interpreter','latex','FontSize', 14);
xlabel('time $\tau$ [h]','Interpreter','latex','FontSize', 14);
xlim([0 24])
ylim([-1 120])
hold off;
% Save figure to .fig and .svg formats
savefig('./fig/03_Prenew.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/03_Prenew.svg');

%% Save initialization results
save('initialization.mat');
clear;
