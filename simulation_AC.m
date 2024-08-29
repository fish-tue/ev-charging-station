%% SIMULATION - Simulate AC incentive mechanism
% Description:
%   Simulation of the daily Nash equilibrium and the decisions of each user
% Outputs:
%   simulation.mat (simulation output)
% Assumptions and limitations:
%   - Urgency distribution is uniform
%   - Initial karma distribution is dicrete uniform with support 
%   {25p_1,25p_1,+1,...,50p_1}
%   - k_ref is 0
% Other m-files required:
%   - n_arcs_individual.m
% MAT-files required:
%   - initialization.mat (generated with initialization.m)
%   - optimization.mat (generated with optimize_x.m)
%   - prices.mat (generated with pricing_r.m)
%   - simulation_baseline.mat (generated with simulation_base.m)
% Toolboxes required: none
% Authors: Lotte Hollander, Leonardo Pedroso, Mauro Salazar, ...
%   Ricardo de Castro, W.P.M.H. (Maurice) Heemels
% Revision history:
%   22/08/2024 - Lotte Hollander
%       * Implementation for MSc thesis [1]
%   05/02/2024 - Leonardo Pedroso
%       * Added final publication reference to [2] 
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Hollander, L. Pedroso, R. de Castro and M. Salazar, "Fair 
%   Incentive Mechanisms for Differentiated Services in a Public Electric 
%   Vehicle Charging Station", Eindhoven University of Technology, 2024.
%   [2] L. Pedroso, W.P.M.H. Heemels and M. Salazar, "Urgency-Aware Routing
%   in Single Origin-Destination Itineraries Through Artificial 
%   Currencies," 2023 62nd IEEE Conference on Decision and Control (CDC), 
%   Singapore, Singapore, 2023, pp. 4142-4149, 
%   doi: 10.1109/CDC49753.2023.10383739.

%% Initialization
clear;

% Load parameters
load('initialization.mat', 'r', 'N', 'M', 'par', 'T', 'u_min', 'u_bar', ...
    'u_max', 'd_xp', 'd_xnp')
% Load xp_star and xnp_star
load('optimization.mat')
% Load prices
load('prices.mat');
% Load baseline test result
load('simulation_baseline.mat')

% Enforce constraints on charging spots
IsConstrained = true;

% Nash iterations
it_nash_max = 50;
it_nash_epsl = 1e-3;

% Simulation
T_sim = 400;                % number of days
k_xp_init_max = 10*max(p_Cp);
k_xnp_init_max = 10*max(p_Cnp);

%% Simulation
% Simulation - Initialization
k = zeros(N,T_sim+1);       % karma/artificial currency
u = zeros(N,T_sim+1);       % urgency
y = zeros(N,T_sim+1);       % user's choice at each t
it_nash = zeros(T_sim+1,1); % number of iterations for the nash equilibrium

% Allocation based on priority
x = zeros(r,T_sim+1);       % aggregate choices
xp = zeros(r,T_sim+1);      % aggregate choices priority users
xnp = zeros(r,T_sim+1);     % aggregate choices non-priority users

% Priority mask
y_p = binornd(1,par.P_p,N,1);

% Karma initialization
% uniform initial karma distribution
% k(:,1) = round(unifrnd(k_init_max/2,k_init_max,N,1));
k(logical(y_p),1) = round(unifrnd(k_xp_init_max*0.75,k_xp_init_max,sum(y_p),1));
k(~y_p,1) = round(unifrnd(k_xnp_init_max*0.75,k_xnp_init_max,sum(~y_p),1));

% k_ref initialization
k_ref = zeros(N,1);

% Simulation loop
for t = 1:T_sim+1
    % Pick daily urgencies
    u(:,t) = unifrnd(0,u_max,N,1);
    % Agents decisions
    y_go = binornd(1,par.P_go,N,1);      % mask of charging agents 
    % Init Nash equilibrium iterations
    if t > 1     
        x(:,t) = x(:,t-1);
        xp(:,t) = xp(:,t-1);
        xnp(:,t) = xnp(:,t-1);
    end
    % Until convergence of Nash equilibrium is reached
    I = eye(r);
    y(:,t) = ones(N,1);
    while true
        % Previous iteration's flows
        x_prev = x(:,t);
        xp_prev = xp(:,t);
        xnp_prev = xnp(:,t);
        
        % Agents decisions
        for i = 1:N
            if ~y_go(i)                 % non-charging agent
                y(i,t) = 0;
                continue; 
            end 

            % With capacity constraint 
            if IsConstrained == true
                x_aux = xp(:,t) + xnp(:,t) - I(:,y(i,t))/N;
                x_cap_indx = x_aux*N > M; 
                x_cap_indx(1) = false;  % exclude r_0 in capacity constraint

                d_aux_xp = d_xp(x(:,t) - I(:,y(i,t))/N + ones(r,1)/N);
                d_aux_xp(x_cap_indx) = 100*max(d_aux_xp);

                d_aux_xnp = d_xnp(x(:,t)-I(:,y(i,t))/N + ones(r,1)/N);
                d_aux_xnp(x_cap_indx) = 100*max(d_aux_xnp);

                if y_p(i)
                    y(i,t) = n_arcs_individual(...
                    d_aux_xp,par.P_go*T,p_Cp,k(i,t),k_ref(i)+(1-par.P_go)*T*p_Cp(1),u(i,t),u_min,u_bar,u_max);
                elseif ~y_p(i)
                    y(i,t) = n_arcs_individual(...
                    d_aux_xnp,par.P_go*T,p_Cnp,k(i,t),k_ref(i)+(1-par.P_go)*T*p_Cnp(1),u(i,t),u_min,u_bar,u_max);
                end
        
            % Without capacity constraint    
            elseif IsConstrained == false
                if y_p(i)
                    y(i,t) = n_arcs_individual(...
                    d_xp(x(:,t)-I(:,y(i,t))/N+ones(r,1)/N),par.P_go*T,p_Cp,k(i,t),k_ref(i)+(1-par.P_go)*T*p_Cp(1),u(i,t),u_min,u_bar,u_max); 
                elseif ~y_p(i)
                    y(i,t) = n_arcs_individual(...
                    d_xnp(x(:,t)-I(:,y(i,t))/N+ones(r,1)/N),par.P_go*T,p_Cnp,k(i,t),k_ref(i)+(1-par.P_go)*T*p_Cnp(1),u(i,t),u_min,u_bar,u_max); 
                end
            end

            % Aggregate behaviour
            for j = 1:r
                xp(j,t) = sum(y(find(y_p==true),t)==j)/N;   % allocation of priority users
                xnp(j,t) = sum(y(find(y_p==false),t)==j)/N; % allocation of non-priority users
            end
            x(:,t) = xp(:,t) + xnp(:,t);                    % allocation of all users
        end
        
        % Catch infeasibility
        if sum(isnan(y(:,t)))
            error("Caught infeasibility!");
        end
        % Catch Nash equilibrium iterations not converging
        it_nash(t) = it_nash(t)+1;
        if norm(x_prev-x(:,t))<1e-3/N, break; end
        if it_nash(t) > it_nash_max
            warning("Nash iterations did not converge.");
            break; 
        end
    end
    % Karma dynamics
    if t < T_sim+1
        for i = 1:N
            if ~y_go(i)
                % k(i,t+1) = k(i,t);
               if y_p(i)
                 k(i,t+1) = k(i,t) - p_Cp(1);
               else
                 k(i,t+1) = k(i,t) - p_Cnp(1);
               end
            % if do travel
            elseif y_p(i)
                k(i,t+1) = k(i,t)-p_Cp(y(i,t));
            else 
                k(i,t+1) = k(i,t)-p_Cnp(y(i,t));
            end
        end
    end  
end

%% Plot simulation evolution
% "Color blind safe" color pallet from IBM Design Library 
% magenta, pink, blue, orange, yellow
newcolors = ["#785EF0";"#DC267F";"#4EC3F9";"#FE6100";"#FFB000"];

% Decision Evolution
figure('Position',4*[0 0 192 144]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
for j = 1:r
    aux_x = (0:T_sim)';
    aux_y = sum(x(j:r,:),1)';
    [aux_x,aux_y] = stairs_vector(aux_x,aux_y);
    area(aux_x,aux_y,'LineWidth',0.5);
end
for j = 1:r
    plot([0 T_sim],[sum(x_star(j:r)) sum(x_star(j:r))],'--','Color','black','LineWidth',3);
end
legend({' Charging Elsewhere',' Timeslot Morning',' Timeslot Afternoon',' Timeslot Evening', ' Timeslot Night',' $\mathrm{x}^\star$'},'Location','northeast','Interpreter','latex');
ylabel('$\mathbf{x}^\mathrm{NE}(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([0 1]);
hold off;
colororder(newcolors)
% Save figure to .fig and .svg formats
savefig('./fig/05_decision_bounded.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_decision_bounded.svg');

% Decision Evolution Xp
figure('Position',4*[0 0 192 144]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
for j = 1:r
    aux_x = (0:T_sim)';
    aux_y = sum(xp(j:r,:),1)';
    [aux_x,aux_y] = stairs_vector(aux_x,aux_y);
    area(aux_x,aux_y,'LineWidth',0.5);
end
for j = 1:r
    plot([0 T_sim],[sum(xp_star(j:r)) sum(xp_star(j:r))],'--','Color','black','LineWidth',3);
end
legend({' Charging Elsewhere',' Timeslot Morning',' Timeslot Afternoon',' Timeslot Evening', ' Timeslot Night',' $\mathrm{x}^{p\star}$'},'Location','northeast','Interpreter','latex');
ylabel('$\mathbf{x}^\mathrm{NE}(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([0 0.2]);
hold off;
colororder(newcolors)
% Save figure to .fig and .svg formats
savefig('./fig/05_decision_xp_bounded.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_decision_xp_bounded.svg');

% Decision Evolution Xnp
figure('Position',4*[0 0 192 144]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
for j = 1:r
    aux_x = (0:T_sim)';
    aux_y = sum(xnp(j:r,:),1)';
    [aux_x,aux_y] = stairs_vector(aux_x,aux_y);
    area(aux_x,aux_y,'LineWidth',0.5);
end
for j = 1:r
    plot([0 T_sim],[sum(xnp_star(j:r)) sum(xnp_star(j:r))],'--','Color','black','LineWidth',3);
end
legend({' Charging Elsewhere',' Timeslot Morning',' Timeslot Afternoon',' Timeslot Evening', ' Timeslot Night',' $\mathrm{x}^{np\star}$'},'Location','northeast','Interpreter','latex');
ylabel('$\mathbf{x}^\mathrm{NE}(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([0 0.5]);
hold off;
colororder(newcolors)
% Save figure to .fig and .svg formats
savefig('./fig/05_decision_xnp_bounded.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_decision_xnp_bounded.svg');

% Karma levels
k_mean = mean(k,1);
k_max = max(k);
k_min = min(k);
k_var = sqrt(var(k));
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
aux_x2 = [0:T_sim fliplr(0:T_sim)]';
aux_y2 = [k_max fliplr(k_min)]';
[aux_x2,aux_y2] = stairs_vector(aux_x2,aux_y2);
fill(aux_x2,aux_y2*1e-2,'k',...
    'LineWidth',2,'FaceColor','black','FaceAlpha',0.2,'EdgeAlpha',0);
aux_x1 = [0:T_sim fliplr(0:T_sim)]';
aux_y1 = [k_mean+k_var fliplr(k_mean-k_var)]';
aux_y1(aux_y1<0) = 0;
[aux_x1,aux_y1] = stairs_vector(aux_x1,aux_y1);
fill(aux_x1,aux_y1*1e-2,'k',...
    'LineWidth',2,'FaceColor','black','FaceAlpha',0.4,'EdgeAlpha',0);
stairs(0:T_sim,k_mean*1e-2,'LineWidth',2,'Color','black');
stairs(aux_x1,aux_y1*1e-2,'LineWidth',1,'Color',[0.4 0.4 0.4]);
stairs(aux_x2,aux_y2*1e-2,'LineWidth',1,'Color',[0.4 0.4 0.4]);
legend({' $\max_i$/$\min_i$ $\{k^i(t)\}$',' $\hat{k}(t)\pm\sigma_k(t)$',' $\hat{k}(t)$'},...
    'Location','northeast','Interpreter','latex');
%ylabel('$k^i(t), i = 1,\ldots,M$','Interpreter','latex');
ylabel('Karma level $\times 10^{-2}$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
hold off;
% Save figure to .fig and .svg formats
savefig('./fig/05_karma.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_karma.svg');

% System's cost
cost_soc = zeros(T_sim+1,1);
cost_soc_rel_opt = zeros(T_sim+1,1);
cost_soc_opt = cost.x(x_star,xp_star);
for t = 1:T_sim+1
    cost_soc(t) = cost.x(x(:,t), xp(:,t));
    cost_soc_rel_opt(t) = (cost_soc(t)-cost_soc_opt)/cost_soc_opt;
end
figure('Position',4*[0 0 192 144/2]);
hold on;
grid on;
box on;
set(gca,'FontSize',16);
set(gca, 'Layer', 'top');
set(gca,'TickLabelInterpreter','latex')
stairs(0:T_sim,cost_soc_rel_opt*100,'LineWidth',1.5,'Color','black');
stairs(0:T_sim, cost_soc_rel_base*100,'LineWidth',1.5,'Color','#FE6100');
legend('With artificial currency', 'Without artificial currency', 'Interpreter','latex');
ylabel('$\Delta$ Societal cost $(\%)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylim([-45 180]);
yticks([-40 0 20 40 60 80 100])
hold off;
% Save figure to .fig and .svg formats
savefig('./fig/05_cost.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_karma.svg');

% Average societal cost
AC_soc_cost = mean(cost_soc_rel_opt(T_sim-20:T_sim,:))*100;
base_soc_cost = mean(cost_soc_rel_base(T_sim-20:T_sim,:))*100;

% Plot converged allocation
x_t = [mean(xnp(:,T_sim-60:T_sim),2), mean(xp(:,T_sim-60:T_sim),2)];
R_names = {' Charging Elsewhere';'Timeslot Morning';'Timeslot Afternoon';'Timeslot Evening';'Timeslot Night'};
figure()
box on
set(gca, 'Layer', 'top');
bar(1:numel(R_names), x_t,'stacked')
grid on
set(gca,'xminorgrid','on','yminorgrid','on')
set(gca, 'XTickLabel',R_names, 'XTick',1:numel(R_names),'FontSize',12,'TickLabelInterpreter','latex');
legend({'Non-Priority Users', 'Priority Users'},'Location','northeast','Interpreter','latex');
ylabel('$\mathbf{x}$','Interpreter','latex');
colororder("sail")
ylim([0 0.3]);
% Save figure to .fig and .svg formats
savefig('./fig/05_simulation_bounded.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_simulationX_bounded.svg');

%% Save results
save('simulation.mat');
clear

%% Auxiliary function
% To user the equivalent of stairs with area and fill functions
function [x,y] = stairs_vector(x,y)
    x = [x';x'];
    y = [y';y'];
    y = y(:);
    x = x([2:end end])';
end