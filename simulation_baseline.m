%% simulation - simulate baseline test
% Description:
%   Simulation of the daily Nash equilibrium and the decisions of each user
% Outputs:
%   simulation_baseline.mat (simulation output)
% Assumptions and limitations:
%   - Urgency distribution is uniform
%   - Initial karma distribution is dicrete uniform with support 
%   {25p_1,25p_1,+1,...,50p_1}
%   - k_ref is 0
% Other m-files required:
%   - n_arcs_individual.m
%   - n_arcs_gamma_individual.m
% MAT-files required:
%   - initialization.mat (generated with initialization.m)
%   - optimization.mat (generated with optimize_x.m)
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

%% Parameters
% Load parameters
load('initialization.mat', 'r', 'N', 'M', 'par', 'T', 'u_min', 'u_bar', ...
    'u_max', 'd_xp', 'd_xnp')
% Load x_star
load('optimization.mat')

% Baseline test prices
p_Cp = [0; 0; 0; 0; 0;];
p_Cnp = [0; 0; 0; 0; 0];

% Enforce constraints on charging spots
IsConstrained = true;  

% Nash iterations
it_nash_max = 50;
it_nash_epsl = 1e-3;

% Simulation
T_sim = 400;    % number of days
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
            if ~y_go(i) % non-charging agent
                y(i,t) = 0;
                continue; 
            end 

            % T = Pgo*T, kref_aux = kref + (1_Pgo)*T*p_r0, p_r0 = p_Phome 
            if IsConstrained == true
                x_aux = xp(:,t) + xnp(:,t) - I(:,y(i,t))/N;
                x_cap_indx = x_aux*N > M; 
                x_cap_indx(1) = false;     % DO NOT INCLUDE R0 IN CAPACITY

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

%% Compute societal cost

% System's cost
cost_soc_base = zeros(T_sim+1,1);
cost_soc_rel_base = zeros(T_sim+1,1);
cost_soc_opt = cost.x(x_star,xp_star);
for t = 1:T_sim+1
    cost_soc_base(t) = cost.x(x(:,t), xp(:,t));
    cost_soc_rel_base(t) = (cost_soc_base(t)-cost_soc_opt)/cost_soc_opt;
end

%% Plot simulation evolution
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
savefig('./fig/05_simX-baseline.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_simX-baseline.svg');

%% Save simulation results
save('simulation_baseline.mat', 'x_t', 'cost_soc_rel_base', 'cost_soc_base');
clear;

%% Auxiliary function
% To user the equivalent of stairs with area and fill functions
function [x,y] = stairs_vector(x,y)
    x = [x';x'];
    y = [y';y'];
    y = y(:);
    x = x([2:end end])';
end