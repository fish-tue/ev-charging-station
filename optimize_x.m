%% OPTIMIZATION - Optimize allocation of priority and non-priority users
% Description:
%   Finds the optimal allocation of users through minimizing grid usage 
%   and underuse of chargers. Used for the optimization problem in [1].
% Outputs:
%   - optimization.mat ('xp_star','xnp_star','x_star', 'd_star_xp', ...
%    'd_star_xnp', 'cost')
% Assumptions and limitations:
%   - The discomfort of priority users is minimized
%   - The allocation must be non-negative
%   - The allocation must not exceed maximum capacity
% Other m-files required:
%   - none
% MAT-files required:
%   - initialization.mat (generated with initialization.m)
% Toolboxes required:
%    - YALMIP [2]
% Authors: Lotte Hollander, Leonardo Pedroso, Mauro Salazar, ...
%    Ricardo de Castro
% Revision history:
%   22/08/2024 - Lotte Hollander
%       * Implementation for MSc thesis [1]
% References: 
%   [1] L. Hollander, L. Pedroso, R. de Castro and M. Salazar, "Fair 
%   Incentive Mechanisms for Differentiated Services in a Public Electric 
%   Vehicle Charging Station", Eindhoven University of Technology, 2024.
%   [2] Lofberg, Johan. "YALMIP: A toolbox for modeling and optimization in
%   MATLAB." In 2004 IEEE international conference on robotics and 
%   automation, pp. 284-289. IEEE, 2004.

%% Initialization
clear;

% Load parameters and renewable energy curve
load('initialization.mat', 'r', 'N', 'M', 'par', 'T', 'u_min', 'u_bar', ...
    'u_max', 'd_xp', 'd_xnp')

%% Compute x_star
yalmip('clear');

% Check feasibility
feasibility = logical(N*par.P_p*par.P_go <= M*(r-1));
if ~feasibility
    warning('Feasibilty constraint violated, change parameter input values')
end

% Define optimization variables
xp_star = sdpvar(r,1);
xnp_star = sdpvar(r,1);

% Constraints
constraints = [];
for i = 2:r
% Do not exceed maximum number of charging spots
constraints = [constraints; (xp_star(i,1)+xnp_star(i,1))*N <= M];     
end
constraints = [constraints;
    0 <= xp_star;                                % xp_star is non-negative
    0 <= xnp_star;                               % xnp_star is non-negative
    ones(1,r)*xp_star == par.P_p*par.P_go;       % number of priority users charging
    ones(1,r)*xnp_star == (1-par.P_p)*par.P_go;  % number of non-priority users charging
];

% Average non-renewable energy from r_0
obj.obj_underuse = (1-par.E0)*(par.P_veh_max*(xp_star(1,1)+xnp_star(1,1))*N); 

% Average non-renewable energy from charging station
obj.P_req = min(par.P_veh_max*(xp_star(2:r,1)+xnp_star(2:r,1))*N, par.P_station_max);
obj.P_renew = min(obj.P_req, par.P_renew_mean);
obj.P_grid = obj.P_req - obj.P_renew;

obj.P_non_renew = (1-par.Egrid)*obj.P_grid*par.delta_r;
obj.obj_griduse = (1\sum(par.delta_r))*sum(obj.P_non_renew);

% Cost function
c_grid = obj.obj_underuse + obj.obj_griduse; 

% Optimize
diagnostics_grid = optimize(constraints,c_grid);

% Get value
xp_star = value(xp_star);
xnp_star = value(xnp_star);
x_star = xp_star + xnp_star;

% Get discomfort
d_star_xp = d_xp(x_star);       % perceived discomfort by xp at optimum
d_star_xnp = d_xnp(x_star);     % perceived discomfort by xnp at optimum

%% Compute xp_star and xnp_star
gamma = 1; 

% Define optimization variables
xp_star = sdpvar(r,1);
xnp_star = sdpvar(r,1);

% Constraints
constraints = [constraints;
    0 <= xp_star;                             % x_star is non-negative
    0 <= xnp_star;                            % x_star is non-negative
    xp_star + xnp_star == x_star;
    ones(1,r)*xp_star == par.P_p*par.P_go;
];

% minimize discomfort for priority users
c_d = gamma * (d_xp(x_star)' * xp_star);

% Optimize
diagnostics_d = optimize(constraints, c_d);

% Get value
xp_star = value(xp_star)
xnp_star = value(xnp_star)

%% Societal Cost function
% Function to compute the societal cost
% Input variables are x and xp, where x = xp + xnp

% Minimize underuse
cost.underuse = @(x) (1-par.E0)*(par.P_veh_max*x*N);

% Minimize grid use
cost.griduse = @(x) (1\sum(par.delta_r))*sum((par.delta_r*(1-par.Egrid)*...
    min(par.P_veh_max*(x)*N, par.P_station_max) - ...
    min(min(par.P_veh_max*(x)*N,par.P_station_max),par.P_renew_mean)));

% Minimize discomfort for xp
cost.d = @(x,xp) 10 * (d_xp(x)' * xp);

% Overall societal cost
cost.x = @(x,xp) cost.underuse(x(1,1)) + cost.griduse(x(2:r,1)) + cost.d(x,xp);

%% Plots

% Plot optimal allocation
x_t = [xnp_star, xp_star];
R_names = {'Charge Elsewhere';'Timeslot Morning';'Timeslot Afternoon';'Timeslot Evening';'Timeslot Night'};

figure()
box on
set(gca, 'Layer', 'top');
bar(1:numel(R_names), x_t,'stacked')
colororder(["#544886";"#A28FF5"])
grid on
ylim([0 0.3])
set(gca,'xminorgrid','on','yminorgrid','on')
set(gca, 'XTickLabel',R_names, 'XTick',1:numel(R_names),'TickLabelInterpreter','latex','FontSize',16);
legend({'Non-Priority Users', 'Priority Users'},'Location','northeast','Interpreter','latex','FontSize', 16);
ylabel('$\mathbf{x}$','Interpreter','latex','FontSize', 16);
% Save figure to .fig and .svg formats
savefig('./fig/05_optimalX.fig');
set(gcf,'renderer','Painters');
saveas(gcf,'./fig/05_optimalX.svg');

%% Save results
save('optimization.mat', 'xp_star', 'xnp_star', 'x_star', 'd_star_xp', 'd_star_xnp', 'cost')
clear