%% R_PRICING - Compute near-optimal resource prices
% Description:
%   - Computes the near-optimal resource prices acc. to Section IV in [1].
%   - Follows the methodology in Section V.A in [2].
% Outputs:
%   - prices.mat (r dimensional vector of prices for xp and xnp, 'p_Cp' and 'p_Cnp')
% Assumptions and limitations:
%   - Urgency distribution is uniform
%   - k_ref is a discrete uniform distribution with support {k_ref \in N : 
%    k_ref = 0 \lor k_ref = p_j, j = 0,...,r}
%   - p_home = p_0
%   - Priority users are not allocated to r_0 and are excluded from 
%    home pricing
% Other m-files required:
%   - n_arcs_gamma_individual.m
% MAT-files required:
%   - initialization.mat (generated with initialization.m)
%   - optimization.mat (generated with optimize_x.m)
% Toolboxes required:
%   - Global Optimization Toolbox
% Authors: Lotte Hollander, Leonardo Pedroso, Mauro Salazar, ...
%   Ricardo de Castro, W.P.M.H. (Maurice) Heemels
% Revision history:
%   22/08/2024 - Lotte Hollander
%       * Implementation for MSc thesis [1]
%       * Includes home-pricing and pricing resource r_0.
%       * Includes decoupled pricing for xp and xnp.
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
% Load x_star
load('optimization.mat')

%% Design parameters
% Bounds on price magnitude
ub_p = 20;             
lb_p = -ub_p;
x_star_quant_exp = 3;
kref = 0;

% Sort by discomfort
[d_star_sorted_xp, d_indx_xp] = sort(d_star_xp);
[d_star_sorted_xnp, d_indx_xnp] = sort(d_star_xnp);

xp_star_sorted = xp_star(d_indx_xp);
xnp_star_sorted = xnp_star(d_indx_xnp);

% Exclude (xp_r)* = 0
check_xp = logical(xp_star_sorted >= 1e-6);       % check for non-zero allocations
H_exclXp = diag(check_xp);                        % create diagonal matrix
H_exclXp = H_exclXp(~all(H_exclXp == 0, 2),:);    % remove all zero rows
nvars_xp = length(find(check_xp));                % number of non-zero variables xp

% Exclude (xnp_r)* = 0
check_xnp = logical(xnp_star_sorted >= 1e-2);     % check for non-zero allocations
H_exclXnp = diag(check_xnp);                      % create diagonal matrix
H_exclXnp = H_exclXnp(~all(H_exclXnp == 0, 2),:); % remove all zero rows
nvars_xnp = length(find(check_xnp));              % number of non-zero variables xnp

% Quantization
xp_star_quant = round((10^x_star_quant_exp)*(H_exclXp*xp_star_sorted))/(10^x_star_quant_exp);
xnp_star_quant = round((10^x_star_quant_exp)*(H_exclXnp*xnp_star_sorted))/(10^x_star_quant_exp);

%% Define pricing problem
% Define cost function

% Cost function without homepricing
F_xp = @(p) norm(n_arcs_stationary_flows(H_exclXp*d_star_sorted_xp,T,p', ...
  u_max,u_bar,u_min,par.P_go,par.P_home)*par.P_p - H_exclXp*xp_star_sorted)/sum(H_exclXp*xp_star_sorted);

% Cost function with homepricing
d_indx_xnp_2 = H_exclXnp*d_indx_xnp;
d_indx_r0 = find(d_indx_xnp_2==1);     % find r0 index in d_sorted

F_xnp = @(p) norm(n_arcs_stationary_flows_kref_phome(H_exclXnp*d_star_sorted_xnp,...
  T*par.P_go,p',kref + (1-par.P_go)*T*p(d_indx_r0),u_max,u_bar,u_min,...
  par.P_go,par.P_home,p(d_indx_r0))*(1-par.P_p)- H_exclXnp*xnp_star_sorted)/sum(H_exclXnp*xnp_star_sorted);

% Linear inequality constraints
% p_1 > 0
% p_n < 0
% p_j > p_j+1 (p_j >= p_j+1 + 1)

% Linear inequality constraints for xp
c_ = zeros(nvars_xp-1,1);
c_(1) = 1;
k = zeros(nvars_xp,1);
k(1) = 1;
k(2) = -1;
s = zeros(nvars_xp,1);
s(1) = -1;
s_ = zeros(nvars_xp,1);
s_(end) = 1;
A_xp = [s'; s_'; -toeplitz(c_,k)];
b_xp = [-1;-1;-ones(nvars_xp-1,1)];

% Linear inequality constraints for xnp
c_ = zeros(nvars_xnp-1,1);
c_(1) = 1;
k = zeros(nvars_xnp,1);
k(1) = 1;
k(2) = -1;
s = zeros(nvars_xnp,1);
s(1) = -1;
s_ = zeros(nvars_xnp,1);
s_(end) = 1;
A_xnp = [s'; s_'; -toeplitz(c_,k)];
b_xnp = [-1;-1;-ones(nvars_xnp-1,1)];

% Linear equality constraints
% Price bounds
lb = lb_p;
ub = ub_p;

% Nonlinear inequality constraints
nonlcon = [];

% Stopping criterion
% - Cost below ga_threshold
% - Or constant cost for ga_stall_it generations

% Hyperparameters
ga_population_size = 1e4;
ga_threshold = 1e-2;
ga_stall_it = 2;

%% Solve pricing problem for xp - Genetic Alg.

% Faster population generation: xp_star perpendicular basis
xp_star_pp_basis = zeros(nvars_xp,nvars_xp-1);
I = eye(nvars_xp);
for i = 1:nvars_xp-1
    % Gram-Schmidt orthogonalization
    for j = 1:i
        if j ~= 1
            aux = aux - xp_star_pp_basis(:,j-1)*(xp_star_pp_basis(:,j-1)'*I(:,i));
        else
            aux = I(:,i) - xp_star_quant*(xp_star_quant'*I(:,i)/(xp_star_quant'*xp_star_quant));
        end
    end
    % normalization
    xp_star_pp_basis(:,i) = aux/sqrt(aux'*aux); 
end

% GA for xp
fprintf("---------------------------------------------------------------------");
tic;
options = optimoptions('ga','Display','diagnose',...
    'CreationFcn', {@CustomCreationFcn, xp_star_pp_basis, ub_p},...
    'PopulationSize',ga_population_size,'FitnessLimit',ga_threshold, ...
    'MaxStallGenerations', ga_stall_it, 'UseParallel', true);
p_opt_xp = ga(F_xp,nvars_xp,A_xp,b_xp,[],[],lb,ub,nonlcon,options)';
p_opt_xp = round(p_opt_xp);
toc;
fprintf("---------------------------------------------------------------------\n");

% Output GA solution
fprintf("Genetic alg. prices:\t\t");
fprintf("%d ", p_opt_xp);
fprintf("\n");

%% Solve pricing problem for xnp - Genetic Alg.

% Faster population generation: xnp_star perpendicular basis
xnp_star_pp_basis = zeros(nvars_xnp,nvars_xnp-1);
I = eye(nvars_xnp);
for i = 1:nvars_xnp-1
    % Gram-Schmidt orthogonalization
    for j = 1:i
        if j ~= 1
            aux = aux - xnp_star_pp_basis(:,j-1)*(xnp_star_pp_basis(:,j-1)'*I(:,i));
        else
            aux = I(:,i) - xnp_star_quant*(xnp_star_quant'*I(:,i)/(xnp_star_quant'*xnp_star_quant));
        end
    end
    % normalization
    xnp_star_pp_basis(:,i) = aux/sqrt(aux'*aux); 
end

% GA for xnp
fprintf("---------------------------------------------------------------------");
tic;
options = optimoptions('ga','Display','diagnose',...
    'CreationFcn', {@CustomCreationFcn, xnp_star_pp_basis, ub_p},...
    'PopulationSize',ga_population_size,'FitnessLimit',ga_threshold, ...
    'MaxStallGenerations', ga_stall_it, 'UseParallel', true);
p_opt_xnp = ga(F_xnp,nvars_xnp,A_xnp,b_xnp,[],[],lb,ub,nonlcon,options)';
p_opt_xnp = round(p_opt_xnp);
toc;
fprintf("---------------------------------------------------------------------\n");

% Output GA solution
fprintf("Genetic alg. prices:\t\t");
fprintf("%d ",p_opt_xnp);
fprintf("\n");

%% Save prices

% Unsort the prices
p_Cp = zeros(1,r);                      % all x_r*=0 are equal to zero
p_Cp(H_exclXp*d_indx_xp) = p_opt_xp; 
p_Cp = p_Cp';

p_Cnp = 20*ub*ones(1,r);                % all x_r*=0 have high prices
p_Cnp(H_exclXnp*d_indx_xnp) = p_opt_xnp;
p_Cnp = p_Cnp';

save('prices.mat','p_Cp','p_Cnp');
clear;

%% Auxiliary functions
% Aggregate stategy over distribution of k_ref
% (k_ref is a discrete uniform distribution with support {k_ref \in N : 
% k_ref = 0 \lor k_ref = p_j, j = 1,...,n})
function x_inf = n_arcs_stationary_flows(d,T,p,s_max,s_bar,s_min,P_go,P_home)
    % Support of theta_p distribution (distribution of k_ref values)
    k_ref = 0; %round([0;p(p>0)]);
    x_inf = zeros(length(p),length(k_ref));
    % Compute x_inf for each value in the support of theta_p
    for j = 1:length(k_ref)
        x_inf(:,j) = n_arcs_stationary_flows_kref(d,T,p,k_ref(j),s_max,s_bar,s_min,P_go,P_home);
    end
    % Expected aggregate (because theta_p is uniform)
    x_inf = mean(x_inf,2);
end

% Aggregate stategy for a particular value of k_ref
% Home resource is now also priced
% (k_ref is a discrete uniform distribution with support {k_ref \in N : 
% k_ref = 0 \lor k_ref = p_j, j = 1,...,n})
% Assumptions on inputs:
%   - d: Sorted
%   - d: No equal discomforts
%   - p: p_1 > p_2 > ... > p_n
function x_inf = n_arcs_stationary_flows_kref_phome(d,T,p,k_ref,s_max,s_bar,s_min,P_go,P_home,p_home)   
    % Parameters 
    epsl_cmp = 1e-10;
    epsl_eta = 1e-5;
    n = length(p);
    % Project p and k_ref into suitable space
    p = round(p);
    p_home = round(p_home);
    k_ref = round(k_ref);
    % Find gcd among all
    p_factor = p(1);
    for t = 1:length(p)
        p_factor = gcd(p_factor,p(t));
    end
    % Karma limits (according to Section III in [1])
    k_min = max([0 k_ref+p(end)*(T+1)]);
    k_max = k_ref+p(1)*(T+1)-p(end);
    % Karma indexing
    % k_i = k_ref + (i+i0-1)*p_factor;
    i0 = ceil((k_min-k_ref)/p_factor);
    i_max = round(((k_max-k_ref)/p_factor) + 1-i0)+5*abs(p_home);
    k_i = @(i) k_ref+(i+i0-1)*p_factor;
    i_p1 = find(k_i((1:i_max)') == p(1));
    % Compute gamma
    gamma = n_arcs_gamma_individual(d,T,p,k_i((1:i_max)'),k_ref,s_min,s_bar,s_max);
    % Probability of choosing j given karma k P_go*P(j|k_i)
    P_j_given_k = zeros(i_max,n);
    for j = 1:n
        if j~= 1
            P_j_given_k(:,j) = P_go*(gamma(j-1,:)-gamma(j,:))'*s_bar/s_max;
        else
            P_j_given_k(:,j) = P_go*(s_max-gamma(j,:)*s_bar)'/s_max;
        end      
    end
    % Build probability transition matrix A
    A = zeros(i_max,i_max);%P_home*eye(i_max);
    for i_2 = 1:i_max
        for i_1 = 1:i_max
            for j = 1:n
                if abs(k_i(i_1)-(k_i(i_2)-p(j))) < epsl_cmp
                    A(i_1,i_2) = A(i_1,i_2) + P_j_given_k(i_2,j);
                end
            end
            % New here
            if abs(k_i(i_1)-(k_i(i_2)-p_home)) < epsl_cmp
                A(i_1,i_2) = A(i_1,i_2) + P_home;
            end
        end
        if k_i(i_max) + epsl_cmp < (k_i(i_2)-p_home)
            A(end,i_2) = A(end,i_2) + P_home;
        end
    end
    % Remove all unreachable components
    reachable = false(i_max,1);
    searched = false(i_max,1);
    reachable(i_p1,1) = true;
    % Search all that are reachable but not searched
    while sum(reachable.*(~searched))
        for i = 1:i_max
            if ~searched(i) && reachable(i) 
                reachable(A(:,i)~=0) = true;
                searched(i) = true;
            end
        end
    end
    A_irreducible = A(reachable,reachable);
    irreducible_idx = (1:i_max)';
    irreducible_idx = irreducible_idx(reachable);

    % Compute irreducible stationary karma distribution
    eta = zeros(i_max,1);
    eta(i_p1,1) = 1;
    eta_irreducible = eta(irreducible_idx);
    count = 0;

    A_irreducible = A_irreducible^50;
    while norm((A_irreducible-eye(size(A_irreducible,1)))*eta_irreducible) > epsl_eta
        eta_irreducible = A_irreducible*eta_irreducible;
        % norm((A_irreducible-eye(size(A_irreducible,1)))*eta_irreducible);
        count = count+1;
        if count >= 100e3
            break
        end    
    end

    % Expand to the whole domain
    eta = zeros(i_max,1);
    eta(reachable) = eta_irreducible;
    
    % Compute stationary flows
    x_inf = (eta'*P_j_given_k)';
end

% Aggregate stategy for a particular value of k_ref
% (k_ref is a discrete uniform distribution with support {k_ref \in N : 
% k_ref = 0 \lor k_ref = p_j, j = 1,...,n})
% Assumptions on inputs:
%   - d: Sorted
%   - d: No equal discomforts
%   - p: p_1 > p_2 > ... > p_n
function x_inf = n_arcs_stationary_flows_kref(d,T,p,k_ref,s_max,s_bar,s_min,P_go,P_home)   
    % Parameters 
    epsl_cmp = 1e-10;
    epsl_eta = 1e-5;
    n = length(p);
    % Project p and k_ref into suitable space
    p = round(p);
    k_ref = round(k_ref);
    % Find gcd among all
    p_factor = p(1);
    for t = 1:length(p)
        p_factor = gcd(p_factor,p(t));
    end
    % Karma limits (according to Section III in [1])
    k_min = max([0 k_ref+p(end)*(T+1)]);
    k_max = k_ref+p(1)*(T+1)-p(end);
    % Karma indexing
    % k_i = k_ref + (i+i0-1)*p_factor;
    i0 = ceil((k_min-k_ref)/p_factor);
    i_max = round(((k_max-k_ref)/p_factor) + 1-i0);
    k_i = @(i) k_ref+(i+i0-1)*p_factor;
    i_p1 = find(k_i((1:i_max)') == p(1));
    % Compute gamma
    gamma = n_arcs_gamma_individual(d,T,p,k_i((1:i_max)'),k_ref,s_min,s_bar,s_max);
    % Probability of choosing j given karma k P_go*P(j|k_i)
    P_j_given_k = zeros(i_max,n);
    for j = 1:n
        if j~= 1
            P_j_given_k(:,j) = P_go*(gamma(j-1,:)-gamma(j,:))'*s_bar/s_max;
        else
            P_j_given_k(:,j) = P_go*(s_max-gamma(j,:)*s_bar)'/s_max;
        end      
    end
    % Build possibly transition matrix A
    A = P_home*eye(i_max);
    for i_1 = 1:i_max
        for i_2 = 1:i_max
            for j = 1:n
                if abs(k_i(i_1)-(k_i(i_2)-p(j))) < epsl_cmp
                    A(i_1,i_2) = A(i_1,i_2) + P_j_given_k(i_2,j);
                end
            end
        end
    end
    % Remove all unreachable components
    reachable = false(i_max,1);
    searched = false(i_max,1);
    reachable(i_p1,1) = true;
    % Search all that are reachable but not searched
    while sum(reachable.*(~searched))
        for i = 1:i_max
            if ~searched(i) && reachable(i) 
                reachable(A(:,i)~=0) = true;
                searched(i) = true;
            end
        end
    end
    A_irreducible = A(reachable,reachable);
    irreducible_idx = (1:i_max)';
    irreducible_idx = irreducible_idx(reachable);

    % Compute irreducible stationary karma distribution
    eta = zeros(i_max,1);
    eta(i_p1,1) = 1;
    eta_irreducible = eta(irreducible_idx);
    while norm((A_irreducible-eye(size(A_irreducible,1)))*eta_irreducible) > epsl_eta
        eta_irreducible = A_irreducible*eta_irreducible;
    end

    % Expand to the whole domain
    eta = zeros(i_max,1);
    eta(reachable) = eta_irreducible;
    
    % Compute stationary flows
    x_inf = (eta'*P_j_given_k)';
end

% Karma thresholds for unitary decisions, i.e., [y_bar]_j = 1 for some j
function k_ij = k_ij(i,j,k_ref,p,T)
    k_ij = k_ref+p(i)+T*p(j);
end

% Faster population generation function (CustomCreationFcn)
function population = CustomCreationFcn(GenomeLength,FitnessFcn,options,...
    x_star_pp_basis,p_norm)
    % Init initial population matrix
    population = zeros(options.PopulationSize,GenomeLength);
    for i = 1:options.PopulationSize
        % Fill with random vectors generated by the basis perpendicular to
        % x_star
        while true         
            aux = x_star_pp_basis*(rand(GenomeLength-1,1)-0.5);
            % Make sure the prices are sorted, the first is positive, and 
            % the last negative
            if issorted(aux,'descend') && aux(1) > 0 && aux(end) < 0
                population(i,:) = round(p_norm*aux'/sqrt(aux'*aux));
                break;
            end
        end
    end
end
