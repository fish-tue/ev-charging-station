function gamma = gamma_n_arcs_individual(d,T,p,k,k_ref,s_min,s_bar,s_max) 
%% GAMMA_N_ARCS_INDIVIDUAL - compute gamma_j, j = 1,...,n 
% Description:
%   Computes user's best response strategy thresholds, gamma_j, j = 1,...,
%   n, according to Theorem III.2 in [1] for the Karma levels k.
% Syntax: gamma = gamma_n_arcs_individual(d,T,p,k,k_ref,s_min,s_bar,s_max);
% Inputs:
%   - d: discomfort vector (nx1)
%   - T: user decision horizon length
%   - p: prices vector (nx1)
%   - k: vector of Karma levels 
%   - k_ref: single value of reference Karma level
%   - s_min: minimum sensitivity
%   - s_bar: expected sensitivity
%   - s_max: maximum sensitivity
% Outputs:
%   - gamma: matrix (n x length(k)) each column corresponding to \gamma_j,
%   j = 1,...,n, evaluated at a different Karma level in k
% Assumptions and limitations:
%   - d: d_1 < d_2 < ... < d_n
%   - p: p_1 > p_2 > ... > p_n
% Other m-files required: none
% MAT-files required: none
% Toolboxes required: none
% Authors: Leonardo Pedroso, W.P.M.H. (Maurice) Heemels, Mauro Salazar
% Revision history:
%   05/02/2024 - Leonardo Pedroso
%       * Added final publication reference to [1] 
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, W.P.M.H. Heemels and M. Salazar, "Urgency-Aware Routing
%   in Single Origin-Destination Itineraries Through Artificial 
%   Currencies," 2023 62nd IEEE Conference on Decision and Control (CDC), 
%   Singapore, Singapore, 2023, pp. 4142-4149, 
%   doi: 10.1109/CDC49753.2023.10383739.
    
%% Init
n = length(d);
k_len = length(k);

%% Compute y_bar_star
% Optimal y_bar assuming the optimal decision is j for karma level k(t)
y_bar_star = zeros(n,n,k_len);
for t = 1:k_len
    for j = 1:n
        % Check if integer decision j is feasible
        if k(t) < k_ij(j,n,k_ref,p,T)
            y_bar_star(:,j,t) = nan(n,1);
            continue;
        end
        
        % Check if y_bar = l_1 is feasible
        if k(t) < k_ij(j,1,k_ref,p,T)
            % Get the solution for each region a 
            y_bar_star_j_aux = zeros(n,n);
            for a = 1:n
                y_bar_star_j_aux(:,a) = y_bar_n_arcs_individual(j,a,d,T,p,k(t),k_ref);
            end
            % Compute a_hat: the region of the best solution
            [~,a_hat] = min(d'*y_bar_star_j_aux);
            % Retrieve the best solution among all regions
            y_bar_star(:,j,t) = y_bar_star_j_aux(:,a_hat);
        elseif k(t) >= k_ij(j,1,k_ref,p,T)
            % If, for decision j, the user can afford to travel only on
            % the fastest link in y_bar, then he/she will do so
            y_bar_star(:,j,t) = zeros(n,1);
            y_bar_star(1,j,t) = 1;
        end
    end
end

%% Compute gamma_ij
% Sensivity threshold for choosing i over j for karma level k(t)
gamma_ij = nan(n,n,k_len);
for t = 1:k_len
    for i = 1:n
        for j = i+1:n  
            % Check if both i and j are feasible
            if k(t) >= max([0, p(i), k_ij(i,n,k_ref,p,T)])
                gamma_ij(i,j,t) = T*d'*(y_bar_star(:,j,t)-y_bar_star(:,i,t))/(d(i)-d(j));
            else
                gamma_ij(i,j,t) = +inf;
            end
        end
    end
end

%% Compute gamma
gamma = zeros(n,k_len);  
for t = 1:k_len
    gamma(n,t) = s_min/s_bar;
    for j = n-1:-1:1
        gamma_j_underbar = max([-inf gamma_ij(j+1,j+2:n,t)]);
        gamma_j_bar = min([+inf; gamma_ij(1:j,j+1,t)]);
        if gamma_j_underbar < gamma_j_bar
            gamma(j,t) = sat(gamma_j_bar,s_min/s_bar,s_max/s_bar);
        else
            gamma(j,t) = gamma(j+1,t);
        end
    end
end    

end
 
%% Auxiliary functions

% y_bar assuming the optimal decision is j in region a for karma level k
function y_bar_j_a = y_bar_n_arcs_individual(j,a,d,T,p,k,k_ref)
    % Compute j_hat
    j_hat = j_hat_n_arcs_individual(j,a,d,T,p,k,k_ref);
    % Compute y_bar_j_a
    I = eye(length(d));
    y_bar_j_a = (1/(T*(p(a)-p(j_hat))))*...
            (I(:,a)*(k-k_ij(j,j_hat,k_ref,p,T)) - ...
            I(:,j_hat)*(k-k_ij(j,a,k_ref,p,T)));
end

% Karma thresholds for unitary decisions, i.e., [y_bar]_j = 1 for some j
function k_ij = k_ij(i,j,k_ref,p,T)
    k_ij = k_ref+p(i)+T*p(j);
end

% Compute j_hat
function j_hat = j_hat_n_arcs_individual(j,a,d,T,p,k,k_ref)
    % Mask for feasible indices of argmin
    n = length(d);
    mask = true(n,1);  
    % Check feasible indices
    for i = 1:n 
        if i == a
            mask(a) = false;
            continue;
        end
        k_thr = [k_ij(j,a,k_ref,p,T);k_ij(j,i,k_ref,p,T)];
        if k < min(k_thr) || k > max(k_thr)
            mask(i) = false;
        end
    end
    % Compute (d_i_d_a)/(p_a-p_i) for feasible indices
    alpha = nan(length(d),1);
    for i = 1:length(d)
        if ~mask(i), continue; end
        alpha(i) = (d(i)-d(a))/(p(a)-p(i));
    end
    % Retrieve argmin (j_hat)
    [~,j_hat] = min(alpha); 
end

% Saturation function
function out = sat(in,out_min,out_max)
    if in >= out_max
        out = out_max;
    elseif in <= out_min
        out = out_min;
    else
        out = in;
    end
end
