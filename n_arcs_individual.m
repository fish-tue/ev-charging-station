function j_star = n_arcs_individual(d,T,p,k,k_ref,s,s_min,s_bar,s_max)
%% N_ARCS_INDIVIDUAL - compute solution to individual user's problem
% Description:
%   Computes user's best response strategy according to Problem II.2 in [1]
% Syntax: j_star = n_arcs_individual(d,T,p,k,k_ref,s,s_min,s_bar,s_max)
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
%   - j_star: a solution to Problem II.2 in [1]
% Assumptions and limitations:
%   - If there are multiple solutions, it outputs the arc with the smaller
%   index
% Other m-files required:
%   - n_arcs_gamma_individual.m
% MAT-files required: none
% Toolboxes required: none
% Authors: Leonardo Pedroso, W.P.M.H. (Maurice) Heemels, Mauro Salazar
% Revision history:
%   05/02/2024 - Leonardo Pedroso
%       * Added final publication reference to [1] 
%   24/11/2023 - Leonardo Pedroso
%       * Bug fixes
%   13/03/2023 - Leonardo Pedroso
%       * Initial implementation
% References: 
%   [1] L. Pedroso, W.P.M.H. Heemels and M. Salazar, "Urgency-Aware Routing
%   in Single Origin-Destination Itineraries Through Artificial 
%   Currencies," 2023 62nd IEEE Conference on Decision and Control (CDC), 
%   Singapore, Singapore, 2023, pp. 4142-4149, 
%   doi: 10.1109/CDC49753.2023.10383739.

%% Check feasibility
% checks if Karma level is greater than or equal to a certain threshold computed on the maximum discomfort and price values
% otherwise it returns NaN as solution

n = length(d);
if k < max([0, k_ij(n,n,k_ref,sort(p,1,'descend'),T)])
    j_star = nan;
    return;
end

%% Reduce to decision that satisfies the following assumptions
% d: d_1 < d_2 < ... < d_n
% p: p_1 > p_2 > ... > p_n
% 1. Sort discomforts and prices in ascending order of discomfort
[d_srt,d_srt_idx] = sort(d);
p_srt = p(d_srt_idx);

% 2. See if there are equal discomforts and select corresponding prices
[d_srt_neq, d_srt_neq_idx, ~] = unique(d_srt,'stable');
p_srt_neq = p_srt(d_srt_neq_idx);

% 3. Take out links p_j such that p_j >= p_i and d_j > d_i
% These will never be solutions
idx_prev_dp = 1;
mask_d_set_neq_idx = false(length(d_srt_neq),1);
mask_d_set_neq_idx(1) = true;
for l = 2:length(d_srt_neq)
    if (p_srt_neq(l) < p_srt_neq(idx_prev_dp))
        idx_prev_dp = l;
        mask_d_set_neq_idx(l) = true;
    end
end
d_srt_neq_dp = d_srt_neq(mask_d_set_neq_idx);
p_srt_neq_dp = p_srt_neq(mask_d_set_neq_idx);

%% Solve reduced individual problem
% Reduced dimension
n_red = length(d_srt_neq_dp);
% Get solution thresholds
gamma = n_arcs_gamma_individual(d_srt_neq_dp,T,p_srt_neq_dp,k,k_ref,s_min,s_bar,s_max);
% Check region of the solution based on the sensitivity  
for j = 1:n_red
    if j ~= 1
        if s/s_bar > gamma(j,1) && s/s_bar <= gamma(j-1,1)
            j_star_srt_neq_dp = j;
        end
    else
        if s/s_bar > gamma(j,1)
            j_star_srt_neq_dp = j;
        end
    end
end

%% Unwind solution
% Maps the solution index back to the original discomfort vector, i.e.,
% before removal of impossible decisions 
% (p_j such that p_j = p_i and d_j > d_i)
j_star_srt_neq = find(d_srt_neq == d_srt_neq_dp(j_star_srt_neq_dp));
% Possibly equal discomfort
% If more than one arc has the same discomfort, then choose the one
% that is typically fastest and still feasible (as a draw criterion)
j_star_srt = find(d_srt == d_srt_neq(j_star_srt_neq));
if length(j_star_srt) > 1
    for l = 1:length(j_star_srt)
        if k >= k_ref+p_srt(j_star_srt(l))+T*min(p)
            j_star_srt = j_star_srt(l);
            break;
        end
    end
end
% Match with original arc labels 
j_star = d_srt_idx(j_star_srt);
if isnan(j_star)
    2
end
end

%% Auxiliary functions
% Karma thresholds for unitary decisions, i.e., [y_bar]_j = 1 for some j
function k_ij = k_ij(i,j,k_ref,p,T)
    k_ij = k_ref+p(i)+T*p(j);
end