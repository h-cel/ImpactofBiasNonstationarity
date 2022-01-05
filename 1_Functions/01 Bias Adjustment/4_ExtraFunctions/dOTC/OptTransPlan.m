function [Pij] = OptTransPlan(cost,r, c, lambda, epsilon)
%   OPTTRANSPLAN This function calculates the optimal transport plan using the
%   Sinkhorn-Knopp algorithm
%
%   This function is launched in the dOTC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   Inputs:
%       cost: the cost associated with transforming the space partition of
%       time series i to time series j  
%       r: optimal transport plan of time series i
%       c: optimal transport plan of time series j
%       lambda: parameter for the calculation
%       epsilon: tolerance
%   Output:
%       Pij: optimal transport plan between time series i and j
%
%   Last update by J. Van de Velde on 26/11/'19

%% Implementation

u = ones(length(r),1);
v = ones(length(c),1);
K = exp(-cost/lambda);

uprev = zeros(size(u));

while sum(abs(u-uprev))> epsilon
    uprev = u;
    u = r./(K*v);
    v = c./(K'*u);
end

Pij = sparse(diag(v)*K*diag(u));
end

