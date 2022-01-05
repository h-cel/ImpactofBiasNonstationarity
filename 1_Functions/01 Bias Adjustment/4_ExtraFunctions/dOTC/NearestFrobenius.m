function [Anew] = NearestFrobenius(A)
%   NEARESTFROBENIUS Estimation of a positive semidefinite matrix
%
%   This function is launched in the dOTC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function calculates the nearest positive semidefinite matrix in
%   the Frobenius norm to an arbitrary real matrix, according to the
%   algorithm given by Higham (1988)
%
%   Input:
%       A: a non-positive-definite matrix
%   Output
%       Anew: an estimated positive semidefinite matrix
%
%   Last update by J. Van de Velde on 26/11/'19

%% Calculation
% Calculation of B

B = (A+A')/2;

% Check for singularity

if sum(isinf(B))> 0 | sum(isnan(B))> 0
    Anew = NaN;
    return
end

if cond(B) == Inf
    Anew = NaN;
    return
end

%Spectral decomposition (as B is symmetric, a singular value decomposition
%should suffice)

[U, S, V] = svd(B);

% Calculation of di

for i =1:size(S,1)
    if S(i,i) <0
        S(i,i) = 0;
    end
end

Anew = U*S*V';

end

