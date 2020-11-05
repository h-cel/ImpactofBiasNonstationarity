function out = T(x,b)
%   T Implements the CDF of the triangular distribution
%
%   This function is launched in the occAdj_TDA.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   Inputs:
%       x: quantile for which the CDF value has to be calculated
%       b: parameter of the triangular distribution
%   Output:
%       out: CDF value of the triangular distribution
%
%   Last update by J. Van de Velde on 12/09/'19

%% Implementation

if x < b %Upper limit
    out = 1 - (b-x)^2 / b^2;
else
    out = 1;
end

end

