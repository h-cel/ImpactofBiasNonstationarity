function [cost] = DistEucl(ci, cj)
%   DISTEUCL This function calculates the distance cost matrix between two matrices
%
%   This function is launched in the dOTC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   Inputs:
%       ci: a variable space partition for time series i
%       cj: a variable space partition for time series j
%   Output:
%       cost: the cost associated with transforming the space partition of
%       time series i to time series j  
%
%   Last update by J. Van de Velde on 26/11/'19

%% Set up

nrowsi = size(ci, 1);
nrowsj = size(cj, 1);
cost = zeros(nrowsi, nrowsj);

%% Calculation

for i= 1:nrowsi
    for j = 1: nrowsj
        cost(i,j)= norm((cj(i,:)-ci(j,:)))^2;
    end
end



end

