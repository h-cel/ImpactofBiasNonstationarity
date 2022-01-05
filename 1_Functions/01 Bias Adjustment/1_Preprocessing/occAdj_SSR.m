function [outho, ouths, outfs, th] = occAdj_SSR(xho, xhs, xfs)
%   occAdj_SSR This function implements a precipitation occurrence bias adjustment 
%   method using the stochastic singularity removal technique
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress) 
%
%   This function implements the stochastic singularity removal as
%   proposed by Vrac, Noël and Vautard (2016) in "Bias correction of
%   precipitation through Singularity Stochastic Removal: Because
%   occurences matter"
%
%   Inputs:
%       xfs: the future simulations, to be corrected, a n x 4-matrix with the respective
%       variable in the last column and Y M D in the first 3 columns
%       xho: the historical observations, idem
%       xhs: the historical simulations, idem
%   Outputs:
%       outho: matrix of the changed historical observations
%       ouths: matrix of the changed historical simulations
%       outfs: matrix of the changed future simulations
%       th: threshold calculated in this function
%
%   Last update by J. Van de Velde on 09/12/'19

%% Implementation

% Set-up
ndays = length(xho);

% Determine th
% th is determined by the smallest wet day precipitation depth

th = min([min(xho(xho(:,4)>0.1, 4)), min(xhs(xhs(:, 4)>0.1, 4)), min(xfs(xfs(:, 4)>0.1, 4))]);

% Loop over the days 

for i= 1:ndays
    if xho(i, 4) < th
        randomnmbr = rand(1);
        xho(i,4) = randomnmbr*(th-0.1)+0.1; %Random selection of a day out of the interval between the threshold and 0.1
    end
    if xhs(i, 4) < th
        randomnmbr = rand(1);
        xhs(i, 4) = randomnmbr*(th-0.1)+0.1;
    end
    if xfs(i,4) < th
        randomnmbr = rand(1);
        xfs(i,4) = randomnmbr*(th-0.1)+0.1;
    end
end

outfs = xfs(:, 4);
ouths = xhs(:, 4);
outho = xho(:, 4);
end

