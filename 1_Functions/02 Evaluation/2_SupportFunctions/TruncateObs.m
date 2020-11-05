function [xobs] = TruncateObs(obs, tComp)
%   TRUNCATEOBS This function truncates the observations to the selected timeframe
%
%   This function is launched in the c_BAEvaluation.m script
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   Inputs:
%       obs: matrix with the observation data. First three columns are YMD,
%       last three columns are ETP
%       tComp: Selected timeframe to evaluate for. Should be the same as the future timeframe
%       for bias correction
%   Outputs:
%       xobs: truncated (in number of rows) observation matrix
%
% Last update by J. Van de Velde on 09/12/'19

%% Truncation

tStart = find(obs(:,1) == tComp(1,1) & obs(:,2) == tComp(1,2) & obs(:,3) == tComp(1,3)); %Index of the start of the control period: year, month and day must be the same for both tControl and OC
tEnd = find(obs(:,1) == tComp(2,1) & obs(:,2) == tComp(2,2) & obs(:,3) == tComp(2,3)); %Index of the end of the control period: year, month and day must be the same for both tControl and OC
xobs = obs(tStart:tEnd,:);

end

