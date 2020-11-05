function [PIndex] = SDII(P)
%   SDII This function calculates the Simple Daily Intensity Index
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function calculates the Simple Daily Intensity Index (SDII), one of
%   the indices by the ETCCDI (Zhang et al. (2011))
%
%   Input:
%       P: precipitation matrix with Y:M:D:Value as columns
%   Output:
%       PIndex: vector with the value for every year
%
%   Last update by J. Van de Velde on 09/12/'19

%% Set-up

Ystart = P(1,1);
Yend = P(end, 1);
nYears= Ystart-Yend+1;
PIndex= zeros(1, nYears);

%% Loop

cnt = 1;

for i = Ystart:Yend
    Pyear = P(P(:,1) == i, end);
    nWetDays = sum(Pyear> 0.1); % Andere threshold dan oorspronkelijke ETCCDI!
    TotalP = sum(Pyear(Pyear > 0.1));
    
    PIndex(cnt) = TotalP/nWetDays;
    cnt = cnt + 1;
end

end

