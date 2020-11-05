function [PIndex] = CDD(P)
%   CDD This function calculates the maximum number of consecutive dry days
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function calculates the maximum number of consecutive dry days, one of
%   the indices by the ETCCDI (Zhang et al. (2011))
%
%   Inputs:
%       P: precipitation vector
%   Outputs:
%       PIndex: vector with the value for every year
%
%   Last update by J. Van de Velde on 09/12/'19

%% Set-up

cntb = 0;
cntN = 0;

%% Making b
%Check whether precip starts rainy or dry

if P(1) < 0.1
    start = find(P >= 0.1, 1);
    cntb = cntb+1;
    b(cntb) = start-1;
    cntN = cntN + 1;
else
    start = 1;
end

%Selection loop

for i = start:length(P)
    if P(i) < 0.1 && P(i-1) >= 0.1
        cntb = cntb+1;
        b(cntb) = 1;
        cntN = cntN + 1;
    elseif P(i) < 0.1 && P(i-1) < 0.1
        b(cntb) = b(cntb)+1;
    end 
end

%% Index calculation

PIndex = max(b);

end

