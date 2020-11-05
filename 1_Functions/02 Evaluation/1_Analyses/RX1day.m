function [PIndex] = RX1day(P)
%   RX1DAY This function calculates the monthly maximum 1-day P
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function calculates the monthly maximum 1-day P, one of
%   the indices by the ETCCDI (Zhang et al. (2011))
%
%   Input:
%       P: precipitation matrix with Y:M:D:Value as columns
%   Output:
%       PIndex: vector with the value for every month
%
%   Last update by J. Van de Velde on 09/12/'19


%% Set-up

Ystart = P(1,1);
Yend = P(end, 1);

nMonths =  (Yend-Ystart+1)*12;

PIndex = zeros(1, nMonths);
cnt = 1;

%% Loop

for i = Ystart:Yend
    for j = 1:12
        PMonth = P(P(:,1) == i & P(:,2) == j,end);
        PIndex(cnt) = max(PMonth);
        cnt = cnt +1;
    end
end
end

