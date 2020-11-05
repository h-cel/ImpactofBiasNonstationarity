function [Spelllengths] = SpellDist(P)
%   SPELLDIST This function calculates the ecdf of the spell length
%   distribution
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   Input:
%       P: vector with precipitation data
%   Output:
%       xecdfSpell: ecdf of the wet spell length
%
%   Last update by J. Van de Velde on 09/12/'19

%% Vector of spell lengths

%Set-up
cnt = 1;
X(cnt) = 0;

%Define start
if P(1) < 0.1
    start = find(P >= 0.1, 1);
else
    start = 1;
end

%Selection loop

for i = start:length(P)
    if P(i) > 0.1 
        X(cnt) = X(cnt) + 1;
    elseif P(i) < 0.1 && P(i-1) > 0.1
        cnt = cnt +1;
        X(cnt) = 0;
    end 
end

%% Result

Spelllengths = X;

end

