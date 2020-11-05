function [results] = TransProb(Precip, threshold)
%   TRANSPROB This function calculates the transition probabilities for a
%   precipitation time series
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function calculates the different transition probabilities for a
%   daily precipitation time series: dry-dry, dry-wet, wet-dry, wet-wet
%
%   Inputs:
%       Precip: daily precipitation time series
%       threshold: threshold to define dry days
%   Output:
%       results: vector with the different transition probabilities
%
%   Last update by J. Van de Velde on 09/12/'19

%% Calculation

%Initiatie

P00 = 0;
P01 = 0;
P10 = 0;
P11 = 0;


%Loop

for i=2: length(Precip)
        if Precip(i-1) < threshold && Precip(i) < threshold 
            P00 = P00 + 1;
        elseif Precip(i-1) >= threshold && Precip(i) < threshold
            P10 = P10 + 1;
        elseif Precip(i-1) < threshold && Precip(i) >= threshold
            P01 = P01 + 1;
        elseif Precip(i-1) >= threshold && Precip(i) >= threshold
            P11 = P11 + 1;
        end
   
end

%Probabilities

p00 = P00/sum(Precip < threshold);
p01 = P01/sum(Precip < threshold);
p10 = P10/sum(Precip >= threshold);
p11 = P11/sum(Precip >= threshold);

results = [p00, p01, p10, p11];

end

