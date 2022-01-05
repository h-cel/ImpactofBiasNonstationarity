function [score] = Perkins(tsmod, tsor, nrepetitions)
%PERKINS Calculation of Perkins' Skill Score (Perkins et al., 2007)
%   This function calculates the Perkins' Skill Score (PSS), a measure for the
%   similarity between PDFs (Perkins et al., 2007). In this code, this
%   measure is used to compare modelled and original time series, but it
%   can be used for other aspects as well.
%   
%   Inputs:
%       tsmod: modelled time series 
%       tsor: original time series
%       nrepetitions: number of repetitions used, the final score will be an
%       average of n repetitions
%   Output:
%       score: the PSS

%% Calculation

scoretmp = zeros(nrepetitions, 1);

for n=1:nrepetitions
    
    [Nho, edgesho] =  histcounts(tsor(:,end), 'Normalization', 'probability');
    
    %expanding edges
    if edgesho(1)>min(tsmod(:,end,n))
        edgesho(1) = min(tsmod(:,end,n));
    end
    if edgesho(end)<max(tsmod(:,end,n))
        edgesho(end) = max(tsmod(:,end,n));
    end
    
    [Nfa] = histcounts(tsmod(:,end,n), edgesho, 'Normalization', 'probability');
    
     
    scoretmp(n) = sum(min(Nho,Nfa));

end

score= mean(scoretmp);

end

