function [val20y] = val20yRP(data)
%   VAL20YRP This function calculates the 20 year return period value
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function calculates the flow value that corresponds to the 20 year
%   return period
%
%   Input:
%       data: discharge time series
%   Output:
%       val20y: 20 year return period value
%
%   Last update by J. Van de Velde on 09/12/'19

%% Yearly maximum

ystart = data(1, 1);
yend = data(end, 1);
nyears = yend-ystart + 1;
yearmax = zeros(nyears, 1);
cnt = 0;

for y = ystart:yend
    cnt = cnt+1;
    datay = data(data(:,1) == y, end);
    yearmax(cnt) = max(datay);
end

%% Ranks

max_sorted = sort(yearmax, 'descend');
[~, max_rnk] = ismember(yearmax,max_sorted);

%% Return periods

T = (nyears + 1)./max_rnk;

%% Reduced variate

yT = -log(log(T./(T-1)));

%% 

p = polyfit(yT, yearmax, 1);

%% Value for 20 years

yT20 = -log(log(20/(20-1)));

val20y = p(2) + p(1)*yT20;

end


