function Dist = EnDist(simulation, target)
%   ENDIST This function calculates the energy distance between to
%   matrices.
%
%   This function is launched in the MBCn.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   The calculation of the energy distance is based on the research by
%   Szekely and Rizzo (2013). In this function, the use is assumed to
%   compare simulated climate variables (a matrix of T1 timesteps and N variables) 
%   with historical observations (a matrix of T2 timesteps and N
%   variables). T1 and T2 are often the same size.
%
%   Inputs:
%       simulation: a matrix of size T1 x N
%       target: a matrix of size T2 x N
%   Outputs:
%       Dist : the energy distance, a scalar
%
%   Last update by J. Van de Velde on 25/11/'19

%% Set-up

[n1, nvars] = size(target);
n2 = size(simulation, 1);
XY = zeros(1, n1*n2);
XX = zeros(1, n1*n1);
YY = zeros(1, n2*n2);

%% Standardization

simulation_st = zeros(n2,nvars);
target_st = zeros(n1, nvars);

simulation_st(:,1) = (simulation(:,1)-nanmean(simulation(:,1)))/nanstd(simulation(:,1));
target_st(:,1) = (target(:,1)-nanmean(target(:,1)))/nanstd(target(:,1));

if nvars > 1
    for i= 2:nvars
        simulation_st(:,i) = (simulation(:,i)-nanmean(simulation(:,i)))/nanstd(simulation(:,i));
        target_st(:,i) = (target(:,i)-nanmean(target(:,i)))/nanstd(target(:,i));
    end
end

%% XY
% Correlation between target and simulation

cnt = 1;
for i=1:n1
    for j = 1:n2
        XY(cnt) = norm(target_st(i,:)-simulation_st(j,:));
        cnt = cnt+1;
    end
end

sumXY = 2/(n1*n2)*nansum(XY);

%% XX
% Correlation within target

cnt = 1;
for i=1:n1
    for j = 1:n1
        XX(cnt) = norm(target_st(i,:)-target_st(j,:));
        cnt = cnt+1;
    end
end

sumXX = 1/(n1^2)*nansum(XX);

%% YY
% Correlation within simulation

cnt = 1;
for i=1:n2
    for j = 1:n2
        XX(cnt) = norm(simulation_st(i,:)-simulation_st(j,:));
        cnt = cnt+1;
    end
end

sumYY = 1/(n2^2)*nansum(YY);

%% Result

Dist = sumXY-sumXX-sumYY;
        

end

