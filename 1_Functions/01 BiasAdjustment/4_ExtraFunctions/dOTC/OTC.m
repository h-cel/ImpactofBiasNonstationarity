function [xfc2] = OTC(xfs,xfc1)
%   OTC Implements optimal transport correction
%
%   This function is launched in the dOTC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function implements the OTC method as proposed by
%   Robin et al. (2019). This is bias adjustment method that is part of the
%   dOTC bias adjustment method.
%
%   Inputs:
%       xfs: original future time series, a n x 6-matrix with the variables 
%       in the last column and Y M D in the first 3 columns 
%       xfc1: first version of the corrected future time series, idem
%   Output:
%       xfc2: final version of the corrected future time series
%
%   Last update by J. Van de Velde on 26/11/'19

%% Set-up

% Ranges
EMax = max([max(xfs(:,1)), max(xfc1(:,1))])+0.01;
EMin = min([min(xfs(:,1)), min(xfc1(:,1))]);
TMax = max([max(xfs(:,2)), max(xfc1(:,2))])+0.01;
TMin = min([min(xfs(:,2)), min(xfc1(:,2))]);
PMax = max([max(xfs(:,3)), max(xfc1(:,3))])+0.01;
PMin = min([min(xfs(:,3)), min(xfc1(:,3))]);

% Making the partitions
numcells = 25; %Number of cells per variable, can also be an input variable

Pcells = linspace(PMin, PMax, numcells +1); % Numcells+1: these are the ranges of the cells, not the centers
Ecells = linspace(EMin, EMax, numcells +1);
Tcells = linspace(TMin, TMax, numcells +1);

% Cell sizes

Ecellsize = (EMax-EMin)/numcells;
Pcellsize = (PMax-PMin)/numcells;
Tcellsize = (TMax-TMin)/numcells;

% Law matrices

Phs_est = zeros(numcells^3, 1);
Pfc1_est = zeros(numcells^3, 1);
chs = zeros(numcells^3, 3);
cfc1 = zeros(numcells^3, 3);

%% Law estimation

% xhs
loc = 0;
nrows = size(xfs,1);
for i = 1:numcells
    for j = 1:numcells
        for k = 1:numcells
            loc = loc+1;
            cmin = [Ecells(i) Tcells(j), Pcells(k)];
            cmax = [Ecells(i+1) Tcells(j+1), Pcells(k+1)];
            cnt = 0;
            for l = 1:nrows
                if (xfs(l,:) >= cmin) == [1 1 1] & (xfs(l,:) < cmax) == [1 1 1]
                    cnt = cnt+1;
                end
            end
            Phs_est(loc) = cnt/nrows;
            chs(loc, :) = (cmin + cmax)./2;
        end
    end
end

% xfs
nrows = size(xfc1,1);
loc = 0;
for i = 1:numcells
    for j = 1:numcells
        for k = 1:numcells
            loc = loc +1;
            cmin = [Ecells(i) Tcells(j), Pcells(k)];
            cmax = [Ecells(i+1) Tcells(j+1), Pcells(k+1)];
            cnt = 0;
            for l = 1:nrows
                if (xfc1(l,:) >= cmin) == [1 1 1] & (xfc1(l,:) < cmax) == [1 1 1]
                    cnt = cnt+1;
                end
            end
            Pfc1_est(loc) = cnt/nrows;
            cfc1(loc, :) = (cmin + cmax)./2;
        end
    end
end

%% Optimal plan gamma

costhist = DistEucl(chs, cfc1);
Pgamma = OptTransPlan(costhist,Phs_est, Pfc1_est, 25, 10^-8);

%% Calculation

xfc2 = zeros(size(xfc1));
for i =1:nrows
    %Find cell of xhs
    if round(xfs(i,1), 3) == 0 && round(xfs(i,3),3) %Some mistakes when comparing numbers which only differ some decimals
        I1 = find(round(chs(:,1) -Ecellsize/2, 2) <= round(xfs(i,1),3) & (chs(:,1) +Ecellsize/2) > xfs(i,1) & (chs(:,2) -Tcellsize/2) <= xfs(i,2) & (chs(:,2) +Tcellsize/2) > xfs(i,2) & round(chs(:,3) -Pcellsize/2, 3) <= round(xfs(i,3),3) & (chs(:,3) +Pcellsize/2) > xfs(i,3));
    elseif round(xfs(i,1), 3) == 0 %Some mistakes when comparing numbers which only differ some decimals
        I1 = find(round(chs(:,1) -Ecellsize/2, 2) <= round(xfs(i,1),3) & (chs(:,1) +Ecellsize/2) > xfs(i,1) & (chs(:,2) -Tcellsize/2) <= xfs(i,2) & (chs(:,2) +Tcellsize/2) > xfs(i,2) & (chs(:,3) -Pcellsize/2) <= xfs(i,3) & (chs(:,3) +Pcellsize/2) > xfs(i,3));
    elseif round(xfs(i,3), 3) == 0
        I1 = find((chs(:,1) -Ecellsize/2) <= xfs(i,1) & (chs(:,1) +Ecellsize/2) > xfs(i,1) & (chs(:,2) -Tcellsize/2) <= xfs(i,2) & (chs(:,2) +Tcellsize/2) > xfs(i,2) & round(chs(:,3) -Pcellsize/2, 3) <= round(xfs(i,3),3) & (chs(:,3) +Pcellsize/2) > xfs(i,3));
    else
        I1 = find((chs(:,1) -Ecellsize/2) <= xfs(i,1) & (chs(:,1) +Ecellsize/2) >= xfs(i,1) & (chs(:,2) -Tcellsize/2) <= xfs(i,2) & (chs(:,2) +Tcellsize/2) >= xfs(i,2) & (chs(:,3) -Pcellsize/2) <= xfs(i,3) & (chs(:,3) +Pcellsize/2) >= xfs(i,3));
    end
    if size(I1, 1) ~= 1
        continue
    end
    condlaw = Pgamma(:, I1)./Phs_est(I1);
    jdef = 0;
    while jdef == 0
        j = randi(length(condlaw), 1);
        z = rand(1);
        if z <= condlaw(j)
            jdef = j;
        else
            jdef = 0;
        end
    end
    
    xfc2(i,:) = [rand(1)*randi([-1 1],1)*Ecellsize/2+cfc1(j,1) rand(1)*randi([-1 1],1)*Tcellsize/2+cfc1(j,2) rand(1)*randi([-1 1],1)*Pcellsize/2+cfc1(j,3)];
    
end

end

