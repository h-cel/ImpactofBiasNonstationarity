function [xfa2] = dOTC(xho, xhs, xfs)
%   DOTC This function implements the dOTC bias adjustment method
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function implements the dOTC bias adjustment method as proposed by
%   Robin et al. (2019). This is a bias adjustment method based on the
%   generalization of the transfer function as an optimal transport
%   problem.
%   
%   Inputs:
%       xho: the observational dataset, a n x 6-matrix with the variables 
%       in the last column and Y M D in the first 3 columns
%       xhs: the historical simulation dataset, idem
%       xfs: the future simulation dataset, idem
%   Outputs:
%       xfa2: final adjusted future dataset
%
%   Last update by J. Van de Velde on 05/11/'20
%% Set-up

% Ranges
% with a slight perturbation to allow values on the boundaries to be used
% in later steps
EMax = max([max(xho(:,1)), max(xhs(:,1)), max(xfs(:,1))])+0.01;
EMin = min([min(xho(:,1)), min(xhs(:,1)), min(xfs(:,1))]);
TMax = max([max(xho(:,2)), max(xhs(:,2)), max(xfs(:,2))])+0.01;
TMin = min([min(xho(:,2)), min(xhs(:,2)), min(xfs(:,2))]);
PMax = max([max(xho(:,3)), max(xhs(:,3)), max(xfs(:,3))])+0.01;
PMin = min([min(xho(:,3)), min(xhs(:,3)), min(xfs(:,3))]);

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

Pho_est = zeros(numcells^3, 1);
Phs_est = zeros(numcells^3, 1);
Pfs_est = zeros(numcells^3, 1);
cho = zeros(numcells^3, 3);
chs = zeros(numcells^3, 3);
cfs = zeros(numcells^3, 3);

% Pre-allocation

xfa2=NaN(size(xfs));

%% Law estimation

% xho
nrows = size(xho,1);
loc = 0;
for i = 1:numcells
    for j = 1:numcells
        for k = 1:numcells
            loc = loc +1;
            cmin = [Ecells(i) Tcells(j), Pcells(k)];
            cmax = [Ecells(i+1) Tcells(j+1), Pcells(k+1)];
            cnt = 0;
            for l = 1:nrows
                if (xho(l,:) >= cmin) == [1 1 1] & (xho(l,:) < cmax) == [1 1 1]
                    cnt = cnt+1;
                end
            end
            Pho_est(loc) = cnt/nrows;
            cho(loc, :) = (cmin + cmax)./2;
        end
    end
end

% xhs
loc = 0;
nrows = size(xhs,1);
for i = 1:numcells
    for j = 1:numcells
        for k = 1:numcells
            loc = loc+1;
            cmin = [Ecells(i) Tcells(j), Pcells(k)];
            cmax = [Ecells(i+1) Tcells(j+1), Pcells(k+1)];
            cnt = 0;
            for l = 1:nrows
                if (xhs(l,:) >= cmin) == [1 1 1] & (xhs(l,:) < cmax) == [1 1 1]
                    cnt = cnt+1;
                end
            end
            Phs_est(loc) = cnt/nrows;
            chs(loc, :) = (cmin + cmax)./2;
        end
    end
end

% xfs
nrows = size(xfs,1);
loc = 0;
for i = 1:numcells
    for j = 1:numcells
        for k = 1:numcells
            loc = loc +1;
            cmin = [Ecells(i) Tcells(j), Pcells(k)];
            cmax = [Ecells(i+1) Tcells(j+1), Pcells(k+1)];
            cnt = 0;
            for l = 1:nrows
                if (xfs(l,:) >= cmin) == [1 1 1] & (xfs(l,:) < cmax) == [1 1 1]
                    cnt = cnt+1;
                end
            end
            Pfs_est(loc) = cnt/nrows;
            cfs(loc, :) = (cmin + cmax)./2;
        end
    end
end

%% Optimal plan gamma (historic time series)

costhist = DistEucl(cho, chs);
Pgamma = OptTransPlan(costhist,Pho_est, Phs_est, 20, 10^-8);

%% Optimal plan phi

costsim = DistEucl(chs, cfs);
Pphi = OptTransPlan(costsim, Phs_est, Pfs_est, 20, 10^-8);

%% Cholesky factor

cxho = cov(xho);
cxhs = cov(xhs);

% Check for positive-definiteness and estimation if not
[~,p1] = chol(cxho);
[~,p2] = chol(cxhs);
if p1 ~= 0
    cxho = NearestFrobenius(cxho);
    if isnan(cxho) == 1
        return
    end
elseif p2 ~= 0
    cxhs = NearestFrobenius(cxhs);
    if isnan(cxhs) == 1
        return
    end
end
% Calculation of D
D = chol(cxho)*inv(chol(cxhs)); 

%% Computation of xfa1

xfa1 = zeros(size(xho));
for i =1:nrows
    %Find a cell of xho
    if round(xho(i,1), 3) == 0 && round(xho(i,3),3) %Some mistakes when comparing numbers which only differ some decimals
        I1 = find(round(cho(:,1) -Ecellsize/2, 2) <= round(xho(i,1),3) & (cho(:,1) +Ecellsize/2) > xho(i,1) & (cho(:,2) -Tcellsize/2) <= xho(i,2) & (cho(:,2) +Tcellsize/2) > xho(i,2) & round(cho(:,3) -Pcellsize/2, 3) <= round(xho(i,3),3) & (cho(:,3) +Pcellsize/2) > xho(i,3));
    elseif round(xho(i,1), 3) == 0 %Some mistakes when comparing numbers which only differ some decimals
        I1 = find(round(cho(:,1) -Ecellsize/2, 2) <= round(xho(i,1),3) & (cho(:,1) +Ecellsize/2) > xho(i,1) & (cho(:,2) -Tcellsize/2) <= xho(i,2) & (cho(:,2) +Tcellsize/2) > xho(i,2) & (cho(:,3) -Pcellsize/2) <= xho(i,3) & (cho(:,3) +Pcellsize/2) > xho(i,3));
    elseif round(xho(i,3), 3) == 0
        I1 = find((cho(:,1) -Ecellsize/2) <= xho(i,1) & (cho(:,1) +Ecellsize/2) > xho(i,1) & (cho(:,2) -Tcellsize/2) <= xho(i,2) & (cho(:,2) +Tcellsize/2) > xho(i,2) & round(cho(:,3) -Pcellsize/2, 3) <= round(xho(i,3),3) & (cho(:,3) +Pcellsize/2) > xho(i,3));
    else
        I1 = find((cho(:,1) -Ecellsize/2) <= xho(i,1) & (cho(:,1) +Ecellsize/2) >= xho(i,1) & (cho(:,2) -Tcellsize/2) <= xho(i,2) & (cho(:,2) +Tcellsize/2) >= xho(i,2) & (cho(:,3) -Pcellsize/2) <= xho(i,3) & (cho(:,3) +Pcellsize/2) >= xho(i,3));
    end
    if size(I1, 1) ~= 1
        continue
    end

    %Find a cell of xhs
    condlaw = Pgamma(:, I1)./Pho_est(I1);
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
    
    
    condlaw = Pphi(:, j)./Phs_est(j);
    kdef = 0;
    while kdef == 0
        k = randi(length(condlaw), 1);
        z= rand(1);
        if z <= condlaw(k) 
            kdef = k;
        else
            kdef = 0;
        end
    end
    
    vect = cfs(k,:)- chs(j,:);
    xfa1(i, :) = xho(i,:)+(D*vect')';
end

%% Applying OTC on xfa1

xfa2 = OTC(xfs,xfa1);
xfa2(xfa2(:,end)<0) = 0;




end

