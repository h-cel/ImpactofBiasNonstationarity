function [xfa_def] = MRQNBC(xho, xhs, xfs)
%   MRQNBC This function implements the MRQNBC bias adjustment method
%
%   This function is launched in the BiasAdjusment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function implements the Multivariate Recursive Quantile Nesting
%   Bias Correction algorithm, as proposed by Mehrotra and Sharma (2016)
%
%   Inputs:
%       xho: the observational dataset, a n x 6-matrix with the variables 
%       in the last column and Y M D in the first 3 columns
%       xhs: the historical simulation dataset, idem
%       xfs: the future simulation dataset, idem
%   Outputs:
%       xfa_def: adjusted future dataset
%
%   Last update by J. Van de Velde on 09/12/'19

%% Set up

xfa_def = xfs; % To enable the use of the final results in the loop

%% Main loop
% Literature suggests 3-5 iterations to have best results
% But I removed this, as this kept increasing mistakes (singular matrices) -> z= 1:1

for z= 1:1

%% Daily correction

% Retrieving columns

time_h = xho(:,1:3);
time_f = xfa_def(:,1:3);

Eho = [time_h, xho(:, 4)];
Tho = [time_h, xho(:, 5)]; 
Pho = [time_h, xho(:, 6)]; 
    
Ehs = [time_h, xhs(:, 4)]; 
Ths = [time_h, xhs(:, 5)]; 
Phs = [time_h, xhs(:, 6)]; 
    
Efs = [time_f, xfa_def(:, 4)]; 
Tfs = [time_f, xfa_def(:, 5)]; 
Pfs = [time_f, xfa_def(:, 6)]; 

% QDM

[~, Pfa] = QDM(Pho,Phs,Pfs,2);
[~, Efa] = QDM(Eho,Ehs,Efs,2);
[~, Tfa] = QDM(Tho,Ths,Tfs,1);

% Standardization

xfa_daily_st = [(Efa-mean(Efa))./std(Efa) (Tfa-mean(Tfa))./std(Tfa) (Pfa-mean(Pfa))./std(Pfa)];
xho_st = [(Eho(:,end)-mean(Eho(:,end)))./std(Eho(:,end)) (Tho(:,end)-mean(Tho(:,end)))./std(Tho(:,end)) (Pho(:,end)-mean(Pho(:,end)))./std(Pho(:,end))];

% Matrices

[C, D] = coeff(xho_st);
[E, F] = coeff(xfa_daily_st);

% Correction

nrows = size(xfa_daily_st, 1);
xfa_daily = xfa_daily_st;

for i = 2:nrows
    xfa_daily(i,:) = C*xfa_daily(i-1,:)' + D/F*xfa_daily_st(i,:)' - D/F*E*xfa_daily_st(i-1,:)';
end

% Adding mean and variance

xfa_daily = [xfa_daily(:,1).*std(Eho(:,end))+mean(Eho(:,end)), xfa_daily(:,2).*std(Tho(:,end))+mean(Tho(:,end)), xfa_daily(:,3).*std(Pho(:,end))+mean(Pho(:,end))];

% Aggregation - future

ystart = time_f(1,1);
yend = time_f(end,1);
nyears = yend-ystart+1;
n = nyears*12;
xfs_monthly = zeros(n,5);
cnt = 0;

for y=ystart:yend
    for m=1:12
        cnt = cnt+1;
        ndays = sum(time_f(:,1) == y & time_f(:,2) == m);
        sumE = sum(xfa_daily(time_f(:,1) == y & time_f(:,2) == m, 1));
        sumT = sum(xfa_daily(time_f(:,1) == y & time_f(:,2) == m, 2))/ndays;
        sumP = sum(xfa_daily(time_f(:,1) == y & time_f(:,2) == m, 3));
        xfs_monthly(cnt, :) = [y, m, sumE, sumT , sumP];
    end
end

% Aggregation - observations

ystart = time_h(1,1);
yend = time_h(end,1);
nyears = yend-ystart+1;
n = nyears*12;
xho_monthly = zeros(n,5);
xhs_monthly = zeros(n,5);
cnt = 0;

for y=ystart:yend
    for m=1:12
        cnt = cnt+1;
        ndays = sum(time_h(:,1) == y & time_h(:,2) == m);
        sumE = sum(xho(time_h(:,1) == y & time_h(:,2) == m, 4));
        sumT = sum(xho(time_h(:,1) == y & time_h(:,2) == m, 5))/ndays;
        sumP = sum(xho(time_h(:,1) == y & time_h(:,2) == m, 6));
        xho_monthly(cnt, :) = [y, m, sumE, sumT , sumP];
        sumE = sum(xhs(time_h(:,1) == y & time_h(:,2) == m, 4));
        sumT = sum(xhs(time_h(:,1) == y & time_h(:,2) == m, 5))/ndays;
        sumP = sum(xhs(time_h(:,1) == y & time_h(:,2) == m, 6));
        xhs_monthly(cnt, :) = [y, m, sumE, sumT , sumP];
    end
end

timef_daily = time_f;

%% Monthly correction

% Retrieving columns

timeh = xho_monthly(:, 1:2);
timef = xfs_monthly(:, 1:2);

Eho = [timeh, xho_monthly(:,3)];
Tho = [timeh, xho_monthly(:,4)]; 
Pho = [timeh, xho_monthly(:,5)]; 
    
Ehs = [timeh, xhs_monthly(:,3)]; 
Ths = [timeh, xhs_monthly(:,4)]; 
Phs = [timeh, xhs_monthly(:,5)]; 
    
Efs = [timef, xfs_monthly(:,3)]; 
Tfs = [timef, xfs_monthly(:,4)]; 
Pfs = [timef, xfs_monthly(:,5)]; 

% QDM

[~, Pfa] = QDM_all(Pho,Phs,Pfs, 1);
[~, Efa] = QDM_all(Eho,Ehs,Efs, 1);
[~, Tfa] = QDM_all(Tho,Ths,Tfs, 0);

% Standardization

xfa_monthly_st = [(Efa-mean(Efa))./std(Efa) (Tfa-mean(Tfa))./std(Tfa) (Pfa-mean(Pfa))./std(Pfa)];
xho_monthly_st = [(Eho(:,end)-mean(Eho(:,end)))./std(Eho(:,end)) (Tho(:,end)-mean(Tho(:,end)))./std(Tho(:,end)) (Pho(:,end)-mean(Pho(:,end)))./std(Pho(:,end))];

% Matrices

xfa_monthly = [timef xfa_monthly_st];
nrows = size(xfa_monthly, 1);

% Periodic parameters
[C, D] = coeffPeriodic(timeh, xho_monthly_st, 12);
[E, F] = coeffPeriodic(timef, xfa_monthly_st, 12);

for i = 2:nrows
    m = xfa_monthly(i,2);
    xfa_monthly(i,3:5) = (C(:,:,m)*xfa_monthly(i-1,3:5)'+D(:,:,m)*inv(F(:,:,m))*xfa_monthly_st(i,:)' -D(:,:,m)*inv(F(:,:,m))*E(:,:,m)*xfa_monthly_st(i-1,:)')';
end

% Adding mean and variance

xfa_monthly(:, 3:5) = [xfa_monthly(:,3)*std(Eho(:,end))+mean(Eho(:,end)), xfa_monthly(:,4)*std(Tho(:,end))+mean(Tho(:,end)), xfa_monthly(:,5)*std(Pho(:,end))+mean(Pho(:,end))];

xfa_monthly(xfa_monthly(:,3)<0.1,3) = 0;
xfa_monthly(xfa_monthly(:,5)<0.1,5) = 0;

% Aggregation - future

ystart = time_f(1,1);
yend = time_f(end,1);
nyears = yend-ystart+1;
n = nyears*4;
xfs_season = zeros(n,5);
cnt =  0;

for y=ystart:yend
    for m=1:4
        cnt= cnt+1;
        nmonths = sum(timef(:,1) == y & timef(:,2) <=(m*3) & timef(:,2)>(m*3)-3);
        sumE = sum(xfa_monthly(timef(:,1) == y & timef(:,2) <=(m*3) & timef(:,2)>(m*3)-3, 3));
        sumT = sum(xfa_monthly(timef(:,1) == y & timef(:,2) <=(m*3) & timef(:,2)>(m*3)-3, 4))/nmonths;
        sumP = sum(xfa_monthly(timef(:,1) == y & timef(:,2) <=(m*3) & timef(:,2)>(m*3)-3, 5));
        xfs_season(cnt, :) = [y, m, sumE, sumT , sumP];
    end
end

% Aggregation - historic

ystart = time_h(1,1);
yend = time_h(end,1);
nyears = yend-ystart+1;
n = nyears*4;
xho_season = zeros(n,5);
xhs_season = zeros(n,5);
cnt = 0;

for y=ystart:yend
    for m=1:4
        cnt = cnt+1;
        nmonths = sum(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3);
        sumE = sum(xho_monthly(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3, 3));
        sumT = sum(xho_monthly(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3, 4))/nmonths;
        sumP = sum(xho_monthly(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3, 5));
        xho_season(cnt, :) = [y, m, sumE, sumT , sumP];
        sumE = sum(xhs_monthly(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3, 3));
        sumT = sum(xhs_monthly(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3, 4))/nmonths;
        sumP = sum(xhs_monthly(timeh(:,1) == y & timeh(:,2) <=(m*3) & timeh(:,2)>(m*3)-3, 5));
        xhs_season(cnt, :) = [y, m, sumE, sumT , sumP];
    end
end

%% Seasonal correction

 
% Retrieving columns 
timeh = xho_season(:, 1:2);
timef = xfs_season(:, 1:2);

Eho = [timeh, xho_season(:,3)];
Tho = [timeh, xho_season(:,4)]; 
Pho = [timeh, xho_season(:,5)]; 
    
Ehs = [timeh, xhs_season(:,3)]; 
Ths = [timeh, xhs_season(:,4)]; 
Phs = [timeh, xhs_season(:,5)]; 
    
Efs = [timef, xfs_season(:,3)]; 
Tfs = [timef, xfs_season(:,4)]; 
Pfs = [timef, xfs_season(:,5)]; 

%QDM

[~, Pfa] = QDM_all(Pho,Phs,Pfs, 1);
[~, Efa] = QDM_all(Eho,Ehs,Efs, 1);
[~, Tfa] = QDM_all(Tho,Ths,Tfs, 0);

%Standardization

xfa_season_st = [(Efa-mean(Efa))./std(Efa) (Tfa-mean(Tfa))./std(Tfa) (Pfa-mean(Pfa))./std(Pfa)];
xho_season_st = [(Eho(:,end)-mean(Eho(:,end)))./std(Eho(:,end)) (Tho(:,end)-mean(Tho(:,end)))./std(Tho(:,end)) (Pho(:,end)-mean(Pho(:,end)))./std(Pho(:,end))];

%Matrices

xfa_season = [timef xfa_season_st];
nrows = size(xfa_season, 1);

%Periodic parameters
[C, D] = coeffPeriodic(timeh, xho_season_st, 4);
[E, F] = coeffPeriodic(timef, xfa_season_st, 4);

for i = 2:nrows
    m = xfa_season(i,2);
    xfa_season(i,3:5) = C(:,:,m)*xfa_season(i-1,3:5)'+D(:,:,m)/(F(:,:,m))*xfa_season_st(i,:)' -D(:,:,m)/F(:,:,m)*E(:,:,m)*xfa_season_st(i-1,:)';
end

%Adding mean and variance

xfa_season(:, 3:5) = [xfa_season(:,3)*std(Eho(:,end))+mean(Eho(:,end)), xfa_season(:,4)*std(Tho(:,end))+mean(Tho(:,end)), xfa_season(:,5)*std(Pho(:,end))+mean(Pho(:,end))];

xfa_season(xfa_season(:,3)<0.1,3) = 0;
xfa_season(xfa_season(:,5)<0.1,5) = 0;

% Aggregation - future

ystart = timef(1,1);
yend = timef(end,1);
nyears = yend-ystart+1;
xfs_year = zeros(nyears,4);
cnt = 0;

for y=ystart:yend
    cnt = cnt+1;
    nmonths = sum(xfa_season(:,1) == y);
    sumE = sum(xfa_season(timef(:,1) == y, 3));
    sumT = sum(xfa_season(timef(:,1) == y, 4))/nmonths;
    sumP = sum(xfa_season(timef(:,1) == y, 5));
    xfs_year(cnt, :) = [y, sumE, sumT , sumP];
end

% Aggregation - historic

ystart = timeh(1,1);
yend = timeh(end,1);
nyears = yend-ystart+1;
n = nyears;
xho_year = zeros(n,4);
xhs_year = zeros(n,4);
cnt = 0;

for y=ystart:yend
    cnt = cnt+1;
    nmonths = sum(xho_season(:,1) == y);
    sumE = sum(xho_season(timeh(:,1) == y, 3));
    sumT = sum(xho_season(timeh(:,1) == y, 4))/nmonths;
    sumP = sum(xho_season(timeh(:,1) == y, 5));
    xho_year(cnt, :) = [y, sumE, sumT , sumP];
    sumE = sum(xhs_season(timeh(:,1) == y, 3));
    sumT = sum(xhs_season(timeh(:,1) == y, 4))/nmonths;
    sumP = sum(xhs_season(timeh(:,1) == y, 5));
    xhs_year(cnt, :) = [y, sumE, sumT , sumP];
end

%% Yearly correction

% Retrieving columns

timeh = xho_year(:, 1);
timef = xfs_year(:, 1);

Eho = [timeh, xho_year(:,2)];
Tho = [timeh, xho_year(:,3)]; 
Pho = [timeh, xho_year(:,4)]; 
    
Ehs = [timeh, xhs_year(:,2)]; 
Ths = [timeh, xhs_year(:,3)]; 
Phs = [timeh, xhs_year(:,4)]; 
    
Efs = [timef, xfs_year(:,2)]; 
Tfs = [timef, xfs_year(:,3)]; 
Pfs = [timef, xfs_year(:,4)]; 

% QDM

[~, Pfa] = QDM_all(Pho,Phs,Pfs, 1);
[~, Efa] = QDM_all(Eho,Ehs,Efs, 1);
[~, Tfa] = QDM_all(Tho,Ths,Tfs, 0);

% Standardization

xfa_year_st = [(Efa-mean(Efa))/std(Efa) (Tfa-mean(Tfa))/std(Tfa) (Pfa-mean(Pfa))/std(Pfa)];
xho_year_st = [(Eho(:,end)-mean(Eho(:,end)))/std(Eho(:,end)) (Tho(:,end)-mean(Tho(:,end)))/std(Tho(:,end)) (Pho(:,end)-mean(Pho(:,end)))/std(Pho(:,end))];

% Matrices

[C, D] = coeff(xho_year_st);
[E, F] = coeff(xfa_year_st);

% Correction

nrows = size(xfa_year_st, 1);
xfa_year = xfa_year_st;

for i = 2:nrows
    xfa_year(i,:) = C*xfa_year(i-1,:)' + D/F*xfa_year_st(i,:)' - D/F*E*xfa_year_st(i-1,:)';
end

% Adding mean and variance

xfa_year = [xfa_year(:,1)*std(Eho(:,end))+mean(Eho(:,end)), xfa_year(:,2)*std(Tho(:,end))+mean(Tho(:,end)), xfa_year(:,3)*std(Pho(:,end))+mean(Pho(:,end))];

%% Weighting factors

nrows = size(xfa_daily, 1);
xfa_def = [timef_daily, xfa_daily];

for i =1:nrows
    y = timef_daily(i,1);
    m = timef_daily(i,2);
    s = ceil(timef_daily(i,2)/3);
    %Weighting factor calculation
    weight_year = xfa_year(xfs_year(:,1) == y, :)/xfs_year(xfs_year(:,1) == y, 2:4);
    weight_season = xfa_season(xfa_season(:,1) == y & xfa_season(:,2) == s, 3:5)/xfs_season(xfs_season(:,1) == y & xfs_season(:,2) == s, 3:5);
    weight_month = xfa_monthly(xfa_monthly(:,1) == y & xfa_monthly(:,2) == m, 3:5)/xfs_monthly(xfs_monthly(:,1) == y & xfs_monthly(:,2) == m, 3:5);
    %Definitive calculation
    xfa_def(i, 4:6) = weight_year*weight_season*weight_month*xfa_daily(i,:);
end

end

end

