function [xha, xfa] = QDM_all(xho,xhs,xfs, rel)
%   QDM_all This function implements the quantile delta mapping method
%
%   This function is launched in the MRQNBC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress). 
%
%   This function implements the following calculations, used to correct the
%   Model output statistics of GCM/RCM-modelled data. The correction used
%   depends on the variable
%
%   xfa = xfs * Fho^1(Fhs(xhs))/Fhs^-1(Fhs(xhs))
%   xfs = xfo +Fho^1(Ffo(xfo))-Fhs^-1(Ffo(xfo))
%
%   With:
%       xfa: the rescaled future period model output
%       xhs: the original future period model output
%       Fho: the cdf of the control dataset
%       Ffs: the cdf of the model dataset during the future period
%       Fhs: the cdf of the model dataset during the control period
%
%   This calculation corrects the modelled data of the future period,
%   according to the ratio or difference between the control period and modelled control
%   period data.
%
%   Inputs:
%       xho: the data of the control period, a matrix with the respective
%       variable in the last column
%       xhs: the modelled data of the control period, a matrix with the respective
%       variable in the last column
%       xfs: the modelled data of the future period, a matrix with the respective
%       variable in the last column
%       type: 0, 1. Indicates whether non-relative changes or changes for all rainy days (2) have to be
%       used.
%   Outputs:
%       xha: adjusted historical time series
%       xfa: adjusted future time series
%
%   Last update by J. Van de Velde on 03/03/'20

%% QDM

if rel == 0 %Non-Relative changes. As non-relative changes are never frequency-corrected, there is no need for wet-day-selection


    %Initialization
    xfa = nan(size(xfs(:,end)));
    xha = nan(size(xhs(:,end)));
    
    %Making ecdf's
    [Fho_ecdf, xho_ecdf] = ecdf(xho(:,end));
    [Fhs_ecdf, xhs_ecdf] = ecdf(xhs(:,end));
    [Ffs_ecdf, xfs_ecdf] = ecdf(xfs(:,end));
    for i = 1:size(xfs,1)
        tmp = max(Ffs_ecdf(xfs_ecdf <= xfs(i,end)));  %Calculation of Ffs_ecdf(xfsm)
        num = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Ffs_ecdf(xfsm))
        denum = interp1(Fhs_ecdf,xhs_ecdf,tmp);  %Calculation of Fhs_ecdf^1(Ffs_ecdf(xfsm))
        xfa(i) = xfs(i,end) + num - denum; %Equidistant calculation
    end
    for i= 1:size(xhs,1)
        tmp = max(Fhs_ecdf(xhs_ecdf <= xhs(i,end))); %Calculation of Fhs_ecdf(xhsm)
        xha(i) = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Fhs_ecdf(xhsm))
    end
    
else
    %Initialization
    xfa = nan(size(xfs(:,end)));
    xha = nan(size(xhs(:,end)));
    
    %Making ecdf's
    [Fho_ecdf, xho_ecdf] = ecdf(xho(:,end));
    [Fhs_ecdf, xhs_ecdf] = ecdf(xhs(:,end));
    [Ffs_ecdf, xfs_ecdf] = ecdf(xfs(:,end));
    for i = 1:size(xfs,1)
        tmp = max(Ffs_ecdf(xfs_ecdf <= xfs(i,end)));  %Calculation of Ffs_ecdf(xfsm)
        num = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Ffs_ecdf(xfsm))
        denum = interp1(Fhs_ecdf,xhs_ecdf,tmp);  %Calculation of Fhs_ecdf^1(Ffs_ecdf(xfsm))
        xfa(i) = xfs(i,end) * num/denum; %Equidistant calculation
    end
    for i= 1:size(xhs,1)
        tmp = max(Fhs_ecdf(xhs_ecdf <= xhs(i,end))); %Calculation of Fhs_ecdf(xhsm)
        xha(i) = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Fhs_ecdf(xhsm))
    end
    

end

