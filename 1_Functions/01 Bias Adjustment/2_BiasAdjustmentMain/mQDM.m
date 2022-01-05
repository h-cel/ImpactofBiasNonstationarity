function xfa = mQDM(xho,xhs,xfs,type)
%   mQDM This function implements the modified Quantile Delta Mapping bias adjustment method
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress). 
%
%   This function implements the following calculation, used to adjust the
%   Model output statistics of GCM/RCM-modelled data:
%
%   xfa = xho * Ffs^1(Fho(xho))/Fhs^-1(Fho(xho))
%   xfa = xho + Ffs^1(Fho(xho)) - Fhs^-1(Fho(xho))
%
%   With:
%       xfa: the rescaled future period model output
%       xho: the control period dataset
%       Fho: the cdf of the control dataset
%       Ffs: the cdf of the future simulations
%       Fhs: the cdf of the historic simulations
%
%   This calculation corrects the control period data according to the
%   climate signal, which is either the difference or the ratio between the future and control
%   period GCM/RCM simulations.
%
%   Inputs:
%       xho: the data of the control period, a n x 4-matrix with the respective
%       variable in the last column and Y M D in the first 3 columns
%       xhs: the modelled data of the control period, idem
%       xfs: the modelled data of the future period, idem
%       type: 1, 2. Indicates whether non-relative changes or changes for all rainy days (2) have to be
%       used.
%   Outputs:
%       xfa: adjusted future time series
%
%   Last update by J. Van de Velde on 09/12/'19

%% mQDM

if type == 1 %Relative changes. As relative changes are never frequency-corrected, there is no need for wet-day-selection
    
    %Initialization
    xfa = nan(size(xho(:,end)));

    %Looping over the days
    nrowsh = size(xho(:,end), 1);
    nrowsf = size(xfs(:,end),1);
    for j = 1:nrowsh
        %Selection of data for a moving window of 91 days, out of each year
        %-> total number of days used = amount of years x 91
        rangewindowh = max(j-45, 1):min(j+45,nrowsh);
        rangewindowf = max(j-45, 1):min(j+45,nrowsf);%Correction for the days at the beginning of the time series
        windstartf = [xfs(rangewindowf(1),2), xfs(rangewindowf(1),3)];
        windstarth = [xhs(rangewindowh(1),2), xhs(rangewindowh(1),3)];
        
        xhom = [];
        xhsm = [];
        xfsm = [];
        
        % Collection of data
        for i = 1:nrowsf
            if xfs(i,2) == windstartf(1) && xfs(i,3) == windstartf(2)
                xfsm = [xfsm; xfs(i:min(i+length(rangewindowf), nrowsf),end)];
            end
        end
        for i=1:nrowsh
            if xhs(i,2) == windstarth(1) && xhs(i,3) == windstarth(2)
                xhsm = [xhsm; xhs(i:min(i+length(rangewindowh), nrowsh),end)];
                xhom = [xhom; xho(i:min(i+length(rangewindowh), nrowsh),end)];
            end
        end
        %Making ecdf's
        [Fho_ecdf, xho_ecdf] = ecdf(xhom(:,end));
        [Fhs_ecdf, xhs_ecdf] = ecdf(xhsm(:,end));
        [Ffs_ecdf, xfs_ecdf] = ecdf(xfsm(:,end));
        tmp = max(Fho_ecdf(xho_ecdf <= xho(j,end)));  %Calculation of Fho_ecdf(xhom)
        num = interp1(Ffs_ecdf,xfs_ecdf,tmp); %Calculation of Ffs_ecdf^1(Fho_ecdf(xhom))
        denum = interp1(Fhs_ecdf,xhs_ecdf,tmp); %Calculation of Fhs_ecdf^1(Fho_ecdf(xhom))
        xfa(j) = xho(j,end) + num - denum; %Equidistant calculation
    end

elseif type == 2
    xfa = xho(:,end); %This ensures that the data of dry days is not disturbed
    
    nrowsh = size(xho(:,end), 1);
    nrowsf = size(xfs(:,end),1);
    for j = 1:nrowsh
        if xho(j, end) > 0.1
        %Selection of data for a moving window of 91 days, out of each year
        %-> total number of days used = amount of years x 91
        rangewindowh = max(j-45, 1):min(j+45,nrowsh);
        rangewindowf = max(j-45, 1):min(j+45,nrowsf);%Correction for the days at the beginning of the time series
        windstartf = [xfs(rangewindowf(1),2), xfs(rangewindowf(1),3)];
        windstarth = [xhs(rangewindowh(1),2), xhs(rangewindowh(1),3)];
        
        xhom = [];
        xhsm = [];
        xfsm = [];
        
        % Collection of data
        for i = 1:nrowsf
            if xfs(i,2) == windstartf(1) && xfs(i,3) == windstartf(2)
                xfsm = [xfsm; xfs(i:min(i+length(rangewindowf), nrowsf),end)];
            end
        end
        for i=1:nrowsh
            if xhs(i,2) == windstarth(1) && xhs(i,3) == windstarth(2)
                xhsm = [xhsm; xhs(i:min(i+length(rangewindowh), nrowsh),end)];
                xhom = [xhom; xho(i:min(i+length(rangewindowh), nrowsh),end)];
            end
        end

        %Making ecdf's
        [Fho_ecdf, xho_ecdf] = ecdf(xhom(:,end));
        [Fhs_ecdf, xhs_ecdf] = ecdf(xhsm(:,end));
        [Ffs_ecdf, xfs_ecdf] = ecdf(xfsm(:,end));

        %Calculation of correction for the future time series
        tmp = max(Fho_ecdf(xho_ecdf <= xho(j,end))); %Calculation of Fho_ecdf(xhom)
        num = interp1(Ffs_ecdf,xfs_ecdf,tmp); %Calculation of Ffs_ecdf^1(Fho_ecdf(xhom))
        denum = interp1(Fhs_ecdf,xhs_ecdf,tmp); %Calculation of Fhs_ecdf^1(Fho_ecdf(xhom))
        xfa(j) = xho(j,end) * num/denum; %Equiratio calculation, but only the wet days are changed.
        end
    end

end

end

