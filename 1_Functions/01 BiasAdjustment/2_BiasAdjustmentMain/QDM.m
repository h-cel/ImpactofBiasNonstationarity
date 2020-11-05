function [xha, xfa] = QDM(xho,xhs,xfs,type)
%   QDM This function implements the Quantile Delta Mapping bias adjustment method
%
%   This function is launched in the Biasadjustment.m function
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
%       xho: the data of the control period, a n x 4-matrix with the respective
%       variable in the last column and Y M D in the first 3 columns
%       xhs: the modelled data of the control period, idem
%       xfs: the modelled data of the future period, idem
%       type: 1, 2. Indicates whether non-relative changes or changes for all rainy days (2) have to be
%       used.
%   Outputs:
%       xha: adjusted historical time series
%       xfa: adjusted future time series
%
%   Last update by J. Van de Velde on 09/12/'19

%% QDM

if type == 1 %Non-Relative changes. 
    
    %Initialization
    xfa = nan(size(xfs(:,end)));
    xha = nan(size(xhs(:,end)));

    %Looping over the days
    nrows = size(xfs(:,end), 1);
    for j = 1:nrows
        %Selection of data for a moving window of 91 days, out of each year
        %-> total number of days used = amount of years x 91
        rangewindow = max(j-45, 1):min(j+45,nrows); %Correction for the days at the beginning of the time series
        windstartf = [xfs(rangewindow(1),2), xfs(rangewindow(1),3)];
        windstarth = [xhs(rangewindow(1),2), xhs(rangewindow(1),3)];
        
        xhom = [];
        xhsm = [];
        xfsm = [];
        
        % Collection of data
        for i = 1:nrows
            if xfs(i,2) == windstartf(1) && xfs(i,3) == windstartf(2)
                xfsm = [xfsm; xfs(i:min(i+length(rangewindow), nrows),end)];
            end
            if xhs(i,2) == windstarth(1) && xhs(i,3) == windstarth(2)
                xhsm = [xhsm; xhs(i:min(i+length(rangewindow), nrows),end)];
                xhom = [xhom; xho(i:min(i+length(rangewindow), nrows),end)];
            end
        end
        %Making ecdf's
        [Fho_ecdf, xho_ecdf] = ecdf(xhom(:,end));
        [Fhs_ecdf, xhs_ecdf] = ecdf(xhsm(:,end));
        [Ffs_ecdf, xfs_ecdf] = ecdf(xfsm(:,end));
        
        %Calculation of correction for the future time series
        tmp = max(Ffs_ecdf(xfs_ecdf <= xfs(j, end)));  %Calculation of Ffs_ecdf(xfsm)
        num = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Ffs_ecdf(xfsm))
        denum = interp1(Fhs_ecdf,xhs_ecdf,tmp);  %Calculation of Fhs_ecdf^1(Ffs_ecdf(xfsm))
        xfa(j) = xfs(j,end) + num - denum; %Equiratio calculation
        
        %Calculation of the correction for the historic time series
        tmp = max(Fhs_ecdf(xhs_ecdf <= xhs(j, end))); %Calculation of Fhs_ecdf(xhsm)
        xha(j) = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Fhs_ecdf(xhsm))
    end
    
elseif type == 2
    xfa = xfs(:,end); %This ensures that the data of dry days is not disturbed
    xha = xhs(:,end); %Idem
    
    nrows = size(xfs(:,end), 1);
    for j = 1:nrows
        %Selection of data for a moving window of 91 days, out of each year
        %-> total number of days used = amount of years x 91
        rangewindow = max(j-45, 1):min(j+45,nrows); %Correction for the days at the beginning of the time series
        windstartf = [xfs(rangewindow(1),2), xfs(rangewindow(1),3)];
        windstarth = [xhs(rangewindow(1),2), xhs(rangewindow(1),3)];
        
        xhom = [];
        xhsm = [];
        xfsm = [];
        
        % Collection of data
        for i = 1:nrows
            if xfs(i,2) == windstartf(1) && xfs(i,3) == windstartf(2)
                xfsm = [xfsm; xfs(i:min(i+length(rangewindow), nrows),end)];
            end
            if xhs(i,2) == windstarth(1) && xhs(i,3) == windstarth(2)
                xhsm = [xhsm; xhs(i:min(i+length(rangewindow), nrows),end)];
                xhom = [xhom; xho(i:min(i+length(rangewindow), nrows),end)];
            end
        end

        %Making ecdf's
        [Fho_ecdf, xho_ecdf] = ecdf(xhom(:,end));
        [Fhs_ecdf, xhs_ecdf] = ecdf(xhsm(:,end));
        [Ffs_ecdf, xfs_ecdf] = ecdf(xfsm(:,end));

        if xfs(j, end) > 0.1 
            %Calculation of correction for the future time series
            tmp = max(Ffs_ecdf(xfs_ecdf <= xfs(j, end))); %Calculation of Ffo_ecdf(xfsm)
            num = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Ffo_ecdf(xfsm))
            denum = interp1(Fhs_ecdf,xhs_ecdf,tmp); %Calculation of Fhs_ecdf^1(Ffo_ecdf(xfsm))
            xfa(j) = xfs(j,end) * num/denum; %Equiratio calculation, but only the wet days are changed
        end

        if xhs(j, end) > 0.1
            %Calculation of the correction for the historic time series
            tmp = max(Fhs_ecdf(xhs_ecdf <= xhs(j, end))); %Calculation of Fhs_ecdf(xhsm)
            xha(j) = interp1(Fho_ecdf,xho_ecdf,tmp); %Calculation of Fho_ecdf^1(Fhs_ecdf(xhsm))
        end
    end

end
    

end

