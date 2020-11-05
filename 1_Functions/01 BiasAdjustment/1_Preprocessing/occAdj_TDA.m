function [outfs, ouths] = occAdj_TDA(xfs, xho, xhs)
%   occAdj_TDA This function implements the precipitation
%   occurrence bias adjustment as proposed by M.T.Pham ('Triangular distribution correction')
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress) 
%
%   This function adjust the precipitation occurrence bias in a given time
%   series.
%   This function is based on the method described in Pham (2016) and is
%   made by merging two functions, drywetCorrectionOC and
%   drywetCorrectionOGF, written by him.
%
%   This occurrence adjustment has several important steps. First, the
%   frequency of dry days in the future has to be calculated. This frequency
%   is used to calculate the number of days deltaN that has to be added or
%   substracted (if deltaN is either positive or negative). This addition and
%   substraction is done using a triangular distribution with an upper limit
%   of 0.9 and an extra check. This check is to make sure that no extremes
%   are altered.
%   
%   Inputs:
%       xfs: the future simulations, to be corrected, a n x 4-matrix with the respective
%       variable in the last column and Y M D in the first 3 columns
%       xho: the historical observations, idem
%       xhs: the historical simulations, idem
%   Outputs:
%       outfs: matrix of the corrected future data
%       ouths: matrix of the corrected historical data
%
%   Last update by J. Van de Velde on 09/12/'19

%% Initialization
outfs = nan(size(xfs,1),1); %outfs is initialized by making it the same size as the last column of xfs
ouths = nan(size(xhs,1),1); %idem

%% Monthly loop

for m = 1:12 %For each month
    idm = find(xfs(:,2) == m); %Gets all the indices of the month m
    
    %Selection of the rows belonging to the month m
    xfsm = xfs(xfs(:,2) == m ,:);
    xhsm = xhs(xhs(:,2) == m,:);
    xhom = xho(xho(:,2) == m,:);
    
    %Calculation of dry day frequencies
    %In the future and control period simulations, < 0.1 mm is counted as a dry day, this is the measurement error of the Uccle timeseries
    fhs = sum(xhsm(:,end) < 0.1)/size(xhsm,1); % Fraction in GCM/RCM simulations for control period 

    fho = sum(xhom(:,end) == 0)/size(xhom,1); %Frequency of dry days in control period
    ffs = sum(xfsm(:,end) < 0.1)/size(xfsm,1); % Frequency of dry days in GCM/RCM simulations for future period 
    fFDF = ffs * fho/fhs; %Calculation of the Future Dry day Frequency, which is a rescaling of ffs according to the equiratio procedure
    deltaN = length(xfsm)*(fFDF-ffs); %Difference in dry days, this needs to be added/substracted if deltaN is either positive or negative
    
    xfs_wet = xfsm((xfsm(:,end) >= 0.1), end); %Selection of wet days
    [Ffs_wet_ecdf, xfs_wet_ecdf] = ecdf(xfs_wet); %ecdf of wet days
    xfs_dry = xfsm((xfsm(:,end) < 0.1), end); %Selection of dry days
    
    %Dry-wet correction
    b = 0.9; %Parameter for the triangular distribution, the upper value
    %Addition of dry days
    if deltaN > 0
        i = 1;
        while i <= deltaN
            t = randi([1 length(xfs_wet)],1); % Random selection of index
            xt = xfs_wet(t); %Selection of a random value: value
            ksi = max(Ffs_wet_ecdf(xfs_wet_ecdf <= xt)); %CDF selection
            
            if T(ksi,b) < rand(1) %Check using the triangular distribution
                xfs_wet(t) = 0.1*rand(1); %Replacement of a wet day by a dry day (<0.1 mm)
                i = i + 1; %Counter for the number of days
            end
        end
    %Substraction of dry days    
    elseif deltaN < 0
        i = 1;
        while i <= abs(deltaN)
            t = randi([1 length(xfs_dry)],1); %Random selection of index
            ksi = b*(1-sqrt(1-rand(1)));
            xfs_dry(t) = max(xfs_wet_ecdf(Ffs_wet_ecdf <= ksi)); %Update of the selected day
            i = i + 1; %Counter for the number of days
        end
    end
      
    % Arrange in original order for month m
    cnt1 = 0;
    cnt2 = 0;
    xfsdw = xfsm;
    for i = 1:size(xfsm,1)
        if (xfsm(i,end) >= 0.1) == true %Wet days
            cnt1 = cnt1 + 1;
            xfsdw(i,end) = xfs_wet(cnt1);
        else %Dry days
            cnt2 = cnt2 + 1;
            xfsdw(i,end) = xfs_dry(cnt2);
        end
    end
    
    % Arrange in original order in total
    for i = 1:length(idm)
        outfs(idm(i)) = xfsdw(i,end);
    end
    
    %% Correction for historic time series
    
    deltaN = length(xhsm)*(fho-fhs); %Difference in dry days, this needs to be added/substracted if deltaN is either positive or negative
    
    xhs_wet = xhsm((xhsm(:,end) >= 0.1), end); %Selection of wet days
    [Fhs_wet_ecdf, xhs_wet_ecdf] = ecdf(xhs_wet); %ecdf of wet days
    xhs_dry = xhsm((xhsm(:,end) < 0.1), end); %Selection of dry days
    
    %Dry-wet correction
    b = 0.9; %Parameter for the triangular distribution, the upper value
    %Addition of dry days
    if deltaN > 0
        i = 1;
        while i <= deltaN
            t = randi([1 length(xhs_wet)],1); % Random selection of index
            xt = xhs_wet(t); %Selection of a random value: value
            ksi = max(Fhs_wet_ecdf(xhs_wet_ecdf <= xt)); %CDF selection
            
            if T(ksi,b) < rand(1) %Check using the triangular distribution
                xhs_wet(t) = 0.1*rand(1); %Replacement of a wet day by a dry day (<0.1 mm)
                i = i + 1; %Counter for the number of days
            end
        end
    %Substraction of dry days    
    elseif deltaN < 0
        i = 1;
        while i <= abs(deltaN)
            t = randi([1 length(xhs_dry)],1); %Random selection of index
            ksi = b*(1-sqrt(1-rand(1)));
            xhs_dry(t) = max(xhs_wet_ecdf(Fhs_wet_ecdf <= ksi)); %Update of the selected day
            i = i + 1; %Counter for the number of days
        end
    end
      
    % Arrange in original order for month m
    cnt1 = 0;
    cnt2 = 0;
    xhsdw = xhsm;
    for i = 1:size(xhsm,1)
        if (xhsm(i,end) >= 0.1) == true %Wet days
            cnt1 = cnt1 + 1;
            xhsdw(i,end) = xhs_wet(cnt1);
        else %Dry days
            cnt2 = cnt2 + 1;
            xhsdw(i,end) = xhs_dry(cnt2);
        end
    end
    
    % Arrange in original order in total
    for i = 1:length(idm)
        ouths(idm(i)) = xhsdw(i,end);
    end
    
end

end

