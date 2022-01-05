function [xho, xhs, xfs] = prepareBiasdata(ho, hs, fs, tControl, tFuture)
% prepareBiasdata truncates the given datasets to the specified periods.
% This function is launched in the b_configurationBiasAdjustment config
% file and is used in the calculations done in Van de Velde et al. (in
% progress)
%
% Input   ho            vector that contains the original observations for the current/control period
%         hs            vector that contains the original GCM output for the current/control period 
%         hf            vector that contains the original GCM output for the future period
%         tControl      3D-row with year, month, day of the current/control period
%         tFuture       3D-row with year, month, day of the current/control period
% Output  xho           vector that contains the truncated original observations for the current/control period
%         xhs           vector that contains the truncated original GCM output for the current/control period 
%         xfs           vector that contains the truncated original GCM output for the future period
%
% Last update: 25/11/'19

%% Truncation
% In case the extended Uccle time series is used, pre, or post-2005 doens't
% matter for historical observations

tStart = find(ho(:,1) == tControl(1,1) & ho(:,2) == tControl(1,2) & ho(:,3) == tControl(1,3)); %Index of the start of the control period: year, month and day must be the same for both tControl and ho
tEnd = find(ho(:,1) == tControl(2,1) & ho(:,2) == tControl(2,2) & ho(:,3) == tControl(2,3)); %Index of the end of the control period: year, month and day must be the same for both tControl and ho
xho = ho(tStart:tEnd,:);

tStart = find(hs(:,1) == tControl(1,1) & hs(:,2) == tControl(1,2) & hs(:,3) == tControl(1,3)); %Same procedure, but for hs and tControl
tEnd = find(hs(:,1) == tControl(2,1) & hs(:,2) == tControl(2,2) & hs(:,3) == tControl(2,3)); %Same procedure, but for hs and tControl
xhs = hs(tStart:tEnd,:);
    
% For simulated series, it is possible that 'future' data that belongs to the 'historical' series has to be used.

if tFuture(1,1) > 2006
    
    %Truncating all data from the future dataset
    tStart = find(fs(:,1) == tFuture(1,1) & fs(:,2) == tFuture(1,2) & fs(:,3) == tFuture(1,3)); %Same procedure, but for fs and tFuture
    tEnd = find(fs(:,1) == tFuture(2,1) & fs(:,2) == tFuture(2,2) & fs(:,3) == tFuture(2,3)); %Same procedure, but for fs and tFuture
    xfs = fs(tStart:tEnd,:);

elseif tFuture(2,1) < 2006
    
    %Truncating data using only historical data
    tStart = find(hs(:,1) == tFuture(1,1) & hs(:,2) == tFuture(1,2) & hs(:,3) == tFuture(1,3)); %Same procedure, but for hs and tFuture
    tEnd = find(hs(:,1) == tFuture(2,1) & hs(:,2) == tFuture(2,2) & hs(:,3) == tFuture(2,3)); %Same procedure, but for hs and tFuture
    xfs = hs(tStart:tEnd,:);
    
else
    
    %Truncating data combined from historical and future datasets
    tStart = find(hs(:,1) == tFuture(1,1) & hs(:,2) == tFuture(1,2) & hs(:,3) == tFuture(1,3)); %Same procedure, but for hs and tFuture
    tEnd = find(fs(:,1) == tFuture(2,1) & fs(:,2) == tFuture(2,2) & fs(:,3) == tFuture(2,3)); %Same procedure, but for fs and tFuture
    xfs = [hs(tStart:end, :); fs(1:tEnd,:)];

    
end

end

