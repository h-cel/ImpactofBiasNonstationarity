function [outfs, ouths] = occAdj_Threshold(xho, xhs, xfs)
%   occAdj_THRESHOLD This function implements a standard precipitation occurrence bias
%   adjustment method using a threshold
%
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress) 
%
%   This function adjust the precipitation occurrence bias in a given time
%   series.
%   Thresholding is based on two assumptions 
%   1)the bias is always wet
%   2)the wet bias is time-invariant: the same bias can be used for control
%   and future time series
%
%   Inputs: 
%       xho: historical observation, a n x 4-matrix with the respective
%       variable in the last column and Y M D in the first 3 columns
%       xhs: historical simulations, idem
%       xfs: future simulations, idem
%   Outputs: 
%       outfs: corrected future time series
%       ouths: corrected historical time series
%
%   Last update by J. Van de Velde on 09/12/'19

%% Initialization
ouths = nan(size(xhs,1),1); %out is initialized by making it the same size as the last column of xhs
outfs = nan(size(xfs,1),1); %out is initialized by making it the same size as the last column of xfs


%% Monthly loop

for m = 1:12
    idmhs = find(xhs(:, 2) == m);
    idmfs = find(xfs(:, 2) == m);
    
    %Selection of the rows belonging to the month m
    xhsm = xhs(xhs(:,2) == m, :);
    xfsm = xfs(xfs(:,2) == m, :);
    xhom = xho(xho(:, 2) == m, :);
    
    %Calculation of dry day frequencies
    fdryho = sum(xhom(:, end) == 0)/size(xhom,1);
    fdryhs = sum(xhsm(:,end) < 0.1)/size(xhsm,1);
    deltaN = length(xhsm)*(fdryho-fdryhs);
    
    %CONTROL
    
    %Adjusting control time series 
    xhswet = xhsm((xhsm(:,end)>=0.1),end);
    xhswet_corr = xhswet;
    xhsdry = xhsm((xhsm(:,end)<0.1),end);
    [xhswet_sort, I] = sort(xhswet);
    for i = 1: deltaN
        xhswet_sort(i) = 0;
    end
    
    for i = 1: length(xhswet)
        xhswet_corr(I(i)) = xhswet_sort(i);    
    end
    
    % Rearranging the control time series
    cnt1 = 0;
    cnt2 = 0;
    xhsdw = xhsm;
    for i = 1: size(xhsm, 1)
        if (xhsm(i,end)>0.1) == true
            cnt1= cnt1 + 1;
            xhsdw(i,end) = xhswet_corr(cnt1);
        else
            cnt2 = cnt2+1;
            xhsdw(i,end) = xhsdry(cnt2);
        end
    end
    
    %FUTURE
    
    %Adjusting future time series 
    xfswet = xfsm((xfsm(:,end)>=0.1),end);
    xfsdry = xfsm((xfsm(:,end)<0.1),end);
    [xfswet_sort, I] = sort(xfswet);
    for i = 1: deltaN
        xfswet_sort(i) = 0;
    end
    
    for i = 1: length(xfswet)
        xfswet(I(i)) = xfswet_sort(i);    
    end
    
    % Rearranging the future time series
    cnt1 = 0;
    cnt2 = 0;
    xfsdw = xfsm; 
    for i = 1: size(xfsm, 1)
        if (xfsm(i,end)>0.1) == true
            cnt1= cnt1 + 1;
            xfsdw(i,end) = xfswet(cnt1);
        else
            cnt2 = cnt2+1;
            xfsdw(i,end) = xfsdry(cnt2);
        end
    end
    
    %FINISHING
    
    % Arrange in original order in total 
    for i = 1:length(idmhs)
        ouths(idmhs(i)) = xhsdw(i,end);  
    end
    
    for i = 1:length(idmfs)
        outfs(idmfs(i)) = xfsdw(i,end);
    end 
    
end

end

