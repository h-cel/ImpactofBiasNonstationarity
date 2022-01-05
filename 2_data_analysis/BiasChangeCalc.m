function [RIndex] = BiasChangeCalc(index, loc1, months, quantnum, loc2, lag)
%BIASCHANGE This function calculates the R Index for a specific result index
%   Calculation of the R Index, proposed by Maurer et al. (2013)
%   
%   Inputs:
%       Index: index for which the R index has to be calculated
%       loc1: location of the variable for which the R index has to be
%       calculated
%       quantnum: in case a quantile index has to be calculated,
%       specification of the quantile
%       loc2: in case a correlation has to be calculated, the location of
%       the second variable
%       lag: in case the crosscorrelation has to be calculated, the lag
%   Output:
%       RIndex: the R Index value
%
%   Last update by J. Van de Velde on 15/02/'21

%% Paths

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\Onderzoek'), genpath('E:\Users\jpvdveld\Onderzoek\Data\0_original')); %Both Code and Data paths need to be added with their subfolders.

%% Load data

xfs = matload('MPI-rcp45_xfs.mat');
xhs = matload('MPI-rcp45_xhs.mat');
xho = matload('Uccle_xho.mat');
xobs = matload('Uccle_xobs.mat');

%% Timeframe selection

xobsm = xobs(ismember(xobs(:,2),months), :);
xhom = xho(ismember(xho(:,2),months), :);
xfsm = xfs(ismember(xfs(:,2),months), :);
xhsm = xhs(ismember(xhs(:,2),months), :);

%% Index selection

switch index
    case 'average'
        ho = mean(xhom(:,loc1), 'omitnan');
        hs = mean(xhsm(:,loc1), 'omitnan'); 
        fs = mean(xfsm(:,loc1), 'omitnan');
        fo = mean(xobsm(:,loc1), 'omitnan');            
    case 'standarddev'
        ho = std(xhom(:,loc1), 'omitnan');
        hs = std(xhsm(:,loc1), 'omitnan');
        fs = std(xfsm(:,loc1), 'omitnan');
        fo = std(xobsm(:,loc1), 'omitnan');
    case 'quant'
        ho = quantile(xhom(:,loc1), quantnum);
        hs = quantile(xhsm(:,loc1), quantnum);
        fo = quantile(xobsm(:,loc1), quantnum);
        fs = quantile(xfsm(:,loc1), quantnum);
    case 'correl'
        ho = corr(xhom(:, loc1),xhom(:, loc2), 'type', 'Spearman', 'rows', 'complete');
        hs = corr(xhsm(:, loc1),xhsm(:, loc2), 'type', 'Spearman', 'rows', 'complete');
        fo = corr(xobsm(:, loc1),xobsm(:, loc2), 'type', 'Spearman', 'rows', 'complete');
        fs = corr(xfsm(:, loc1),xfsm(:, loc2), 'type', 'Spearman', 'rows', 'complete');       
    case 'crosscor'
        ho_temp = xcorr(xhom(~isnan(xhom(:,loc1))&~isnan(xhom(:,loc2)), loc1),xhom(~isnan(xhom(:,loc1))&~isnan(xhom(:,loc2)),loc2), 1, 'coeff');
        hs_temp = xcorr(xhsm(~isnan(xhsm(:,loc1))&~isnan(xhsm(:,loc2)), loc1),xhsm(~isnan(xhsm(:,loc1))&~isnan(xhsm(:,loc2)),loc2), 1, 'coeff');
        fo_temp = xcorr(xobsm(~isnan(xobsm(:,loc1))&~isnan(xobsm(:,loc2)), loc1),xobsm(~isnan(xobsm(:,loc1))&~isnan(xobsm(:,loc2)),loc2), 1, 'coeff');
        fs_temp = xcorr(xfsm(~isnan(xfsm(:,loc1))&~isnan(xfsm(:,loc2)), loc1),xfsm(~isnan(xfsm(:,loc1))&~isnan(xfsm(:,loc2)),loc2), 1, 'coeff');
        
        ho = ho_temp(lag+2);
        hs = hs_temp(lag+2);
        fo = fo_temp(lag+2);
        fs = fs_temp(lag+2);
    case 'autocorr'
        ho_temp = xcorr(xhom(:,loc1), 1, 'coeff');
        hs_temp = xcorr(xhsm(:,loc1), 1, 'coeff');
        fo_temp = xcorr(xobsm(:,loc1), 1, 'coeff');
        fs_temp = xcorr(xfsm(:,loc1), 1, 'coeff');
        ho = ho_temp(3, 1);
        hs = hs_temp(3, 1);
        fo = fo_temp(3, 1);
        fs = fs_temp(3, 1);
    case 'trans00'
        transho = TransProb(xhom(:, loc1), 0.1);
        transhs = TransProb(xhsm(:, loc1), 0.1);
        transfo = TransProb(xobsm(:, loc1), 0.1);
        transfs = TransProb(xfsm(:, loc1), 0.1);
        ho = transho(1);
        hs = transhs(1);
        fo = transfo(1);
        fs = transfs(1);
    case 'trans10'
        transho = TransProb(xhom(:, loc1), 0.1);
        transhs = TransProb(xhsm(:, loc1), 0.1);
        transfo = TransProb(xobsm(:, loc1), 0.1);
        transfs = TransProb(xfsm(:, loc1), 0.1);
        ho = transho(3);
        hs = transhs(3);
        fo = transfo(3);
        fs = transfs(3);
    case 'ndry'
        ho =sum(xhom(:, loc1)<0.1);
        hs =sum(xhsm(:, loc1)<0.1);
        fo =sum(xobsm(:, loc1)<0.1);
        fs =sum(xfsm(:, loc1)<0.1);
end

bias_h = hs-ho;
bias_f = fs-fo;
RIndex = abs(bias_f-bias_h)/((abs(bias_f)+abs(bias_h))/2);
end

