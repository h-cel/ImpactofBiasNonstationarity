function [RIndex] = BiasChangeCalc(index, loc1, quantnum, loc2, lag)
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
%   Last update by J. Van de Velde on 26/03/'20

%% Paths

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\Onderzoek'), genpath('E:\Users\jpvdveld\Onderzoek\Data\0_original')); %Both Code and Data paths need to be added with their subfolders.

%% Load data

xfs = matload('MPI-rcp45_xfs.mat');
xhs = matload('MPI-rcp45_xhs.mat');
xho = matload('Uccle_xho.mat');
xobs = matload('Uccle_xobs.mat');

%% Index selection

switch index
    case 'average'
        ho = mean(xho(:,loc1), 'omitnan');
        hs = mean(xhs(:,loc1), 'omitnan'); 
        fs = mean(xfs(:,loc1), 'omitnan');
        fo = mean(xobs(:,loc1), 'omitnan');            
    case 'standarddev'
        ho = std(xho(:,loc1), 'omitnan');
        hs = std(xhs(:,loc1), 'omitnan');
        fs = std(xfs(:,loc1), 'omitnan');
        fo = std(xobs(:,loc1), 'omitnan');
    case 'quant'
        ho = quantile(xho(:,loc1), quantnum);
        hs = quantile(xhs(:,loc1), quantnum);
        fo = quantile(xobs(:,loc1), quantnum);
        fs = quantile(xfs(:,loc1), quantnum);
    case 'correl'
        ho = corr(xho(:, loc1),xho(:, loc2), 'type', 'Spearman', 'rows', 'complete');
        hs = corr(xhs(:, loc1),xhs(:, loc2), 'type', 'Spearman', 'rows', 'complete');
        fo = corr(xobs(:, loc1),xobs(:, loc2), 'type', 'Spearman', 'rows', 'complete');
        fs = corr(xfs(:, loc1),xfs(:, loc2), 'type', 'Spearman', 'rows', 'complete');       
    case 'crosscor'
        ho_temp = xcorr(xho(~isnan(xho(:,loc1))&~isnan(xho(:,loc2)), loc1),xho(~isnan(xho(:,loc1))&~isnan(xho(:,loc2)),loc2), 1, 'coeff');
        hs_temp = xcorr(xhs(~isnan(xhs(:,loc1))&~isnan(xhs(:,loc2)), loc1),xhs(~isnan(xhs(:,loc1))&~isnan(xhs(:,loc2)),loc2), 1, 'coeff');
        fo_temp = xcorr(xobs(~isnan(xobs(:,loc1))&~isnan(xobs(:,loc2)), loc1),xobs(~isnan(xobs(:,loc1))&~isnan(xobs(:,loc2)),loc2), 1, 'coeff');
        fs_temp = xcorr(xfs(~isnan(xfs(:,loc1))&~isnan(xfs(:,loc2)), loc1),xfs(~isnan(xfs(:,loc1))&~isnan(xfs(:,loc2)),loc2), 1, 'coeff');
        
        ho = ho_temp(lag+2);
        hs = hs_temp(lag+2);
        fo = fo_temp(lag+2);
        fs = fs_temp(lag+2);
    case 'autocorr'
        ho_temp = xcorr(xho(:,loc1), 1, 'coeff');
        hs_temp = xcorr(xhs(:,loc1), 1, 'coeff');
        fo_temp = xcorr(xobs(:,loc1), 1, 'coeff');
        fs_temp = xcorr(xfs(:,loc1), 1, 'coeff');
        ho = ho_temp(3, 1);
        hs = hs_temp(3, 1);
        fo = fo_temp(3, 1);
        fs = fs_temp(3, 1);
    case 'trans00'
        transho = TransProb(xho(:, loc1), 0.1);
        transhs = TransProb(xhs(:, loc1), 0.1);
        transfo = TransProb(xobs(:, loc1), 0.1);
        transfs = TransProb(xfs(:, loc1), 0.1);
        ho = transho(1);
        hs = transhs(1);
        fo = transfo(1);
        fs = transfs(1);
    case 'trans10'
        transho = TransProb(xho(:, loc1), 0.1);
        transhs = TransProb(xhs(:, loc1), 0.1);
        transfo = TransProb(xobs(:, loc1), 0.1);
        transfs = TransProb(xfs(:, loc1), 0.1);
        ho = transho(3);
        hs = transhs(3);
        fo = transfo(3);
        fs = transfs(3);
    case 'ndry'
        ho =sum(xho(:, loc1)<0.1);
        hs =sum(xhs(:, loc1)<0.1);
        fo =sum(xobs(:, loc1)<0.1);
        fs =sum(xfs(:, loc1)<0.1);
end

bias_h = ho-hs;
bias_f = fo-fs;
RIndex = abs(bias_f-bias_h)/((abs(bias_f)+abs(bias_h))/2);
end

