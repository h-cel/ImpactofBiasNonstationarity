%% Configuration file to launch bias adjustment
% Configuration script for the calculations done in Van de Velde et al. (in
% pub.)

% This files launches the bias adjustment of the GCM/RCM input.
% This was written by M.T. Pham (2016) and adapted by K. De Roos (2018) and
% J. Van de Velde (2018 and onwards)

% Last update by J. Van de Velde on 09/12/'19

close all; clear all;
clc;

%% 1. Add bias correction code and its subfolders to Matlab's path

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\Onderzoek'), genpath('E:\Users\jpvdveld\Onderzoek\Data')); %Both Code and Data paths need to be added with their subfolders.
save_loc = 'E:\Users\jpvdveld\Onderzoek\Data\0_original\';

%% 2. Specify datasets, datatype, control and future period

% Format of dataset:
% year month day var1 (...) varN
% Remark: Fill datagaps with NaN. If data is loaded with the script
% a_loadClimateData, this should already have been done.

% ho  Original observations for the historical period
% hs  Original GCM output (simulations) for the historical period
% fs  Original GCM output (simulations) for the future period

tmp = matload('ETP_117y.mat'); %Historical observations: full file.
ho = tmp;
hs = matload('MPI_hist_or.mat'); %Depending on the file name implemented in a_loadClimateData
fs = matload('MPI_rcp45_or.mat'); %Depending on the file name implemented in a_loadClimateData

% Datatype specifies the input variable:
% "P" for Precipitation
% "E" for Evapotranspiration
% "T" for Temperature

datatype = {'E', 'T', 'P'};

% Format of both control and future period:
% [year month day]

% The datasets have the following timeframes:
% hs: 1/1/1956 - 31/12/2005
% fs: 1/1/2006 - 31/12/2100
% ho:  2/1/1901 - 31/12/2017
% However, a limited timeframe needs to be chosen, to be sure that the future period really is the 
% future and the control period really represents the history.
tControl =  [1970 1 1; 1989 12 31]; %Changed to 1970-1989/1998-2017 for the MPI-RCA dataset
tFuture = [1970 1 1; 1989 12 31]; %Comparison calibration period
%tFuture = [1998 1 1; 2017 12 31]; %Standard future period: [2071 1 1; 2100 12 31]; BC Comparison period: [1998 1 1; 2017 12 31]
%tControl =  [1971 1 1; 2000 12 31]; %Control period: [1971-1-1; 2000 12 31]
%tFuture = [2071 1 1; 2100 12 31]; %Standard future period: [2071 1 1; 2100 12 31]
%% 3. Prepare bias Data
% Data preparation

[xho, xhs, xfs] = prepareBiasdata(ho, hs, fs, tControl, tFuture); %Truncation of the datasets to the given timeframes
%save(strcat(save_loc, 'MPI-rcp45_calSim.mat'), 'xfs');
%save(strcat(save_loc, 'MPI-rcp45_xhs.mat'), 'xhs');
%save(strcat(save_loc, 'Uccle_xho.mat') , 'xho');
% TransformToNetCDF(xho, 'E:\Users\jpvdveld\Onderzoek\Data\0_original\', 'xho');
% TransformToNetCDF(xhs, 'E:\Users\jpvdveld\Onderzoek\Data\0_original\', 'xhs');
% TransformToNetCDF(xfs, 'E:\Users\jpvdveld\Onderzoek\Data\0_original\', 'xfs');

%% 4. Specify bias adjustment method

% A combination of different bias correction methods can be given. This is
% done via two Boolean vectors, in which each element represents a method

occ_methods = [0 0 0 1]; %Occurrence
int_methods = [0 0 0 0 1]; %Intensity

names_occ = {'none', 'ssr', 'tda', 'threshold'}; 
names_int = {'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc'};

%% 5. Specify ensemble members
% This is implemented as a correction for the enlargement of uncertainty by
% the stochastic steps. It is used as a base repetitions number, i.e. in
% each stochastic step n repetitions are used.

n = 10; %Standard: 20

%% 6. Start the procedure
% Make sure the selection is set as wanted!

for i=1:length(occ_methods)
    if occ_methods(i) == 1 %Checks which occurrence methods are being calculated
        name_occ = names_occ{i};
    for j = 1:length(int_methods)
        if int_methods(j) == 1 % Checks which intensity methods are being calculated
            name_int = names_int{j};
            tic
            BiasAdjustment(xho, xhs, xfs, datatype, name_occ , name_int, n, 'MPI-rcp45calibration')
            toc
        end
    end
    end
end

