%% Load climate data
% Configuration script for the calculations done in Van de Velde et al. (in
% pub.)
%
% This script loads the following climate datasets: observation data from
% the Uccle dataset, and climate modeled for both historic and
% future time spans. For this setup, RCA4 (Strandberg et al., 2015) forced with MPI-ESM-LR (Popke et al., 2013) was used. 
%
% The Climate model data is ordered by the RCP's: the Representative
% Concentration Pathways. These resemble the possible increases for net
% radiation in comparison to pre-industrial values according to different scenarios (e.g. an increase in
% 2100 of 4.5 W/m^2 is used in RCP 4.5). (IPCC, 2013)
%
% The climate data is used in the netcdf-format (.nc). The functions ncdisp
% and nc-read are used to work with this format.
% The ncdisp-function gives information about the used format, dimensions 
% and attributes of the specified variable e.g. the timeframe of the 
% .nc-dataset. The ncread-function extracts the specified variable of 
% .nc-dataset in Matlab format. 
% For more info about the netcdf-format, see e.g. https://en.wikipedia.org/wiki/NetCDF
%
% Code originally written by M.T. Pham (2016), adapted by J. Van de Velde
% (Last update: 25/11/'19)

clear all; 
close all; 
clc;

%% Set paths
%Extra paths

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

%Specific paths to save to

save_or = 'E:\Users\jpvdveld\Onderzoek\Data\0_original\'; %Original or expanded original data path

%% Historical data

% MPI-RCA4
ncdisp('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'time');
tmp = ncread('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'time');
test = datevec(datenum([1949 12 1]) + tmp); %Returns vector with Y M D H M S
time = test(:,1:3); %Only the first three columns are relevant, H M S aren't used

ncdisp('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'evspsblpot');
vardata = ncread('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'evspsblpot');
E = reshape(vardata(2,2,:),1,13149)*24*60*60; %Transformation to amount of evaporation per day

ncdisp('pr__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'pr');
vardata = ncread('pr__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'pr');
P = reshape(vardata(2,2,:),1,13149)*24*60*60;

ncdisp('tas__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'tas');
vardata = ncread('tas__EUR-11_MPI-M-MPI-ESM-LR_historical_r1i1p1_SMHI-RCA4_v1a_day_19700101-20051231_npix.nc', 'tas');
T = reshape(vardata(2,2,:),1,13149)-273.15; %Transformation from Kelvin to °C

MPI_hist_or = [time E' T' P']; %Placing all data together in one dataset
save(strcat(save_or, 'MPI_hist_or.mat'), 'MPI_hist_or'); %Saving in the right folder!

%% RCP 4.5

% MPI-RCA4
%For an explanation of the different steps, see MPI-RCA4, historical data
ncdisp('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','time')
tmp = ncread('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','time');
test = datevec(datenum([1949 12 1]) + tmp);
time = test(:,1:3);

ncdisp('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','evspsblpot');
vardata = ncread('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','evspsblpot');
E = reshape(vardata(2,2,:),1,34698)*24*60*60;

ncdisp('pr__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','pr');
vardata = ncread('pr__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','pr');
P = reshape(vardata(2,2,:),1,34698)*24*60*60;

ncdisp('tas__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','tas');
vardata = ncread('tas__EUR-11_MPI-M-MPI-ESM-LR_rcp45_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','tas');
T = reshape(vardata(2,2,:),1,34698)-273.15;

MPI_rcp45 = [time E' T' P']; 
save(strcat(save_or, 'MPI_rcp45_or.mat'), 'MPI_rcp45');

%% RCP 8.5

% MPI
%For an explanation of the different steps, see MPI-RCA4, historical data
ncdisp('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','time')
tmp = ncread('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','time');
test = datevec(datenum([1949 12 1]) + tmp);
time = test(:,1:3);

ncdisp('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','evspsblpot');
vardata = ncread('evspsblpot__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','evspsblpot');
E = reshape(vardata(2,2,:),1,34698)*24*60*60;

ncdisp('pr__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','pr');
vardata = ncread('pr__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','pr');
P = reshape(vardata(2,2,:),1,34698)*24*60*60;

ncdisp('tas__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','tas');
vardata = ncread('tas__EUR-11_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_SMHI-RCA4_v1a_day_20060101-21001231_npix.nc','tas');
T = reshape(vardata(2,2,:),1,34698)-273.15;

MPI_rcp85 = [time E' T' P'];
save(strcat(save_or, 'MPI_rcp85_or.mat'), 'MPI_rcp85');

%%

clear tmp time E P T test vardata i days dyVec mnts mntVec time2 yrs yrVec