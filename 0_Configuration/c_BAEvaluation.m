%% Configuration file to evaluate bias adjustment methods
%   Configuration script for the evaluation of the calculations done in Van de Velde et al. (in
%   pub.)
%
%   Last update by J. Van de Velde on 09/12/'19

%% Clearing files

clc
clear all

%% Path selection

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\Onderzoek'), genpath('E:\Users\jpvdveld\Onderzoek\Data')); %Both Code and Data paths need to be added with their subfolders.

%% Set-up

% Window
% same as the tFuture used to run the bias adjustment to be evaluated

tComp = [1998 1 1; 2017 12 31];

% Methods
% Selection the methods that have to be included in the current evaluation
% calculations

occ_methods = [0 0 0 1]; %Occurrence 
int_methods = [1 1 1 1 1 0]; %Intensity

% Model name (can be looped to evaluate different dataset!)
% Same name as used to save the bias adjustment methods

name = 'MPI-rcp45';

%% Time series

% Observations

tmp = matload('ETP_117y.mat'); %Historical observations: full file
ho = tmp;
xobs = TruncateObs(ho, tComp);
save('E:\Users\jpvdveld\Onderzoek\Data\0_original\Uccle_xobs.mat', 'xobs');

% Original timeseries

data= load(strcat(name, '_xfs'));
xfs = data.xfs;

%% Function

[MPI_rcp45_indices, MPI_rcp45_bias, MPI_rcp45_RB_O, MPI_rcp45_RB_MB] = BA_Evaluation(xfs, xobs, name, '_threshOnly', occ_methods, int_methods, 0);