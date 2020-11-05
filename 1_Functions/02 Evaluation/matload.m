function X=matload(matfile)
% MATLOAD Loads .mat-file

%% Implementation

A=load(matfile);
vars=whos('-file',matfile);
X=A.(vars.name);

