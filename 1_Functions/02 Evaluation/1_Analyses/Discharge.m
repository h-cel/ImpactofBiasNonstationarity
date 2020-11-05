function [Qsim] = Discharge(data, name)
%   DISCHARGE Runs PDM for a certain dataset
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function runs the Probability-Distributed Model for runoff
%   generation.
%
%   Inputs:
%       data: the dataset for which the PDM will run
%       name: the name of the dataset, to be able to save this clearly
%   Output:
%       Qsim: simulated discharge
%
%   Last update by J. Van de Velde on 09/12/'19

%% Inputs

time = data(:,1:3); %Used in the case
E = data(:, 4);
E(isnan(E)) = 0;
P = data(:, 6);
P(isnan(P)) = 0;
PDM_par = matload('paramPDM.mat');
names = {name, '', '', '', 'Qsim', ''};

run('case0.m')


end

