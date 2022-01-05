%% Case 0
%   This is a script that runs the PDM model using measurements
%
%   This script is launched in the Discharge.m function
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   Last update by J. Van de Velde on 27/11/'19

%% Simulations of Q using original P and original E 
%Initialization
ndays = size(time,1);
Qsim = nan(ndays,4);
Qsim(:,1:3) = time;
%Input preparation
inputs.P = reshape(kron(P/24, ones(1,24))', ndays*24,1); %Hourly precipitation input data
inputs.A = 385;
inputs.E = reshape(kron(E/24, ones(1,24))', ndays*24,1); %Hourly evaporation input data
%Simulation
Qsim(:,4) = PDMPieter(inputs, PDM_par);
disp('Q case 0 simulated.')

%% Save output

save(['E:\Users\jpvdveld\Onderzoek\Data\3_cases\' names{1} '_' names{5} '_c0.mat'], 'Qsim');
