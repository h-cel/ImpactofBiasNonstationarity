function [indices, bias, RB_O, RB_MB] = BA_Evaluation(xfs, xobs, name, save_name, occ_methods, int_methods, wet)
%   BA_EVALUATION Evaluation of bias adjustment methods
%
%   This function is launched in the c_BAEvaluation.m script
%   file and is used in the evaluation done in Van de Velde et al. (in
%   progress) 
%
%   This function calculates indices of selected bias adjustment methods by comparing them
%   to the original simulations and to the observations. The methods compared
%   here are always a combination occurence and intensity correcting methods.
%
%   Inputs:
%       xfs: the RCM simulation data
%       xobs: observational data
%       name: the name of the RCM used
%       save_name: name that has to be included in the savefile as a marker
%       for what is to be evaluated.
%       occ_methods: Boolean vector [1x4] indicating which occurence methods
%       have been used. 1: 'none', 2: 'ssr', 3: 'tdc', 4: 'threshold'.
%       int_methods: Boolean vector [1x5] indicating which intensity methods
%       have been used. 1: 'qdm', 2: 'mqdm', 3: 'mbcn', 4: 'mrqnbc', 5:
%       'dotc', 6: 'other' -> 'other' is in case a special case has to be
%       evaluated
%       wet: Boolean indicating whether the P quantiles have to be
%       calculated on wet days only (1) or on all days (0)
%   Outputs:
%       indices: cell with indices calculated for observations, the original RCM simulations and every method
%       bias: cell with the biases (= difference) between 1) the observations and 2) the original RCM simulations and the methods used
%       RB_O: Residual Bias relative to the Observations
%       RB_MB: Residual Bias relative to the Model Bias
%
%   Last update by J. Van de Velde on 27/02/'20
 
%% Set-up

names_occ = {'none', 'ssr', 'tdc', 'threshold'};
names_int = {'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc'};
cnt_methods = sum(occ_methods)*sum(int_methods);
names_methods = cell(1, cnt_methods);
name_figs = cell(1, cnt_methods);
barQ = zeros(25,cnt_methods+1);

num_indices = 52;
indices = cell(num_indices+1, cnt_methods+3);
indices(2:end, 1) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
indices(1, 2:3) = {'Observations', 'Original RCM simulations'};
bias = cell(num_indices+1, cnt_methods+2);
bias(2:end,1) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
bias(1,2) = {'Original RCM simulations'};
RB_O = cell(num_indices+1, cnt_methods+1);
RB_O(2:end,1) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
RB_MB = cell(num_indices+1, cnt_methods+1);
RB_MB(2:end,1) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
save_loc = 'E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\';
save_loc_fig = 'D:\Users\jpvdveld\Documents\PhD\FigurenResultaten\Evaluation\';


%% OBSERVATIONS
%Calculation of the index values for the observations

% Discharge

    Qsim_obs = Discharge(xobs, 'Obs'); %Misschien meer mogelijkheden voor de naam doen als ik verschillende observatieperiodes wil gebruiken
    % I keep all the data for debugging/checking purposes)

    % Statistics (= general overview)
    
    qnt_standard = quantile(Qsim_obs(:, end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95 0.99, 0.995]);
    indices{2,2} = qnt_standard(1);
    indices{3,2} = qnt_standard(2);
    indices{4,2} = qnt_standard(3);
    indices{5,2} = qnt_standard(4);
    indices{6,2} = qnt_standard(5);
    indices{7,2} = qnt_standard(6);
    indices{8,2} = qnt_standard(7);
    indices{9,2} = qnt_standard(8);

    % Indices
    
    val20yr = val20yRP(Qsim_obs);
    
    indices{10,2} = val20yr;

    % ECDFs
    
    [F_Q_obs_ecdf, x_Q_obs_ecdf] = ecdf(Qsim_obs(:, end));
    [F_Q_obs95_ecdf, x_Q_obs95_ecdf] = ecdf(Qsim_obs(Qsim_obs(:, end) > qnt_standard(6), end));
    
% Precipitation

P = xobs(:, [1:3,6]);
P(isnan(P)) = 0;
Pobs = P;

    % Statistics (= general overview)
    
    if wet == 1
        qnt_standard = quantile(P(P(:,end)>0.1,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    else
        qnt_standard = quantile(P(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    end
    indices{11,2} = qnt_standard(1); 
    indices{12,2} = qnt_standard(2);
    indices{13,2} = qnt_standard(3);
    indices{14,2} = qnt_standard(4);
    indices{15,2} = qnt_standard(5);
    indices{16,2} = qnt_standard(6);
    indices{17,2} = qnt_standard(7);
    indices{18,2} = qnt_standard(8);
    
    % Indices
    
    trans = TransProb(P(:, end), 0.1);
    Spelllengthobs = SpellDist(P(:, end));
    rx1day = RX1day(P);
    rx5day = RX5day(P);
    sdii =  SDII(P);
    r10 = R10(P);
    r20 = R20(P);
    prcptot = PRCPTOT(P);
    ndry = sum(P(:,end)<0.1);
    
    indices{19,2} = trans(1);
    indices{20,2} = trans(3);
    indices{21,2} = rx1day;
    indices{22,2} = rx5day;
    indices{23,2} = sdii;
    indices{24,2} = r10;
    indices{25,2} = r20;
    indices{26,2} = ndry;
    indices{27,2} = prcptot;
    
    % ECDFs
    
    [F_P_obs_ecdf, x_P_obs_ecdf] = ecdf(P(:, end));
    [F_P_obs95_ecdf, x_P_obs95_ecdf] = ecdf(P(P(:, end) > qnt_standard(6),end));

% Evaporation

E = xobs(:, 1:4);
E(isnan(E)) = 0;
Eobs = E;

    % Statistics (= general overview)
    
    qnt_standard = quantile(E(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    indices{28,2} = qnt_standard(1); 
    indices{29,2} = qnt_standard(2);
    indices{30,2} = qnt_standard(3);
    indices{31,2} = qnt_standard(4);
    indices{32,2} = qnt_standard(5);
    indices{33,2} = qnt_standard(6);
    indices{34,2} = qnt_standard(7);
    indices{35,2} = qnt_standard(8);
    
    % ECDFs
    
    [F_E_obs_ecdf, x_E_obs_ecdf] = ecdf(E(:,end));
    
% Temperature

T = xobs(:, [1:3, 5]);
T(isnan(T)) = 0;
Tobs = T;

        % Statistics (= general overview)
    
        qnt_standard = quantile(T(:, end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
        indices{36,2} = qnt_standard(1); 
        indices{37,2} = qnt_standard(2);
        indices{38,2} = qnt_standard(3);
        indices{39,2} = qnt_standard(4);
        indices{40,2} = qnt_standard(5);
        indices{41,2} = qnt_standard(6);
        indices{42,2} = qnt_standard(7);
        indices{43,2} = qnt_standard(8);
        
        % ECDFs
    
        [F_T_obs_ecdf, x_T_obs_ecdf] = ecdf(T(:, end));

% General

corrPE = corr(P(:, end),E(:,end), 'type', 'Spearman');
corrPT = corr(P(:, end),T(:, end), 'type', 'Spearman');
corrET = corr(E(:,end),T(:, end), 'type', 'Spearman');
autocorrP = xcorr(P(:,end), 1, 'coeff');
crosscorrPE = xcorr(P(:, end),E(:,end), 1, 'coeff');
crosscorrPT = xcorr(P(:, end),T(:, end), 1, 'coeff');
crosscorrET = xcorr(E(:,end),T(:, end), 1, 'coeff');

indices{44,2} = corrPE;
indices{45,2} = corrPT;
indices{46,2} = corrET;
indices{47,2} = autocorrP(3);
indices{48,2} = crosscorrPE(2);
indices{49,2} = crosscorrPE(3);
indices{50,2} = crosscorrPT(2);
indices{51,2} = crosscorrPT(3);
indices{52,2} = crosscorrET(2);
indices{53,2} = crosscorrET(3);

%% ORIGINAL
%Calculation of the index values for the raw simulations

% Discharge

Qsim_xfs = Discharge(xfs, strcat(name, '_xfs'));

    % Statistics (= general overview)
    
    qnt_standard = quantile(Qsim_xfs(:, end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    indices{2,3} = qnt_standard(1); 
    indices{3,3} = qnt_standard(2);
    indices{4,3} = qnt_standard(3);
    indices{5,3} = qnt_standard(4);
    indices{6,3} = qnt_standard(5);
    indices{7,3} = qnt_standard(6);
    indices{8,3} = qnt_standard(7);
    indices{9,3} = qnt_standard(8);

    % Indices
    
    val20yr = val20yRP(Qsim_xfs);
    
    indices{10,3} = val20yr;
    
    % ECDFs
    
    [F_Q_xfs_ecdf, x_Q_xfs_ecdf] = ecdf(Qsim_xfs(:, end));
    [F_Q_xfs95_ecdf, x_Q_xfs95_ecdf] = ecdf(Qsim_xfs(xfs(:, end) > qnt_standard(6), end)); %Dit zal waarschijnlijk een ander quantile hebben? Wat zijn de gevolgen?
    
    % Bar plot
    
    barQ(:,1) = Qsim_xfs(60:84,4);
    
% Precipitation

P = xfs(:, [1:3,6]);
P(isnan(P)) = 0;
Pxfs = P;

    % Statistics (= general overview)
    
    if wet == 1
        qnt_standard = quantile(P(P(:,end)>0.1,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    else
        qnt_standard = quantile(P(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    end
    indices{11,3} = qnt_standard(1); 
    indices{12,3} = qnt_standard(2);
    indices{13,3} = qnt_standard(3);
    indices{14,3} = qnt_standard(4);
    indices{15,3} = qnt_standard(5);
    indices{16,3} = qnt_standard(6);
    indices{17,3} = qnt_standard(7);
    indices{18,3} = qnt_standard(8);
    
    % Indices
    
    trans = TransProb(P(:, end), 0.1);
    Spelllengthxfs = SpellDist(P(:, end));
    rx1day = RX1day(P);
    rx5day = RX5day(P);
    sdii =  SDII(P);
    r10 = R10(P);
    r20 = R20(P);
    prcptot = PRCPTOT(P);
    ndry = sum(P(:,end)<0.1);
    
    indices{19,3} = trans(1);
    indices{20,3} = trans(3);
    indices{21,3} = rx1day;
    indices{22,3} = rx5day;
    indices{23,3} = sdii;
    indices{24,3} = r10;
    indices{25,3} = r20;
    indices{26,3} = ndry;
    indices{27,3} = prcptot;
    
    % ECDFs
    
    [F_P_xfs_ecdf, x_P_xfs_ecdf] = ecdf(P(:, end));
    [F_P_xfs95_ecdf, x_P_xfs95_ecdf] = ecdf(P(P(:, end) > qnt_standard(6), end));
    
    % Bar plot
            
    barP(:,1) = P(60:84,end);
    
% Evaporation

E = xfs(:, 1:4);
E(isnan(E)) = 0;
Exfs = E;

    % Statistics (= general overview)
    
    qnt_standard = quantile(E(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    indices{28,3} = qnt_standard(1); 
    indices{29,3} = qnt_standard(2);
    indices{30,3} = qnt_standard(3);
    indices{31,3} = qnt_standard(4);
    indices{32,3} = qnt_standard(5);
    indices{33,3} = qnt_standard(6);
    indices{34,3} = qnt_standard(7);
    indices{35,3} = qnt_standard(8);
    
    % ECDFs
    
    [F_E_xfs_ecdf, x_E_xfs_ecdf] = ecdf(E(:, end));
    
% Temperature

T = xfs(:, [1:3,5]);
T(isnan(T)) = 0;
Txfs = T;

    % Statistics (= general overview)
    
    qnt_standard = quantile(T(:, end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    indices{36,3} = qnt_standard(1); 
    indices{37,3} = qnt_standard(2);
    indices{38,3} = qnt_standard(3);
    indices{39,3} = qnt_standard(4);
    indices{40,3} = qnt_standard(5);
    indices{41,3} = qnt_standard(6);
    indices{42,3} = qnt_standard(7);
    indices{43,3} = qnt_standard(8);
    
    % ECDFs
    
    [F_T_xfs_ecdf, x_T_xfs_ecdf] = ecdf(T(:, end));

% General (correlation structure etc.)

corrPE = corr(P(:, end),E(:,end), 'type', 'Spearman');
corrPT = corr(P(:, end),T(:, end), 'type', 'Spearman');
corrET = corr(E(:,end),T(:, end), 'type', 'Spearman');
autocorrP = xcorr(P(:,end), 1, 'coeff');
crosscorrPE = xcorr(P(:, end),E(:,end), 1, 'coeff');
crosscorrPT = xcorr(P(:, end),T(:, end), 1, 'coeff');
crosscorrET = xcorr(E(:,end),T(:, end), 1, 'coeff');

indices{44,3} = corrPE;
indices{45,3} = corrPT;
indices{46,3} = corrET;
indices{47,3} = autocorrP(3);
indices{48,3} = crosscorrPE(2);
indices{49,3} = crosscorrPE(3);
indices{50,3} = crosscorrPT(2);
indices{51,3} = crosscorrPT(3);
indices{52,3} = crosscorrET(2);
indices{53,3} = crosscorrET(3);

%% BA METHODS
%Calculation of the index values for the corrected simulations for each selected method 

% Method selection

cnt = 0;
check = 0;

for i=1:length(occ_methods)
    if occ_methods(i) == 1
        name_occ = names_occ{i};
        for j = 1:length(int_methods)
            if int_methods(j) == 1
                if j == 6
                    cnt = cnt + 1;
                    names_methods{cnt} = input('Give the name of the file you want to test: ');
                    check = cnt; %For later purposes
                else
                    name_int = names_int{j};
                    cnt = cnt + 1;
                    names_methods{cnt} = strcat(name, '_', name_occ, '_', name_int, '_results.mat');
                    name_figs{cnt} = strcat(name, '-', name_occ, '-', name_int);
                end
            end
            
        end
    end
end

% Loop over methods

for i = 1:cnt_methods
    
    if i == check
        name_method = input('What should the internal output name of this method be? '); %Name for headers
        datacell = matload(names_methods{i});
        name_figs{check} = name_method;
    else
        name_method = names_methods{i};
        datacell = matload(name_method);
    end
    
    indices(1, i+3)= {name_method};
    bias(1, i+2) = {name_method};
    RB_O(1,i+1) = {name_method};
    RB_MB(1,i+1) = {name_method};

    time = datacell{1,1}(:,1:3);
    
    nrepeats = size(datacell,2);
    
    % Pre-allocation
    
    qnt_standardQ = nan(nrepeats, 8);
    val20yr = nan(nrepeats, 1);
    qnt_standardP = nan(nrepeats, 8);
    trans = nan(nrepeats, 4);
    rx1day = nan(nrepeats, 240); %Juiste instelling?
    rx5day = nan(nrepeats, 240);
    sdii = nan(nrepeats, 20);
    r10 = nan(nrepeats, 20);
    r20 = nan(nrepeats, 20);
    prcptot = nan(nrepeats, 20);
    ndry = nan(nrepeats, 1);
    qnt_standardE = nan(nrepeats, 8);
    qnt_standardT = nan(nrepeats, 8);
    corrPE = nan(nrepeats, 1);
    corrPT = nan(nrepeats, 1);
    corrET = nan(nrepeats, 1);
    autocorrP = nan(nrepeats, 3);
    crosscorrPE = nan(nrepeats, 3);
    crosscorrPT = nan(nrepeats, 3);
    crosscorrET = nan(nrepeats, 3);

    
    for coln = 1:nrepeats
    
        E = [time, datacell{1,coln}(:,4)];
        T = [time, datacell{1,coln}(:,5)];
        P = [time, datacell{1,coln}(:,6)];
        
        if strcmp(name_method, 'MPI-rcp45_threshold_dotc_results.mat') == 1 % Post-processing step for dOTC: negative values should be removed
            P(P(:,end)<0,end) = 0;
            E(E(:,end)<0,end) = 0;
        end
            
        data = [time, E(:, end), T(:, end), P(:, end)]; %Standard method
        
        %For testing: is it better to use only adjusted P?
        %E = xfs(:, 4);
        %E(isnan(E)) = 0;
        %T = xfs(:, 5);
        %T(isnan(T)) = 0;
        %data = [time, E(:,end), T(:, end), P(:, end)];
    
        % Discharge

            Qsim_fc = Discharge(data, name_method);

            % Statistics (= general overview)

            qnt_standardQ(coln, :) = quantile(Qsim_fc(:, end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
        
            % Indices

            val20yr(coln) = val20yRP(Qsim_fc);
            
            % Bar plot
            
            barQ(:,1+i) = Qsim_fc(60:84,4);
    
        % Precipitation
    
            % Statistics (= general overview)
            
            if wet == 1
                qnt_standardP(coln, :) = quantile(P(P(:,end)>0.1,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
            else
                qnt_standardP(coln, :) = quantile(P(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
            end
    
            qnt_standardP(coln, :) = quantile(P(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    
            % Indices
    
            trans(coln, :) = TransProb(P(:, end), 0.1);

            rx1day(coln, :) = RX1day(P);
            rx5day(coln, :) = RX5day(P);
            sdii(coln, :) =  SDII(P);
            r10(coln, :) = R10(P);
            r20(coln, :) = R20(P);
            prcptot(coln, :) = PRCPTOT(P);
            ndry(coln) = sum(P(:,end)<0.1);
            
            % Bar plot
            
            barP(:,1+i) = P(60:84,end);

        % Evaporation
    
            % Statistics (= general overview)
    
            qnt_standardE(coln, :) = quantile(E(:,end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);

        % Temperature
    
            % Statistics (= general overview)
    
            qnt_standardT(coln, :) = quantile(T(:, end), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);

        % General
    
        corrPE(coln) = corr(P(:, end),E(:,end), 'type', 'Spearman', 'rows', 'complete');
        corrPT(coln) = corr(P(:, end),T(:, end), 'type', 'Spearman', 'rows', 'complete');
        corrET(coln) = corr(E(:,end),T(:, end), 'type', 'Spearman', 'rows', 'complete');
        autocorrP(coln, :) = xcorr(P(:, end), 1, 'coeff');
        crosscorrPE(coln, :) = xcorr(P(~isnan(P(:,end))&~isnan(E(:,end)), end),E(~isnan(P(:,end))&~isnan(E(:,end)),end), 1, 'coeff');
        crosscorrPT(coln, :) = xcorr(P(~isnan(P(:,end))&~isnan(T(:,end)), end),T(~isnan(P(:,end))&~isnan(T(:,end)), end), 1, 'coeff');
        crosscorrET(coln, :) = xcorr(E(~isnan(E(:,end))&~isnan(T(:,end)), end),T(~isnan(E(:,end))&~isnan(T(:,end)), end), 1, 'coeff');
    
    end
    
    % Mean index calculation
    
    indices{2,3+i} = mean(qnt_standardQ(:,1));
    indices{3,3+i} = mean(qnt_standardQ(:,2));
    indices{4,3+i} = mean(qnt_standardQ(:,3));
    indices{5,3+i} = mean(qnt_standardQ(:,4));
    indices{6,3+i} = mean(qnt_standardQ(:,5));
    indices{7,3+i} = mean(qnt_standardQ(:,6));
    indices{8,3+i} = mean(qnt_standardQ(:,7));
    indices{9,3+i} = mean(qnt_standardQ(:,8));
    indices{10,3+i} = mean(val20yr);
    indices{11,3+i} = mean(qnt_standardP(:,1));
    indices{12,3+i} = mean(qnt_standardP(:,2));
    indices{13,3+i} = mean(qnt_standardP(:,3));
    indices{14,3+i} = mean(qnt_standardP(:,4));
    indices{15,3+i} = mean(qnt_standardP(:,5));
    indices{16,3+i} = mean(qnt_standardP(:,6));
    indices{17,3+i} = mean(qnt_standardP(:,7));
    indices{18,3+i} = mean(qnt_standardP(:,8));
    indices{19,3+i} = mean(trans(:, 1));
    indices{20,3+i} = mean(trans(:, 3));
    indices{21,3+i} = mean(rx1day, 1);
    indices{22,3+i} = mean(rx5day, 1);
    indices{23,3+i} = mean(sdii, 1);
    indices{24,3+i} = mean(r10, 1);
    indices{25,3+i} = mean(r20, 1);
    indices{26,3+i} = mean(ndry);
    indices{27,3+i} = mean(prcptot, 1);
    indices{28,3+i} = mean(qnt_standardE(:,1));
    indices{29,3+i} = mean(qnt_standardE(:,2));
    indices{30,3+i} = mean(qnt_standardE(:,3));
    indices{31,3+i} = mean(qnt_standardE(:,4));
    indices{32,3+i} = mean(qnt_standardE(:,5));
    indices{33,3+i} = mean(qnt_standardE(:,6));
    indices{34,3+i} = mean(qnt_standardE(:,7));
    indices{35,3+i} = mean(qnt_standardE(:,8));
    indices{36,3+i} = mean(qnt_standardT(:,1));
    indices{37,3+i} = mean(qnt_standardT(:,2));
    indices{38,3+i} = mean(qnt_standardT(:,3));
    indices{39,3+i} = mean(qnt_standardT(:,4));
    indices{40,3+i} = mean(qnt_standardT(:,5));
    indices{41,3+i} = mean(qnt_standardT(:,6));
    indices{42,3+i} = mean(qnt_standardT(:,7));
    indices{43,3+i} = mean(qnt_standardT(:,8));
    indices{44,3+i} = mean(corrPE);
    indices{45,3+i} = mean(corrPT);
    indices{46,3+i} = mean(corrET);
    indices{47,3+i} = mean(autocorrP(:, 3));
    indices{48,3+i} = mean(crosscorrPE(:, 2));
    indices{49,3+i} = mean(crosscorrPE(:, 3));
    indices{50,3+i} = mean(crosscorrPT(:, 2));
    indices{51,3+i} = mean(crosscorrPT(:, 3));
    indices{52,3+i} = mean(crosscorrET(:, 2));
    indices{53,3+i} = mean(crosscorrET(:, 3));
    
    % Discharge ECDFs
    
    [F_Q_fc_ecdf, x_Q_fc_ecdf] = ecdf(Qsim_fc(:, end));
    [F_Q_fc95_ecdf, x_Q_fc95_ecdf] = ecdf(Qsim_fc(Qsim_fc(:, end) > indices{7,3+i}, end));
        
    f1 = figure(1);
    plot(x_Q_fc_ecdf, F_Q_fc_ecdf, 'LineWidth', 2)
    title('Empirical CDF for discharge')
    xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    hold on
    
    f2 = figure(2);
    plot(x_Q_fc95_ecdf, F_Q_fc95_ecdf, 'LineWidth', 2)
    title('Empirical CDF for discharge larger than 95th quantile')
    xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    hold on
    
    % Discharge PDFs
    
    f16 = figure(16);
    histogram(Qsim_fc(:, end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
    f17 = figure(17);
    histogram(Qsim_fc(Qsim_fc(:, end) > indices{7,3+i}, end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
    % Precipitation ECDFs
    
    [F_P_fc_ecdf, x_P_fc_ecdf] = ecdf(P(:, end));
    [F_P_fc95_ecdf, x_P_fc95_ecdf] = ecdf(P(P(:, end) > indices{16,3+i}, end));
    
    f3 = figure(3);
    plot(x_P_fc_ecdf, F_P_fc_ecdf, 'LineWidth', 2)
    title('Empirical CDF for precipitation', 'Interpreter', 'Latex')
    xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    hold on
    
    f4 = figure(4);
    plot(x_P_fc95_ecdf, F_P_fc95_ecdf, 'LineWidth', 2)
    title('Empirical CDF for precipitation larger than 95th quantile', 'Interpreter', 'Latex')
    xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    hold on
    
    % Precipitation PDFs
    
    Spelllength = SpellDist(P(:, end));
    
    f8 = figure(8);
    histogram(P(:, end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
    f9 = figure(9);
    histogram(P(P(:, end) > indices{16,3+i}, end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
    f7 = figure(7);
    histogram(Spelllength, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
    % Evaporation ECDFs
    
    [F_E_fc_ecdf, x_E_fc_ecdf] = ecdf(E(:, end));
    
    f5 = figure(5);
    plot(x_E_fc_ecdf, F_E_fc_ecdf, 'LineWidth', 2)
    title('Empirical CDF of evaporation', 'Interpreter', 'Latex')
    xlabel('Evaporation (mm)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    hold on
    
    % Evaporation PDFs
    
    f10 = figure(10);
    histogram(E(:, end), 'BinWidth', 0.25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
    % Temperature ECDFs
    
    [F_T_fc_ecdf, x_T_fc_ecdf] = ecdf(T(:, end));
    
    f6 = figure(6);
    plot(x_T_fc_ecdf, F_T_fc_ecdf, 'LineWidth', 2)
    title('Empirical CDF for temperature', 'Interpreter', 'Latex')
    xlabel('Temperature (°C)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    hold on
    
    % Temperature PDFs
    
    f11 = figure(11);
    histogram(T(:, end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
    hold on
    
end

%% COMBINATION

% Pictures

% Finishing ECDF plots

lgnd1 = [name_figs, 'Raw RCM simulations', 'Observations'];
lgndbar = ['Raw RCM simulations', name_figs];
lgndbox = ['Observations', 'Raw simulations',  name_figs];

figure(1)
plot(x_Q_xfs_ecdf, F_Q_xfs_ecdf, 'LineWidth', 4) 
plot(x_Q_obs_ecdf, F_Q_obs_ecdf, 'LineWidth', 4)
legend(lgnd1, 'Location', 'southeast', 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f1);
saveas(f1, strcat(save_loc_fig, 'Qcdf.png'))
hold off

figure(2)
plot(x_Q_xfs95_ecdf, F_Q_xfs95_ecdf, 'LineWidth', 4) 
plot(x_Q_obs95_ecdf, F_Q_obs95_ecdf, 'LineWidth', 4)
legend(lgnd1, 'Location', 'southeast', 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f2);
saveas(f2, strcat(save_loc_fig, 'Q95cdf.png'))
hold off

figure(3)
plot(x_P_xfs_ecdf, F_P_xfs_ecdf, 'LineWidth', 4) 
plot(x_P_obs_ecdf, F_P_obs_ecdf, 'LineWidth', 4)
legend(lgnd1, 'Location', 'southeast', 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f3);
saveas(f3, strcat(save_loc_fig, 'Pcdf.png'))
hold off

figure(4)
plot(x_P_xfs95_ecdf, F_P_xfs95_ecdf, 'LineWidth', 4) 
plot(x_P_obs95_ecdf, F_P_obs95_ecdf, 'LineWidth', 4)
legend(lgnd1, 'Location', 'southeast', 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f4);
saveas(f4, strcat(save_loc_fig, 'P95cdf.png'))
hold off

figure(5)
plot(x_E_xfs_ecdf, F_E_xfs_ecdf, 'LineWidth', 4) 
plot(x_E_obs_ecdf, F_E_obs_ecdf, 'LineWidth', 4)
legend(lgnd1, 'Location', 'southeast', 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f5);
saveas(f5, strcat(save_loc_fig, 'Ecdf.png'))
hold off

figure(6)
plot(x_T_xfs_ecdf, F_T_xfs_ecdf, 'LineWidth', 4) 
plot(x_T_obs_ecdf, F_T_obs_ecdf, 'LineWidth', 4)
legend(lgnd1, 'Location', 'southeast', 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f6);
saveas(f6, strcat(save_loc_fig, 'Tcdf.png'))
hold off

% Finishing other pictures

figure(7)
histogram(Spelllengthxfs, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Spelllengthobs, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Wet spell length (days)', 'Interpreter', 'Latex')
ylabel('Probability density function', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
%set(gca,'FontName', 'UGent Panno Text')
fullfig(f7)
saveas(f7, strcat(save_loc_fig, 'Spelllength.png'))
hold off

figure(8)
histogram(Pxfs(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Pobs(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
ylabel('Probability density', 'Interpreter', 'Latex')
title('PDF of precipitation', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f8)
saveas(f8, strcat(save_loc_fig, 'Ppdf.png'))
hold off

figure(9)
histogram(Pxfs(Pxfs(:, end) > indices{16,3},end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Pobs(Pobs(:, end) > indices{16,2},end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
ylabel('Probability density', 'Interpreter', 'Latex')
title('PDF of precipitation > 95th quantile', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f9)
saveas(f9, strcat(save_loc_fig, 'P95pdf.png'))
hold off

figure(10)
histogram(Exfs(:,end), 'BinWidth', 0.25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Eobs(:,end), 'BinWidth', 0.25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Evaporation (mm)', 'Interpreter', 'Latex')
ylabel('Probability density', 'Interpreter', 'Latex')
title('PDF of evapotranspiration', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f10)
saveas(f10, strcat(save_loc_fig, 'Epdf.png'))
hold off

figure(11)
histogram(Txfs(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Tobs(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Temperature (°C)', 'Interpreter', 'Latex')
ylabel('Probability density', 'Interpreter', 'Latex')
title('PDF of Temperature', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f11)
saveas(f11, strcat(save_loc_fig, 'Tpdf.png'))
hold off

figure(16)
histogram(Qsim_xfs(:,end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Qsim_obs(:,end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
ylabel('Probability density', 'Interpreter', 'Latex')
title('PDF of discharge', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f16)
saveas(f16, strcat(save_loc_fig, 'Qpdf.png'))
hold off

figure(17)
histogram(Qsim_xfs(Qsim_xfs(:, end) > indices{7,3},end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
histogram(Qsim_obs(Qsim_obs(:, end) > indices{7,2},end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
ylabel('Probability density', 'Interpreter', 'Latex')
title('PDF of discharge > 95th quantile', 'Interpreter', 'Latex')
legend(lgnd1, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f17)
saveas(f17, strcat(save_loc_fig, 'Q95pdf.png'))
hold off

% Box plots

RX1daycomb = zeros(240, cnt_methods+2);
RX5daycomb = zeros(240, cnt_methods+2);
SDIIcomb = zeros(20, cnt_methods+2);
PRCPTOTcomb = zeros(20, cnt_methods+2);
R10comb = zeros(20, cnt_methods+2);
R20comb = zeros(20, cnt_methods+2);
for i = 1:cnt_methods+2
    RX1daycomb(:, i) = indices{21,1+i};
    RX5daycomb(:, i) = indices{22,1+i};
    SDIIcomb(:,i) = indices{23,1+i};
    PRCPTOTcomb(:,i) = indices{27,1+i};
    R10comb(:,i) = indices{24,1+i};
    R20comb(:,i) = indices{25,1+i};
end

f12 = figure(12);
boxplot(RX1daycomb, 'Labels', lgndbox, 'LabelOrientation', 'inline')
ylabel('Precipitation (mm)', 'Interpreter', 'Latex')
xlabel('Methods', 'Interpreter', 'Latex')
title('Box plots of monthly 1-day maximum precipitation amounts for each BC method', 'Interpreter', 'Latex')
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f12)
saveas(f12, strcat(save_loc_fig, 'BoxRX1day.png'))

f13 = figure(13);
boxplot(RX5daycomb, 'Labels', lgndbox, 'LabelOrientation', 'inline')
ylabel('Precipitation (mm)', 'Interpreter', 'Latex')
xlabel('Methods', 'Interpreter', 'Latex')
title('Box plots of monthly consecutive 5-day maximum precipitation amounts for each BC method', 'Interpreter', 'Latex')
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f13)
saveas(f13, strcat(save_loc_fig, 'BoxRX5day.png'))

f14 = figure(14);
boxplot(SDIIcomb, 'Labels', lgndbox, 'LabelOrientation', 'inline')
ylabel('Intensity index (mm/day)', 'Interpreter', 'Latex')
xlabel('Methods', 'Interpreter', 'Latex')
title('Box plots of yearly SDII indices for each BC method', 'Interpreter', 'Latex')
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f14)
saveas(f14, strcat(save_loc_fig, 'BoxSDII.png'))

f15 = figure(15);
boxplot(PRCPTOTcomb, 'Labels', lgndbox, 'LabelOrientation', 'inline')
ylabel('Annual precipitation (mm)', 'Interpreter', 'Latex')
xlabel('Methods', 'Interpreter', 'Latex')
title('Box plots of annual precipitation for each BC method', 'Interpreter', 'Latex')
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f15)
saveas(f15, strcat(save_loc_fig, 'BoxPRCPTOT.png'))

f18 = figure(18);
boxplot(R10comb, 'Labels', lgndbox, 'LabelOrientation', 'inline')
ylabel('Number of days', 'Interpreter', 'Latex')
xlabel('Methods', 'Interpreter', 'Latex')
title('Box plots of annual number of heavy precipitation days', 'Interpreter', 'Latex')
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f18)
saveas(f18, strcat(save_loc_fig, 'BoxR10.png'))

f19 = figure(19);
boxplot(R20comb, 'Labels', lgndbox, 'LabelOrientation', 'inline')
ylabel('Number of days', 'Interpreter', 'Latex')
xlabel('Methods', 'Interpreter', 'Latex')
title('Box plots of annual number of very heavy precipitation days', 'Interpreter', 'Latex')
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f19)
saveas(f19, strcat(save_loc_fig, 'BoxR20.png'))

f20 = figure(20);
x = 1:1:25;
bar(x, barQ)
xlabel('Days', 'Interpreter', 'Latex')
ylabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
legend(lgndbar, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f20)
saveas(f20, strcat(save_loc_fig, 'BarQ.png'))

f21 = figure(21);
x = 1:1:25;
bar(x, barP)
xlabel('Days', 'Interpreter', 'Latex')
ylabel('Precipitation (mm)', 'Interpreter', 'Latex')
legend(lgndbar, 'FontSize', 30)
set(gca,'FontSize',30)
set(gca,'FontName', 'UGent Panno Text')
fullfig(f21)
saveas(f21, strcat(save_loc_fig, 'BarP.png'))
    

% Bias

for i = 3:cnt_methods + 3
    for j = 2:num_indices+1
        if j ~= (21 && 22 && 23 && 24 && 25 && 27) 
            %Difference
            bias{j,i-1} = indices{j,i}-indices{j,2};
        end
    end
end

for i = 3:cnt_methods + 2
    for j = 2:num_indices+1
        if j ~= (21 && 22 && 23 && 24 && 25 && 27)
            RB_O{j,i-1} = 1 - (abs(bias{j,2})-abs(bias{j,i}))/abs(indices{j,2});
            RB_MB{j,i-1} = 1 - (abs(bias{j,2})-abs(bias{j,i}))/abs(bias{j,2});
        end
    end
end

%% Opslaan

save(strcat(save_loc, name, save_name, '_indices.mat'), 'indices')
save(strcat(save_loc, name, save_name, '_biases.mat'), 'bias')
save(strcat(save_loc, name, save_name, '_RB_O.mat'), 'RB_O')
save(strcat(save_loc, name, save_name, '_RB_MB.mat'), 'RB_MB')



end



