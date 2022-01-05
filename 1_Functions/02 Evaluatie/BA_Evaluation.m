function [indices, bias, RB_O, RB_MB] = BA_Evaluation(xfs, xobs, name, save_name, occ_methods, int_methods, wet, months, NameTimes)
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
%       'dotc', 6: 'r2d2', 7: 'other' -> 'other' is for evaluation of
%       special cases
%       wet: Boolean indicating whether the P quantiles have to be
%       calculated on wet days only (1) or on all days (0)
%       months: indicates which months have to be pooled for the
%       evaluation. A [ixj] matrix, with with i the number of timeframes
%       and j the number of months in the time frame. Usual combinations
%       are 3x4 (seasonal evaluation), 1x12 (yearly evaluation), 2x6 (false
%       winter vs summer evaluation) and 12x1 (monthly evaluation)
%       NameTimes: cell with the name for each Period, e.g. {'Summer',
%       'Winter'). Order should correspond to order in months matrix.
%   Outputs:
%       indices: cell with indices calculated for observations, the original RCM simulations and every method
%       bias: cell with the biases (= difference) between 1) the observations and 2) the original RCM simulations and the methods used
%       RB_O: Residual Bias relative to the Observations
%       RB_MB: Residual Bias relative to the Model Bias
%
%   Last update by J. Van de Velde on 11/02/'21
 
%% Set-up

names_occ = {'none', 'ssr', 'tdc', 'threshold'};
names_int = {'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc', 'r2d2'};
cnt_methods = sum(occ_methods)*sum(int_methods);
names_methods = cell(1, cnt_methods);
name_figs = cell(1, cnt_methods);
nTimes = length(NameTimes);

num_indices = 52;
indices = cell(num_indices+1, cnt_methods+3, nTimes);
bias = cell(num_indices+1, cnt_methods+2, nTimes);
RB_O = cell(num_indices+1, cnt_methods+1, nTimes);
RB_MB = cell(num_indices+1, cnt_methods+1, nTimes);
barQ = cell(cnt_methods +1,nTimes);
barP = cell(cnt_methods +1,nTimes);
F_Q_obs_ecdf = cell(1,nTimes);
x_Q_obs_ecdf = cell(1,nTimes);
F_Q_obs95_ecdf = cell(1,nTimes);
x_Q_obs95_ecdf = cell(1,nTimes); 
F_P_obs_ecdf = cell(1,nTimes); 
x_P_obs_ecdf = cell(1,nTimes); 
F_P_obs95_ecdf = cell(1,nTimes); 
x_P_obs95_ecdf = cell(1,nTimes);
F_E_obs_ecdf = cell(1,nTimes);
x_E_obs_ecdf = cell(1,nTimes); 
F_T_obs_ecdf = cell(1,nTimes); 
x_T_obs_ecdf = cell(1,nTimes); 
Spelllengthobs = cell(1,nTimes); 
Pobs = cell(1,nTimes); 
Qsim_obs = cell(1,nTimes); 
Eobs = cell(1,nTimes); 
Tobs = cell(1,nTimes);
F_Q_xfs_ecdf = cell(1,nTimes);
x_Q_xfs_ecdf = cell(1,nTimes);
F_Q_xfs95_ecdf = cell(1,nTimes);
x_Q_xfs95_ecdf = cell(1,nTimes); 
F_P_xfs_ecdf = cell(1,nTimes); 
x_P_xfs_ecdf = cell(1,nTimes); 
F_P_xfs95_ecdf = cell(1,nTimes); 
x_P_xfs95_ecdf = cell(1,nTimes);
F_E_xfs_ecdf = cell(1,nTimes);
x_E_xfs_ecdf = cell(1,nTimes); 
F_T_xfs_ecdf = cell(1,nTimes); 
x_T_xfs_ecdf = cell(1,nTimes); 
Spelllengthxfs = cell(1,nTimes); 
Pxfs = cell(1,nTimes); 
Qsim_xfs = cell(1,nTimes); 
Exfs = cell(1,nTimes); 
Txfs = cell(1,nTimes);
F_Q_fc_ecdf = cell(cnt_methods,nTimes);
x_Q_fc_ecdf = cell(cnt_methods,nTimes);
F_Q_fc95_ecdf = cell(cnt_methods,nTimes);
x_Q_fc95_ecdf = cell(cnt_methods,nTimes); 
F_P_fc_ecdf = cell(cnt_methods,nTimes); 
x_P_fc_ecdf = cell(cnt_methods,nTimes); 
F_P_fc95_ecdf = cell(cnt_methods,nTimes); 
x_P_fc95_ecdf = cell(cnt_methods,nTimes);
F_E_fc_ecdf = cell(cnt_methods,nTimes);
x_E_fc_ecdf = cell(cnt_methods,nTimes); 
F_T_fc_ecdf = cell(cnt_methods,nTimes); 
x_T_fc_ecdf = cell(cnt_methods,nTimes); 
Spelllength_fc = cell(1,nTimes); 
Qsim_fc = cell(cnt_methods,nTimes); 
Qsimfc95 = cell(cnt_methods,nTimes); 
Phist = cell(cnt_methods,nTimes); 
Phist95 = cell(cnt_methods,nTimes); 
Ehist = cell(cnt_methods,nTimes); 
Thist = cell(cnt_methods,nTimes);
PSS_P = cell(cnt_methods,nTimes);
PSS_T = cell(cnt_methods,nTimes);
PSS_E = cell(cnt_methods,nTimes);
PSS_Q = cell(cnt_methods,nTimes);

for m=1:nTimes
    indices(2:end, 1, m) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
    indices(1, 2:3, m) = {'Observations', 'Original RCM simulations'};
    bias(2:end,1, m) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
    bias(1,2,m) = {'Original RCM simulations'};
    RB_O(2:end,1, m) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
    RB_MB(2:end,1,m) = {'Q-qnt5'; 'Q-qnt25';'Q-qnt50';'Q-qnt75';'Q-qnt90';'Q-qnt95';'Q-qnt99';'Q-qnt995'; 'Q-val20yr'; 'P-qnt5'; 'P-qnt25';'P-qnt50';'P-qnt75';'P-qnt90';'P-qnt95';'P-qnt99';'P-qnt995'; 'P-trans00'; 'P-trans10'; 'RX1day'; 'RX5day'; 'SDII'; 'R10'; 'R20'; 'ndry'; 'PRCPTOT'; 'E-qnt5'; 'E-qnt25';'E-qnt50';'E-qnt75';'E-qnt90';'E-qnt95';'E-qnt99';'E-qnt995'; 'T-qnt5'; 'T-qnt25';'T-qnt50';'T-qnt75';'T-qnt90';'T-qnt95';'T-qnt99';'T-qnt995'; 'corrPE'; 'corrPT'; 'corrET'; 'autocorrP1'; 'crosscorrPE0'; 'crosscorrPE1'; 'crosscorrPT0'; 'crosscorrPT1'; 'crosscorrET0'; 'crosscorrET1'};
end

save_loc = 'E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\';
save_loc_fig = 'D:\Users\jpvdveld\Documents\PhD\FigurenResultaten\Evaluation\';


%% OBSERVATIONS
%Calculation of the index values for the observations

for m = 1:nTimes
    xobsm = xobs(ismember(xobs(:,2),months(m,:)),:);
    
    [indices(2:end,2,m), F_Q_obs_ecdf{m}, x_Q_obs_ecdf{m}, F_Q_obs95_ecdf{m}, ...
        x_Q_obs95_ecdf{m}, F_P_obs_ecdf{m}, x_P_obs_ecdf{m}, F_P_obs95_ecdf{m}, x_P_obs95_ecdf{m}, ...
        F_E_obs_ecdf{m}, x_E_obs_ecdf{m}, F_T_obs_ecdf{m}, x_T_obs_ecdf{m}, Spelllengthobs{m}, Pobs{m}, Qsim_obs{m}, Eobs{m}, Tobs{m}, ~, ~] = IndexCalcStandard(xobsm, wet, num_indices, strcat(name, '_xfs'), 1, months(m,:));
    
end

%% ORIGINAL
%Calculation of the index values for the raw simulations

for m = 1:nTimes
    xfsm = xfs(ismember(xfs(:,2),months(m,:)),:);
    
    [indices(2:end,3,m), F_Q_xfs_ecdf{m}, x_Q_xfs_ecdf{m}, F_Q_xfs95_ecdf{m}, ...
        x_Q_xfs95_ecdf{m}, F_P_xfs_ecdf{m}, x_P_xfs_ecdf{m}, F_P_xfs95_ecdf{m}, x_P_xfs95_ecdf{m}, ...
        F_E_xfs_ecdf{m}, x_E_xfs_ecdf{m}, F_T_xfs_ecdf{m}, x_T_xfs_ecdf{m}, Spelllengthxfs{m}, Pxfs{m}, Qsim_xfs{m}, Exfs{m}, Txfs{m}, barQ{1,m}, barP{1,m}] = IndexCalcStandard(xfsm, wet, num_indices, strcat(name, '_xfs'), 1, months(m,:));
    
end
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
                if j == 7
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
    
    for m = 1:nTimes
        indices(1, i+3, m)= {name_method};
        bias(1, i+2, m) = {name_method};
        RB_O(1,i+1,m) = {name_method};
        RB_MB(1,i+1,m) = {name_method};
    end

    if strcmp(name_method, [name, '_threshold_r2d2_results.mat']) == 1 % R2D2 is run in R, thus data loading is different
        time = xfs(:,1:3);
        nrepeats =  numel(fieldnames(datacell));
        fields = fieldnames(datacell); % extract names of features
    else
        time = datacell{1,1}(:,1:3);
        nrepeats = size(datacell,2);
    end
    
    data = nan(length(time), 6, nrepeats);
    E = nan(length(time), 4, nrepeats);
    T = nan(length(time), 4, nrepeats);
    P = nan(length(time), 4, nrepeats);

    
    for coln = 1:nrepeats
        if strcmp(name_method, [name, '_threshold_r2d2_results.mat']) == 1
            fieldn = getfield(datacell, fields{coln});
            E(:,:, coln) = [time, fieldn(:,1)];
            T(:,:, coln) = [time, fieldn(:,2)];
            P(:,:, coln) = [time, fieldn(:,3)];
        else
            E(:,:, coln) = [time, datacell{1,coln}(:,4)];
            T(:,:, coln) = [time, datacell{1,coln}(:,5)];
            P(:,:, coln) = [time, datacell{1,coln}(:,6)];
        end
        
        if strcmp(name_method, [name, '_threshold_dotc_results.mat']) == 1 % Post-processing step for dOTC: negative values should be removed
            P(P(:,end, coln)<0,end, coln) = 0;
            E(E(:,end, coln)<0,end, coln) = 0;
        end

        if strcmp(name_method, [name, '_threshold_mrqnbc_results.mat']) == 1 % Post-processing step for MRQNBC: negative values should be removed
            P(P(:,end, coln)<0,end, coln) = 0;
            E(E(:,end, coln)<0,end, coln) = 0;
        end
            
        data(:,:, coln) = [time, E(:, end, coln), T(:, end, coln), P(:, end, coln)]; %Standard method
        
    end
        
        %For testing: is it better to use only adjusted P?
        %E = xfs(:, 4);
        %E(isnan(E)) = 0;
        %T = xfs(:, 5);
        %T(isnan(T)) = 0;
        %data = [time, E(:,end), T(:, end), P(:, end)];
        
        for m = 1:nTimes
            
            datam = data(ismember(data(:,2), months(m,:)),:,:);
            Em = data(ismember(E(:,2),months(m,:)),1:4,:);
            Tm = data(ismember(T(:,2),months(m,:)),[1:3,5],:);
            Pm = data(ismember(P(:,2),months(m,:)),[1:3,6],:);
            
            %Index calculation
            [indices(2:end,3+i,m), F_Q_fc_ecdf{i, m}, x_Q_fc_ecdf{i,m}, F_Q_fc95_ecdf{i,m}, ...
                x_Q_fc95_ecdf{i,m}, F_P_fc_ecdf{i,m}, x_P_fc_ecdf{i,m}, F_P_fc95_ecdf{i,m}, x_P_fc95_ecdf{i,m}, ...
                F_E_fc_ecdf{i,m}, x_E_fc_ecdf{i,m}, F_T_fc_ecdf{i,m}, x_T_fc_ecdf{i,m}, ~, ~, Qsim_fc{i,m}, ~, ~, barQ{1+i,m}, barP{1+i,m}] = IndexCalcStandard(datam, wet, num_indices, strcat(name, '_xfc'), nrepeats, months(m,:), Em, Tm, Pm);
            
            PSS_P{i,m} =  Perkins(Pm, Pobs{m}, nrepeats);
            PSS_T{i,m} =  Perkins(Tm, Tobs{m}, nrepeats);
            PSS_E{i,m} =  Perkins(Em, Eobs{m}, nrepeats);
            PSS_Q{i,m} =  Perkins(Qsim_fc{i,m}, Qsim_obs{m}, nrepeats);
            
            %Extra calculations
            
            Qsimfc95{i,m} = Qsim_fc{i,m}(Qsim_fc{i,m}(:, end) > indices{7,3+i,m}, end);
            Spelllength_fc{i,m} = SpellDist(Pm(:, end));
            Phist{i,m} = Pm(:, end, coln);
            Phist95{i,m} = Pm(Pm(:, end, coln) > indices{16,3+i,m},end, coln);
            Ehist{i,m} = Em(:, end, coln);
            Thist{i,m} = Tm(:, end, coln);
            
        end
    
end

%% COMBINATION

% Pictures

% Finishing ECDF plots

lgnd1 = [name_figs, 'Raw RCM simulations', 'Observations'];
lgndbar = ['Raw RCM simulations', name_figs];
lgndbox = {'(a)', '(b)',  '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'};

for m = 1: nTimes
    
    for i = 1:cnt_methods
        f1 = figure(1);
        plot(x_Q_fc_ecdf{i,m}, F_Q_fc_ecdf{i,m}, 'LineWidth', 2)
        hold on
        
        f2 = figure(2);
        plot(x_Q_fc95_ecdf{i,m}, F_Q_fc95_ecdf{i,m}, 'LineWidth', 2)
        hold on
        
        f3 = figure(3);
        histogram(Qsim_fc{i,m}(:, end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        f4 = figure(4);
        histogram(Qsimfc95{i,m}, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        f5 = figure(5);
        plot(x_P_fc_ecdf{i,m}, F_P_fc_ecdf{i,m}, 'LineWidth', 2)
        hold on
        
        f6 = figure(6);
        plot(x_P_fc95_ecdf{i,m}, F_P_fc95_ecdf{i,m}, 'LineWidth', 2)
        hold on
        
        f7 = figure(7);
        histogram(Phist{i,m}, 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        f8 = figure(8);
        histogram(Phist95{i,m}, 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        f9 = figure(9);
        histogram(Spelllength_fc{i,m}, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        f10 = figure(10);
        plot(x_E_fc_ecdf{i,m}, F_E_fc_ecdf{i,m}, 'LineWidth', 2)
        hold on
        
        f11 = figure(11);
        histogram(Ehist{i,m}, 'BinWidth', 0.25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        f12 = figure(12);
        plot(x_T_fc_ecdf{i,m}, F_T_fc_ecdf{i,m}, 'LineWidth', 2)
        hold on
        
        f13 = figure(13);
        histogram(Thist{i,m}, 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2)
        hold on
        
        barQplot(i+1,:) = barQ{i+1,m};
%         f20 = figure(20);
%         x = 1:1:25;
%         bar(x, barQ{i+1,m})
%         hold on

        barPplot(i+1,:) = barP{i+1,m};
%         f21 = figure(21);
%         x = 1:1:25;
%         bar(x, barP{i+1,m})
%         hold on
    end
    
    figure(1)
    plot(x_Q_xfs_ecdf{m}, F_Q_xfs_ecdf{m}, 'LineWidth', 4)
    plot(x_Q_obs_ecdf{m}, F_Q_obs_ecdf{m}, 'LineWidth', 4)
        title(['Empirical CDF for discharge for ' NameTimes(m)],'Interpreter', 'Latex')
    xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f1);
    saveas(f1, strcat(save_loc_fig, ['Qcdf' NameTimes{m} '.png']))
    hold off
    
    figure(2)
    hold on
    plot(x_Q_xfs95_ecdf{m}, F_Q_xfs95_ecdf{m}, 'LineWidth', 4)
    plot(x_Q_obs95_ecdf{m}, F_Q_obs95_ecdf{m}, 'LineWidth', 4)
    %title(['Empirical CDF for discharge larger than 95th quantile for ' NameTimes(m)], 'Interpreter', 'Latex')
    xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f2);
    saveas(f2, strcat(save_loc_fig, ['Q95cdf' NameTimes{m} '.png']))
    hold off
    
    figure(3)
    histogram(Qsim_xfs{m}(:,end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Qsim_obs{m}(:,end), 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    ylabel('Probability density', 'Interpreter', 'Latex')
    %title(['PDF of discharge for ' NameTimes(m)], 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f3)
    saveas(f3, strcat(save_loc_fig, ['Qpdf' NameTimes{m} '.png']))
    hold off
    
    Qsim_xfshist = Qsim_xfs{m}(Qsim_xfs{m}(:, end) > indices{7,3},end);
    Qsim_obshist = Qsim_obs{m}(Qsim_obs{m}(:, end) > indices{7,2},end);
    
    figure(4)
    histogram(Qsim_xfshist, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Qsim_obshist, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    ylabel('Probability density', 'Interpreter', 'Latex')
    %title(['PDF of discharge $>$ 95th quantile for ' NameTimes(m)], 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f4)
    saveas(f4, strcat(save_loc_fig, ['Q95pdf' NameTimes{m} '.png']))
    hold off
    
    figure(5)
    plot(x_P_xfs_ecdf{m}, F_P_xfs_ecdf{m}, 'LineWidth', 4)
    plot(x_P_obs_ecdf{m}, F_P_obs_ecdf{m}, 'LineWidth', 4)
    %title(['Empirical CDF for precipitation for ' NameTimes(m)], 'Interpreter', 'Latex')
    xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'Location', 'southeast', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f5);
    saveas(f5, strcat(save_loc_fig, ['Pcdf' NameTimes{m} '.png']))
    hold off
    
    figure(6)
    plot(x_P_xfs95_ecdf{m}, F_P_xfs95_ecdf{m}, 'LineWidth', 4)
    plot(x_P_obs95_ecdf{m}, F_P_obs95_ecdf{m}, 'LineWidth', 4)
    %title(['Empirical CDF for precipitation $>$ 95th quantile for ' NameTimes(m)], 'Interpreter', 'Latex')
    xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f6);
    saveas(f6, strcat(save_loc_fig, ['P95cdf' NameTimes{m} '.png']))
    hold off
    
    figure(7)
    histogram(Pxfs{m}(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Pobs{m}(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
    ylabel('Probability density', 'Interpreter', 'Latex')
    %title(['PDF of precipitation for ' NameTimes(m)], 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f7)
    saveas(f7, strcat(save_loc_fig, ['Ppdf' NameTimes{m} '.png']))
    hold off
    
    Pxfs_hist = Pxfs{m}(Pxfs{m}(:, end) > indices{16,3},end);
    Pobs_hist = Pobs{m}(Pobs{m}(:, end) > indices{16,2},end);
    
    figure(8)
    histogram(Pxfs_hist, 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Pobs_hist, 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Precipitation (mm)', 'Interpreter', 'Latex')
    ylabel('Probability density', 'Interpreter', 'Latex')
    %title(['PDF of precipitation $>$ 95th quantile for ' NameTimes(m)], 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f8)
    saveas(f8, strcat(save_loc_fig, ['P95pdf' NameTimes{m} '.png']))
    hold off
    
    figure(9)
    histogram(Spelllengthxfs{m}, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Spelllengthobs{m}, 'BinWidth', 1, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Wet spell length (days)', 'Interpreter', 'Latex')
    ylabel('Probability density function', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f9)
    saveas(f9, strcat(save_loc_fig, ['Spelllength' NameTimes{m} '.png']))
    hold off
    
    figure(10)
    plot(x_E_xfs_ecdf{m}, F_E_xfs_ecdf{m}, 'LineWidth', 4)
    plot(x_E_obs_ecdf{m}, F_E_obs_ecdf{m}, 'LineWidth', 4)
    %title(['Empirical CDF of evaporation for ' NameTimes(m)], 'Interpreter', 'Latex')
    xlabel('Evaporation (mm)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f10);
    saveas(f10, strcat(save_loc_fig, ['Ecdf' NameTimes{m} '.png']))
    hold off
    
    figure(11)
    histogram(Exfs{m}(:,end), 'BinWidth', 0.25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Eobs{m}(:,end), 'BinWidth', 0.25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Evaporation (mm)', 'Interpreter', 'Latex')
    ylabel('Probability density', 'Interpreter', 'Latex')
    %title(['PDF of evapotranspiration for ' NameTimes(m)], 'Interpreter', 'Latex')
    legend(lgnd1,'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f11)
    saveas(f11, strcat(save_loc_fig, ['Epdf' NameTimes{m} '.png']))
    hold off
    
    figure(12)
    plot(x_T_xfs_ecdf{m}, F_T_xfs_ecdf{m}, 'LineWidth', 4)
    plot(x_T_obs_ecdf{m}, F_T_obs_ecdf{m}, 'LineWidth', 4)
    title(['Empirical CDF for temperature for ' NameTimes(m)], 'Interpreter', 'Latex')
    xlabel('Temperature (°C)', 'Interpreter', 'Latex')
    ylabel('Cumulative probability', 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f12);
    saveas(f12, strcat(save_loc_fig, ['Tcdf' NameTimes{m} '.png']))
    hold off
    
    figure(13)
    histogram(Txfs{m}(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    histogram(Tobs{m}(:,end), 'BinWidth', 2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 4)
    xlabel('Temperature (°C)', 'Interpreter', 'Latex')
    ylabel('Probability density', 'Interpreter', 'Latex')
    %title(['PDF of Temperature for ' NameTimes(m)], 'Interpreter', 'Latex')
    legend(lgnd1, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f13)
    saveas(f13, strcat(save_loc_fig, ['Tpdf' NameTimes{m} '.png']))
    hold off
    
    % Box plots
    
    %RX1daycomb = zeros(240, cnt_methods+2); %Pre-allocation depends on
    %number of years and number months/period -> how to efficiently code
    %this?
    %RX5daycomb = zeros(240, cnt_methods+2);
    %SDIIcomb = zeros(20, cnt_methods+2);
    %PRCPTOTcomb = zeros(20, cnt_methods+2);
    %R10comb = zeros(20, cnt_methods+2);
    %R20comb = zeros(20, cnt_methods+2);
    for i = 1:cnt_methods+2
        RX1daycomb(:, i) = indices{21,1+i, m};
        RX5daycomb(:, i) = indices{22,1+i,m};
        SDIIcomb(:,i) = indices{23,1+i,m};
        PRCPTOTcomb(:,i) = indices{27,1+i,m};
        R10comb(:,i) = indices{24,1+i,m};
        R20comb(:,i) = indices{25,1+i,m};
    end
    
    f14 = figure(14);
    boxplot(RX1daycomb, 'Labels', lgndbox(1:cnt_methods+2))
    ylabel('Precipitation (mm)', 'Interpreter', 'Latex')
    xlabel('Methods', 'Interpreter', 'Latex')
    title(['Box plots of monthly 1-day max P amounts for each BA method for ' NameTimes(m)], 'Interpreter', 'Latex')
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f14)
    saveas(f14, strcat(save_loc_fig, ['BoxRX1day' NameTimes{m} '.png']))
    
    f15 = figure(15);
    boxplot(RX5daycomb, 'Labels', lgndbox(1:cnt_methods+2))
    ylabel('Precipitation (mm)', 'Interpreter', 'Latex')
    xlabel('Methods', 'Interpreter', 'Latex')
    %title(['Box plots of monthly consecutive 5-day max P amounts for each BA method for ' NameTimes(m)], 'Interpreter', 'Latex')
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f15)
    saveas(f15, strcat(save_loc_fig, ['BoxRX5day' NameTimes{m} '.png']))
    
    f16 = figure(16);
    boxplot(SDIIcomb, 'Labels', lgndbox(1:cnt_methods+2))
    ylabel('Intensity index (mm/day)', 'Interpreter', 'Latex')
    xlabel('Methods', 'Interpreter', 'Latex')
    %title(['Box plots of yearly SDII indices for each BA method for ' NameTimes(m)], 'Interpreter', 'Latex')
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f16)
    saveas(f16, strcat(save_loc_fig, ['BoxSDII' NameTimes{m} '.png']))
    
    f17 = figure(17);
    boxplot(PRCPTOTcomb, 'Labels', lgndbox(1:cnt_methods+2), 'LabelOrientation', 'inline')
    ylabel('Annual precipitation (mm)', 'Interpreter', 'Latex')
    xlabel('Methods', 'Interpreter', 'Latex')
    %title(['Box plots of annual precipitation for each BA method for ' NameTimes(m)], 'Interpreter', 'Latex')
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f17)
    saveas(f17, strcat(save_loc_fig, ['BoxPRCPTOT' NameTimes{m} '.png']))
    
    f18 = figure(18);
    boxplot(R10comb, 'Labels', lgndbox(1:cnt_methods+2))
    ylabel('Number of days', 'Interpreter', 'Latex')
    xlabel('Methods', 'Interpreter', 'Latex')
    %title(['Box plots of annual number of heavy precipitation days for ' NameTimes(m)], 'Interpreter', 'Latex')
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f18)
    saveas(f18, strcat(save_loc_fig, ['BoxR10' NameTimes{m} '.png']))
    
    f19 = figure(19);
    boxplot(R20comb, 'Labels', lgndbox(1:cnt_methods+2), 'LabelOrientation', 'inline')
    ylabel('Number of days', 'Interpreter', 'Latex')
    xlabel('Methods', 'Interpreter', 'Latex')
    %title(['Box plots of annual number of very heavy precipitation days for ' NameTimes(m)], 'Interpreter', 'Latex')
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f19)
    saveas(f19, strcat(save_loc_fig, ['BoxR20' NameTimes{m} '.png']))
    
    barQplot(1,:) = barQ{1,m};
    f20 = figure(20);
    x = 1:1:25;
    bar(x, barQplot)
    xlabel('Days', 'Interpreter', 'Latex')
    ylabel('Discharge ($\mathrm{m}^3/\mathrm{s}$)', 'Interpreter', 'Latex')
    legend(lgndbar, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f20)
    saveas(f20, strcat(save_loc_fig, ['BarQ' NameTimes{m} '.png']))
    
    barPplot(1,:) = barP{1,m};
    f21 = figure(21); %Take x from previous bar plot
    bar(x, barPplot)
    xlabel('Days', 'Interpreter', 'Latex')
    ylabel('Precipitation (mm)', 'Interpreter', 'Latex')
    legend(lgndbar, 'Location', 'northeastoutside', 'FontSize', 30)
    set(gca,'FontSize',30)
    set(gca,'FontName', 'UGent Panno Text')
    fullfig(f21)
    saveas(f21, strcat(save_loc_fig, ['BarP' NameTimes{m} '.png']))
    
    close all
end
    

% Bias
for m = 1: nTimes
    for i = 3:cnt_methods + 3
        for j = 2:num_indices+1
            if j ~= (21 && 22 && 23 && 24 && 25 && 27)
                %Difference
                bias{j,i-1,m} = indices{j,i,m}-indices{j,2,m};
            end
        end
    end
    
    for i = 3:cnt_methods + 2
        for j = 2:num_indices+1
            if j ~= (21 && 22 && 23 && 24 && 25 && 27)
                RB_O{j,i-1,m} = 1 - (abs(bias{j,2,m})-abs(bias{j,i,m}))/abs(indices{j,2,m});
                RB_MB{j,i-1,m} = 1 - (abs(bias{j,2,m})-abs(bias{j,i,m}))/abs(bias{j,2,m});
            end
        end
    end
end

%% Opslaan

save(strcat(save_loc, name, save_name, '_indices.mat'), 'indices')
save(strcat(save_loc, name, save_name, '_biases.mat'), 'bias')
save(strcat(save_loc, name, save_name, '_RB_O.mat'), 'RB_O')
save(strcat(save_loc, name, save_name, '_RB_MB.mat'), 'RB_MB')
save(strcat(save_loc, name, save_name, '_PSS_P.mat'), 'PSS_P')
save(strcat(save_loc, name, save_name, '_PSS_T.mat'), 'PSS_T')
save(strcat(save_loc, name, save_name, '_PSS_E.mat'), 'PSS_E')
save(strcat(save_loc, name, save_name, '_PSS_Q.mat'), 'PSS_Q')

end



