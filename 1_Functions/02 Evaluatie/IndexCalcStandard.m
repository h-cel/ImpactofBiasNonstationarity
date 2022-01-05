function [indices, F_Q_ecdf, x_Q_ecdf, F_Q_95_ecdf, x_Q_95_ecdf, F_P_ecdf, x_P_ecdf, F_P_95_ecdf, x_P_95_ecdf, F_E_ecdf, x_E_ecdf, F_T_ecdf, x_T_ecdf, Spelllength, P, Qsim, E, T, barQ, barP] = IndexCalcStandard(x, wet, num_indices, name, nrepeats, months, E, T, P)
%INDEXCALCSTANDARD Calculation of BA evaluation indices
%
%   This function is launched in the BA_Evaluation.m function
%   file and is used in the evaluation done in Van de Velde et al. (2020) 
%
%   This function calculates indices of selected bias adjustment methods by comparing them
%   to the original simulations and to the observations. The methods compared
%   here are always a combination occurence and intensity correcting methods.
%
%   Inputs:
%       x: the data of which the indices have to be calculated [nx6]
%       matrix, with n timestep and Y:M:D:E:T:P as columns
%       wet: Boolean indicating whether the P quantiles have to be
%       calculated on wet days only (1) or on all days (0)
%       num_indices: number of indices in the final result
%       name: base name of the time period studied (string), used for saving the
%       discharge
%       nrepeats: number of repeats of the bias adjustment procedure.
%       Depends on the input data
%       months: [1xn] vector indicating which n months have to be selected
%       E, T, P: (not necessary): preprocessed E, T and P matrices. If
%       these are not given as an input, they will be calculated based on
%       x.
%   Outputs:
%       indices: cell with indices calculated for observations, the
%       original RCM simulations and every method for all months in the
%       months-vector
%       *_ecdf: empirical CDF data, used for example plotting
%       Spellength: spell length PDF for the months in the months-vector
%       E, T, P: E, T and P matrices as calculated in this function. Will
%       be the same as the input if these are given as input.
%       Qsim: simulated discharge for the months given in the months
%       vector.
%       bias: cell with the biases (= difference) between 1) the observations and 2) the original RCM simulations and the methods used
%       barQ, barP: examples of the P and Q timeseries. Based on a
%       subjective selection of days
%
%   Last update by J. Van de Velde on 11/02/'21

%% Pre-allocation

indices = cell(num_indices,1);

qnt_standardQ = nan(nrepeats, 8);
Qsim = nan(length(x),4, nrepeats);
val20yr = nan(nrepeats, 1);
qnt_standardP = nan(nrepeats, 8);
trans = nan(nrepeats, 4);
%rx1day = nan(nrepeats, 240); %pre-allocation depends on the number of
%years and the number of months/period
%rx5day = nan(nrepeats, 240);
%sdii = nan(nrepeats, 20);
%r10 = nan(nrepeats, 20);
%r20 = nan(nrepeats, 20);
%prcptot = nan(nrepeats, 20);
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

%% Calculations

for coln = 1:nrepeats

% Discharge

    Qsim(:,:,coln) = Discharge(x(:,:, coln), name);
    % I keep all the data for debugging/checking purposes)

    % Statistics (= general overview)
    
    qnt_standardQ(coln, :) = quantile(Qsim(:, end,coln), [0.05 0.25, 0.50, 0.75, 0.90, 0.95 0.99, 0.995]);

    % Indices
    
    val20yr(coln) = val20yRP(Qsim(:,:,coln));
    
    % Bar plot
    if coln == nrepeats
        barQ(:) = Qsim(60:84,4, coln);
    end
    
% Precipitation
if nargin <7
    P = x(:, [1:3,6]);
    P(isnan(P)) = 0;
end

    % Statistics (= general overview)
    
    if wet == 1
        qnt_standardP(coln, :) = quantile(P(P(:,end, coln)>0.1,end, coln), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    else
        qnt_standardP(coln, :) = quantile(P(:,end, coln), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    end
    
    % Indices
    
    trans(coln, :) = TransProb(P(:, end, coln), 0.1);
    rx1day(coln, :) = RX1day(P(:,:, coln), months);
    rx5day(coln, :) = RX5day(P(:,:, coln), months);
    sdii(coln, :) =  SDII(P(:,:, coln));
    r10(coln, :) = R10(P(:,:, coln));
    r20(coln, :) = R20(P(:,:, coln));
    prcptot(coln, :) = PRCPTOT(P(:,:, coln));
    ndry(coln) = sum(P(:,end, coln)<0.1);
    
    % Bar plot
    if coln == nrepeats
        barP(:,1) = P(60:84,end, coln);
    end

% Evaporation
if nargin <7
    E = x(:, 1:4);
    E(isnan(E)) = 0;
end

    % Statistics (= general overview)
    
    qnt_standardE(coln, :) = quantile(E(:,end, coln), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);
    
% Temperature
if nargin <7
    T = x(:, [1:3, 5]);
    T(isnan(T)) = 0;
end

        % Statistics (= general overview)
    
        qnt_standardT(coln, :) = quantile(T(:, end, coln), [0.05 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.995]);


% General

corrPE(coln) = corr(P(:, end, coln),E(:,end, coln), 'type', 'Spearman', 'rows', 'complete');
corrPT(coln) = corr(P(:, end, coln),T(:, end, coln), 'type', 'Spearman', 'rows', 'complete');
corrET(coln) = corr(E(:,end, coln),T(:, end, coln), 'type', 'Spearman', 'rows', 'complete');
autocorrP(coln, :) = xcorr(P(:, end, coln), 1, 'coeff');
crosscorrPE(coln, :) = xcorr(P(~isnan(P(:,end, coln))&~isnan(E(:,end, coln)), end, coln),E(~isnan(P(:,end, coln))&~isnan(E(:,end, coln)),end, coln), 1, 'coeff');
crosscorrPT(coln, :) = xcorr(P(~isnan(P(:,end, coln))&~isnan(T(:,end, coln)), end),T(~isnan(P(:,end, coln))&~isnan(T(:,end, coln)), end, coln), 1, 'coeff');
crosscorrET(coln, :) = xcorr(E(~isnan(E(:,end,coln))&~isnan(T(:,end, coln)), end, coln),T(~isnan(E(:,end,coln))&~isnan(T(:,end,coln)), end, coln), 1, 'coeff');

end

%% Indices

indices{1} = mean(qnt_standardQ(:,1));
indices{2} = mean(qnt_standardQ(:,2));
indices{3} = mean(qnt_standardQ(:,3));
indices{4} = mean(qnt_standardQ(:,4));
indices{5} = mean(qnt_standardQ(:,5));
indices{6} = mean(qnt_standardQ(:,6));
indices{7} = mean(qnt_standardQ(:,7));
indices{8} = mean(qnt_standardQ(:,8));
indices{9} = mean(val20yr);
indices{10} = mean(qnt_standardP(:,1));
indices{11} = mean(qnt_standardP(:,2));
indices{12} = mean(qnt_standardP(:,3));
indices{13} = mean(qnt_standardP(:,4));
indices{14} = mean(qnt_standardP(:,5));
indices{15} = mean(qnt_standardP(:,6));
indices{16} = mean(qnt_standardP(:,7));
indices{17} = mean(qnt_standardP(:,8));
indices{18} = mean(trans(:, 1));
indices{19} = mean(trans(:, 3));
indices{20} = mean(rx1day, 1);
indices{21} = mean(rx5day, 1);
indices{22} = mean(sdii, 1);
indices{23} = mean(r10, 1);
indices{24} = mean(r20, 1);
indices{25} = mean(ndry);
indices{26} = mean(prcptot, 1);
indices{27} = mean(qnt_standardE(:,1));
indices{28} = mean(qnt_standardE(:,2));
indices{29} = mean(qnt_standardE(:,3));
indices{30} = mean(qnt_standardE(:,4));
indices{31} = mean(qnt_standardE(:,5));
indices{32} = mean(qnt_standardE(:,6));
indices{33} = mean(qnt_standardE(:,7));
indices{34} = mean(qnt_standardE(:,8));
indices{35} = mean(qnt_standardT(:,1));
indices{36} = mean(qnt_standardT(:,2));
indices{37} = mean(qnt_standardT(:,3));
indices{38} = mean(qnt_standardT(:,4));
indices{39} = mean(qnt_standardT(:,5));
indices{40} = mean(qnt_standardT(:,6));
indices{41} = mean(qnt_standardT(:,7));
indices{42} = mean(qnt_standardT(:,8));
indices{43} = mean(corrPE);
indices{44} = mean(corrPT);
indices{45} = mean(corrET);
indices{46} = mean(autocorrP(:, 3));
indices{47} = mean(crosscorrPE(:, 2));
indices{48} = mean(crosscorrPE(:, 3));
indices{49} = mean(crosscorrPT(:, 2));
indices{50} = mean(crosscorrPT(:, 3));
indices{51} = mean(crosscorrET(:, 2));
indices{52} = mean(crosscorrET(:, 3));

%% CDFs

coln = 1;

% Making plots of the last column as an example
% This should be improved: not representative in case of R2D2 and dOTC

[F_Q_ecdf, x_Q_ecdf] = ecdf(Qsim(:, end));
[F_Q_95_ecdf, x_Q_95_ecdf] = ecdf(Qsim(Qsim(:, end) > indices{6}, end));

[F_P_ecdf, x_P_ecdf] = ecdf(P(:, end, coln));
[F_P_95_ecdf, x_P_95_ecdf] = ecdf(P(P(:, end, coln) > indices{15}, end, coln));

[F_T_ecdf, x_T_ecdf] = ecdf(T(:, end, coln));

[F_E_ecdf, x_E_ecdf] = ecdf(E(:, end, coln));

Spelllength = SpellDist(P(:, end, coln));


end

