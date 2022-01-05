%   BiasChange
%   This script calculates some biases to check how these have changed
%   between the calibration and validation period.
%
%   This file uses xho, xhs, xfs as loaded by b_configurationBiasAdjustment
%   and xobs (the validation period observations) as loaded in
%   c_BAEvaluation.

%   Last update by J. Van de Velde on 15/02/'21

%% Setup

months = [12, 1, 2; 3, 4, 5; 6, 7, 8; 9, 10,11];
nTimes = size(months, 1);

for m=1:nTimes
    
    %% Precipitation
    
    % Average daily precipitation
    
    R_Pav(m) = BiasChangeCalc('average', 6, months(m, :));
    
    % Standard deviation daily precipitation
    R_Pstd(m) = BiasChangeCalc('standarddev', 6, months(m, :));
    
    % Quantiles
    
    R_P5(m) = BiasChangeCalc('quant', 6, months(m, :), 0.05);
    R_P25(m) = BiasChangeCalc('quant', 6, months(m, :), 0.25);
    R_P50(m) = BiasChangeCalc('quant', 6, months(m, :), 0.50);
    R_P75(m) = BiasChangeCalc('quant', 6, months(m, :), 0.75);
    R_P90(m) = BiasChangeCalc('quant', 6, months(m, :), 0.90);
    R_P95(m) = BiasChangeCalc('quant', 6, months(m, :), 0.95);
    R_P99(m) = BiasChangeCalc('quant', 6, months(m, :), 0.99);
    R_P995(m) = BiasChangeCalc('quant', 6, months(m, :), 0.995);
    
    % Precipitation occurrence
    
    R_autocorrP(m) = BiasChangeCalc('autocorr', 6, months(m, :));
    R_ndry(m) = BiasChangeCalc('ndry', 6, months(m, :));
    R_P00(m) = BiasChangeCalc('trans00', 6, months(m, :));
    R_P10(m) = BiasChangeCalc('trans10', 6, months(m, :));
    
    %% Temperature
    
    % Average daily average temperature
    
    R_Tav(m) = BiasChangeCalc('average', 5, months(m, :));
    
    % Standard deviation daily precipitation
    
    R_Tstd(m) = BiasChangeCalc('standarddev', 5, months(m, :));
    
    % Quantiles
    
    R_T5(m) = BiasChangeCalc('quant', 5, months(m, :), 0.05);
    R_T25(m) = BiasChangeCalc('quant', 5, months(m, :), 0.25);
    R_T50(m) = BiasChangeCalc('quant', 5, months(m, :), 0.50);
    R_T75(m) = BiasChangeCalc('quant', 5, months(m, :), 0.75);
    R_T90(m) = BiasChangeCalc('quant', 5, months(m, :), 0.90);
    R_T95(m) = BiasChangeCalc('quant', 5, months(m, :), 0.95);
    R_T99(m) = BiasChangeCalc('quant', 5, months(m, :), 0.99);
    R_T995(m) = BiasChangeCalc('quant', 5, months(m, :), 0.995);
    
    %% Potential evaporation
    
    % Average daily average potential evapotranspiration
    
    R_Eav(m) = BiasChangeCalc('average', 4, months(m, :));
    
    % Standard deviation daily precipitation
    
    R_Estd(m) = BiasChangeCalc('standarddev', 4, months(m, :));
    
    % Quantiles
    
    R_E5(m) = BiasChangeCalc('quant', 4, months(m, :), 0.05);
    R_E25(m) = BiasChangeCalc('quant', 4, months(m, :), 0.25);
    R_E50(m) = BiasChangeCalc('quant', 4, months(m, :), 0.50);
    R_E75(m) = BiasChangeCalc('quant', 4, months(m, :), 0.75);
    R_E90(m) = BiasChangeCalc('quant', 4, months(m, :), 0.90);
    R_E95(m) = BiasChangeCalc('quant', 4, months(m, :), 0.95);
    R_E99(m) = BiasChangeCalc('quant', 4, months(m, :), 0.99);
    R_E995(m) = BiasChangeCalc('quant', 4,months(m, :), 0.995);
    
    %% General
    
    %Correlations
    
    R_corrPE(m) = BiasChangeCalc('correl', 6, months(m, :), [], 4);
    R_corrPT(m) = BiasChangeCalc('correl', 6, months(m, :), [], 5);
    R_corrET(m) = BiasChangeCalc('correl', 5, months(m, :), [], 4);
    
    %Crosscorrelations
    
    R_crosscorrPE0(m) = BiasChangeCalc('crosscor', 6, months(m, :), [], 4, 0);
    R_crosscorrPT0(m) = BiasChangeCalc('crosscor', 6, months(m, :), [], 5, 0);
    R_crosscorrET0(m) = BiasChangeCalc('crosscor', 4, months(m, :), [], 5, 0);
    
    R_crosscorrPE1(m) = BiasChangeCalc('crosscor', 6, months(m, :), [], 4, 1);
    R_crosscorrPT1(m) = BiasChangeCalc('crosscor', 6, months(m, :), [], 5, 1);
    R_crosscorrET1(m) = BiasChangeCalc('crosscor', 4, months(m, :), [], 5, 1);
    
end

