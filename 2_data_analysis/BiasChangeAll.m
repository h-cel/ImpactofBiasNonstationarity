%   BiasChange
%   This script calculates some biases to check how these have changed
%   between the calibration and validation period.
%
%   This file uses xho, xhs, xfs as loaded by b_configurationBiasAdjustment
%   and xobs (the validation period observations) as loaded in
%   c_BAEvaluation.

%   Last update by J. Van de Velde on 17/01/'20

%% Precipitation

% Average daily precipitation

R_Pav = BiasChangeCalc('average', 6);

% Standard deviation daily precipitation
R_Pstd = BiasChangeCalc('standarddev', 6);

% Quantiles

R_P5 = BiasChangeCalc('quant', 6, 0.05);
R_P25 = BiasChangeCalc('quant', 6, 0.25);
R_P50 = BiasChangeCalc('quant', 6, 0.50);
R_P75 = BiasChangeCalc('quant', 6, 0.75);
R_P90 = BiasChangeCalc('quant', 6, 0.90);
R_P95 = BiasChangeCalc('quant', 6, 0.95);
R_P99 = BiasChangeCalc('quant', 6, 0.99);
R_P995 = BiasChangeCalc('quant', 6, 0.995);

% Precipitation occurrence

R_autocorrP = BiasChangeCalc('autocorr', 6);
R_ndry = BiasChangeCalc('ndry',6);
R_P00 = BiasChangeCalc('trans00', 6);
R_P10 = BiasChangeCalc('trans10', 6);

%% Temperature

% Average daily average temperature

R_Tav = BiasChangeCalc('average', 5);

% Standard deviation daily precipitation

R_Tstd = BiasChangeCalc('standarddev', 5);

% Quantiles

R_T5 = BiasChangeCalc('quant', 5, 0.05);
R_T25 = BiasChangeCalc('quant', 5, 0.25);
R_T50 = BiasChangeCalc('quant', 5, 0.50);
R_T75 = BiasChangeCalc('quant', 5, 0.75);
R_T90 = BiasChangeCalc('quant', 5, 0.90);
R_T95 = BiasChangeCalc('quant', 5, 0.95);
R_T99 = BiasChangeCalc('quant', 5, 0.99);
R_T995 = BiasChangeCalc('quant', 5, 0.995);

%% Potential evapotranspiration

% Average daily average potential evapotranspiration

R_Eav = BiasChangeCalc('average', 4);

% Standard deviation daily precipitation

R_Estd = BiasChangeCalc('standarddev', 4);

% Quantiles

R_E5 = BiasChangeCalc('quant', 4, 0.05);
R_E25 = BiasChangeCalc('quant', 4, 0.25);
R_E50 = BiasChangeCalc('quant', 4, 0.50);
R_E75 = BiasChangeCalc('quant', 4, 0.75);
R_E90 = BiasChangeCalc('quant', 4, 0.90);
R_E95 = BiasChangeCalc('quant', 4, 0.95);
R_E99 = BiasChangeCalc('quant', 4, 0.99);
R_E995 = BiasChangeCalc('quant', 4, 0.995);

%% General

%Correlations

R_corrPE = BiasChangeCalc('correl', 6, [], 4);
R_corrPT = BiasChangeCalc('correl', 6, [], 5);
R_corrET = BiasChangeCalc('correl', 5, [], 4);

%Crosscorrelations

R_crosscorrPE0 = BiasChangeCalc('crosscor', 6, [], 4, 0);
R_crosscorrPT0 = BiasChangeCalc('crosscor', 6, [], 5, 0);
R_crosscorrET0 = BiasChangeCalc('crosscor', 4, [], 5, 0);

R_crosscorrPE1 = BiasChangeCalc('crosscor', 6, [], 4, 1);
R_crosscorrPT1 = BiasChangeCalc('crosscor', 6, [], 5, 1);
R_crosscorrET1 = BiasChangeCalc('crosscor', 4, [], 5, 1);

