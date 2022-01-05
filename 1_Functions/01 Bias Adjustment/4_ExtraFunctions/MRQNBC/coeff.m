function [C,D] = coeff(Z)
%	COEFF This function calculates the coefficient matrices for a timeseries
%
%   This function is launched in the MRQNBC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function calculates the lag-0 and lag-1 cross-correlation matrices
%   of a time series and uses these to construct the coefficient matrices.
%   The calculations are based on the papers by Mehrotra and Sharma, esp.
%   "A multivariate quantile-matching bias correction approach with auto- and 
%   cross-dependence across multiple time scales: Implications for
%   downscaling" (2016) and the book referenced therein, Salas: "Applied
%   modeling of hydrologic time series" (1980)
%   
%   Input:
%       Z: a n x 3 standardized timeseries of which the coefficients have to be calculated
%   Output:
%       C and D: Coefficient matrices
%
%   Last update by J. Van de Velde on 26/11/'19

%% M1 and M0 matrices

[nrows, nvars] = size(Z);
M0 = zeros(nvars);
M1 = zeros(nvars);

% M0
cnt0 = zeros(nrows,1);

for i = 1:nvars
    for j = 1: nvars
        for N = 1: nrows
            cnt0(N) = Z(N,i)*Z(N,j);
        end
        M0(i,j) = sum(cnt0)/nrows;
    end
end

%M1

cnt1 = zeros(nrows,1);

for i = 1:nvars
    for j = 1: nvars
        for N = 1: nrows-1
            cnt1(N) = Z(N,i)*Z(N+1,j);
        end
        M1(i,j) = sum(cnt1)/(nrows-1);
    end
end

%% C

C= M1/M0;

%% D

DDtrans = M0-C*M1';
[V, S] = eig(DDtrans);

D = V*sqrt(S)*V';


end

