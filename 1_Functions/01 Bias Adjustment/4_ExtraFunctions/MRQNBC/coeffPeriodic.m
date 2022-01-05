function [C,D] = coeffPeriodic(time, Z, nPeriods)
%   COEFF This function calculates the periodic coefficient matrices for a timeseries
%
%   This function is launched in the MRQNBC.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function calculates the lag-0 and lag-1 cross-correlation matrices
%   of a periodic time series and uses these to construct the periodic coefficient matrices
%   The calculations are based on the papers by Mehrotra and Sharma, esp.
%   "A multivariate quantile-matching bias correction approach with auto- and 
%   cross-dependence across multiple time scales: Implications for
%   downscaling" (2016) and the book referenced therein, Salas: "Applied
%   modeling of hydrologic time series" (1980)
%
%   Inputs:
%       time: a n x 3 Y M D time series
%       Z: a n x 3 standardized time series of which the coefficients have to be calculated
%       nPeriods: number of periods, 12 if months are used and 4 if seasons are used 
%   Output:
%       C and D: Periodic coefficient matrices
%
%   Last update by J. Van de Velde on 26/11/'19

%Changes needed: remove all mentions of Z, outdated and unclear

%% Set up

C = zeros(3,3,nPeriods);
D = zeros(3,3,nPeriods);
periods = [1:1:nPeriods, 1];

%% Calculation of M0 and M1

for p= 1:nPeriods+1
    m = periods(p);
    Zm = Z(time(:,2) == m, :);
    if m ~=1
        Zmprev = Z(time(:,2) == m-1, :);
    else
        Zmprev = Z(time(:,2) == nPeriods, :);
    end
    [nrows, nvars] = size(Zm);
    
    M0 = zeros(nvars);
    M1 = zeros(nvars);
    cnt0 = zeros(nrows,1);
    
    % M0
    for i = 1:nvars
        for j = 1: nvars
            for N=1:nrows
                cnt0(N) = (Zm(N,i)-mean(Zm(:,i)))*(Zm(N,j)-mean(Zm(:,j)));
            end
            M0(i,j) = sum(cnt0)/(nrows*std(Zm(:,i))*std(Zm(:,j)));
        end
    end
    
    if p >1
    
    %M1
    
    cnt1 = zeros(nrows,1);
    
    for i = 1:nvars
        for j = 1: nvars
            if m ~= 1
                for N=1:nrows
                    cnt1(N) = (Zm(N,i)-mean(Zm(:,i)))*(Zmprev(N,j)-mean(Zmprev(:,j)));
                end
                M1(i,j) = sum(cnt1)/(nrows*std(Zm(:,i))*std(Zmprev(:,j)));
            else
                for N=2:nrows
                    cnt1(N) = (Zm(N,i)-mean(Zm(:,i)))*(Zmprev(N-1,j)-mean(Zmprev(:,j)));
                end
                M1(i,j) = sum(cnt1)/((nrows-1)*std(Zm(:,i))*std(Zmprev(:,j)));
            end
        end
    end
    
    % C
    
    C(:,:,m)= M1/M0_prev;
    
    % D
    
    DDtrans = M0-C(:,:,m)*M1';
    
    [V, S] = eig(DDtrans);
    D(:,:,m) = V*sqrt(S)*V';

% Extra computation in case the eigenvalues method does not work:
    
 %   for j = 1:nvars
 %       for i = 1:nvars
 %           DDdiag = DDtrans(j,j);
 %           if j == 1
 %               sumD = 0;
 %           else
 %                sommation1 = zeros(j-1,1);
 %                for k=1:j-1
 %                       sommation1(k) = D(j,k,m)^2;
 %                end
 %                sumD = sum(sommation1(k));
 %           end
 %           poscheck = DDdiag - sumD; 
 %           if i < j
 %               D(i,j,m) = 0;
 %           elseif poscheck <= 0
 %               D(i,j,m) = 0;
 %           elseif poscheck > 0
 %               if j == 1
 %                   num = DDtrans(i,j);
 %               else
 %                   sommation2 = zeros(j-1,1);
 %                   for k=1:j-1
 %                       sommation2(k) = D(j,k,m)*D(i,k,m);
 %                   end
 %                   num = DDtrans(i,j)-sum(sommation2(k));
 %               end
 %               denum = sqrt(poscheck);
 %               D(i,j,m) = num/denum;
 %           end
 %       end
 %   end
    
    end
    
    M0_prev = M0;

end


end

