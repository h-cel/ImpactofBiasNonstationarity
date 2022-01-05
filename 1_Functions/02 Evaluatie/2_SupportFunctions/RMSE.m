function r = RMSE(x1,x2)
%   RMSE calculates RMSE
%
%   Last update by J. Van de Velde on 26/11/'19: documentation

%% Calculation

MSE = mean((x1-x2).^2);
r = sqrt(MSE);

end

