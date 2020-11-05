%% Scriptje voor kleine berekening: hoeveel schieten de laatste jaren uit t.o.v. 20e eeuw?

%% Set paths
%Extra paths

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\Onderzoek'), genpath('E:\Users\jpvdveld\Onderzoek\Data')) %Both Code and Data paths need to be added with their subfolders.

%% Loading data

tmp = load('ETP_117y');
ETP = tmp.ETP;

%% Calculations

%Yearly timeseries

E_year = zeros(117,1);
T_year = zeros(117,1);
P_year = zeros(117,1);

for i = 1901:2017
    E_year(i-1900) = nansum(ETP(ETP(:,1) == i, 4));
    T_year(i-1900) = nanmean(ETP(ETP(:,1) == i, 5));
    P_year(i-1900) = nansum(ETP(ETP(:,1) == i, 6));
end

%Means (from 1920-1990)

E_mean = mean(E_year(20:90));
T_mean = mean(T_year(20:90));
P_mean = mean(P_year(20:90));

%Differences

E_diff = (E_year-E_mean)/E_mean*100;
T_diff = T_year-T_mean;
P_diff = (P_year-P_mean)/P_mean*100;

%% Plots

figure(1)
plot(T_diff(1:100), P_diff(1:100), 'o')
xlim([-2.5 2.5])
ylim([-50, 50])
hold on
plot(T_diff(101:117), P_diff(101:117), 'ro')
plot(T_diff(1:100), E_diff(1:100), 'b*')
plot(T_diff(101:117), E_diff(101:117), 'r*')
plot([-2.5,2.5], [0,0], 'r--')
plot([0,0], [-50,50], 'r--')
xlabel('$\mathrm{Temperature}\; \mathrm{anomaly} \left(\textordmasculine \mathrm{C}\right)$', 'Interpreter', 'Latex')
ylabel('$\mathrm{Relative}\; \mathrm{anomaly} \left(\%\right)$', 'Interpreter', 'Latex')
%title({'T and P changes from 1901-2017'; 'compared with mean values from 1920-1990'})
%set(gca,'FontName', 'UGent Panno Text')
set(gca,'FontSize',20)
legend({'Precipitation 1901-2000', 'Precipitation 2001-2017', 'Pot. Evaporation 1901-2000', 'Pot. Evaporation 2001-2017'}, 'Location', 'Southwest')

%figure(2)
%plot(T_diff(1:100), E_diff(1:100), '*')
%xlim([-2.5 2.5])
%ylim([-40, 40])
%hold on
%plot(T_diff(101:117), E_diff(101:117), 'r*')
%plot([-2.5,2.5], [0,0], 'r--')
%plot([0,0], [-50,50], 'r--')
%xlabel('$\mathrm{Temperature} \;\mathrm{difference} \left(\textordmasculine \mathrm{C}\right)$', 'Interpreter', 'Latex')
%ylabel('$\mathrm{Relative}\; \mathrm{potential}\; \mathrm{evaporation} \;\mathrm{difference} \left(\%\right)$', 'Interpreter', 'Latex')
%title({'T and E changes from 1901-2017'; 'compared with mean values from 1920-1990'})
%set(gca,'FontName', 'UGent Panno Text')
%set(gca,'FontSize',20)