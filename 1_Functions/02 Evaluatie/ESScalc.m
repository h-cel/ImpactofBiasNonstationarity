%% Energy distance Skill Score
%   This is a script to calculate the Energy distance skill score of multiple
%   methods, based on the formula in Cannon, 2016
%
%   Configuration script for the evaluation calculations done in Van de Velde et al. (in
%   pub.)

%% Path selection

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\StochasticModelling'), genpath('E:\Users\jpvdveld\Onderzoek\Data')); %Both Code and Data paths need to be added with their subfolders.

%% Method selection

occ_methods = [0 0 0 1]; %Occurrence 
int_methods = [1 1 1 1 1]; %Intensity
%Model used
name = 'MPI-rcp45';
%Time frame used for the comparison
tComp = [1998 1 1; 2017 12 31]; 

%% Set-up

cnt_methods = sum(occ_methods)*sum(int_methods);
names_occ = {'none', 'ssr', 'tdc', 'threshold'};
names_int = {'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc'};
names_methods = cell(1, cnt_methods);
ESS_def = nan(1, cnt_methods);

%% Timeseries Selection

% Observations

tmp = matload('ETP_117y.mat'); %Historical observations: full file
%ho = tmp{1}(:,1:end-1); % D (fraction of day that is dry) is ignored
ho = tmp;
xobs = TruncateObs(ho, tComp);

% Original timeseries

data= load(strcat(name, '_xfs'));
xfs = data.xfs;

% Methods

cnt = 0;

for i=1:length(occ_methods)
    if occ_methods(i) == 1
        name_occ = names_occ{i};
        for j = 1:length(int_methods)
            if int_methods(j) == 1
                name_int = names_int{j};
                cnt = cnt + 1;
                names_methods{cnt} = strcat(name, '_', name_occ, '_', name_int, '_results.mat');
            end
            
        end
    end
end

%% Score calculation

% Original score

or_score = EnDist(xfs(:, 4:6),xobs(:, 4:6));

% Methods

for i = 1:cnt_methods
    
    
    name_method = names_methods{i};
    datacell = matload(name_method);

    time = datacell{1,1}(:,1:3);
    
    nrepeats = size(datacell,2);
    
    ESS = nan(nrepeats, 1);
    
    for coln = 1:nrepeats
    
        E = [time, datacell{1,coln}(:,4)];
        T = [time, datacell{1,coln}(:,5)];
        P = [time, datacell{1,coln}(:,6)];
    
        data = [time, E(:, end), T(:, end), P(:, end)];
        
        ESS(coln) = 1-EnDist(data(:, 4:6),xobs(:, 4:6))/or_score;
        
    end
    
    ESS_def(i) = mean(ESS);
    
end

save(strcat('E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\', name, '_threshOnly', 'ESS.mat'), 'ESS_def') 


