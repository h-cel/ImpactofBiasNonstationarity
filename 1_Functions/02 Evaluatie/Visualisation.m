function bad_performers = Visualisation(figure_name, name, RB_O_name, RB_MB_name, NameTimes)
%   Visualisation Function to visualise a choice of error and bias index combinations
%
%   This function is used in the evaluation done in Van de Velde et al. (2020) 
%
%   This function compares certain methods and indices of bias adjustment methods based on their
%   relative performance on a bias scale (how much bias is removed) and
%   relative added value (is the bias removal worth it?). To do this,
%   the uses is asked for input, so a selection of indices can be made.
%
%   Occurrence and intensity methods are asked each after another, so a
%   very specific selection can be made.
%   For the indices, a number has to be given. These number correspond to
%   the rows of the error and bias files. If the standard files are used,
%   these are the numbers:
%
%   Inputs:
%       figure_name: a name for the figure you want to make
%       RB_O: a file with all values relative to the error (not necessary)
%       RB_MB: a file with all values relative to the bias (not necessary)
%       NameTimes: cell with strings. These are the names of the different
%       seasons/months as also used in c_BAEvaluation
%   Outputs:
%       bad_performers: a file with the names of the method-index combinations
%       that fall out of normal error or bias removal ranges
%
%   Last update by J. Van de Velde on 12/02/'21

%% Path selection

addpath(genpath('D:\Users\jpvdveld\Documents\PhD\Code\Onderzoek'), genpath('E:\Users\jpvdveld\Onderzoek\Data')); %Both Code and Data paths need to be added with their subfolders.

%% Preparation

clf

if nargin == 1
    RBO_tmp = load('MPI-rcp45all_RB_O.mat');
    RBMB_tmp = load('MPI-rcp45all_RB_MB.mat');
    RB_O = RBO_tmp.RB_O;
    RB_MB = RBMB_tmp.RB_MB;
else
    RBO_tmp = load(RB_O_name);
    RBMB_tmp = load(RB_MB_name);
    RB_O = RBO_tmp.RB_O;
    RB_MB = RBMB_tmp.RB_MB;
end

save_loc_fig = 'D:\Users\jpvdveld\Documents\PhD\FigurenResultaten\Evaluation\';
save_loc = 'E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\';
methods_occurrence = {'none', 'ssr', 'tda', 'threshold'};
methods_intensity = {'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc', 'r2d2'};
indices = [];
icons = {'o', '+', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'}; % Possible icons. Let's assume for now that more methods than there are icons will never be compared 

%% Selection

str_more = 'Y';
str_more2 = 'Y';
methods = cell(1,1);
name_figs = cell(1,1);
cnt = 0;

while strcmp(str_more,'Y') ==  1
    x1 = input('Which occurrence method do you want to compare?: ');
    if sum(strcmp(x1, methods_occurrence(:))) == 0
        disp('Please check your answer for typos')
        continue
    end
    x2 = input('Which intensity method do you want to compare?: ');
    if strcmp(x2, 'other')
        x2 = input('Specify another method: ');
    elseif sum(strcmp(x2, methods_intensity(:))) == 0
        disp('Please check your answer for typos')
        continue
    end
    methods{cnt+1} = strcat(name, '_', x1, '_', x2, '_results.mat');
    name_figs{cnt+1} = strcat(name, '-', x1, '-', x2);
    cnt = cnt + 1;
    str_more = input('Do you want to add a method? [Y/N]: ');
end
   
while strcmp(str_more2,'Y') ==  1
    x3 = input('Which index number do you want to use in the comparison?: '); %It's possible to give a vector instead of giving every index as a seperate input!
    indices = [indices x3];
    str_more2 = input('Do you want to add an index? [Y/N]: ');
end

num_indices = length(indices);


%% Interpretation plots

titles = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'}; %Only used for publication
legend_labels = {' ', 'Q$_5$', 'Q$_{25}$', 'Q$_{50}$', 'Q$_{75}$', 'Q$_{90}$', 'Q$_{95}$', 'Q$_{99}$', 'Q$_{99.5}$', 'Q$_{\mathrm{T}20}$', 'P$_5$', 'P$_{25}$', ...
'P$_{50}$', 'P$_{75}$', 'P$_{90}$', 'P$_{95}$', 'P$_{99}$', 'P$_{99.5}$', ...
'P$_{\mathrm{P}00}$', 'P$_{\mathrm{P}10}$', '','','','','', 'N$_{\mathrm{dry}}$', '',...
'E$_5$', 'E$_{25}$', 'E$_{50}$', 'E$_{75}$', 'E$_{90}$', 'E$_{95}$', 'E$_{99}$', 'E$_{99.5}$', ...
'T$_5$', 'T$_{25}$', 'T$_{50}$', 'T$_{75}$', 'T$_{90}$', 'T$_{95}$', 'T$_{99}$', 'T$_{99.5}$',...
'corr$_{\mathrm{P,E}}$', 'corr$_{\mathrm{P,T}}$', 'corr$_{\mathrm{E,T}}$', 'P$_{\mathrm{lag1}}$', 'crosscorr$_{\mathrm{P,E,0}}$', 'crosscorr$_{\mathrm{P,E,1}}$', 'crosscorr$_{\mathrm{P,T,0}}$', 'crosscorr$_{\mathrm{P,T,1}}$','crosscorr$_{\mathrm{E,T,0}}$', 'crosscorr$_{\mathrm{E,T,1}}$'};
colours_times = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980] , [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
figure_handle = figure(1);
t = tiledlayout('flow');
for i = 1:length(methods)
    nexttile
    data_RBMB = zeros(length(indices), 1);
    data_RBO = zeros(length(indices), 1);
    p = zeros(length(indices), 1);
    bad_performers = cell(1,3);
    cnt_perf = 1;
    for m = 1: length(NameTimes)
        for j=1:num_indices
            data_RBMB(j) = RB_MB{indices(j),strcmp(methods(i), RB_MB(1,:,m)),m};
            data_RBO(j) = RB_O{indices(j),strcmp(methods(i), RB_MB(1,:,m)),m};
            % Keeping track of the values that perform too bad to be on the
            % plot
            if data_RBMB(j) > 2 || data_RBO(j) > 2 || data_RBMB(j) <0 || data_RBO(j) < 0 || isinf(data_RBMB(j)) || isinf(data_RBO(j))
                bad_performers{cnt_perf, 1} = [methods{i}, RB_MB{indices(j),1,m}];
                bad_performers{cnt_perf, 2} = data_RBMB(j);
                bad_performers{cnt_perf, 3} = data_RBO(j);
                cnt_perf = cnt_perf + 1;
            end
            p(j) = scatter(data_RBMB(j), data_RBO(j), 80, colours_times{m}, icons{j}, 'LineWidth', 2);
            hold on
        end
        xlim([0, 1.4])
        ylim([0, 1.4])
        plot([1,1], [0,2], 'r--')
        plot([0,2], [1,1], 'r--')
        xlabel('$\mathrm{RB}_{\mathrm{MB}}$', 'interpreter', 'latex', 'FontName', 'Helvetica') %Set off for presentation styling
        ylabel('$\mathrm{RB}_{\mathrm{O}}$', 'interpreter', 'latex', 'FontName', 'Helvetica') %Set off for presentation styling
        grid on
        set(gca,'FontSize',20)
        title(titles{i}) %Uses 'titles' in publication, but this can also be set to use 'name_figs'
        axis square
        
    end
    
end

save(strcat(save_loc ,figure_name, NameTimes{m}, '_bad_performers.mat'), 'bad_performers')

lg = legend(p(:), legend_labels{indices}, 'FontSize', 15, 'interpreter', 'latex');
lg.Layout.Tile = 'East';
%title(t, ['Indices for ', NameTimes{m}])
fullfig(figure_handle)
saveas(figure_handle, strcat(save_loc_fig, figure_name, '.png'))  


close all

end

