function bad_performers = Visualisation(figure_name, RB_O_name, RB_MB_name)
%   Visualisation Function to visualise a choice of error and bias index combinations
%
%   This function is used in the evaluation done in Van de Velde et al. (in
%   progress) 
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
%   Outputs:
%       bad_performers: a file with the names of the method-index combinations
%       that fall out of normal error or bias removal ranges
%
%   Last update by J. Van de Velde on 26/11/'19

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
methods_occurrence = {'none', 'ssr', 'tdc', 'threshold'};
methods_intensity = {'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc'};
indices = [];
icons = {'o', '+', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'}; % Possible icons. Let's assume for now that more methods than there are icons will never be compared 
name = 'MPI-rcp45';

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

titles = {'(a)', '(b)', '(c)', '(d)', '(e)'}; %Only used for publication
figure_handle = figure(1);
data_RBMB = zeros(length(indices), 1);
data_RBO = zeros(length(indices), 1);
p = zeros(length(indices), 1);
bad_performers = cell(1,3);
cnt_perf = 1;
size = ceil(length(methods)/2);
for i = 1:cnt+1
    subplot(2,size,i) %size+1 for presentation use
hold on
    if i<3
        for j=1:num_indices
            data_RBMB(j) = RB_MB{indices(j),strcmp(methods(i), RB_MB(1,:))};
            data_RBO(j) = RB_O{indices(j),strcmp(methods(i), RB_MB(1,:))};
            % Keeping track of the values that perform too bad to be on the
            % plot
            if data_RBMB(j) > 2 || data_RBO(j) > 2 || data_RBMB(j) <0 || data_RBO(j) < 0 || isinf(data_RBMB(j)) || isinf(data_RBO(j))
                bad_performers{cnt_perf, 1} = [methods{i}, RB_MB{j+1,1}];
                bad_performers{cnt_perf, 2} = data_RBMB(j);
                bad_performers{cnt_perf, 3} = data_RBO(j);
                cnt_perf = cnt_perf + 1;
            end
            p(j) = scatter(data_RBMB(j), data_RBO(j), 80, icons{j}, 'LineWidth', 2);
        end
        xlim([0, 1.4])
        ylim([0, 1.4])
        plot([1,1], [0,2], 'r--')
        plot([0,2], [1,1], 'r--')
        xlabel('$\mathrm{RB}_{\mathrm{MB}}$', 'interpreter', 'latex', 'FontName', 'Helvetica') %Set off for presentation styling
        ylabel('$\mathrm{RB}_{\mathrm{O}}$', 'interpreter', 'latex', 'FontName', 'Helvetica') %Set off for presentation styling
        grid on
        %legend(p(:), RB_O(indices,1), 'FontSize', 15, 'Location',
        %'Southeast') %Set of for the introduction of a general legend
        set(gca,'FontSize',20)
        title(titles{i}) %Uses 'titles' in publication, but this can also be set to use 'name_figs'
        axis square
    elseif i == 3 %Option to enable a separate legend at a good location
        data_RBMB = NaN(length(indices), 1);
        data_RBO = NaN(length(indices), 1);
        for j=1:num_indices
            p(j) = scatter(data_RBMB(j), data_RBO(j), 80, icons{j}, 'LineWidth', 2);
        end
        xlim([0, 1.4])
        ylim([0, 1.4])
        axis off
        legend(p(:), RB_O(indices,1), 'FontSize', 25, 'Location','Northwest');
    else %Continuation of the main plots, but with i-1
         for j=1:num_indices
            data_RBMB(j) = RB_MB{indices(j),strcmp(methods(i-1), RB_MB(1,:))};
            data_RBO(j) = RB_O{indices(j),strcmp(methods(i-1), RB_MB(1,:))};
            % Keeping track of the values that perform too bad to be on the
            % plot
            if data_RBMB(j) > 2 || data_RBO(j) > 2 || data_RBMB(j) <0 || data_RBO(j) < 0 || isinf(data_RBMB(j)) || isinf(data_RBO(j))
                bad_performers{cnt_perf, 1} = [methods{i-1}, RB_MB{j+1,1}];
                bad_performers{cnt_perf, 2} = data_RBMB(j);
                bad_performers{cnt_perf, 3} = data_RBO(j);
                cnt_perf = cnt_perf + 1;
            end
            p(j) = scatter(data_RBMB(j), data_RBO(j), 80, icons{j}, 'LineWidth', 2);
        end
        xlim([0, 1.4])
        ylim([0, 1.4])
        plot([1,1], [0,2], 'r--')
        plot([0,2], [1,1], 'r--')
        xlabel('$\mathrm{RB}_{\mathrm{MB}}$', 'interpreter', 'latex', 'FontName', 'Helvetica') %Set off for presentation styling
        ylabel('$\mathrm{RB}_{\mathrm{O}}$', 'interpreter', 'latex', 'FontName', 'Helvetica') %Set off for presentation styling
        grid on
        %legend(p(:), RB_O(indices,1), 'FontSize', 15, 'Location',
        %'Southeast') %Set of for the introduction of a general legend
        set(gca,'FontSize',20)
        title(titles{i-1}) %Will be set off when making up the figure for publication
        axis square
    end
end

fullfig(figure_handle)
%saveas(figure_handle, strcat(save_loc_fig, figure_name, '.png'))

save(strcat(save_loc ,figure_name, '_bad_performers.mat'), 'bad_performers')

end

