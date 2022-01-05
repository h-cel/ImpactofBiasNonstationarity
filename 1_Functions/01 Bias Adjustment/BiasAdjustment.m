function BiasAdjustment(xho, xhs, xfs, var, methodocc, methodint, n, modeltype)
%   BiasAdjustment This function implements the bias adjustment according
%   to the input given.
% 
%   This function is launched in the b_configurationBiasAdjustment config
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress) 
%
%   Inputs:   
%       xho: vector that contains the original observations for the current/control period
%       xhs: vector that contains the original GCM output for the current/control period 
%       xfs: vector that contains the original GCM output for the future period
%       var: specifies the input variable: "P" for Precipitation, "E" for Evapotranspiration or "T" for Temperature
%       methodocc: 'none', 'tda', 'threshold', 'ssr'
%       methodint: 'qdm', 'mqdm', 'mbcn', 'mrqnbc', 'dotc'
%       n: ensemble members or repetitions base number
%       modeltype: string with the name of the model run used
%
%   Last updated by Jorn Van de Velde on 25/02/'20

%% Set-up

% Adding save location
save_bias = strcat('E:\Users\jpvdveld\Onderzoek\Data\1_biascorrection\',modeltype,'_', methodocc,'_',methodint); %Bias adjustment path, file named is created by combining the specifications

% Data preparation
ndaysh = length(xho);
ndaysf = length(xfs);
timef = xfs(:,1:3);
timeh = xho(:,1:3);

% Making selections of the variables
% The datasets are structured as follows: Y M D VAR VAR VAR, so to select
% the variable, its location in 'var' + 3 is needed.
% This is based on the structure of the Uccle dataset, and might need
% changes if another dataset is used.

for i = 1:length(var)
    if strcmp(var{i},'T') == 1
        t = i+3;
        tho = xho(:,[1:3,t]); %All rows, columns 1:3 (Y M D) + the variable 
        ths = xhs(:,[1:3,t]);
        tfs = xfs(:,[1:3,t]);
    elseif strcmp(var{i},'E') == 1
        e = i+3;
        eho = xho(:,[1:3,e]);
        ehs = xhs(:,[1:3,e]);
        efs = xfs(:,[1:3,e]);
    else
        p = i+3;
        pho = xho(:,[1:3,p]);
        phs = xhs(:,[1:3,p]);
        pfs = xfs(:,[1:3,p]);
    end
end

%% Preprocessing: Occurrence bias adjustment

% In this part, occurrence bias adjustment takes place

if strcmp(methodocc, 'none') == 0
    switch methodocc
        case 'tda'
            
            %Initialize dry-wet corrected vectors
            pfsdw = nan(ndaysf,n+3); %20 repetitions, extra column for each member
            pfsdw(:,1:3) = timef;
            phsdw = nan(ndaysh,n+3); 
            phsdw(:,1:3) = timeh;
            phodw = nan(ndaysh, n+3);
            phodw(:, 1:3) = timeh;
            
            for i= 1:n
                tic
                %Adjustment
                [pfsdw(:,i+3), phsdw(:,i+3)] = occAdj_TDA(pfs,pho,phs);
                %Repetition of unchanged variables
                phodw(:,i+3) = pho(:,4); % Unchanged, so repetition of the same column
                disp([methodocc, ': Run ', num2str(i), ' of ',num2str(n), ' completed in: ', num2str(toc), 's (Combination:', methodocc, ' + ', methodint, ')'])
            end

        case 'threshold'
            tic
            %Initialize dry-wet corrected vectors
            pfsdw = nan(ndaysf,4); 
            pfsdw(:,1:3) = timef;
            phsdw = nan(ndaysh,4); 
            phsdw(:,1:3) = timeh;
            phodw = nan(ndaysh,4);
            phodw(:, 1:3) = timeh;
            
            %Adjustment
            [pfsdw(:,end), phsdw(:,end)] = occAdj_Threshold(pho,phs,pfs);
            
            %Repetition of unchanged variables
            phodw(:,end) = pho(:,4); % Unchanged, so repetition of the same column
            
            disp([methodocc, ': completed in: ', num2str(toc), 's (Combination: ', methodocc, ' + ', methodint, ')'])

           
        case 'ssr'
            tic
            %Initialize dry-wet corrected vectors
            pfsdw = nan(ndaysf,n+3); %20 repetitions, extra column for each member
            pfsdw(:,1:3) = timef;
            phsdw = nan(ndaysh,n+3); %20 repetitions, extra column for each member
            phsdw(:,1:3) = timeh;
            phodw = nan(ndaysh, n+3); %20 repetitions, extra column for each member
            phodw(:, 1:3) = timeh;
            
            for i= 1:n
            [phodw(:,i+3), phsdw(:,i+3), pfsdw(:,i+3), pth] = occAdj_SSR(pho,phs,pfs);
            disp([methodocc, ': Run ', num2str(i), ' of ',num2str(n), ' completed in: ', num2str(toc), 's (Combination:', methodocc, ' + ', methodint, ')'])
            end
    end
else
    %'No explicit adjustment' of precipitation
    pfsdw = [timef pfs(:, 4)];
    phsdw = [timeh phs(:, 4)];
    phodw = [timeh pho(:, 4)];
end

%% Intensity Bias Adjustment

% In this part, the intensity bias adjustment takes place

%Initialize rescaled vectors
pfa = nan(ndaysf,n+3);
pfa(:,1:3) = timef;

switch methodint
    case 'qdm'
        tic
        n2 = (size(phodw,2)-3); %Repetition number based on occurrence adjustment
        selectionBA = GenerateSelectionBA({timeh, timef});
        efa = NaN(length(timef),1);
        tfa = NaN(length(timef),1);
        for season = 1:4
            for i = 1:n2 % For each of the repetitions made
                [~, pfa(selectionBA{2,season},i+3)] = QDM(phodw(selectionBA{1,season},[1:3,i+3]),phsdw(selectionBA{1,season},[1:3,i+3]),pfsdw(selectionBA{2,season},[1:3,i+3]), 2);
            end
            [~, efa(selectionBA{2,season},:)] = QDM(eho(selectionBA{1,season},:),ehs(selectionBA{1,season},:),efs(selectionBA{2,season},:), 2);
            [~, tfa(selectionBA{2,season},:)] = QDM(tho(selectionBA{1,season},:),ths(selectionBA{1,season},:),tfs(selectionBA{1,season},:), 1); %QDM and occurrence adjustment aren't needed for temperature data
        end
        
        %Saving
        out = cell(1,n2); %Cell size depending on repetitions needed
        for i = 1:n2 %Save everything, because some variables will be overwritten
            out{i} = nan(size(xfs));
            out{i}(:,1:3) = timef;
            out{i}(:,t) = tfa(:,end);
            out{i}(:,p) = pfa(:,i+3);
            out{i}(:,e) = efa(:,end);
        end
        save(strcat(save_bias, '_results'),'out')
        disp([methodint, ': completed in: ', num2str(toc), 's (Combination:', methodocc, ' + ', methodint, ')'])
    case 'mqdm'
        tic
        pfa = nan(ndaysh,n+3);
        pfa(:,1:3) = timeh;
        efa = NaN(length(timeh),1);
        tfa = NaN(length(timeh),1);
        selectionBA = GenerateSelectionBA({timeh, timef});
        n2 = (size(phodw,2)-3); %Repetition number based on occurrence adjustment
        for season= 1:4
            for i = 1:n2 % For each of the repetitions made
                pfa(selectionBA{1,season},i+3) = mQDM(phodw(selectionBA{1,season},[1:3,i+3]),phsdw(selectionBA{1,season},[1:3,i+3]),pfsdw(selectionBA{2,season},[1:3,i+3]), 2);
            end
            efa(selectionBA{1,season}) = mQDM(eho(selectionBA{1,season},:),ehs(selectionBA{1,season},:),efs(selectionBA{2,season},:), 2);
            tfa(selectionBA{1,season}) = mQDM(tho(selectionBA{1,season},:),ths(selectionBA{1,season},:),tfs(selectionBA{2,season},:), 1); %QDM and occurrence adjustment aren't needed for temperature data
        end
        
        %Save
        out = cell(1,n2); %Cell size depending on repetitions needed
        for i = 1:n2
            out{i} = nan(size(xho));
            out{i}(:,1:3) = timeh;
            out{i}(:,t) = tfa(:,end);
            out{i}(:,p) = pfa(:,i+3);
            out{i}(:,e) = efa(:,end);
        end
        save(strcat(save_bias, '_results'),'out')
        disp([methodint, ': completed in: ', num2str(toc), 's (Combination:', methodocc, ' + ', methodint, ')'])
    case 'mbcn' 
        tic
        n2 = (size(phodw,2)-3); %Repetition number based on occurrence adjustment
        tfa = nan(ndaysf,n2+3);
        tfa(:,1:3) = timef;
        efa = nan(ndaysf,n2+3);
        efa(:,1:3) = timef;
        selectionBA = GenerateSelectionBA({timeh, timef});
        for i = 1:n2 % For each of the repetitions made
            for season = 1:4
                % Recombination of the variables
                Xho = [eho(selectionBA{1,season},:), tho(selectionBA{1,season},4), phodw(selectionBA{1,season}, i+3)];
                Xhs = [ehs(selectionBA{1,season},:), ths(selectionBA{1,season},4), phsdw(selectionBA{1,season}, i+3)];
                Xfs = [efs(selectionBA{2,season},:), tfs(selectionBA{2,season},4), pfsdw(selectionBA{2,season}, i+3)];
                [Xfa] = MBCn(Xho,Xhs, Xfs, 0.0001);
                efa(selectionBA{2,season},i+3) = Xfa(:,1);
                tfa(selectionBA{2,season},i+3) = Xfa(:,2);
                pfa(selectionBA{2,season},i+3) = Xfa(:,3);
            end
        end
        %Save
        out = cell(1,n2); %Cell size depending on repetitions needed
        for i = 1:n2
            out{i} = nan(size(xfs));
            out{i}(:,1:3) = timef;
            out{i}(:,t) = tfa(:,i+3);
            out{i}(:,p) = pfa(:,i+3);
            out{i}(:,e) = efa(:,i+3);
        end
        save(strcat(save_bias, '_results'),'out')
        disp([methodint, ': completed in: ', num2str(toc), 's (Combination:', methodocc, '+', methodint, ')'])
    case 'mrqnbc'
        tic
        n2 = (size(phodw,2)-3); %Repetition number based on occurrence adjustment
        tfa = nan(ndaysf,n2+3);
        tfa(:,1:3) = timef;
        efa = nan(ndaysf,n2+3);
        efa(:,1:3) = timef;
        for i = 1:n2 % For each of the repetitions made
            % Recombination of the variables
            Xho = [eho, tho(:,4), phodw(:, i+3)];
            Xhs = [ehs, ths(:,4), phsdw(:, i+3)];
            Xfs = [efs, tfs(:,4), pfsdw(:, i+3)];
            [Xfa] = MRQNBC(Xho,Xhs, Xfs);
            % Updating variables
            efa(:,i+3) = Xfa(:,4);
            tfa(:,i+3) = Xfa(:,5);
            pfa(:,i+3) = Xfa(:,6);          
        end
        %Save
        out = cell(1,n2);
        for i = 1:n2
            out{i} = nan(size(xfs));
            out{i}(:,1:3) = timef;
            out{i}(:,t) = tfa(:,i+3);
            out{i}(:,p) = pfa(:,i+3);
            out{i}(:,e) = efa(:,i+3);
        end
        save(strcat(save_bias, '_results'),'out')
        disp([methodint, ': completed in: ', num2str(toc), 's (Combination:', methodocc, ' + ', methodint, ')'])
    case 'dotc'
        n2 = (size(phodw,2)-3)*n;
        tfa = nan(ndaysf,n2+3);
        tfa(:,1:3) = timef;
        efa = nan(ndaysf,n2+3);
        efa(:,1:3) = timef;
        pfa = nan(ndaysf,n2+3);
        pfa(:,1:3) = timef;
        cnt = 1;
        selectionBA = GenerateSelectionBA({timeh, timef});
        for i = 1:size(phodw,2)-3
            for j = 1:n
                tic
                for season = 1:4
                % Recombination of the variables
                Xho = [eho(selectionBA{1,season},4), tho(selectionBA{1,season},4), phodw(selectionBA{1,season}, i+3)];
                Xhs = [ehs(selectionBA{1,season},4), ths(selectionBA{1,season},4), phsdw(selectionBA{1,season}, i+3)];
                Xfs = [efs(selectionBA{2,season},4), tfs(selectionBA{2,season},4), pfsdw(selectionBA{2,season}, i+3)];
                [Xfa] = dOTC(Xho,Xhs, Xfs);
                if sum(isnan(Xfa)) > 700 % Escape option if dOTC fails
                    cnt = cnt +1;
                    continue
                end
                % Updating variables
                efa(selectionBA{2,season},cnt+3) = Xfa(:,1);
                tfa(selectionBA{2,season},cnt+3) = Xfa(:,2);
                pfa(selectionBA{2,season},cnt+3) = Xfa(:,3);
                end
                cnt = cnt +1;
                disp([methodint, ': Run ', num2str(j), ' of ',num2str(n), ' completed in: ', num2str(toc), 's (Combination:', methodocc, ' + ', methodint, ')'])
            end
        end
        %Save
        out = cell(1,n2);
        for i = 1:n2
            out{i} = nan(size(xfs));
            out{i}(:,1:3) = timef;
            out{i}(:,t) = tfa(:,i+3);
            out{i}(:,p) = pfa(:,i+3);
            out{i}(:,e) = efa(:,i+3);
        end
        save(strcat(save_bias, '_results'),'out')
end

%% Post-processing

% Post-processing steps, like in SSR, end up here
switch methodocc
    case 'ssr'
         [out] = postprocessingSSR(out, pth, p);
        save(strcat(save_bias, '_results'),'out') %Overwrites previous version, though it might be interesting to compare before and after postprocessing
end

end
