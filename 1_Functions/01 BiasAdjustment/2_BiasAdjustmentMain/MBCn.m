function [Xfa_def] = MBCn(Xho_or,Xhs_or, Xfs_or, tolerance)
%   MBCn This function implements the MBCn bias adjustment method
%   
%   This function is launched in the BiasAdjustment.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function implements the Multivariate Bias Correction in
%   N-dimensions bias adjustment method, as proposed by Cannon (2017). This
%   method combines random rotations with Quantile Delta Mapping to increase
%   persistence attributes of the corrected data.
%    
%   Inputs:
%       Xho_or: the observational dataset, a n x 6-matrix with the variables 
%       in the last column and Y M D in the first 3 columns
%       Xhs_or: the historical simulation dataset, idem
%       Xfs_or: the future simulation dataset, idem
%       tolerance: the tolerance for difference in distances between two
%       following corrections
%   Outputs:
%       Xfa_def: finally adjusted values of the future time series
%
%   Last update by J. Van de Velde on 09/12/'19
   
%% Set-up

% Tolerance set-up
prevdist = 1;
diffdist = 1;
% Selection of variable data
Xho = Xho_or(:,4:6); 
Xhs = Xhs_or(:,4:6);
Xfs = Xfs_or(:,4:6);
% Selection of time
time_h = Xho_or(:,1:3);
time_f = Xfs_or(:,1:3);

%% Rotation loop

while diffdist > tolerance 
    % Rotation
    
    RotMatrix = RandMatrix();
    
    XhoRot = Xho*RotMatrix; % Selection of variable data only, time data needs to be kept for the correction
    XhsRot = Xhs*RotMatrix;
    XfsRot = Xfs*RotMatrix;
    
    % Retrieving the columns
    
    Eho = [time_h, XhoRot(:, 1)];
    Tho = [time_h, XhoRot(:, 2)]; 
    Pho = [time_h, XhoRot(:, 3)]; 
    
    Ehs = [time_h, XhsRot(:, 1)]; 
    Ths = [time_h, XhsRot(:, 2)]; 
    Phs = [time_h, XhsRot(:, 3)]; 
    
    Efs = [time_f, XfsRot(:, 1)]; 
    Tfs = [time_f, XfsRot(:, 2)]; 
    Pfs = [time_f, XfsRot(:, 3)]; 
    
    % QDM
    %Absolute change form is used here!
    
    [Pha, Pfa] = QDM(Pho,Phs,Pfs,1); 
    [Eha, Efa] = QDM(Eho,Ehs,Efs,1);
    [Tha, Tfa] = QDM(Tho,Ths,Tfs,1);
    
    %Rebuilding the resulting matrices
    
    Xha = [Eha, Tha, Pha];
    Xfa = [Efa, Tfa, Pfa];
    
    % Rotation back
    
    Xhs = Xha/RotMatrix; 
    Xfs = Xfa/RotMatrix;
    
    % Energy distance calculation
    
    dist = EnDist(Xhs, Xho);
    diffdist = abs(prevdist-dist);
    prevdist = dist;
    
end

%% Final QDM

%Retrieving the variables

Eho = Xho_or(:, [1:3,4]); 
Tho = Xho_or(:, [1:3,5]); 
Pho = Xho_or(:, [1:3,6]); 
    
Ehs = Xhs_or(:, [1:3,4]); 
Ths = Xhs_or(:, [1:3,5]); 
Phs = Xhs_or(:, [1:3,6]); 
    
Efs = Xfs_or(:, [1:3,4]); 
Tfs = Xfs_or(:, [1:3,5]); 
Pfs = Xfs_or(:, [1:3,6]); 

[~, Pfa] = QDM(Pho,Phs,Pfs,2);
[~, Efa] = QDM(Eho,Ehs,Efs,2);
[~, Tfa] = QDM(Tho,Ths,Tfs,1);


%% Reordening (Schaake shuffle method)

% Removing ties
Xfs = Xfs + rand(size(Xfs,1), 3)*0.001;

% Checking for the rank structure

E_sorted_multi = sort(Xfs(:,1));
[~, E_rnk] = ismember(E_sorted_multi, Xfs(:,1));
T_sorted_multi = sort(Xfs(:,2));
[~, T_rnk] = ismember(T_sorted_multi, Xfs(:,2));
P_sorted_multi = sort(Xfs(:,3));
[~, P_rnk] = ismember(P_sorted_multi, Xfs(:,3));

E_sort_uni = sort(Efa);
T_sort_uni = sort(Tfa);
P_sort_uni = sort(Pfa);

Xfa_def =[E_sort_uni(E_rnk) T_sort_uni(T_rnk) P_sort_uni(P_rnk)];


end

