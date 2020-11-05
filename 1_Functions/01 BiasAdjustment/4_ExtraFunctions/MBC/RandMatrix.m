function RotMatrix = RandMatrix()
%   RANDMATRIX This function generates a random rotation matrix
%
%   This function is launched in the MBCn.m function
%   file and is used in the calculations done in Van de Velde et al. (in
%   progress).
%
%   This function generates a random orthogonal rotation matrix according to the
%   algorithm given by Mezzadri (2007)
%
%   Output:
%       RotMatrix: randomly generated rotation matrix
%
%   Last update by J. Van de Velde on 25/11/'19

%% Creating Z matrix

Z = randn(3,3);

%% QR-decompositie

[Q,R] = qr(Z);

%% Lambda

num = diag(R);
denum = abs(num);

lambda = diag(num./denum);

%% Rotation Matrix

RotMatrix =  Q*lambda;

end

