function TransformToNetCDF(data, save_loc, filename)
%   TRANSFORMTONETCDF This function transforms a dataset to a .nc-file  
%   
%   This function transforms a .mat time series (result from the
%   PrepareBiasData function) to an .nc-file, so it can be used in in other
%   tools, such as Python packages
%
%   Last adjusted for ISIMIP3, which requires the variables to be in
%   separate .nc-files
%
%   Last update by J. Van de Velde on 05/06/'20

%% Set-up

nametas = strcat(save_loc, filename, '_tas.nc');
namepr = strcat(save_loc, filename, '_pr.nc');

%% Dimensions

mySchema.Name   = '/';
mySchema.Format = 'classic';
mySchema.Dimensions(1).Name   = 'time';
mySchema.Dimensions(1).Length = size(data, 1);
%mySchema.Dimensions(2).Name   = 'rlat';
%mySchema.Dimensions(2).Length = 1;
%mySchema.Dimensions(3).Name   = 'rlon';
%mySchema.Dimensions(3).Length = 1;
ncwriteschema(nametas, mySchema);
ncwriteschema(namepr, mySchema);

%% Netcdf creation & time

nccreate(nametas ,'time', 'Dimensions', {'time', size(data, 1)});
ncwriteatt(nametas,'time', 'units', 'days since 1950-01-01 12:00:00');
ncwriteatt(nametas,'time', 'calendar', 'proleptic_gregorian');
%ncwriteatt(nametas,'time', 'axis', 'T');
ncwriteatt(nametas,'time', 'standard_name', 'time');
ncwriteatt(nametas,'time', 'long_name', 'time');
ncwrite(nametas,'time', (datenum(data(:,1:3))-datenum(1950, 1, 1)));

nccreate(namepr ,'time', 'Dimensions', {'time', size(data, 1)});
ncwriteatt(namepr,'time', 'units', 'days since 1950-01-01 12:00:00');
ncwriteatt(namepr,'time', 'calendar', 'proleptic_gregorian');
%ncwriteatt(namepr,'time', 'axis', 'T');
ncwriteatt(namepr,'time', 'standard_name', 'time');
ncwriteatt(namepr,'time', 'long_name', 'time');
ncwrite(namepr,'time', (datenum(data(:,1:3))-datenum(1950, 1, 1)));

%% Latitudes 
% Data from original .nc-files

%nccreate(nametas ,'rlat', 'Dimensions', {'rlat', 1});
%ncwrite(nametas,'rlat', 0.9350);
%ncwriteatt(nametas,'rlat', 'axis', 'Y');
%nccreate(namepr ,'rlat', 'Dimensions', {'rlat', 1});
%ncwrite(namepr,'rlat', 0.9350);
%ncwriteatt(namepr,'rlat', 'axis', 'Y');


%nccreate(nametas ,'rlon', 'Dimensions', {'rlon', 1});
%ncwrite(nametas,'rlon', -8.5750);
%ncwriteatt(nametas,'rlon', 'axis', 'X');
%nccreate(namepr ,'rlon', 'Dimensions', {'rlon', 1});
%ncwrite(namepr,'rlon', -8.5750);
%ncwriteatt(namepr,'rlon', 'axis', 'X');


%% P & T

nccreate(nametas ,'tas', 'Dimensions', {'time', size(data, 1)});
ncwrite(nametas, 'tas', data(:, 5))

nccreate(namepr ,'pr', 'Dimensions', {'time', size(data, 1)});
ncwrite(namepr, 'pr', data(:, 6))

%% Global variables

ncwriteatt(nametas,'/','Conventions','CF-1.4');
ncwriteatt(namepr,'/','Conventions','CF-1.4');

end

