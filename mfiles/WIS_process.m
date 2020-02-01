%-- WIS_process.m
%-- 

clear all;
yyyymm = '201610';
src_dir = '\\h.cdip.ucsd.edu\data09\cdip\WIS\data\';
ndbc_id = '41112';

fname = ['WIS-ocean_waves_ST' ndbc_id '_' yyyymm '.nc'];
nc_file = [src_dir ndbc_id '/' fname];

%% Get general information
wis_info = ncinfo(nc_file);

%% Get a list of variables to process
wis_vars = {};
for vid = 1:length(wis_info.Variables)
    wis_vars{vid} = wis_info.Variables(vid).Name;
end

%% Load all variables into a matlab structure
disp(['Loading data from ' fname]);
wis_data = read_nc(nc_file,wis_vars);

%% Get dimensions of direction, freq, time
[ND,NF,NT] = size(wis_data.directionalWaveEnergyDensity);

%% Make some shortcut var names
wisdir = wis_data.waveDirectionBins;
wisfreq = wis_data.waveFrequency;
wis2dspec = wis_data.directionalWaveEnergyDensity;


%% Calculate freq bandwidths
bw = diff(wis_data.waveFrequency);
bw(end+1) = bw(end);
wis_data.waveBandwidth = bw;

%%  Assign WIS directions to correct bins 
%-- (direction coming from, starting with 5 deg bin)
deg_per_bin = 360 / ND;
rot2dspec = zeros(ND,NF,NT);
bin_dirs = zeros(size(wisdir));
%-- directions start with 2.5, so need to subtract this to find idx
for i = 1:ND
    angle = wisdir(i) + 180.0;
    if (angle >= 360.0) 
        angle = angle - 360.0;
    end
    dir_idx(i) = (angle-min(wisdir))/deg_per_bin + 1;
    rot2dspec(dir_idx(i),:,:) = wis2dspec(i,:,:);
    bin_dirs(dir_idx(i)) = single(angle);    
end
rot2dspec = rot2dspec .* (2.0*3.14159/ND);

%% Call function to compute A,B coeffs using rotated spectra and directions
ab_data = WIS_calc_ab(bin_dirs,wisfreq,rot2dspec);


%% Function: read_nc
%-- Function to read remote netcdf file and cell arrays of variables and
%-- create a matlab struct with given vars
function [data] = read_nc(url,var_list)
    data = {};
    %% Iterate over vars of interest and get data
    for vid = 1:length(var_list)
        the_var = var_list{vid};
        %-- Find the index of this variableda
        vindex = find(strcmp(var_list,the_var)) - 1;
        %disp(['--> Loading ' the_var]);
        try
            temp = ncread(url,the_var);
        catch
            disp(['Cant read variable ' the_var])
            return
        end
        eval(['data.' the_var '=temp;']);
    end
end