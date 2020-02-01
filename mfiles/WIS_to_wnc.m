%% WIS_to_wnc.m
%-------------------------------------------------------------------------
%- Aggregate WIS matlab files into matlab
%- Calculate additional parameters/variables
%- and write output to CDIP netcdf format suitable for CDIP THREDDS server.
%- The WIS model data is originally stored in NetCDF format on following server:
%- https://chlthredds.erdc.dren.mil/thredds/catalog/wis/catalog.html
%- Files are stored as individual months in specific regions.
%-------------------------------------------------------------------------

%function [ data ] = WIS_to_wnc(cdip_id,region,yyyymm)
clear all;
data = {};

%% Initialize variables
cdip_id = '132';
start_date = '201601';
end_date = '201612';
data_dir = '../data/';

%% Get CDIP template file
try
    src_url = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/test/WW3/41112_WW3_realtime.nc';
    S0 = ncinfo(src_url);
catch
    disp(['-->Problem loading: ' src_url]);
    return
end

%% Aggregate data from monthly files
mat_fname = [data_dir cdip_id '_' start_date '-' end_date '_WIS.mat'];
disp(['Aggregating ' mat_fname]);
if(~ isfile(mat_fname))
    wis_data = WIS_aggregate(cdip_id,start_date,end_date);
else
    load(mat_fname);
end
ndbc_id = num2str(wis_data.station_name);
mat_time = time_correct(wis_data.time);

%% Compute BW,EnergyDensity,A1,B1,A2,B2
disp(['Compute additional parameters ... ']);
[ND,NF,NT] = size(wis_data.directionalWaveEnergyDensity);
energy2d = wis_data.directionalWaveEnergyDensity;     
ab_data = calc_a_b_wis(wis_data.waveFrequency, ...
    wis_data.waveDirectionBins,energy2d);
wis_data.waveBandwidth = ab_data.bw;
wis_data.waveA1Value = ab_data.a1;
wis_data.waveA2Value = ab_data.a2;
wis_data.waveB1Value = ab_data.b1;
wis_data.waveB2Value = ab_data.b2;
wis_data.waveEnergyDensity = ab_data.energy;
wis_data.station_name = num2str(wis_data.station_name);

%% Compute time, freq, dir bounds and wave flags
%-- Add time_bounds 
time_bounds = [];
delta_time = wis_data.time(2)-wis_data.time(1);
for i = 1:length(wis_data.time)-1
    time_bounds(:,i) = [wis_data.time(i); wis_data.time(i+1);];
end
time_bounds(:,i+1) = [wis_data.time(end); wis_data.time(end)+delta_time;];
wis_data.time_bounds = time_bounds;

%-- Add frequency bounds 
%wis_data.f = wis_data.f';
waveFrequencyBounds = [];
for i = 1:NF-1
    waveFrequencyBounds(i,:) = [wis_data.waveFrequency(i); ...
        wis_data.waveFrequency(i+1);];
end
waveFrequencyBounds(i+1,:) = [wis_data.waveFrequency(end); ...
    wis_data.waveFrequency(end)+wis_data.waveBandwidth(end);];
wis_data.waveFrequencyBounds = waveFrequencyBounds;
%-- Add in waveFlag parameters
wis_data.waveFlagPrimary = ones(1,length(wis_data.time));
wis_data.waveFlagSecondary = ones(1,length(wis_data.time));
%-- Add in waveDirectionBounds
 waveDirectionBounds = [];
for i = 1:ND-1
    waveDirectionBounds(i,:) = [wis_data.waveDirectionBins(i); ...
        wis_data.waveDirectionBins(i+1);];
end
delta_dir = wis_data.waveDirectionBins(2)-wis_data.waveDirectionBins(1);
waveDirectionBounds(i+1,:) = [wis_data.waveDirectionBins(end); ...
    wis_data.waveDirectionBins(end)+delta_dir;];
wis_data.waveDirectionBounds = waveDirectionBounds;   

%% Get variable key name matches between cdip and matlab-wis files
M = regexp(fileread('../cdip-wis-wnckey.csv'),'[\n\r]+','split');
cdip_key = [];wis_key = [];dim_key=[];
for i = 1:length(M)-1
    var_names = strsplit(char(M(i)),',');
    cdip_key{i} = char(var_names(1));
    wis_key{i} = char(var_names(2));
    dim_key{i} = char(var_names(3));
end 

%% Create a variable list from both template and key file
var_index_list = {};
idx = 1;
for vid = 1:length(cdip_key)
    cdip_var = cdip_key{vid};
    wis_var = wis_key{vid};
    if (~ strcmp(wis_var,'nan'))
        wis_var_list{idx} = wis_var;
        cdip_var_list{idx} = cdip_var;
        idx = idx + 1;
    end
end



%% Create a netcdf file in same fashion as cdip schema template
out_dir = '../data/';
fname = [ndbc_id '_WIS.nc'];
nc_file = [out_dir fname];
if (isfile(nc_file))
    try
        delete(nc_file);
    catch
        disp(['File exists: ' nc_file])
        return
    end
end
try
    nccreate(nc_file,'waveTime','Dimensions', ...
        {'waveTime' NT}, ...
        'DataType','int32', 'Format','classic');
    nccreate(nc_file,'waveTimeBounds','Dimensions', ...
        {'metaBoundsCount' 2 'waveTime' NT}, ...
        'DataType','int32');
    nccreate(nc_file,'waveFlagPrimary','Dimensions', ...
        {'waveTime' NT}, ...
        'DataType','int8');
    nccreate(nc_file,'waveFlagSecondary','Dimensions', ...
        {'waveTime' NT}, ...
        'DataType','int8');    
    nccreate(nc_file,'waveHs','Dimensions', ...
        {'waveTime' NT},'DataType','single');
    nccreate(nc_file,'waveTp','Dimensions', ...
        {'waveTime' NT},'DataType','single');    
    nccreate(nc_file,'waveTa','Dimensions', ...
        {'waveTime' NT},'DataType','single'); 
    nccreate(nc_file,'waveDp','Dimensions', ...
        {'waveTime' NT},'DataType','single'); 
    nccreate(nc_file,'waveModelInputSource','Dimensions', ...
        {'waveTime' NT}); 
    nccreate(nc_file,'waveFrequency','Dimensions', ...
        {'waveFrequency' NF},'DataType','single'); 
    nccreate(nc_file,'waveFrequencyBounds','Dimensions', ...
        {'waveFrequency' NF 'metaBoundsCount' 2}, ...
        'DataType','single');
    nccreate(nc_file,'waveBandwidth','Dimensions', ...
        {'waveFrequency' NF},'DataType','single');
    nccreate(nc_file,'waveEnergyDensity','Dimensions', ...
        {'waveTime' NT 'waveFrequency' ...
        NF},'DataType','single');
    nccreate(nc_file,'waveA1Value','Dimensions', ...
        {'waveTime' NT 'waveFrequency' ...
        NF},'DataType','single');
    nccreate(nc_file,'waveB1Value','Dimensions', ...
        {'waveTime' NT 'waveFrequency' ...
        NF},'DataType','single');
    nccreate(nc_file,'waveA2Value','Dimensions', ...
        {'waveTime' NT 'waveFrequency' ...
        NF},'DataType','single');
    nccreate(nc_file,'waveB2Value','Dimensions', ...
        {'waveTime' NT 'waveFrequency' ...
        NF},'DataType','single');
    nccreate(nc_file,'waveDirection','Dimensions', ...
        {'waveDirection' ND},'DataType','single');
    nccreate(nc_file,'waveDirectionBounds','Dimensions', ...
        {'waveDirection' ND 'metaBoundsCount' 2},...
        'DataType','single');
%     nccreate(nc_file,'waveDirectionSpectrum','Dimensions', ...
%         {'waveTime' NT ...
%         'waveFrequency' NF ...
%         'waveDirection' ND});
    nccreate(nc_file,'metaSiteLabel','Dimensions', ...
        {'maxStrlen64' 64},'DataType','char');
    nccreate(nc_file,'metaLatitude','DataType','single');  
    nccreate(nc_file,'metaLongitude','DataType','single');  
    nccreate(nc_file,'metaWaterDepth','DataType','single');  
    nccreate(nc_file,'metaShoreNormal','DataType','single');  
catch
    disp('Could not create netcdf variables')
    return
end

%% Write Global Attributes
atts = wis_data.ncinfo.Attributes;
for i = 1:length(atts)
    ncwriteatt(nc_file,'/',atts(i).Name,atts(i).Value);
end

%% Write variable data
for i = 1:length(cdip_var_list)
    var_name = cdip_var_list{i};
    wis_var_name = wis_var_list{i};
    %-- Figure which index this is in netcdf template
    for j=1:length(S0.Variables)
        if(strcmp(S0.Variables(j).Name,var_name))
            idx = j;
        end
    end
    var_type = S0.Variables(idx).Datatype;
    %-- Get original attributes
    W = wis_data.ncinfo.Variables;
    for j=1:length(W)   
        if(strcmp(W(j).Name,var_name))
            idx = j;
        end
    end   
    var_atts = W(idx).Attributes;   %-- Struct with name->value
    
    %-- Get data
    eval(['data=wis_data.' wis_var_name ';']);
    %-- Cast variable correctly
    eval([var_type '(data);']);
    %-- Get attributes
    %-- Write data to netcdf variable
    eval(['ncwrite(nc_file,var_name,wis_data.' wis_var_name ')'])
    %- Write variable attributres;
    for j = 1:length(var_atts)
        att = var_atts(j);
        if (~ strcmp(att.Name,'_FillValue'))
            ncwriteatt(nc_file,var_name,att.Name,att.Value);
        end       
    end
    
end

S = ncinfo(nc_file);




%% Function: read_nc
%-- Function to read remote netcdf file and cell arrays of variables and
%-- create a matlab struct with given vars
function [data] = read_nc(url,var_list)
    %% Iterate over vars of interest and get data
    for vid = 1:length(var_list)
        the_var = var_list{vid};
        %-- Find the index of this variableda
        vindex = find(strcmp(var_list,the_var)) - 1;
        disp(['--> Loading ' the_var]);
        temp = ncread(url,the_var);
        eval(['data.' the_var '=temp;']);
    end
end


%end

%% Function to correct time epoch from Allie H cdipxww3 
function [ mat_time ] = time_correct(ts)
    toff = datenum(1970,1,1,0,0,0);
    mat_time = ts./(24*60*60) + toff;
end