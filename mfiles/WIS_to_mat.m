%% WIS_to_mat.m
%-------------------------------------------------------------------------
%- Load monthly individual WIS netcdf file into matlab and write output
%- to CDIP netcdf format suitable for use on CDIP THREDDS server.
%- The WIS model data is stored in NetCDF format on following server:
%- https://chlthredds.erdc.dren.mil/thredds/catalog/wis/catalog.html
%- Files are stored as individual months in specific regions.
%-------------------------------------------------------------------------

%function [ data ] = WIS_to_mat(cdip_id,region,yyyymm)
    region = 'Atlantic';         %- Need to figure out which region buoy is in.
    cdip_id = '132';
    yyyymm = '201610';

    %% Initialize variables
    data = {};
    out_dir = '../data/';
    savefile = [cdip_id '_',yyyymm '_' 'WIS.mat'];

    
    %% Find ndbc_id using table: ../ndbc_id_table.csv
    M = csvread('../ndbc_id_table.csv');
    index = find(M(:,1) == str2num(cdip_id));
    ndbc_id = num2str(M(index,2));       %'46219';

    %% Check to see if file was already processed,
    if(isfile([out_dir savefile]))
        load([out_dir savefile]);
    else
        
        %% Here is information from USACE file
        url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/';
        url = [url region '/ST' ndbc_id '/'];                  %Pacific/ST46219/
        yyyy = yyyymm(1:4);                 %- '2009'
        src_dir = [url,yyyy,'/'];
        wis_fname = ['WIS-ocean_waves_ST',ndbc_id,'_',yyyymm,'.nc'];
        wis_url = [src_dir wis_fname];
        try
            disp(['-->Reading: ' wis_url]);
            wis_info = ncinfo(wis_url);
        catch
            disp(['-->Problem loading: ' wis_url]);
            return;
        end

        %% Create a variable list from both template and key file
        wis_vars = {};
        for vid = 1:length(wis_info.Variables)
            wis_vars{vid} = wis_info.Variables(vid).Name;
        end

        %% Get the data into matlab structure
        wis_data = read_nc(wis_url,wis_vars);
        wis_data.ncinfo = wis_info;

        %% Save to temporary outfile
        out_dir = '../data/';
        savefile = [cdip_id '_',yyyymm '_' 'WIS.mat'];
        save([out_dir savefile],['wis_data']);
    end
    
    data = wis_data;

%end


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
        disp(['--> Loading ' the_var]);
        try
            temp = ncread(url,the_var);
        catch
            disp(['Cant read variable ' the_var])
            return
        end
        eval(['data.' the_var '=temp;']);
    end
end


%% Function to correct time epoch from Allie H cdipxww3 
function [ mat_time ] = time_correct(ts)
    toff = datenum(1970,1,1,0,0,0);
    mat_time = ts./(24*60*60) + toff;
end