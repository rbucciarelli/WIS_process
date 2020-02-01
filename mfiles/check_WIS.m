%% check_WIS.m
%-------------------------------------------------------------------------
%- Check WIS output for duplicate freqs and start/end dates.
%- The WIS model data is stored in NetCDF format on following server:
%- https://chlthredds.erdc.dren.mil/thredds/catalog/wis/catalog.html
%- https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/Pacific/ST46221/2003/WIS-ocean_waves_ST46221_200303.nc
%- Files are stored as individual months in specific regions.
%-------------------------------------------------------------------------

clear all;

data = {};

%% Initialize variables
src_url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/';
region_list = {'Atlantic', 'Pacific', 'GulfOfMexico'};
region_list = {'Atlantic'};

for ri = 1:length(region_list)
    region = string(region_list(ri,:));
    dup_file = strcat('../',region,'_duplicate_freqs.txt');
    missing_file = strcat('../',region,'_missing.txt');
    %- Get a list of stations in region
    region_file = strcat('../cdip-',region,'.csv');
    M = csvread(region_file);
    stn_list = M(:,1);
    stn_list = num2str(stn_list, '%03g');
    for si = 1:length(stn_list)
        cdip_id = stn_list(si,:);
        %- Find ndbc_id using table: ../cdip-Pacific.csv
        index = find(M(:,1) == str2num(cdip_id));
        ndbc_id = num2str(M(index,2));       %'46219';
        %- Get/Set start and end dates
        start_date = num2str((M(index,3)));
        start_num = datenum(start_date,'yyyymm');
        start_date = '201301';
        end_date = '201712';    %- Pacific
        end_date = '201703';   %- Atlantic
        
        end_num = datenum(end_date,'yyyymm');
        %-- Set data source URL
        url = strcat(src_url,region,'/ST',ndbc_id,'/');                  %Pacific/ST46219/
        %-- Set a list of dates to check netcdf file information      
        if(start_num < end_num)
            fprintf('%s/%s\n',cdip_id,ndbc_id);
            date_list = get_dates(start_date,end_date);
            %-- Call a function to generate a filename list to process
            file_list = gen_file_list(url,ndbc_id,date_list);
            %-- Iterate over dates checking year first and last month 
            bad_list = [];
            no_list = [];
            for di = 1:length(date_list)
                fname = string(file_list(di));
                the_date = num2str(date_list(di));
                yr = the_date(1:4);
                mo = the_date(5:6);
                if (mo == "01") || (mo == "12")  
                    %fprintf('-->%s\n',the_date);
                    try
                        finfo = ncinfo(fname);
                        fvars = finfo.Variables;
                        freqs = ncread(fname,'waveFrequency');
                        %-- Check if there are duplicate freqs
                        idx = find(diff(freqs) == 0);
                        if(length(idx) > 0)
                            disp(strcat("--> ",the_date,': DUPLICATE FREQS'));
                            bad_list = [bad_list str2num(the_date)];
                        else
                            disp(strcat("--> ",the_date,': GOOD'));
                        end
                    catch
                        % code to handle the error
                        disp(strcat("--> ",the_date,': NODATA'));
                        no_list = [no_list str2num(the_date)];
                    end                       
                end
            end
            if(length(bad_list) > 1)
                fid = fopen(dup_file,'a');
                fprintf(fid,'%s\t%s\t%g-%g\n',cdip_id,ndbc_id,bad_list(1),bad_list(end));
                fclose(fid);
            end
            if(length(no_list) > 1)
                fid = fopen(missing_file,'a');
                fprintf(fid,'%s\t%s\t%g-%g\n',cdip_id,ndbc_id,no_list(1),no_list(end));
                fclose(fid);
            end            
                
           
        end
        
    end

end

%-- Generate a list of files to process given date_list
function [ file_list ] = gen_file_list(url,ndbc_id,date_list)
    file_list = {};
    for i = 1:length(date_list)
        yyyymm = num2str(date_list(i));
        yyyy = yyyymm(1:4);                 %- '2009'
        src_dir = strcat(url,yyyy,'/');
        fname = ['WIS-ocean_waves_ST',ndbc_id,'_',yyyymm,'.nc'];
        file_list{i} = strcat(src_dir,fname);
    end
end

function [ data ] = check_nc_file(fname)
    data = {};
    disp(['Checking ' fname]);
    
    %ncid = netcdf.open(strcat(src_dir,fname])),'NC_NOWRITE');
    
    
end

% region = 'Pacific';         %- Need to figure out which region buoy is in.
% cdip_id = '067';
% start_time = '199801';      % YYYYMM
% end_time = '201312';
% 
% 
% %- List of variables to extract from NetCDF file
% var_list = {'time', 'waveFrequency', 'directionalWaveEnergyDensity', ...
%     'waveDirectionBins', 'waveHs', 'waveTp', 'waveTm', ...
%     'latitude', 'longitude'};
% 
% %- Find ndbc_id using table: ../ndbc_id_table.csv
% M = csvread('../ndbc_id_table.csv');
% index = find(M(:,1) == str2num(cdip_id));
% ndbc_id = num2str(M(index,2));       %'46219';
% 
% 
% %- WW3 files are by month, create an array of dates
% date_list = get_dates(start_time,end_time); 
% 
% si = 1;     %- start_index
% 
% url = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/';
% url = [url region '/ST' ndbc_id '/'];                  %Pacific/ST46219/
% 
% %- Iterate over dates and load netcdf files from usace thredds server
% for i = 1:length(date_list)
%     yyyymm = num2str(date_list(i));
%     yyyy = yyyymm(1:4);                 %- '2009'
%     src_dir = [url,yyyy,'/'];
%     fname = ['WIS-ocean_waves_ST',ndbc_id,'_',yyyymm,'.nc'];
%     disp(['Loading ' fname]);
%     ncid = netcdf.open([src_dir fname],'NC_NOWRITE');
%     [ndims,nvars,natts,unlimdimID] = netcdf.inq(ncid);
%     finfo = ncinfo([src_dir fname]);
%     %-- Read depth from global attributes
%     depth = netcdf.getAtt(ncid, netcdf.getConstant('NC_GLOBAL'),'depth');
%     
%     %-- Create a list of variables from netcdf
%     fvars = {};
%     for vid = 1:length(finfo.Variables)
%         fvars{vid} = finfo.Variables(vid).Name;
%     end
%     
%     %% Iterate over vars of interest and get data
%     for vid = 1:length(var_list)
%         the_var = var_list{vid};
%         %-- Find the index of this variableda
%         vindex = find(strcmp(fvars,the_var)) - 1;
%         disp(['--> Loading ' the_var]);
%         eval([the_var '=netcdf.getVar(ncid,' num2str(vindex) ');']);
%     end
%     netcdf.close(ncid);
%     
%     %% Aggregate monthly data based on start/end times
%     mat_time = time_correct(double(time));
%     %-- Figure out end index
%     ei = si + length(time) - 2;
%     
%     data.lat = latitude;
%     data.lon = longitude;
%     data.depth = depth;
%     data.f = waveFrequency;
%     data.dir = waveDirectionBins;
%    
%     %-- Calculate a1,b1,a2,b2
%     [ND,NF,NT] = size(directionalWaveEnergyDensity);
%     %-- Function requires 2d spectra to be [time x freq x dir]
%     energy2d = directionalWaveEnergyDensity; 
%     %-- Calculate 1D energy
% %     rdir = deg2rad(waveDirectionBins);
% %     dtheta=abs(rdir(2)-drir(1));
% %     icnt = 1;
% %     for i = 1:length(data.f)
% %         for j = 1:length
% %     energy = b.sp1d{index+1} = sum(sp2d')*dtheta;
% %    
%     
%     ab_data = calc_a_b_wis(waveFrequency,waveDirectionBins,energy2d);
%     data.bw = ab_data.bw;
%     data.time(si:ei) = mat_time(1:end-1);
%     %energy = squeeze(energy2d(1,:,1:end-1));
%     data.energy(:,si:ei) = ab_data.energy(:,1:end-1);
%     data.a0(:,si:ei) = ab_data.a0(:,1:end-1);
%     data.a1(:,si:ei) = ab_data.a1(:,1:end-1);
%     data.b1(:,si:ei) = ab_data.b1(:,1:end-1);
%     data.a2(:,si:ei) = ab_data.a2(:,1:end-1);
%     data.b2(:,si:ei) = ab_data.b2(:,1:end-1);
%     data.hs(si:ei) = waveHs(1:end-1)';
%  
%     si = ei + 1;
%     
%     %-- Save annual output in case FRF THREDDS server closes connection
%     
%     
% end
% 
% %% Save data to .mat file
% eval(['A' cdip_id '=data;']);   
% out_dir = '../data/';
% savefile = ['A',cdip_id,'.mat'];
% save([out_dir savefile],['A' cdip_id])
% 
% 
% 
% %% Function to correct time epoch from Allie H cdipxww3 
% function [ mat_time ] = time_correct(ts)
%     toff = datenum(1970,1,1,0,0,0);
%     mat_time = ts./(24*60*60) + toff;
% end