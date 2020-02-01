%% WIS_aggregate.m
%-------------------------------------------------------------------------
%- Aggregate monthly individual WIS mat files between start and end dates.
%- Will download any missing files from FRF THREDDS server if not found 
%- locally.
%-------------------------------------------------------------------------

function [ big_data ] = WIS_aggregate(cdip_id,start_date,end_date)

    %% Initialize variables
    big_data = {};
    src_dir = '../data/';

    %% Get WIS region
    [region,stn_start] = WIS_region(cdip_id);

    %% Call get_dates.m to generate an array of dates from start to end date
    date_list = get_dates(start_date,end_date);
    
    
    %% Iterate over dates
    %-- Define start index
    si = 1;
    for i = 1:length(date_list)
        yyyymm = num2str(date_list(i));
        %-- WIS_to_mat.m creates files: 132_201607_WIS.mat, contains 
        %-- struct wis_data
        fname = [cdip_id '_' yyyymm '_WIS.mat'];
        if(isfile([src_dir fname]))
            load([src_dir fname]);
        else
            %-- run script to load data and save to disk
            try
                wis_data = WIS_to_mat(cdip_id,region,yyyymm);
                load([src_dir fname]);
            catch
                disp(['WIS_to_mat Failed ' yyyymm]);
                return
            end
        end
        NT = length(wis_data.time);
        NF = length(wis_data.waveFrequency);
        ND = length(wis_data.waveDirectionBins);
        
        ei = si + NT - 2;
        %-- Initialize if start date
        if(i == 1)
            big_data = wis_data;
            big_data.directionalWaveEnergyDensity = zeros(ND,NF,NT-1);
        end
        
        %-- Find which vars not to increment
        dim_list = {'waveFrequency';'ncinfo';'waveDirectionBins'};
        var_list = fieldnames(wis_data);
        for vi = 1:length(var_list)
            var_name = var_list{vi};
            var_size = size([wis_data.(var_name)]);
            input_data = [wis_data.(var_name)];
            output_data = [big_data.(var_name)];
            %-- first increment variables with unlimited time-bounds
            if (~isempty(find(var_size == NT) > 0) && (var_size(end) == 1)) 
                output_data(si:ei) = input_data(1:end-1);
            elseif (strcmp(var_name,'directionalWaveEnergyDensity'))
                idx = si-1;
                for dd=1:ND
                    for ff=1:NF
                        for tt=1:NT-1
                            output_data(dd,ff,idx+tt) = input_data(dd,ff,tt);
                        end
                    end
                end
            end
            [big_data.(var_name)] = output_data;
        end
        si = ei + 1;

        
    end
    
    %% Save to temporary outfile
    out_dir = '../data/';
    savefile = [cdip_id '_',start_date '-' end_date '_' 'WIS.mat'];
    wis_data = big_data;
    save([out_dir savefile],['wis_data']);
    
end

