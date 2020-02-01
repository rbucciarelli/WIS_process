%% hurricane_plots.m
%-- Reproduce CDIP hurricane plots

%% Initialize
clear all;
close all;
ww3 = true;                %-- Set to 'true' if ww3 to be plotted
data_dir = '../data/';
hname = 'Matthew';
yyyymm = '201610';
start_date = '20161006';
end_date = '20161012';
start_dn = datenum(start_date,'yyyymmdd');
end_dn = datenum(end_date,'yyyymmdd');

%% Get a list of stations to process
stn_list = dlmread('../matthew_stns.txt');
stn_list = num2str(stn_list,'%03d\n');

%% Iterate over stations
for ii = 1:length(stn_list)
    cdip_id = stn_list(ii,:);               %-- '132'
    disp(['Station ' cdip_id]);
    out_file = [hname '_' cdip_id '.png'];
    %cdip_name = 'FERNANDINA BEACH, FL';


    %% Load in CDIP data
    disp('Loading CDIP Data')
    cdip_fname = [data_dir,'C',cdip_id,'_',yyyymm,'.mat'];
    if (isfile(cdip_fname))
        load(cdip_fname);       %-- C132
        eval(['cdip_data=','C',cdip_id,';']);
    else
        bulk = 0;
        cdip_data = load_CDIP(cdip_id,start_date(1:6),end_date(1:6),bulk);
    end
    cdip_name = char(cdip_data.name);   %-- 'FORT PIERCE, FL BUOY - 134p1'
    cdip_name = cdip_name(1:strfind(cdip_name,' BUOY')-1);  %-- 'FORT PIERCE, FL

    %% Load in WIS data
    disp('Loading WIS Data')
    wis_fname = [data_dir,'A',cdip_id,'_',yyyymm,'.mat'];
    if (isfile(wis_fname))
        load(wis_fname);
        eval(['wis_data=','A',cdip_id,';']);
    else
        wis_data = load_WIS(cdip_id,'Atlantic',yyyymm);
    end

    %eval(['data=' 'A' cdip_id ';']);
    ndbc_id = get_ndbc(cdip_id);

    %% Load in WW3 data
    if(ww3)
        disp('Loading WW3 Data')
        ww3_fname = [data_dir,'W',cdip_id,'_',yyyymm,'.mat'];
        if (isfile(ww3_fname))
            load(ww3_fname);
            eval(['ww3_data=','W',cdip_id,';']);
        else
            %-- Load ww3 data from CDIP thredds server
            bulk = 0;
            ww3_data = load_WW3_CDIP(cdip_id,start_date(1:6),end_date(1:6),bulk);
        end
    end
    
    %% Get Hmax values. For now get them from text file
    maxhs_fname = '../matthew_hmax_up.txt';
    M = importdata(maxhs_fname,'\t',1);
    %-- Find station Hmax 
    si = find(strcmp(cdip_id,M.textdata(:,1)));
    if (si > 0)
        h_max = M.data(si-1);           %-- 8.93 (meters)
        date_max = char(M.textdata(si,2));    %-- '2016/10/07'
        time_max = char(M.textdata(si,3));    %-- '08:39:10'
        t_max = datenum([date_max ' ' time_max],'yyyy/mm/dd HH:MM:SS');
    else
        h_max = 0;
    end

    %-- Set variables for plotting
    the_title = [cdip_id '-' cdip_name];
    m2ft = 3.2804;

    %% Set up figure and plot data
    figure;
    fig = gcf;
    
    %-- Plot cdip data
    idx = get_indices(cdip_data.time,start_dn,end_dn);
    plot(cdip_data.time(idx),cdip_data.hs(idx),'k','LineWidth',0.8);
    hold on;
    %-- Plot WIS data
    idx = get_indices(wis_data.time,start_dn,end_dn);
    plot(wis_data.time(idx),wis_data.hs(idx),'-','Color',[0.5 0.5 0.5]);
    %plot(wis_data.time(idx),wis_data.hs(idx),'Color',[0.25 0.25 0.25]);
    if (ww3)
    %-- Plot WW3 data
        idx = get_indices(ww3_data.time,start_dn,end_dn);
        plot(ww3_data.time(idx),ww3_data.hs(idx),'--','Color',[0.5 0.5 0.5]);
    end

    xlims = [start_dn end_dn];
    ylims = [0 8];

    %-- Plot Hmax
    if (h_max ~= 0)
        %-- plot a vertical line at t_max,h_max
        plot([t_max t_max],ylims,'k')
        %-- Plot a textbox with 
        dist = 0.325;
        str = [num2str(h_max) ' m'];
        %P = patch(xdim,ydim,[0.9 0.9 0.9]);
        rectangle('Position',[t_max-dist 0.175 dist*1.75 dist*1.75],'Curvature',0.75, ...
            'FaceColor',[0.95 0.95 0.95]);
        text(t_max-0.215,0.525,str,'FontSize',9);
        %annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end

    font_label = 9;
    yyaxis left;
    xlabel('Day (UTC)','FontSize',font_label);
    ylabel('Hs (m)','FontSize',font_label);
    set(gca,'XLim',xlims);
    set(gca,'YLim',ylims);
    set(gca,'XTickLabel',datestr(get(gca,'XTick'),'dd'))
    set(gca,'YTick',ylims(1):2:ylims(2));
    yTicks = get(gca,'YTick');
    grid on;
    title(the_title, 'FontSize',font_label);
    %-- 2nd yaxis in feet
    yyaxis right;
    set(gca,'YLim',ylims * m2ft);
    set(gca,'YTick',round(yTicks * m2ft))
    ylabel('Hs (ft)','FontSize',font_label);
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8.33 2.78];
    print(['../export/' out_file],'-dpng','-r0');
    %print('5by3DimensionsFigure','-dpng','-r0');
    close all;
end


%% Function to return the indices of array between start/end values
function [ idx ] = get_indices(A,start,finish)
    idx = find((A >= start) & (A <= finish));
end

%% Function to get ndbc id given cdip_id
function [ ndbc_id ] = get_ndbc(cdip_id)
    %- Find ndbc_id using table: ../ndbc_id_table.csv
    M = csvread('../ndbc_id_table.csv');
    index = find(M(:,1) == str2num(cdip_id));
    ndbc_id = num2str(M(index,2));       %'46219';
end

%% Function to correct time epoch 
function [ mat_time ] = time_correct(ts)
    toff = datenum(1970,1,1,0,0,0);
    mat_time = ts./(24*60*60) + toff;
end
