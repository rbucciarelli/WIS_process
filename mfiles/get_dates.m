%% Function to return char array of dates YYYYMM 
function [ dates ] = get_dates(start_time,end_time)
    syear = start_time(1:4);
    smonth = start_time(5:6);
    eyear = end_time(1:4);
    emonth = end_time(5:6);
    years = str2num(syear):str2num(eyear);
    dates = [];
    N = 0;
    index = 1;
    for i = 1:length(years)
        yyyy = num2str(years(i),'%04i');
        if (yyyy == syear)
            months = str2num(smonth):str2num(emonth);
        elseif (yyyy == eyear)
            months = 1:str2num(emonth);
        else
            months = 1:12;
        end
        for j = 1:length(months)
            mm = num2str(months(j),'%02i');
            dates(index) = str2num([yyyy,mm]);
            index = index + 1;
        end
    end
end