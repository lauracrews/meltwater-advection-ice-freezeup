%Loads MODIS sea surface temperature data from a specified modis file
%within given geographic bounds. Specifying geographic range helps with
%speed, especially when iterating through multiple files.

%Returns arrays of sea surface temperature and corresponding
%latitude/longitude coordinates, as well as the starting and ending time
%for data included in the image. 

%modisTitle is a handy description of the data, good for using on colorbars
%or figure titles (or just printing what data is used)

function [lon, lat, sst, startTime, endTime, modisTitle] = loadMODISsst(modisFile, minlon, maxlon, minlat, maxlat)
    
    %Determine if this is a 1-day or 8-day file. Extract the time of the
    %image from the file title
    if contains(modisFile, 'DAY') %Daily file
        startTime = datenum([modisFile(17:18), ' ', modisFile(19:20), ' ', modisFile(13:16)]);
        modisTitle = ['MODIS-Terra NSST ', datestr(startTime, 'mmm dd')];
        endTime = startTime + 1;
    elseif contains(modisFile, '8D') %8-day avereage
        modisStart = modisFile(6:8); modisStart = [datestr(str2num(modisStart) + 1, 'mmm dd'), ' 2018']; %Start day of the average
        modisEnd = modisFile(13:15); modisEnd = [datestr(str2num(modisEnd) + 1, 'mmm dd'), ' 2018']; %End day of the average  
        modisTitle = ['MODIS-Terra NSST ', datestr(modisStart, 'mmm dd'), ' to ', datestr(modisEnd, 'mmm dd')];
        startTime = datenum(modisStart); endTime = datenum(modisEnd);
    else
        disp('Unknown file format')
        return
    end

    sst_lon = ncread(modisFile, 'lon'); sst_lat = ncread(modisFile, 'lat'); %In descending order from 90 to -90
    %Identify MODIS cells to load and plot
    ind1_lon = find(sst_lon >= minlon, 1, 'first');
    ind2_lon = find(sst_lon <= maxlon, 1, 'last');
    ind1_lat = find(sst_lat <= maxlat, 1, 'first');
    ind2_lat = find(sst_lat >= minlat, 1, 'last');

    lon = sst_lon(ind1_lon:ind2_lon);
    lat = sst_lat(ind1_lat:ind2_lat);
    [lon, lat] = meshgrid(lon, lat);
    clear sst_lat sst_lon

    %Read in sea surface temperature data
    sst = ncread(modisFile, 'sst', [ind1_lon, ind1_lat], [ind2_lon - ind1_lon + 1, ind2_lat - ind1_lat + 1]);
    sst = sst';

    %Read quality flag
    qual_sst = ncread(modisFile, 'qual_sst', [ind1_lon, ind1_lat], [ind2_lon - ind1_lon + 1, ind2_lat - ind1_lat + 1]);
    qual_sst = qual_sst';

    %Reject SST values for which the quality flag exceeds 1
    sst(qual_sst > 1) = nan;
end