%Averages ERA5 data within the geographic coordinates specified in
%pts. Pts is an array in which the first row contains the
%longitudes and the second row contains the corresponding latitudes.
%Returns timeseries between startTime and endTime of the spatially-averaged 
%parameters given as strings in params.  
%
%Optionally make the function return [avgTimeseries, numERA5cells] to see
%how many ERA5 grid cells go into the average
function avgTimeseries = makeERA5timeseries(pts, startTime, endTime, params)
    filename = 'ERA5.nc';
    
    %Read in ERA5 time, convert to Matlab datetime format. Identify
    %earliest and latest ERA5 measurements to be included in the timeseries
    eraTime = double(ncread(filename, 'time'))/24 + datenum('1900-01-01 00:00:00');    
    [~, startInd] = min(abs(eraTime - startTime));
    [~, endInd] = min(abs(eraTime - endTime)); 
    time = num2cell(eraTime(startInd:endInd));

    %Read in the entire longitude/latitude arrays for the data in ERA5.nc
    eraLon = double(ncread(filename, 'longitude'));
    eraLat = double(ncread(filename, 'latitude'));
    
    %Identify the indices of ERA5 data points within the desired geographic
    %bounds. These indices are used later to only read in the subset when
    %loading all the parameters
    minlon = min(pts(1, :)); maxlon = max(pts(1,:));
    minlat = min(pts(2, :)); maxlat = max(pts(2,:));    
    [~, startLonInd] = min(abs(eraLon - minlon));
    [~, endLonInd] = min(abs(eraLon - maxlon));
    [~, startLatInd] = min(abs(eraLat - maxlat));
    [~, endLatInd] = min(abs(eraLat - minlat));

    %Reload longitude and latitude to only include the desired geographic range
    eraLon = double(ncread(filename, 'longitude', startLonInd, [endLonInd - startLonInd + 1]));
    eraLat = double(ncread(filename, 'latitude', startLatInd, [endLatInd - startLatInd + 1]));
    
    %Convert longitude and latitude to arrays
    [eraLon, eraLat] = meshgrid(eraLon, eraLat);
    eraLon = eraLon'; eraLat = eraLat';

    %Identify points that are actually in the bounding box given in poins
    %(not all points within the min/max lat/lon bounds will be inside the
    %polygon if the polygon is irregularly shaped)
    [in, on] = inpolygon(eraLon, eraLat, pts(1, :), pts(2, :));
    inBoxMask = in; 
%     inBoxMask(on) = 0; %Do not include cells on the edge

%     if sum(double(inBoxMask(:))) > 4
%         disp('More than four ERA5 cells used')
%     end
    
    %Iterate through the atmospheric variables, 
    timeseries = cell(length(time), length(params)); 
    for i = 1:length(params)
        curParam = params{i};
        
        %Load the data for the current atmospheric variable at alltimesteps
        data = double(ncread(filename, curParam, [startLonInd startLatInd startInd],...
            [endLonInd - startLonInd + 1,  endLatInd - startLatInd + 1, endInd - startInd + 1]));
        
        %Iterate through all timesteps and average within the geographic
        %boudnaries at each timestep
        for curTime = 1:length(time)
            curData = data(:, :, curTime);
            timeseries{curTime, i} = mean(curData(inBoxMask));
        end
    end
    
    avgTimeseries = cell2struct([time, timeseries], ['time', params], 2);
    numERA5cells = sum(double(inBoxMask(:)));
end