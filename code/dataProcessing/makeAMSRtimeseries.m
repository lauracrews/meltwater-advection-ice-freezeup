%Averages AMSR2 sea ice concentration within the geographic bounds
%specified in pts. Pts is an array in which the first row contains the
%longitudes and the second row contains the corresponding latitudes.
%Returns a timeseries of mean ice concentration on each day between the
%input startTime and endTime

function [times, meanIceConcentration] = makeAMSRtimeseries(pts, startTime, endTime)
    
    load AMSR2_2018.mat
    
%Identify AMRS2 data within the time range between startTime and endTime
    [~, startInd] = min(abs(AMSR2.mattime - startTime));
    [~, endInd] = min(abs(AMSR2.mattime - endTime)); 
    times = AMSR2.mattime(startInd:endInd);
    
    %Identify which AMSR2 data points are within the geographic region
    %specified by pts. 
    [in, on] = inpolygon(AMSR2.lon, AMSR2.lat, pts(1, :), pts(2, :));
    inBoxMask = in; 
    %inBoxMask(on) = 0; %Do not include cells on the edge

    %Iterate through each day and average ice concentration in the region
    count = 1;
    for curTime = startInd:endInd
        curData = AMSR2.SIC(:, :, curTime);
        meanIceConcentration(count) = mean(curData(inBoxMask));
        count = count + 1;
    end
end

    