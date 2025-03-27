%Function to identify Seaglider and uCTD profiles in specified geographic
%bounds and time period
function profNums = identifyProfileNumbers(profiles, minlon, minlat, maxlon, maxlat, startTime, endTime)
        
    inTimeMask = zeros(size(profiles.times));
    inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;

    inRegionMask = zeros(size(profiles.times));    
    inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
        & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;

    profNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1);     
end