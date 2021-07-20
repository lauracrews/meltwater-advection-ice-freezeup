%Reads CyroSat-2 .asc files into Matlab and extracts lat/lon and dynamic ocean topography (DOT) 
%Calculates time of each data point (which are referenced to the satelite's equator crossing
%Files were downloaded from http://rads.tudelft.nl/rads/data/authentication.cgi

clearvars -except AMSR2 profiles wvdata metData 
close all;
rawDataDir = [userpath, '/meltwaterAdvection/data/DOT/raw/'];
files = [dir([rawDataDir, '*109.asc']); dir([rawDataDir, '*110.asc'])];
if isempty(files)
    disp('Dynamic ocean topography data not found.', newline, 'See instructions for data access at https://github.com/lauracrews/meltwaterAdvection')
    return
end

saveFiles = true;
saveMatDir = [userpath, '/meltwaterAdvection/data/DOT/mat/'];
if ~exist(saveMatDir, 'dir'); mkdir(saveMatDir); end
removeNans = true;
%%

%Empty arrays to combine data from all files into 
allLats = [];
allLons = [];
allTimes = [];
allDOT = [];

%Iterate through all files
for i = 1:length(files)
    curFile = importdata([rawDataDir, files(i).name], '\t');
    disp(['Processing DOT file ', num2str(i), ' of ', num2str(length(files))])
    
    %string manipulation to extract the time the satellite crossed the equator
    eqTime = curFile{6}; 
    ind = strfind(eqTime, '(');
    eqTime = datenum(eqTime(ind:end-1));
    
    startLine = 1; %Find first line of data
    while contains(curFile(startLine), '#')
        startLine = startLine + 1; 
    end
      
    %Read in each line of data
    times = nan.*ones(size(curFile(startLine:end)));
    lats = nan.*ones(size(curFile(startLine:end)));
    lons = nan.*ones(size(curFile(startLine:end)));
    DOT = nan.*ones(size(curFile(startLine:end)));    
    lineCount = 1;
    for lineNum = startLine:length(curFile)
        
        curLine = curFile{lineNum};
        while strcmp(curLine(1), ' ') %Remove leading spaces
            curLine = curLine(2:end);
        end
        
        inds = strfind(curLine, ' '); %Data divided by spaces
        indlon = strfind(curLine, '-'); %Dash (negative sign) indicates start of the lon measurement
        times(lineCount) = str2double(curLine(1:inds(1)-1));
        lats(lineCount) = str2double(curLine(inds(2)+1:inds(3)-1));
        lons(lineCount) = str2double(curLine(indlon:inds(4)-1));
        DOT(lineCount) = str2double(curLine(inds(end)+1:end));       
        lineCount = lineCount + 1;
    end
    
    %Calcuate time of measurement, referenced to when the satelite crosses
    %the equator
    times = (times ./ (60*60*24)) + eqTime; %Convert seconds to days
    
    if removeNans
        times(isnan(DOT)) = [];
        lons(isnan(DOT)) = [];
        lats(isnan(DOT)) = [];
        DOT(isnan(DOT)) = [];
    end
    
    %Add in the data from this file
    allTimes = [allTimes; times];
    allLats = [allLats; lats];
    allLons = [allLons; lons];
    allDOT = [allDOT; DOT];
    
    if saveFiles
        cd(saveMatDir)
        save([files(i).name(1:end-4), '.mat'], 'times', 'lats', 'lons', 'DOT')
    end
end
    
lats = allLats;
times = allTimes;
lons = allLons;
DOT = allDOT;


if removeNans
    nanInds = [find(isnan(times)); find(isnan(lats)); find(isnan(lons)); find(isnan(DOT))];
    times(nanInds) = [];
    lons(nanInds) = [];
    lats(nanInds) = [];
    DOT(nanInds) = [];
end

%Sort to be in order of collection time
[~, idx] = sort(times);
lats = lats(idx);
lons = lons(idx);
DOT = DOT(idx);
times = times(idx);

cd([userpath, '/meltwaterAdvection/data/'])
save('allDOT.mat', 'times', 'lats', 'lons', 'DOT')
       