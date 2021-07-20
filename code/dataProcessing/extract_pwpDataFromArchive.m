function [modelOutput, forc] = extract_pwpDataFromArchive(profNum)
    fileName = [userpath, '/meltwaterAdvection/data/pwpResults/pwpResults_profile', num2str(profNum), '.nc'];
    modelOutput = []; forc = [];
    if ~exist(fileName, 'file')
        disp(['PWP model results for this profile not found.', newline,...
            'Download and unzip the results from http://hdl.handle.net/1773/47135', newline, ...
            'to the directory ~/meltwaterAdvection/data/pwpResults/'])
        return
    end
    
    forc.time = ncread(fileName, 'time');
    forc.heatFlux = ncread(fileName, 'heatFlux');
    
    modelOutput.z = ncread(fileName, 'z');
    modelOutput.config.lat = ncreadatt(fileName, '/', 'lat_observedProfile');
    modelOutput.time = ncread(fileName, 'time');
    modelOutput.temp = ncread(fileName, 'CT');
    modelOutput.sal = ncread(fileName, 'SA');
end