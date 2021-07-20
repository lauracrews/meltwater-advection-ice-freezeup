dataDir = [userpath, '/meltwaterAdvection/data/'];
if exist([dataDir, 'healyUnderway.csv'], 'file') == 0 || ...
        exist('waveglider.nc', 'file') == 0 || ...
        exist('profiles.nc', 'file') == 0
    
    disp('Observational data not found. Download data from http://hdl.handle.net/1773/47135 to the directory ~/meltwaterAdvection/data/')
    return
end

%%

metData = table2struct(readtable([dataDir, 'healyUnderway.csv']), 'toScalar', true);

%%
curFile = [dataDir, 'waveglider.nc'];
wvdata.vehicle = ncread(curFile, 'vehicle');
wvdata.lats = ncread(curFile, 'lats');
wvdata.lons = ncread(curFile, 'lons');
wvdata.times = ncread(curFile, 'times');
wvdata.depths = ncread(curFile, 'z');
wvdata.temps = ncread(curFile, 'temps');
wvdata.salts = ncread(curFile, 'salts');
wvdata.CT = ncread(curFile, 'CT');
wvdata.SA = ncread(curFile, 'SA');
wvdata.sigthes = ncread(curFile, 'sigthes');

%%
curFile = [dataDir, 'profiles.nc'];

dataset = cellstr(num2str(ncread(curFile, 'vehicle')));
inds = find(strcmp('  1', dataset) == 1);
for i = 1:length(inds)
    dataset{inds(i)} = 'uCTD';
end

inds = find(~strcmp('  1', dataset) == 1);
for i = 1:length(inds)
    dataset{inds(i)} = ['sg', dataset{inds(i)}];
end

profiles.dataset = dataset;
profiles.qualFlag = ncread(curFile, 'qualFlag');
profiles.z = ncread(curFile, 'z');
profiles.lats = ncread(curFile, 'lats');
profiles.lons = ncread(curFile, 'lons');
profiles.times = ncread(curFile, 'times');
profiles.temps = ncread(curFile, 'temps');
profiles.CT = ncread(curFile, 'CT');
profiles.salts = ncread(curFile, 'salts');
profiles.SA = ncread(curFile, 'SA');
profiles.sigthes = ncread(curFile, 'sigthes');
profiles.pwpFreezeTime = ncread(curFile, 'pwpFreezeTime');
profiles.pwpz = ncread(curFile, 'pwpz');
profiles.pwpFreezeCTprof = ncread(curFile, 'pwpCT');
profiles.pwpFreezeSAprof = ncread(curFile, 'pwpSA');
profiles.pwpFreezePdenProf = ncread(curFile, 'pwpSigthes');


%If surface data is missing, fill with the shallowest available measurement
for k = 1:length(profiles.lats)
    firstTemp = find(~isnan(profiles.temps(:, k)), 1, 'first');
    profiles.temps(1:firstTemp, k) = profiles.temps(firstTemp, k);

    firstSalt = find(~isnan(profiles.salts(:, k)), 1, 'first');
    profiles.salts(1:firstSalt, k) = profiles.salts(firstSalt, k);
    
    firstCT = find(~isnan(profiles.CT(:, k)), 1, 'first');
    profiles.CT(1:firstCT, k) = profiles.CT(firstCT, k);

    firstSA = find(~isnan(profiles.SA(:, k)), 1, 'first');
    profiles.SA(1:firstSA, k) = profiles.SA(firstSA, k);

    firstDense = find(~isnan(profiles.sigthes(:, k)), 1, 'first');
    profiles.sigthes(1:firstDense, k) = profiles.sigthes(firstDense, k);  
end