%Combines data from all the Wave Gliders into one data structure
function wvdata = loadWaveglider
    %Empty arrays  
    times = []; lons = []; lats = []; temps = []; salts = []; depths = []; vehicles = [];   
    waveglider_files = dir([userpath, '/meltwaterAdvection/data/waveglider/*.mat']);
    
    %Iterate through each Wave Glider file 
    for i = 1:length(waveglider_files)
        load(waveglider_files(i).name)
         
        %Extract coordinates and times for all measurements
        lons = [lons; [SV3.lon]']; lats = [lats; [SV3.lat]'];
        times = [times; [SV3.time]'];

        curVehicle = str2num(waveglider_files(i).name(5:7)); %Extract vehicle number from file name. 
        vehicles = [vehicles; curVehicle .* ones(size([SV3.lon]'))];
        
        %Due to the format of the Wave Glider files, need to iterate
        %through and extract the temperature and salinity at each CTD
        %separately
        for j = 1:length(SV3)
            depths = [depths; SV3(j).CTdepth(1, 1:3)];
            salts = [salts; SV3(j).salinity(1, 1:3)];
            temps = [temps; SV3(j).watertemp(1, 1:3)];
        end
    end
    
    %Create data structure    
    wvdata.vehicle = vehicles;
    wvdata.lats = lats;
    wvdata.lons = lons;
    wvdata.times = times;
    wvdata.depths = depths;
    wvdata.temps = temps;
    wvdata.salts = salts;

    %Rename field to be consistent with the other data
    wvdata = renameStructField(wvdata, 'depths', 'z');
    
    %Calculate conservative temperature, absolute salinity, potential density
    [SA, ~, ~] = gsw_SA_Sstar_from_SP(wvdata.salts, sw_pres(wvdata.z, repmat(wvdata.lats, [1, 3])), wvdata.lons, wvdata.lats);
    CT = gsw_CT_from_t(SA, wvdata.temps, sw_pres(wvdata.z,repmat(wvdata.lats, [1, 3])));
    sigthes = gsw_rho(SA, CT, 0); %Reference sea pressure of 0 (the surface)
    
    wvdata.CT = CT; wvdata.SA = SA; wvdata.sigthes = sigthes;

    
end