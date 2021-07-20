%Makes a data structure of temperature and salinity profiles collected by
%all platforms (all Seagliders and uCTD) in a specified geographic region
%and date range

%Calls calculateHeatContent to depth-integrate the heat in the temperature
%profiles
function profiles = profilesInRegion_timeseries(pts, startDay, endDay, integrationType, integrationLimit)
%%    
    minDepth = 1;
    maxDepth = 300;
    depthProf = minDepth:1:maxDepth;
    profiles.z = depthProf';

    %Arrays to hold unknown number of profiles
    allTemps = nan .* ones(length(depthProf), 1);
    allSalts = nan .* ones(length(depthProf), 1);
    allSigthes = nan .* ones(length(depthProf), 1);
    allIntegrationDepths = nan;
    allHeatContents = nan;
    
    %Iterate through each of the possible sources of data, find all
    %profiles from that source that fall within current lat/lon and
    %time boundaries
    profCount = 1;
    for dataset = 1:8
        switch dataset
            case 1
                load uCTD_SODA_2018.mat
                load goodProfiles.mat
                txt = 'uCTD'; 
                goodProfiles(find(goodProfiles == 129, 1)) = [];
            case 2
                load sg198_data.mat
                load goodProfiles_sg198.mat
                txt = 'sg198';
            case 3
                load sg199_data.mat
               load goodProfiles_sg199.mat
                txt = 'sg199';
            case 4
                load sg227_data.mat
                load goodProfiles_sg227.mat
                txt = 'sg227';
            case 5
                load sg228_data.mat
                load goodProfiles_sg228.mat
                txt = 'sg228';
            case 6 
                load sg229_data.mat
                load goodProfiles_sg229.mat
                txt = 'sg229';
            case 7
                load sg230_data.mat
                load goodProfiles_sg230.mat
                txt = 'sg230';
            case 8
                load sg226_data.mat
                load goodProfiles_sg226.mat
                txt = 'sg226';
        end

        if dataset == 1 %uCTD data
            lats = uCTD.lat';
            lons = uCTD.lon';
            times = uCTD.time';
        else % glider data - need to average lon, lats, times over dive due to sawtooth           
            minDepthInd = find(sg_data.z == minDepth);
            maxDepthInd = find(sg_data.z == maxDepth);
            lats = nanmean(sg_data.lat(minDepthInd:maxDepthInd, :), 1)'; %Average lat, lon, time over upper part of dive     
            lons = nanmean(sg_data.lon(minDepthInd:maxDepthInd, :), 1)';
            times = nanmean(sg_data.time(minDepthInd:maxDepthInd, :), 1)';
            
%Optional code to plot lat/lon vs salinity (a proxy for depth) to see that the 
%upgoing and downgoing dives are being used
%             if dataset == 1
%                 figure
%                 scatter3(sg_data.lon(:, 53), sg_data.lat(:, 53), sg_data.S(:, 53),'filled', 'b')
%                 hold on
%                 scatter3(sg_data.lon(:, 54), sg_data.lat(:, 54), sg_data.S(:, 54),'filled', 'r')
%                 scatter3(sg_data.lon(:, 55), sg_data.lat(:, 55), sg_data.S(:, 55),'filled', 'b')
%                 scatter3(sg_data.lon(:, 56), sg_data.lat(:, 56), sg_data.S(:, 56),'filled', 'r')
%                 scatter3(sg_data.lon(:, 57), sg_data.lat(:, 57), sg_data.S(:, 57),'filled', 'b')
%                 scatter3(sg_data.lon(:, 58), sg_data.lat(:, 58), sg_data.S(:, 58),'filled', 'r')
%                 set(gca, 'zdir', 'reverse')

                %Optional code to see horizontal distance between profiles
%                 delta = m_lldist(sg_data.lon(1, :), sg_data.lat(1, :)); %Convert km to m
%                 figure; hist(delta(:), [0:.5:10])
%             end

                
        end

        %Find points within the lat/lon box
        inBoxMask = inpolygon(lons, lats, pts(1, :), pts(2, :));

        %Only use the good profiles
        goodProfMask = ismember(1:1:length(lats), goodProfiles)';
        
        %Combine all of the above filters to determine which profiles
        %should be included in the data structure 
        profNums = zeros(size(times));
        profNums(inBoxMask == 1 & goodProfMask == 1) = 1;
        profNums = find(profNums == 1);

        if isempty(profNums) %No profiles from this platform meet the criteria
            continue             
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Add profiles to the matrix of profile data
        if dataset == 1 %uCTD data
            temps = uCTD.temperature;
            salts = uCTD.salinity;
            firstMeasurementIndex = find(uCTD.z == minDepth);
            depths = uCTD.z(firstMeasurementIndex:end);
            for i = 1:length(profNums)
                
                profNum = profNums(i);

                if times(profNum) < startDay | times(profNum) > endDay %Profile is outside of the time range
                    continue
                end
                
                profiles.dataset{profCount} = txt;
                profiles.lats(profCount) = uCTD.lat(profNum);
                profiles.lons(profCount) = uCTD.lon(profNum);
                profiles.times(profCount) = uCTD.time(profNum);
                
                saltProf = salts(firstMeasurementIndex:end, profNum);
                tempProf = temps(firstMeasurementIndex:end, profNum);
                
                allTemps(:, profCount) = interp1(depths(isfinite(tempProf)), tempProf(isfinite(tempProf)), depthProf);
                allSalts(:, profCount) = interp1(depths(isfinite(saltProf)), saltProf(isfinite(saltProf)), depthProf);

                %Calculate conservative temperature, absolute salinity,
                %potential density
                [SA, ~, ~] = gsw_SA_Sstar_from_SP(saltProf, uCTD.pressure(firstMeasurementIndex:end, profNum), uCTD.lon(profNum), uCTD.lat(profNum));
                CT = gsw_CT_from_t(SA, tempProf, uCTD.pressure(firstMeasurementIndex:end, profNum));
                densProf = gsw_rho(SA, CT, 0);
                
                allCTs(:, profCount) = interp1(depths(isfinite(CT)), CT(isfinite(CT)), depthProf);
                allSAs(:, profCount) = interp1(depths(isfinite(SA)), SA(isfinite(SA)), depthProf);
                allSigthes(:, profCount) = interp1(depths(isfinite(densProf)), densProf(isfinite(densProf)), depthProf);
                
                %Depth integrate the heat content
                [heatContent, integrationDepth] = calculateHeatContent(allCTs(:, profCount), allSAs(:, profCount), allSigthes(:, profCount), depthProf, profiles.lats(profCount), integrationType, integrationLimit);
                allHeatContents(profCount) = heatContent;
                allIntegrationDepths(profCount) = integrationDepth;
                                  
                profCount = profCount + 1;

            end

        else %glider data 
            firstMeasurementIndex = find(sg_data.z == minDepth);
            depths = sg_data.z(firstMeasurementIndex:end);
            for i = 1:length(profNums)

                profNum = profNums(i);
                
                if times(profNum) < startDay | times(profNum) > endDay
                    continue
                end
                
                profiles.dataset{profCount} = txt;
                profiles.lats(profCount) = lats(profNum);
                profiles.lons(profCount) = lons(profNum);
                profiles.times(profCount) = times(profNum);
                
                tempProf = sg_data.T(firstMeasurementIndex:end, profNum);
                saltProf = sg_data.S(firstMeasurementIndex:end, profNum);
                allTemps(:, profCount) = interp1(depths(isfinite(tempProf)), tempProf(isfinite(tempProf)), depthProf);
                allSalts(:, profCount) = interp1(depths(isfinite(saltProf)), saltProf(isfinite(saltProf)), depthProf);
                 
                %Calculate conservative temperature, absolute salinity,
                %potential density
                [SA, ~, ~] = gsw_SA_Sstar_from_SP(saltProf, gsw_p_from_z(-depths, lats(profNum)), lons(profNum), lats(profNum));
                CT = gsw_CT_from_t(SA, tempProf, gsw_p_from_z(-depths, lats(profNum)));
                densProf = gsw_rho(SA, CT, 0);

                allCTs(:, profCount) = interp1(depths(isfinite(CT)), CT(isfinite(CT)), depthProf);
                allSAs(:, profCount) = interp1(depths(isfinite(SA)), SA(isfinite(SA)), depthProf);
                allSigthes(:, profCount) = interp1(depths(isfinite(densProf)), densProf(isfinite(densProf)), depthProf);
                  
                %Depth integrate the heat content
                [heatContent, integrationDepth] = calculateHeatContent(allTemps(:, profCount), allSalts(:, profCount), allSigthes(:, profCount), ...
                    depthProf, profiles.lats(profCount), integrationType, integrationLimit);
                
                allHeatContents(profCount) = heatContent;
                allIntegrationDepths(profCount) = integrationDepth;

                profCount = profCount + 1;

            end               
        end
    end %End of loop through each dataset
       
    %Create the data structure containing all data for the profiles
    if profCount ~= 1

        profiles.temps = allTemps;
        profiles.CT = allCTs;
        profiles.salts = allSalts;
        profiles.SA = allSAs;
        profiles.sigthes = allSigthes;
        profiles.integrationDepths = allIntegrationDepths;
        profiles.heatContents = allHeatContents;

        %Update arrays to be sorted chronologically
        [~, idx] = sort(profiles.times);
        profiles.dataset = profiles.dataset(idx);
        profiles.lats = profiles.lats(idx);
        profiles.lons = profiles.lons(idx);
        profiles.times = profiles.times(idx);
        profiles.temps = profiles.temps(:, idx);
        profiles.CT = profiles.CT(:, idx);
        profiles.salts = profiles.salts(:, idx);
        profiles.SA = profiles.SA(:, idx);
        profiles.sigthes = profiles.sigthes(:, idx);
        profiles.integrationDepths = profiles.integrationDepths(idx);
        profiles.heatContents = profiles.heatContents(idx);

      else
        profiles = [];
        disp('No Profiles Found');
    end
end

            