%Calculates the heat content relative to freezing of a T-S profiles. Uses
%salinity data to calculate freezing point. If surface data is missing,
%fills with the shallowest measement. Then if additional data is missing
%fills with linear interpolation. 
%
%Input integrationDepthType:
% 1 = Integrate to constant depth
% 2 = Integrate to constant isopycnal
% 3 = Integrate to mixed layer depth
%
%Input integrationLimit should be set according to the chosen integrationDepthType
% 1 = depth to integrate to
% 2 = isopycnal to integrate to
% 3 = density change relative ot surface to define the mixed layer
function [heatContent, integrationDepth] = calculateHeatContent(tempProf, saltProf, densProf, depthProf, lat, integrationDepthType, integrationLimit)
    heatContent = nan; integrationDepth = nan;
                
    %Fill in nan values in the middle of the profile
    tempProf = interp1(depthProf(~isnan(tempProf)), tempProf(~isnan(tempProf)), depthProf);
    saltProf = interp1(depthProf(~isnan(saltProf)), saltProf(~isnan(saltProf)), depthProf);
    densProf = interp1(depthProf(~isnan(densProf)), densProf(~isnan(densProf)), depthProf);

    %If surface data is missing, fill with the shallowest available
    %measurement up to the sea surface.
    firstTemp = find(~isnan(tempProf), 1, 'first');
    if isnan(firstTemp); return; end
    tempProf(1:firstTemp) = tempProf(firstTemp);

    firstSalt = find(~isnan(saltProf), 1, 'first');
    if isnan(firstSalt); return; end
    saltProf(1:firstSalt) = saltProf(firstSalt);

    firstDensity = find(~isnan(densProf), 1, 'first');
    densProf(1:firstDensity) = densProf(firstDensity);
       
    %Calculate heat at each depth
    presProf = gsw_p_from_z(-depthProf, lat); 
    freezingPtProf = gsw_CT_freezing(saltProf, presProf);
    cp = cp_t_exact(saltProf, tempProf, presProf);
    heatContentProf = cp .* densProf .* (tempProf - freezingPtProf);

    %If the heat content prof has Nans at some subsurface depths,
    %fill them with linear interpolation
    if numel(find(~isnan(heatContentProf))) > 8 %There are at least 8 values in the profile (it is not NaN throughout). This is an arbitrary minimum number of data points

        %Fill with linear interpolation if there are some NaNs in the profile that need to be filled
        heatContentProf = interp1(depthProf(~isnan(heatContentProf)), heatContentProf(~isnan(heatContentProf)), depthProf);
        
        if integrationDepthType == 1 %Constant depth
            shallowInd = 1;
            [~, deepInd] = min(abs(depthProf - integrationLimit));
        elseif integrationDepthType == 2 %Constant isopycnal
            shallowInd = 1;
            [~, deepInd] = min(abs(densProf - integrationLimit));
        elseif integrationDepthType == 3 %Mixd layer depth (surface density change) 
            shallowInd = 1;
            deepInd = find(densProf - densProf(1) >= integrationLimit, 1, 'first');
        end 
        
        integrationDepth = depthProf(deepInd);
    
        %Depth integrate the profile's heat
        heatContent = trapz(depthProf(shallowInd:deepInd), heatContentProf(shallowInd:deepInd)) / 10^6; %Convert J to MJ
    
    else %The whole profile is NaN - don't calculate heat content
        heatContent = nan;
    end

end