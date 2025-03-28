%pwp_StaticAdj: Subroutine to do static instability adjustment
%
%	This routine relies on the fact that density increases
%	are always forced from the top (there is no in situ
%	cooling mechanism) so any instability will result in a 
%	downward deepening from the top. The algorithm starts
%	every time with a mixed layer depth of 1 cell and goes
%	downward until static stability is achieved. 

%Laura Note: sw_pden returns less dense values for temperatures below freezing.
%Sometimes if a lot of heat is extracted from the top cell, it falls below
%freezing so the density profile is stable since the top cell, which is very
%cold, is calculated to be less dense. To see this, run (as an example): 
% plot([-5:.01:4], sw_pden(28 .* ones(1, 901), [-5:.01:4], 0, 0) - 1000)

%If the top cell is below freezing, distribute the cold water over the
%mixed layer from the previous timestep since the heat loss will not
%stratify the ocean, so the new mixed layer will be at least as deep as the
%previous mixed layer. 

%If there are stratifying freshwater fluxes at the same time as cooling,
%will need to revisit this. At the moment (cooling until freezing) there
%are no freshwater fluxes
% if T(1) < sw_fp(S(1), 0) %Top cell is below freezing
%         
%     Sig = sw_pden(S,T,z,0);
% 
%     %Since the cold water is on top and salinity is homogenous within the 
%     %previous mixed layer, need to mix to at least the depth of at least
%     %the previous mld
%     mld_prev = mld;
%     mld = mld - 1; %Not sure why the -1 is needed, but otherwise it deepens the mixed layer from the previous step
%     pwp_MixMixedLayer; %Mix the cold water throughout the previous step's mixed layer
%     
%     
%     %Now proceed with static adjustment (the same as in the original code and
%     %included below) starting from the base of the
%     %previous mixed layer and seeing if further mixing is necessary
%     Sig = sw_pden(S,T,z,0);           % Compute potential density profile
%     
%     for i=mld_prev:nz-1			% Starting from the surface and working down
%         if Sig(i+1) >= Sig(i)       % Check for static instability 
%             mld=i;                  % If stable, set mixed layer depth at current depth
%             break
%         else
%         mld=i+1;                % If unstable, set mixed layer depth at next lowest depth
%         Sigm = mean(Sig(1:mld));  % Calculate the average density of the mixed layer
%         Sig(1:mld) = Sigm*ones(mld,1);  % Set the density of the mixed layer to the mean
%         nstin = nstin+1;                % Keep track of the activity
%         end
%     end
%     
%     if mld > mld_prev %The layer remains statically unstable, need to mix again
%         pwp_MixMixedLayer;
%     end
%     
% else %Surface cell is above freezing, complete the static adjustment as usual
%Now proceed with static adjustment
    Sig = sw_pden(S,T,z,0);           % Compute potential density profile
    for i=1:nz-1			% Starting from the surface and working down
        if Sig(i+1) >= Sig(i)       % Check for static instability 
            mld=i;                  % If stable, set mixed layer depth at current depth
            break
        else
        mld=i+1;                % If unstable, set mixed layer depth at next lowest depth
        Sigm = mean(Sig(1:mld));  % Calculate the average density of the mixed layer
        Sig(1:mld) = Sigm*ones(mld,1);  % Set the density of the mixed layer to the mean
        nstin = nstin+1;                % Keep track of the activity
        end
    end
    
    if mld > 1 
        pwp_MixMixedLayer; % Now mix all other properties %(temperature, salinity, density, and velocity)
    end  
% end

zmld_static(nn) = z(mld); %Added by Laura, depth to which static adjustment penetrated

    
