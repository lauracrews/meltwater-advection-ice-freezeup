%pwp_MixMixedLayer: Subroutine that mixes up the mixed layer



Tm=mean(T(1:mld));                 % Calculate the mean temperature of the mixed layer

T(1:mld)=Tm*ones(mld,1);           % Set the temperature of the mixed layer to the mean temperature

Sm=mean(S(1:mld));                 % Calculate the mean salinity of the mixed layer

S(1:mld)=Sm*ones(mld,1);           % Set the salinity of the mixed layer to the mean salinity


%Laura note: this was the original code. I don't think it's correct to just
%average the mixed layer density due to effect of cabelling. 
% Sigm=mean(Sig(1:mld));             % Calculate the mean density of the mixed layer
% Sig(1:mld)=Sigm*ones(mld,1);       % Set the density of the mixed layer to the mean density

%Instead, recalculate density from the mixed layer averaged temp and
%salinity
Sig = sw_pden(S,T,z,0);

UVm=mean(UV(1:mld,:),1);           % Calculate the mean velocities of the mixed layer

UV(1:mld,1)=UVm(1)*ones(mld,1);    % Set the u velocity to the mean u velocity

UV(1:mld,2)=UVm(2)*ones(mld,1);    % Set the v velocity to the mean v velocity

