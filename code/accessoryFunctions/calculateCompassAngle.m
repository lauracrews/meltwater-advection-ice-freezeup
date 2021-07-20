%Takes in zonal and meridional velocity vector components and calculates 
%the compass angle of that velocity (i.e. north = 0°, east = 90°). 
%Used in plot_frontTracking.m
function degN = calculateCompassAngle(u, v)

        theta = atand(abs(v/u));
        degN = 0;
        
        %Convert to compass coordinates
        if u < 0 & v < 0
            degN = 180 + (90-theta);
        elseif u < 0 & v > 0
            degN = 270 + theta;
        elseif u > 0 & v < 0 
            degN = 90 + theta;
        elseif u > 0 & v > 0
            degN = 90 - theta;
        elseif v == 0
            if u > 0; degN = 90;
            elseif u < 0; degN = 270; end
        elseif u == 0
            if v > 0; degN = 0;
            elseif v < 0; degN = 180;
            end              
        end
                   
end

%Alternative method (which I think works)
% deg = atan2d(v, u);
% degN = 90 - deg; 
% degN(degN < 0) = degN + 360;