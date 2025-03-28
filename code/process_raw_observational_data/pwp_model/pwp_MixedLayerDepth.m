%pwp_MixedLayerDepth PWP subroutine that finds the mixed layer depth index

D = find((Sig - Sig(1)) >=  forc.mldThreshold, 1);

% If the mixed layer criterion is not found, 

% set the mixed layer depth to 1 (the first grid point)

if isempty(D)
    D = 1; 
end

% Otherwise, choose first value greater than 10^-4 for mixed layer depth

mld = D;  % Note mld is an index, not a depth.                    

