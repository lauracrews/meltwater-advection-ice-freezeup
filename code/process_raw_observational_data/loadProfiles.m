%Loads the datastructure of all profiles in the study. Fills missing
%surface data and eliminates profiles with density inversions
function profiles = loadProfiles
    useFreshwaterForcing = true;
    
    cd([userpath, '/meltwaterAdvection/data/'])
    if useFreshwaterForcing
        load('profiles_withFreshwater.mat')
    else
        load('profiles_noFreshwater.mat')
    end
    cd([userpath, '/meltwaterAdvection/'])
    
    %If surface data is missing, fill with the shallowest available measurement
    badProfs = [];
    for k = 1:length(profiles.lats)
        firstTemp = find(~isnan(profiles.CT(:, k)), 1, 'first');
        profiles.CT(1:firstTemp, k) = profiles.CT(firstTemp, k);

        firstSalt = find(~isnan(profiles.SA(:, k)), 1, 'first');
        profiles.SA(1:firstSalt, k) = profiles.SA(firstSalt, k);

        firstDense = find(~isnan(profiles.sigthes(:, k)), 1, 'first');
        profiles.sigthes(1:firstDense, k) = profiles.sigthes(firstDense, k);

        %Eliminate if profile has density inversions > 0.05 kg/m^3 in the upper 40 m 
        %(since we aren't using data below 40 m) 
        deltaDensity = diff(profiles.sigthes(:, k));
        densityInversions = find(deltaDensity < -0.05);
        if min(densityInversions) < 40;  badProfs = [badProfs; k]; end
    end
    
    goodProfsMask = ~ismember(1:length(profiles.lats), badProfs);
    profiles.qualFlag = goodProfsMask;
end