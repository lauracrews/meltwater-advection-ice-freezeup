%Calculates the average heat content in the upper ocea at modeled freeze up for profiles inside
%and outside the meltwater - these values are used in the discussion
%section of the paper

%% Differentiate profiles inside and outside the meltwater, average heat content
close all; 
clearvars -except AMSR2 profiles wvdata metData

defineSODAconstants;
integrationDepth = 40; %Integrate to a constant depth of 40 m

%Load data structure of all profiles
% profiles = loadProfiles;

%Geographic bounds for profiles on the southern Seaglider transect
minlat = 73 + 25/60; maxlat = 74.33; minlon = -148; maxlon = -146;

startTime = datenum('sept 26 2018') + .5;
endTime = datenum('oct 15 2018'); 

%Identify profiles in the geographic range of the Seaglider transect
inRegionMask = zeros(size(profiles.times));    
inRegionMask(profiles.lats >= minlat & profiles.lats <= maxlat ...
    & profiles.lons >= minlon & profiles.lons <= maxlon) = 1;

%Indentify profiles in the time range
inTimeMask = zeros(size(profiles.times));
inTimeMask(profiles.times >= startTime & profiles.times < endTime) = 1;

%Average surface salinity in the upper 5 m, then identify profiles
%inside/outside meltwater using a salinity threshold of 26 g/kg
surfaceSalinity = nanmean(profiles.SA(1:5, :));
inMeltwater = zeros(size(profiles.times));
inMeltwater(surfaceSalinity < 26) = 1;

%Identify the profiles in this geographic region in this time period
meltwaterProfNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1 & inMeltwater == 1);
noMeltwaterProfNums = find(profiles.qualFlag == 1 & inTimeMask == 1 & inRegionMask == 1 & inMeltwater == 0);

%Calculate heat content in PWP model output at freeze up for profiles inside the meltwater
meltwaterPWPupperHeatContent_freezeTime = nan .* ones(size(meltwaterProfNums));
for i = 1:length(meltwaterProfNums)
    profNum = meltwaterProfNums(i); 
    densProfPWP = gsw_rho(profiles.pwpFreezeSAprof(:, profNum), profiles.pwpFreezeCTprof(:, profNum), 0);
    [upperHeatContent_pwp, ~] = calculateHeatContent(profiles.pwpFreezeCTprof(:, profNum), profiles.pwpFreezeSAprof(:, profNum), densProfPWP, profiles.pwpz, profiles.lats(profNum),...
        1, integrationDepth);

    meltwaterPWPupperHeatContent_freezeTime(i) = upperHeatContent_pwp;    
end

%Calculate heat content in PWP model output at freeze up for profiles outside the meltwater
noMeltwaterPWPupperHeatContent_freezeTime = nan .* ones(size(meltwaterProfNums));
for i = 1:length(noMeltwaterProfNums)
    profNum = noMeltwaterProfNums(i); 
    densProfPWP = gsw_rho(profiles.pwpFreezeSAprof(:, profNum), profiles.pwpFreezeCTprof(:, profNum), 0);
    [upperHeatContent_pwp, ~] = calculateHeatContent(profiles.pwpFreezeCTprof(:, profNum), profiles.pwpFreezeSAprof(:, profNum), densProfPWP, profiles.pwpz, profiles.lats(profNum),...
        1, integrationDepth);

    noMeltwaterPWPupperHeatContent_freezeTime(i) = upperHeatContent_pwp;    
end

%Average heat content at freeze up time
mean(meltwaterPWPupperHeatContent_freezeTime)
mean(noMeltwaterPWPupperHeatContent_freezeTime)

%% Plot average profiles inside and outside the meltwater

%Axis limits for plotting
xlimits_temp = [-1.7, 0];
xlimits_salt = [25, 30];
ylimits_profiles = [0, 50];
linstl = '-';

%Temperature profiles
subplot(1, 2, 1) 
plot(nanmean(profiles.pwpFreezeCTprof(:, meltwaterProfNums), 2), profiles.pwpz, 'linestyle', linstl, 'LineWidth', 2, 'color', blue);
hold on
plot(nanmean(profiles.pwpFreezeCTprof(:, noMeltwaterProfNums), 2), profiles.pwpz, 'linestyle', linstl, 'LineWidth', 2, 'color', purple);

freezingPtProf = gsw_CT_freezing(nanmean(profiles.pwpFreezeSAprof(:, meltwaterProfNums), 2), gsw_p_from_z(-profiles.pwpz, 74));
plot(freezingPtProf, profiles.pwpz, 'k', 'linestyle', '--', 'linewidth', 1)

freezingPtProf = gsw_CT_freezing(nanmean(profiles.pwpFreezeSAprof(:, noMeltwaterProfNums), 2), gsw_p_from_z(-profiles.pwpz, 74));
plot(freezingPtProf, profiles.pwpz, 'k', 'linestyle', '--', 'linewidth', 1)

text(-1.6, 25, 'Freezing temperature', 'Rotation', 90, 'color', 'k', 'fontsize', 12)
xlim(xlimits_temp);
ylim(ylimits_profiles)
xlabel(['Conservative temperature (', sprintf(char(176)), 'C)'], 'fontsize', 12, 'FontWeight', 'bold')
ylabel('Depth (m)', 'fontsize', 12, 'FontWeight', 'bold')
set(gca, 'ydir', 'reverse', 'fontsize', 12)
grid on

%Salinity profiles
subplot(1, 2, 2)
plot(nanmean(profiles.pwpFreezeSAprof(:, meltwaterProfNums), 2), profiles.pwpz, 'linestyle', linstl, 'LineWidth', 2, 'color', blue);
hold on
plot(nanmean(profiles.pwpFreezeSAprof(:, noMeltwaterProfNums), 2), profiles.pwpz, 'linestyle', linstl, 'LineWidth', 2, 'color', purple);
xlim(xlimits_salt);
ylim(ylimits_profiles)
grid on
set(gca, 'ydir', 'reverse', 'fontsize', 12)
xlabel('Absolute salinity (g/kg)', 'fontsize', 12, 'FontWeight', 'bold')
ylabel('Depth (m)', 'fontsize', 12, 'FontWeight', 'bold')



