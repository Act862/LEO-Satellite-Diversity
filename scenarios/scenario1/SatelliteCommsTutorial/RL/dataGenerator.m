%% ===================== Generate Training Data from TLE =====================
clear; clc; close all;

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";  % your TLE file
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
gsAlt     = 0;
startTime = datetime(2020,1,11,14,0,0);
stopTime  = startTime + days(1);   % simulation span
sampleSec = 60;                         % sample interval in seconds

%% ===================== SCENARIO SETUP =====================
sc = satelliteScenario(startTime, stopTime, sampleSec);
sats = satellite(sc, tlePath);
nSats = numel(sats);

gs = groundStation(sc, Name="GS", Latitude=gsLat, Longitude=gsLon, Altitude=gsAlt);

% Precompute time vector
T = sc.StartTime : seconds(sc.SampleTime) : sc.StopTime;
NT = numel(T);

%% ===================== FEATURE & TARGET ARRAYS =====================
% Features: time-of-day (sin/cos), day-of-year (sin/cos), #visible satellites
features = zeros(NT,5);
targets  = zeros(NT,1);  % max slant range to visible satellites

visibleMask = false(NT,nSats);
slantAll_m  = NaN(NT,nSats);

wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gsLat, gsLon, gsAlt);

%% ===================== COMPUTE SLANT RANGES & VISIBILITY =====================
for k = 1:nSats
    % Access object (optional)
    acc = access(gs, sats(k));
    
    for ii = 1:NT
        [posECEF, ~] = states(sats(k), T(ii), "CoordinateFrame","ecef");
        slantAll_m(ii,k) = norm(posECEF - [gsX; gsY; gsZ]);
        visibleMask(ii,k) = ~isnan(slantAll_m(ii,k)); % mark visible if slant is valid
    end
end

%% ===================== GENERATE FEATURES =====================
tod = seconds(timeofday(T)) / (24*3600); % fractional time-of-day
doy = day(T, "dayofyear") / 366;        % fractional day-of-year (approx)
fracVis = sum(visibleMask,2) / max(1,nSats);

features(:,1) = sin(2*pi*tod(:));
features(:,2) = cos(2*pi*tod(:));
features(:,3) = sin(2*pi*doy(:));
features(:,4) = cos(2*pi*doy(:));
features(:,5) = fracVis(:);

%% ===================== GENERATE TARGET =====================
% Target: maximum slant range among visible satellites
for ii = 1:NT
    visIdx = find(visibleMask(ii,:));
    if ~isempty(visIdx)
        targets(ii) = max(slantAll_m(ii,visIdx));
    else
        targets(ii) = NaN;
    end
end

%% ===================== SAVE TO CSV =====================
data = [features, targets];
header = {'sin_tod','cos_tod','sin_doy','cos_doy','frac_visible','max_slant_m'};

% Create table
TBL = array2table(data,'VariableNames',header);

% Save CSV
csvFileName = "satellite_training_data.csv";
writetable(TBL, csvFileName);

disp("Training data saved to " + csvFileName);
