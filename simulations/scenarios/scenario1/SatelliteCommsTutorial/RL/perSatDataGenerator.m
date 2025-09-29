%% ===================== LEO Constellation Dataset Generator =====================
clear; clc;

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";  % your TLE file
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
gsAlt     = 0;  % meters
startTime = datetime(2020,1,11,14,0,0);
stopTime  = startTime + minutes(10);
sampleSec = 60;

%% ===================== SCENARIO =====================
sc = satelliteScenario(startTime, stopTime, sampleSec);
sats = satellite(sc, tlePath);
nSats = numel(sats);
gs = groundStation(sc, Name="GS", Latitude=gsLat, Longitude=gsLon, Altitude=gsAlt);

%% Single GS receiver
gGs = gimbal(gs, "MountingAngles",[0;180;0],"MountingLocation",[0;0;-5]);
gsRx = receiver(gGs, Name="GS_RX", MountingLocation=[0;0;1], ...
    GainToNoiseTemperatureRatio=20, RequiredEbNo=1);

%% Scenario times
T = sc.StartTime : seconds(sc.SampleTime) : sc.StopTime;
NT = numel(T);

%% ===================== TRUE SLANT RANGES & VISIBILITY =====================
wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gsLat, gsLon, gsAlt);

slantAll_m = NaN(NT, nSats);
visibleMask = false(NT, nSats);
trueMaxSlant_m = zeros(NT,1);

for k = 1:nSats
    acc = access(gs, sats(k));
    ints = accessIntervals(acc); %#ok<NASGU>
    for ii = 1:NT
        [posECEF, ~] = states(sats(k), T(ii), "CoordinateFrame","ecef");
        slant = norm(posECEF - [gsX; gsY; gsZ]);
        slantAll_m(ii,k) = slant;
        % Visibility proxy: if slant is finite
        visibleMask(ii,k) = isfinite(slant);
    end
end

for ii = 1:NT
    visIdx = find(visibleMask(ii,:));
    if ~isempty(visIdx)
        trueMaxSlant_m(ii) = max(slantAll_m(ii,visIdx));
    else
        trueMaxSlant_m(ii) = NaN;
    end
end

%% ===================== FEATURE ENGINEERING =====================
tod = seconds(timeofday(T)) / (24*3600); % fractional time-of-day
doy = day(T, "dayofyear") / 366;        % fractional day-of-year
fracVis = sum(visibleMask,2) / max(1,nSats);

% Base features
baseFeatures = [sin(2*pi*tod(:)), cos(2*pi*tod(:)), ...
                sin(2*pi*doy(:)), cos(2*pi*doy(:)), fracVis(:)];

% Add per-satellite slant ranges (NaN where not visible)
perSatSlants = slantAll_m; % NT x nSats

% Combine all features
dataset = [baseFeatures, perSatSlants, trueMaxSlant_m];

%% ===================== SAVE TO CSV =====================
% Column headers
colHeaders = ["sin_tod","cos_tod","sin_doy","cos_doy","frac_visible"];
for k = 1:nSats
    colHeaders = [colHeaders, "slant_sat"+k];
end
colHeaders = [colHeaders, "max_slant_m"];

% Convert to table
datasetTable = array2table(dataset, 'VariableNames', colHeaders);

% Save CSV
csvFileName = "leo_dataset.csv";
writetable(datasetTable, csvFileName);

disp("Dataset saved to: " + csvFileName);
