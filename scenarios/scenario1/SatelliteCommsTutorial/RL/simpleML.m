%% ===================== LEO Constellation + Simple Predictor =====================
clear; clc; close all;

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";  % your TLE file
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
startTime = datetime(2020,1,11,14,0,0);
stopTime  = startTime + minutes(10);
sampleSec = 60;

% RF / link setup
txFrequencyHz = 20e9;
txPowerW      = 15;
satDish_m     = 0.5;
gsDish_m      = 2.0;
gsG_T_dB      = 20;
reqEbNo_dB    = 1;
bitRate       = 20e6;
c             = physconst('LightSpeed');

%% ===================== SCENARIO =====================
sc = satelliteScenario(startTime, stopTime, sampleSec);
sats = satellite(sc, tlePath);
gs = groundStation(sc, Name="GS", Latitude=gsLat, Longitude=gsLon, Altitude=0);

% Single GS receiver
gGs = gimbal(gs, "MountingAngles",[0;180;0],"MountingLocation",[0;0;-5]);
gsRx = receiver(gGs, Name="GS_RX", MountingLocation=[0;0;1], ...
    GainToNoiseTemperatureRatio=gsG_T_dB, RequiredEbNo=reqEbNo_dB);
gaussianAntenna(gsRx, "DishDiameter", gsDish_m);

% Satellite TX and gimbals
nSats = numel(sats);
satTx = cell(1,nSats);
gSats = cell(1,nSats);
links = cell(1,nSats);
for k = 1:nSats
    gSats{k} = gimbal(sats(k), "MountingLocation",[0;0;2]);
    satTx{k} = transmitter(gSats{k}, "Name","SatTx"+k, ...
        Frequency=txFrequencyHz, Power=txPowerW, BitRate=bitRate);
    gaussianAntenna(satTx{k}, "DishDiameter", satDish_m);
    % point both ways
    pointAt(gSats{k}, gs);
    pointAt(gGs, sats(k));
    links{k} = link(satTx{k}, gsRx);
end

T = sc.StartTime : seconds(sc.SampleTime) : sc.StopTime;
NT = numel(T);

%% ===================== TRUE SLANT RANGES & VISIBILITY =====================
wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gsLat, gsLon, 0);

trueMaxSlant_m = zeros(NT,1);
visibleMask    = false(NT,nSats);
slantAll_m     = NaN(NT,nSats);

for k = 1:nSats
    [ebno_k, t_k] = ebno(links{k});
    [~, ia, ib] = intersect(T,t_k);
    
    % Access visibility
    acc = access(gs, sats(k));
    ints = accessIntervals(acc);
    
    for ii = 1:NT
        [posECEF, ~] = states(sats(k), T(ii), "CoordinateFrame","ecef");
        slant = norm(posECEF - [gsX; gsY; gsZ]);
        slantAll_m(ii,k) = slant;
        if isfinite(ebno_k(min(end,ii)))
            visibleMask(ii,k) = true;
        end
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

%% ===================== SIMPLE PREDICTOR: SHIFTING BUFFER =====================
bufferLen = 3; 
predMaxSlant_m = NaN(NT,1);

% Fill initial values (first bufferLen steps) with first valid truth
firstValid = find(~isnan(trueMaxSlant_m), 1, 'first');
predMaxSlant_m(1:bufferLen) = repmat(trueMaxSlant_m(firstValid), bufferLen, 1);

% Apply shifting window predictor
for ii = bufferLen+1:NT
    buf = trueMaxSlant_m(ii-bufferLen:ii-1);
    if all(isnan(buf))
        predMaxSlant_m(ii) = predMaxSlant_m(ii-1); % hold last value
    else
        predMaxSlant_m(ii) = max(buf, [], 'omitnan');
    end
end

%% ===================== PREDICTION ERROR =====================
err_m = predMaxSlant_m - trueMaxSlant_m;

%% ===================== PLOT RESULTS =====================
validIdx = ~isnan(predMaxSlant_m);
figure;
plot(T(validIdx), trueMaxSlant_m(validIdx)/1e3, '.', 'DisplayName','Truth');
hold on;
plot(T(validIdx), predMaxSlant_m(validIdx)/1e3, '-', 'DisplayName','Shifting Window');
xlabel('Time'); ylabel('Max Slant Range (km)');
grid on; legend('Location','best');
title('Max Slant Range: Truth vs Simple Predictor');

figure;
plot(T(validIdx), err_m(validIdx)/1e3, 'k.-'); 
grid on; xlabel('Time'); ylabel('Error (km)');
title('Prediction Error of Shifting Buffer Predictor');

disp('Simulation completed.');
