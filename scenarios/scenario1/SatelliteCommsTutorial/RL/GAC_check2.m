%% ===================== Full GAC Simulation + Prediction =====================
clear; clc; close all;

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
startTime = datetime(2025,7,8,0,0,0);
stopTime  = startTime + minutes(8);  % adjust for testing
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
nSats = numel(sats);

gs = groundStation(sc, Name="GS", Latitude=gsLat, Longitude=gsLon, Altitude=0);
gGs = gimbal(gs, "MountingAngles",[0;180;0],"MountingLocation",[0;0;-5]);
gsRx = receiver(gGs, Name="GS_RX", MountingLocation=[0;0;1], ...
    GainToNoiseTemperatureRatio=gsG_T_dB, RequiredEbNo=reqEbNo_dB);
gaussianAntenna(gsRx, "DishDiameter", gsDish_m);

satTx = cell(1,nSats);
gSats = cell(1,nSats);
links = cell(1,nSats);
for k = 1:nSats
    gSats{k} = gimbal(sats(k), "MountingLocation",[0;0;2]);
    satTx{k} = transmitter(gSats{k}, "Name","SatTx"+k, ...
        Frequency=txFrequencyHz, Power=txPowerW, BitRate=bitRate);
    gaussianAntenna(satTx{k}, "DishDiameter", satDish_m);
    pointAt(gSats{k}, gs);
    pointAt(gGs, sats(k));
    links{k} = link(satTx{k}, gsRx);
end

T = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
NT = numel(T);

%% ===================== TRUE SLANT RANGES & VISIBILITY =====================
wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gsLat, gsLon, 0);

trueMaxSlant_m = NaN(NT,1);
visibleMask    = false(NT,nSats);
slantAll_m     = NaN(NT,nSats);

for k = 1:nSats
    [posECEF, ~] = states(sats(k), T, "CoordinateFrame","ecef");
    for ii = 1:NT
        slantAll_m(ii,k) = norm(posECEF(:,ii) - [gsX; gsY; gsZ]);
        visibleMask(ii,k) = isfinite(slantAll_m(ii,k));  % simple visibility
    end
end

for ii = 1:NT
    visIdx = find(visibleMask(ii,:));
    if ~isempty(visIdx)
        trueMaxSlant_m(ii) = max(slantAll_m(ii,visIdx));
    end
end

%% ===================== SIMPLE PREDICTOR: SHIFTING BUFFER =====================
bufferLen = 3;
predMaxSlant_m = NaN(NT,1);
firstValid = find(~isnan(trueMaxSlant_m),1,'first');
predMaxSlant_m(1:bufferLen) = repmat(trueMaxSlant_m(firstValid), bufferLen,1);

for ii = bufferLen+1:NT
    predMaxSlant_m(ii) = max(trueMaxSlant_m(ii-bufferLen:ii-1), [], 'omitnan');
end

%% ===================== PREDICTION ERROR =====================
err_m = predMaxSlant_m - trueMaxSlant_m;
fprintf('Shifting Buffer Predictor Errors:\n');
fprintf('MAE  = %.2f m\n', mean(abs(err_m),'omitnan'));
fprintf('MSE  = %.2f m^2\n', mean(err_m.^2,'omitnan'));
fprintf('RMSE = %.2f m\n', sqrt(mean(err_m.^2,'omitnan')));

%% ===================== PLOTS =====================
validIdx = ~isnan(predMaxSlant_m);

% Max slant range prediction
figure;
plot(T(validIdx), trueMaxSlant_m(validIdx)/1e3,'b.-','DisplayName','Truth'); hold on;
plot(T(validIdx), predMaxSlant_m(validIdx)/1e3,'r.-','DisplayName','Shifting Buffer');
xlabel('Time'); ylabel('Max Slant Range (km)');
grid on; legend('Location','best'); title('Max Slant Range: Truth vs Predictor');

% Prediction error
figure;
plot(T(validIdx), err_m(validIdx)/1e3,'k.-'); grid on;
xlabel('Time'); ylabel('Error (km)'); title('Prediction Error of Shifting Buffer Predictor');

%% ===================== GAC, SC, MRC COMPARISON =====================
% Example: use first 3 satellites
satIds = [1, 2, 3]; 
numSatsAvailable = numel(satIds);

% SC: select max Eb/N0 among available
% MRC: sum SNR of all available
% GAC: use shifting predictor for slant-range weighting
fc = 30e9; % carrier freq
SNR_ul = zeros(NT,numSatsAvailable);

[alt,elevation,slantRangeAll] = aer(sats(satIds),gs);

for k = 1:numSatsAvailable
    % Free-space path loss
    FSPL = 20*log10(slantRangeAll(:,k)) + 20*log10(fc) - 147.55;
    SNR_ul(:,k) = 20 - 0.5 + 6.9 - FSPL;  % simplified link budget
end

% Selection combining (SC)
SC = max(SNR_ul,[],2);

% Max ratio combining (MRC)
MRC = sum(SNR_ul,2);

% GAC (slant-range weighted)
predSlant = predMaxSlant_m(1:NT);
GAC = zeros(NT,1);
for ii = 1:NT
    weights = zeros(numSatsAvailable,1);
    for k = 1:numSatsAvailable
        weights(k) = SNR_ul(ii,k)/slantRangeAll(ii,k);
    end
    GAC(ii) = sum(weights);
end

% Plot SC, MRC, GAC
figure;
plot(T,SC,'b.-'); hold on;
plot(T,MRC,'r.-');
plot(T,GAC,'g.-'); hold off;
xlabel('Time'); ylabel('Metric (arb. units)');
legend('SC','MRC','GAC'); title('SC vs MRC vs GAC');
grid on;
