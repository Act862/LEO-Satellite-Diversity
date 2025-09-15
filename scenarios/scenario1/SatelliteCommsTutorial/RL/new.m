%% ===================== LEO Constellation + NN with Quantization =====================
clear; clc; close all;

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";  % TLE file
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
startTime = datetime(2020,1,11,14,0,0);
stopTime  = startTime + minutes(10);
sampleSec = 60;

% RF / link setup
c = physconst('LightSpeed');

%% ===================== SCENARIO =====================
sc = satelliteScenario(startTime, stopTime, sampleSec);
sats = satellite(sc, tlePath);
nSats = numel(sats);
gs = groundStation(sc, Name="GS", Latitude=gsLat, Longitude=gsLon, Altitude=0);

T = sc.StartTime : seconds(sc.SampleTime) : sc.StopTime;
NT = numel(T);

%% ===================== PER-SATELLITE SLANT RANGES & VISIBILITY =====================
wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gsLat, gsLon, 0);

slantAll_m   = NaN(NT,nSats);
visibleMask  = false(NT,nSats);
trueMaxSlant_m = zeros(NT,1);

for k = 1:nSats
    for ii = 1:NT
        [posECEF, ~] = states(sats(k), T(ii), "CoordinateFrame","ecef");
        slant = norm(posECEF - [gsX; gsY; gsZ]);
        slantAll_m(ii,k) = slant;
        if slant < 1e7 % assume visible if slant reasonable (approx)
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

%% ===================== FEATURE ENGINEERING =====================
tod = seconds(timeofday(T)) / (24*3600); % fractional time-of-day
doy = day(T, "dayofyear") / 366;        % fractional day-of-year
fracVis = sum(visibleMask,2) / max(1,nSats);

% Assemble dataset: [sin_tod, cos_tod, sin_doy, cos_doy, frac_visible, per-sat slants..., max_slant_m]
features = [sin(2*pi*tod(:)), cos(2*pi*tod(:)), sin(2*pi*doy(:)), cos(2*pi*doy(:)), fracVis(:), slantAll_m];
targets  = trueMaxSlant_m;

% Remove rows with NaN
validRows = all(isfinite(features),2) & isfinite(targets);
features = features(validRows,:);
targets  = targets(validRows);

%% ===================== TRAIN NEURAL NETWORK =====================
inputSize = size(features,2);
hiddenSize = 32;

net = feedforwardnet(hiddenSize,'trainlm'); % Levenberg-Marquardt
net.layers{1}.transferFcn = 'relu';
net.layers{2}.transferFcn = 'purelin';

net = configure(net, features', targets');
net.trainParam.epochs = 200;

[net, tr] = train(net, features', targets');

%% ===================== MANUAL QUANTIZATION =====================
% Simulate reduced-size network by reducing weight precision
for i = 1:numel(net.IW)
    net.IW{i} = single(round(net.IW{i}*127)/127);
end
for i = 1:numel(net.LW)
    net.LW{i} = single(round(net.LW{i}*127)/127);
end
for i = 1:numel(net.b)
    net.b{i}  = single(round(net.b{i}*127)/127);
end

%% ===================== PREDICTION =====================
predMaxSlant_m = net(features');

%% ===================== ERROR CALCULATION =====================
err_m = predMaxSlant_m' - targets;
MAE  = mean(abs(err_m));
MSE  = mean(err_m.^2);
RMSE = sqrt(MSE);

fprintf('Prediction Error Metrics:\n');
fprintf('MAE  = %.2f m\n', MAE);
fprintf('MSE  = %.2f m^2\n', MSE);
fprintf('RMSE = %.2f m\n', RMSE);

%% ===================== PLOT RESULTS =====================
figure;
plot(T(validRows), targets/1e3,'b.-','DisplayName','Truth max slant (km)'); hold on;
plot(T(validRows), predMaxSlant_m'/1e3,'r.-','DisplayName','Predicted max slant (km)');
xlabel('Time'); ylabel('Max Slant Range (km)'); grid on; legend('Location','best');
title('Neural Network Prediction vs Truth');

figure;
plot(T(validRows), err_m/1e3,'k.-');
xlabel('Time'); ylabel('Error (km)'); grid on;
title('Prediction Error');

disp('Simulation and prediction completed.');
