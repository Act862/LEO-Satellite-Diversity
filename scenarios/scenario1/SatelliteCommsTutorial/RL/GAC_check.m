%% ===================== FULL GAC SCRIPT WITH PRE-TRAINED NN =====================
clear; clc; close all;

%% ===================== USER INPUTS =====================
tlePath = "leoSatelliteConstellation.tle";  % TLE file
sTime   = datetime(2025,7,8,0,0,0);         % start time
eTime   = sTime + hours(3);                 % simulation duration
sampleTime = 60;                             % seconds

gsLat = 40.43139;
gsLon = -4.24806;
gsName = "Madrid Deep Space Communication Complex";

fc = 30e9;      % carrier frequency
Pt = 20;        % Tx power (dB)
Rb = 10;        % bit rate (bps)
c = physconst('LightSpeed');

%% ===================== CREATE SCENARIO =====================
sc = satelliteScenario(sTime,eTime,sampleTime);
sat = satellite(sc,tlePath);
gs  = groundStation(sc, Name=gsName, Latitude=gsLat, Longitude=gsLon);

T = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
NT = numel(T);

%% ===================== LATENCY =====================
[delay, ~] = latency(sat, gs);  % seconds
latencyRate = diff(delay,1,2)/sampleTime;  % rate of change

%% ===================== DOPPLER =====================
[fshift,timeD,dopplerInfo] = dopplershift(sat, gs, Frequency=fc);

%% ===================== AER AND SLANT RANGE =====================
[alt,elevation,slantRange] = aer(sat, gs);  % meters

% Interpolate slantRange to match NT if needed
numSats = size(slantRange,2);
oldTime = linspace(0,1,size(slantRange,1));
newTime = linspace(0,1,NT);
slantRangeInterp = zeros(NT,numSats);
for k = 1:numSats
    slantRangeInterp(:,k) = interp1(oldTime, slantRange(:,k), newTime, 'linear', 'extrap');
end
slantRange = slantRangeInterp;

% Fraction of visible satellites
frac_visible = sum(~isnan(slantRange),2)/numSats;

%% ===================== TIME-BASED FEATURES =====================
tod = hour(T)*3600 + minute(T)*60 + second(T);
doy = day(T,'dayofyear');

sin_tod = sin(2*pi*tod/86400); sin_tod = sin_tod(:);
cos_tod = cos(2*pi*tod/86400); cos_tod = cos_tod(:);
sin_doy = sin(2*pi*doy/365);   sin_doy = sin_doy(:);
cos_doy = cos(2*pi*doy/365);   cos_doy = cos_doy(:);

%% ===================== PREPARE NN INPUT =====================
Xnn = [sin_tod, cos_tod, sin_doy, cos_doy, frac_visible];  % NT x 5

%% ===================== LOAD PRE-TRAINED NN =====================
nnFile = "trainedMaxSlantNN.mat";  % your saved network
load(nnFile,'net');  % assumes variable 'net'

yPredMaxSlant = net(Xnn')';  % transpose to match feedforward input

%% ===================== FSPL AND SNR =====================
FSPL = 20*log10(slantRange) + 20*log10(fc) - 147.55;  % dB
SNR_ul = Pt - FSPL + 6.9 - 10*log10(c);               % simplified

% GAC computation: max SNR weighted by predicted slant
GAC = (SNR_ul(:,1).*yPredMaxSlant + ...
       SNR_ul(:,2).*yPredMaxSlant + ...
       SNR_ul(:,3).*yPredMaxSlant) ./ Rb;  % example

%% ===================== PLOT RESULTS =====================
figure;
plot(T, delay(1,:)*1000, 'b.-', 'DisplayName','Sat1 Latency'); hold on;
xlabel('Time'); ylabel('Latency (ms)'); grid on;
title('Satellite Latency');

figure;
plot(T, fshift(1,:)/1e3, 'r.-', 'DisplayName','Sat1 Doppler'); hold on;
xlabel('Time'); ylabel('Doppler Shift (kHz)'); grid on;
title('Satellite Doppler Shift');

figure;
plot(T, GAC, 'k.-', 'DisplayName','GAC'); hold on;
xlabel('Time'); ylabel('GAC (dB)'); grid on;
title('GAC vs Time');
legend('Location','best');

disp('GAC calculation completed.');
