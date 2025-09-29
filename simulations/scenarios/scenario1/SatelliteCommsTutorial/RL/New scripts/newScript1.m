%% ===================== Ephemeris-Free Dataset for Max-Slant Prediction =====================
% This script uses satelliteScenario geometry only to produce the LABEL
% (true max slant). Inputs mimic *receiver-observable* features:
%  - Time-of-day / day-of-year (cyclic)
%  - Per-epoch SNR statistics over detected signals (no sat IDs)
%  - Per-epoch Doppler statistics (magnitude only, no sat IDs)
%  - Short-history dynamics (slope of top-1 SNR, handover rate, Doppler rate)
%
% You will train a model to predict max_slant_m from these *ephemeris-free* inputs.
%
% NOTE: In reality, SNR & Doppler are measured by the receiver. Here we
% synthesize them from geometry (for simulation), but we DO NOT feed true slants as inputs.

clear; clc; rng(42);  % reproducibility

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";  % your TLE file
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
gsAlt     = 0;                      % meters
startTime = datetime(2020,1,11,14,0,0);
stopTime  = startTime + minutes(10);
sampleSec = 10;                     % finer sampling helps dynamics

% "RF-ish" parameters for synthetic observables
fc_Hz          = 2.0e9;            % carrier freq (set to your band of interest)
c              = 299792458;        % speed of light (m/s)
lambda         = c / fc_Hz;
SNRref_dB      = 15;               % nominal SNR at r_ref
r_ref_m        = 1.0e6;            % reference range for SNRref (1,000 km)
shadowStd_dB   = 4;                % lognormal shadowing (dB)
snrDetThresh_dB= -5;               % detection threshold (dB)

% Short-history (sliding window) for dynamics features
winSec   = 60;                                    % history window (seconds)
W        = max(2, round(winSec / sampleSec));     % at least 2 samples

% Doppler rate uses simple finite differences; no smoothing by default
topK     = 3;    % for "top-K" summary features (SNR & Doppler)

%% ===================== SCENARIO =====================
sc   = satelliteScenario(startTime, stopTime, sampleSec);
sats = satellite(sc, tlePath);
nSats = numel(sats);
gs   = groundStation(sc, Name="GS", Latitude=gsLat, Longitude=gsLon, Altitude=gsAlt);

T  = sc.StartTime : seconds(sc.SampleTime) : sc.StopTime;
NT = numel(T);

%% ===================== GEOMETRY for LABEL + VISIBILITY =====================
% We use geometry only to compute:
%   - slant ranges (for label)
%   - elevation (for "visible" mask)
%   - Doppler (derived from LOS and relative velocity) -> becomes an observable
%
% No per-satellite slants are exported as INPUT features.

wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gsLat, gsLon, gsAlt);

slantAll_m   = NaN(NT, nSats);
elevAll_deg  = NaN(NT, nSats);
fdAll_Hz     = NaN(NT, nSats);          % signed Doppler (Hz)
visibleMask  = false(NT, nSats);
trueMaxSlant_m = NaN(NT,1);

for k = 1:nSats
    for ii = 1:NT
        % Satellite position & velocity in ECEF
        [posECEF, velECEF] = states(sats(k), T(ii), "CoordinateFrame","ecef"); % 3x1 each
        satPos = posECEF(:).';
        satVel = velECEF(:).';

        % Slant range
        gsPos  = [gsX, gsY, gsZ];
        rVec   = satPos - gsPos;
        r_m    = norm(rVec);
        slantAll_m(ii,k) = r_m;

        % Elevation (via ENU)
        [gsLatDeg, gsLonDeg, gsAltM] = deal(gsLat, gsLon, gsAlt);
        [e, n, u] = ecef2enu(satPos(1),satPos(2),satPos(3), gsLatDeg, gsLonDeg, gsAltM, wgs);
        % Elevation angle in degrees
        horiz  = hypot(e, n);
        elev   = atan2d(u, horiz);
        elevAll_deg(ii,k) = elev;

        % Visibility (simple horizon check)
        visibleMask(ii,k) = elev > 0;

        % Doppler (LOS projection of relative velocity)
        if isfinite(r_m) && r_m > 0
            u_los = rVec / r_m;                 % LOS unit vector (ECEF)
            v_rel = satVel;                     % GS is stationary in ECEF
            fd    = - dot(v_rel, u_los) / lambda;   % Hz
            fdAll_Hz(ii,k) = fd;
        end
    end
end

% Label: true maximum slant among *visible* satellites
for ii = 1:NT
    visIdx = find(visibleMask(ii,:));
    if ~isempty(visIdx)
        trueMaxSlant_m(ii) = max(slantAll_m(ii, visIdx));
    else
        trueMaxSlant_m(ii) = NaN;
    end
end

%% ===================== SYNTHETIC "MEASURED" SNR (ephemeris-free at inference) =====================
% You will NOT feed slant ranges to the model.
% Here SNR is generated from geometry solely to emulate what a receiver would measure.
% Path loss ~ 1/r^2  =>  -20*log10(r/r_ref), plus lognormal shadowing.

SNR_dB = SNRref_dB - 20*log10(max(slantAll_m, 1) ./ r_ref_m) + shadowStd_dB*randn(NT, nSats);

% Consider only visible satellites for detection
detMask = visibleMask & (SNR_dB >= snrDetThresh_dB);

%% ===================== BUILD RECEIVER-REALISTIC FEATURES (NO SAT IDs) =====================
% Per-epoch stats over detected links.
% Also build short-history dynamics of the *strongest* detected link.

% Preallocate feature arrays
n_detected           = zeros(NT,1);
snr_top1             = NaN(NT,1);
snr_top3_mean        = NaN(NT,1);
snr_median           = NaN(NT,1);
snr_var              = NaN(NT,1);

fd_abs_top1          = NaN(NT,1);
fd_abs_top3_mean     = NaN(NT,1);
fd_abs_std           = NaN(NT,1);

% For dynamics
snr_top1_slope       = NaN(NT,1);   % d(SNR_top1)/dt over window
dfd_dt_top1          = NaN(NT,1);   % Doppler rate for top1 over window (Hz/s)
handover_rate_permin = NaN(NT,1);   % number of top1 changes per minute over window

% Track index of strongest link each epoch (for handover rate)
top1_idx = NaN(NT,1);

for ii = 1:NT
    det = detMask(ii,:);
    if any(det)
        snrs = SNR_dB(ii, det);
        fds  = fdAll_Hz(ii, det);

        % Sort by SNR descending
        [snrs_sorted, order] = sort(snrs, 'descend', 'MissingPlacement','last');
        fds_sorted           = fds(order);

        kEff = min(topK, numel(snrs_sorted));

        n_detected(ii)     = numel(snrs_sorted);
        snr_top1(ii)       = snrs_sorted(1);
        snr_top3_mean(ii)  = mean(snrs_sorted(1:kEff), 'omitnan');
        snr_median(ii)     = median(snrs_sorted, 'omitnan');
        snr_var(ii)        = var(snrs_sorted, 0, 'omitnan');

        fd_abs_sorted      = abs(fds_sorted);
        fd_abs_top1(ii)    = fd_abs_sorted(1);
        fd_abs_top3_mean(ii)= mean(fd_abs_sorted(1:kEff), 'omitnan');
        fd_abs_std(ii)     = std(fd_abs_sorted, 0, 'omitnan');

        % Store the absolute sat index of top1 (for handover rate)
        absIdx = find(det);              % absolute sat indices visible/detected
        top1_idx(ii) = absIdx(order(1));
    else
        n_detected(ii) = 0;
    end
end

% Short-history dynamics (slopes & handover rate)
dt = sampleSec;
for ii = 1:NT
    i0 = max(1, ii - W + 1);
    winLen = ii - i0 + 1;
    if winLen >= 2
        % SNR top1 slope (simple finite diff across window)
        if all(isfinite(snr_top1(i0:ii)))
            snr_top1_slope(ii) = (snr_top1(ii) - snr_top1(i0)) / ((winLen-1)*dt); % dB/s
        end

        % Doppler rate for top1:
        % align Doppler of the *current* top1 when available; fallback: diff of fd_abs_top1
        if all(isfinite(fd_abs_top1(i0:ii)))
            dfd = fd_abs_top1(ii) - fd_abs_top1(i0);
            dfd_dt_top1(ii) = dfd / ((winLen-1)*dt); % Hz/s
        end

        % Handover rate: count changes in top1_idx within window
        idxSeq = top1_idx(i0:ii);
        changes = sum(abs(diff(idxSeq(~isnan(idxSeq)))) > 0);
        winMin  = (winLen*dt)/60;
        if winMin > 0
            handover_rate_permin(ii) = changes / winMin;
        end
    end
end

%% ===================== TIME
