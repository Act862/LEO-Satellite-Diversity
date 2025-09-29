%% ============================================================
%  Comparison of Selection Combining, MRC, and GAC
%  Author: (Your Name)
%  ============================================================
clear; clc; rng(1);

%% PARAMETERS
N_sats   = 5;        % number of satellites
N_trials = 1000;     % number of Monte Carlo trials

% Average SNRs per satellite (dB) -> represents link quality differences
avgSNRdB = [10 8 6 4 2];  

% Geometry: assign some random slant ranges (km)
slantRanges = [800, 1000, 1200, 1500, 2000];  
Rmax = max(slantRanges);

%% STORAGE
snr_SC  = zeros(N_trials,1);
snr_MRC = zeros(N_trials,1);
snr_GAC = zeros(N_trials,1);

%% SIMULATION
for t = 1:N_trials
    % --- Step 1: Generate random instantaneous SNRs (Rayleigh fading)
    instSNR_lin = zeros(1,N_sats);
    for k = 1:N_sats
        avgSNR_lin = 10^(avgSNRdB(k)/10);
        h = (randn+1i*randn)/sqrt(2);    % Rayleigh fading channel
        instSNR_lin(k) = avgSNR_lin * abs(h)^2;
    end

    % --- Step 2: Selection Combining (pick best satellite)
    snr_SC(t) = max(instSNR_lin);

    % --- Step 3: Maximal Ratio Combining (weight ∝ SNR)
    w_MRC = instSNR_lin ./ sum(instSNR_lin);  % normalized weights
    snr_MRC(t) = sum(instSNR_lin);            % equivalent post-combining SNR

    % --- Step 4: Geometry Assisted Combining
    geoFactor = (Rmax ./ slantRanges);        % geometry weights
    w_GAC = (instSNR_lin .* geoFactor) ./ sum(instSNR_lin .* geoFactor);
    snr_GAC(t) = sum(instSNR_lin .* geoFactor); % geometry-biased combining
end

%% RESULTS
meanSC  = mean(10*log10(snr_SC));
meanMRC = mean(10*log10(snr_MRC));
meanGAC = mean(10*log10(snr_GAC));

fprintf("Average post-combining SNR (dB):\n");
fprintf("  Selection Combining (SC): %.2f dB\n", meanSC);
fprintf("  Maximal Ratio Combining (MRC): %.2f dB\n", meanMRC);
fprintf("  Geometry-Assisted Combining (GAC): %.2f dB\n", meanGAC);

%% PLOTS
figure;
boxplot([10*log10(snr_SC), 10*log10(snr_MRC), 10*log10(snr_GAC)], ...
    'Labels',{'SC','MRC','GAC'});
ylabel('Post-Combining SNR (dB)');
title('Comparison of Combining Techniques');
grid on;
