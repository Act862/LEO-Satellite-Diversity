%% ==================== Simple Combining Simulation ====================
clear; clc; rng(1);

% Parameters
N_sats   = 5;        % number of satellites
N_trials = 1e4;      % Monte Carlo trials
SNR0_dB  = 10;       % baseline SNR per link (dB)
SNR0     = 10^(SNR0_dB/10);

% Geometry: assign random slant ranges (km)
slantRanges = randi([500 2000], [1, N_sats]); % fixed for all trials
maxSlant = max(slantRanges);

% Storage
SNR_SC  = zeros(N_trials,1);
SNR_MRC = zeros(N_trials,1);
SNR_GAC = zeros(N_trials,1);

for t = 1:N_trials
    % Random Rayleigh fading (per trial, per sat)
    h = (randn(1,N_sats) + 1i*randn(1,N_sats))/sqrt(2);
    
    % Per-satellite SNR
    SNR_i = SNR0 * abs(h).^2 ./ (slantRanges.^2);
    
    % --- Selection Combining ---
    SNR_SC(t) = max(SNR_i);
    
    % --- Maximal Ratio Combining ---
    SNR_MRC(t) = sum(SNR_i);
    
    % --- Geometry-Assisted Combining ---
    weights = maxSlant ./ slantRanges; % geometry weights
    SNR_GAC(t) = sum(SNR_i .* weights);
end

%% ==================== Results ====================
avgSNR_SC  = mean(SNR_SC);
avgSNR_MRC = mean(SNR_MRC);
avgSNR_GAC = mean(SNR_GAC);

fprintf("Average SNR (dB):\n");
fprintf("  SC  = %.2f dB\n", 10*log10(avgSNR_SC));
fprintf("  MRC = %.2f dB\n", 10*log10(avgSNR_MRC));
fprintf("  GAC = %.2f dB\n", 10*log10(avgSNR_GAC));

%% ==================== Plot CDFs ====================
edges = linspace(0, 25, 200);

figure; hold on; grid on;
cdfplot(10*log10(SNR_SC));  
cdfplot(10*log10(SNR_MRC));
cdfplot(10*log10(SNR_GAC));
legend("Selection Combining","MRC","GAC","Location","southeast");
xlabel("Output SNR (dB)"); ylabel("CDF");
title("SC vs MRC vs GAC");
