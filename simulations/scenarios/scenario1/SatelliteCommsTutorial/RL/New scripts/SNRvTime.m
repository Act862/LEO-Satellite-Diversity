%% ============================================================
% SNR vs Time Comparison: SC, MRC, GAC
clear; clc; rng(10);

%% PARAMETERS
N_branches   = 4;                       % number of satellites/branches
N_time      = 200;                       % number of time samples
avgSNRdB    = [10 8 6 4];               % average SNR per branch (dB)
avgSNR_lin  = 10.^(avgSNRdB/10);

% Slant ranges for GAC (geometry-assisted weighting)
slantRanges = [800 1000 1200 1500];  % in km
Rmax = max(slantRanges);
geoFactor = Rmax ./ slantRanges;

%% STORAGE
SNR_SC_time  = zeros(N_time,1);
SNR_MRC_time = zeros(N_time,1);
SNR_GAC_time = zeros(N_time,1);
SNR_branches = zeros(N_time,N_branches);

%% SIMULATION LOOP
for t = 1:N_time
    % Generate Rician fading per branch
    K = 5; % moderate Rician K-factor
    h = sqrt(K/(K+1))*exp(1j*2*pi*rand(1,N_branches)) + ...
        sqrt(1/(K+1))*(randn(1,N_branches)+1j*randn(1,N_branches))/sqrt(2);
    
    % Instantaneous SNR per branch
    gamma = avgSNR_lin .* abs(h).^2;
    SNR_branches(t,:) = 10*log10(gamma);
    
    % --- SC
    gamma_SC = max(gamma);
    SNR_SC_time(t) = 10*log10(gamma_SC);
    
    % --- MRC
    gamma_MRC = sum(gamma);
    SNR_MRC_time(t) = 10*log10(gamma_MRC);
    
    % --- GAC
    gamma_GAC = sum(gamma .* geoFactor);
    SNR_GAC_time(t) = 10*log10(gamma_GAC);
end

%% PLOT
figure; hold on; grid on;
timeVec = 1:N_time;

% Per-branch SNR
plot(timeVec, SNR_branches,'--','Color',[0.7 0.7 0.7]);

% Combined SNR
plot(timeVec, SNR_SC_time,'r','LineWidth',2);
plot(timeVec, SNR_MRC_time,'b','LineWidth',2);
plot(timeVec, SNR_GAC_time,'g','LineWidth',2);

xlabel('Time index'); ylabel('SNR (dB)');
title('SNR vs Time: SC, MRC, GAC');
legend([arrayfun(@(x) sprintf('Branch %d',x),1:N_branches,'UniformOutput',false),...
        'SC','MRC','GAC']);