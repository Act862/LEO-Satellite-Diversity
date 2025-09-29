%% ============================================================
% SC vs MRC vs GAC - BER and Outage Probability
clear; clc; rng(2);

N_branches = 4;               % number of satellites/branches
N_trials   = 1e5;             % Monte Carlo trials
avgSNRdB  = [10 8 6 4];       % avg SNR per branch
avgSNR_lin = 10.^(avgSNRdB/10);

% Slant ranges for geometry-assisted MRC (km)
slantRanges = [800 1000 1200 1500];
Rmax = max(slantRanges);
geoFactor = Rmax ./ slantRanges;

% SNR threshold for outage
SNR_thresh_dB = 5;
SNR_thresh_lin = 10^(SNR_thresh_dB/10);

% Storage
SNR_SC  = zeros(N_trials,1);
SNR_MRC = zeros(N_trials,1);
SNR_GAC = zeros(N_trials,1);
BER_SC = zeros(N_trials,1);
BER_MRC = zeros(N_trials,1);
BER_GAC = zeros(N_trials,1);

for t = 1:N_trials
    % Generate Rayleigh fading per branch
    h = (randn(1,N_branches)+1j*randn(1,N_branches))/sqrt(2);
    
    % Instantaneous SNR per branch
    gamma = avgSNR_lin .* abs(h).^2;
    
    % --- SC
    gamma_SC = max(gamma);
    SNR_SC(t) = gamma_SC;
    BER_SC(t)  = 0.5*erfc(sqrt(gamma_SC));  % BPSK
    
    % --- MRC
    gamma_MRC = sum(gamma);
    SNR_MRC(t) = gamma_MRC;
    BER_MRC(t) = 0.5*erfc(sqrt(gamma_MRC));
    
    % --- GAC
    gamma_GAC = sum(gamma .* geoFactor);
    SNR_GAC(t) = gamma_GAC;
    BER_GAC(t) = 0.5*erfc(sqrt(gamma_GAC));
end

%% Outage Probability
outage_SC  = mean(SNR_SC < SNR_thresh_lin);
outage_MRC = mean(SNR_MRC < SNR_thresh_lin);
outage_GAC = mean(SNR_GAC < SNR_thresh_lin);

fprintf('Outage probability @ %d dB:\n',SNR_thresh_dB)
fprintf('SC  : %.4f\n',outage_SC)
fprintf('MRC : %.4f\n',outage_MRC)
fprintf('GAC : %.4f\n',outage_GAC)

%% Plots
figure;
semilogy(avgSNRdB(1),BER_SC(1),'ro'); hold on;
histogram(10*log10(SNR_SC),50,'Normalization','probability'); grid on;
xlabel('Post-Combining SNR (dB)'); ylabel('Probability / PDF');
title('Post-Combining SNR Distribution');

figure;
semilogy(sort(BER_SC),'r'); hold on;
semilogy(sort(BER_MRC),'b');
semilogy(sort(BER_GAC),'g');
grid on; ylabel('BER'); xlabel('Trial index');
legend('SC','MRC','GAC'); title('BER across trials');

% Outage vs SNR (smooth)
figure;
SNRdB_vec = 0:1:20;
outage_SC_vec = arrayfun(@(snr) mean(SNR_SC < 10^(snr/10)), SNRdB_vec);
outage_MRC_vec = arrayfun(@(snr) mean(SNR_MRC < 10^(snr/10)), SNRdB_vec);
outage_GAC_vec = arrayfun(@(snr) mean(SNR_GAC < 10^(snr/10)), SNRdB_vec);
semilogy(SNRdB_vec,outage_SC_vec,'r','LineWidth',2); hold on;
semilogy(SNRdB_vec,outage_MRC_vec,'b','LineWidth',2);
semilogy(SNRdB_vec,outage_GAC_vec,'g','LineWidth',2); grid on;
xlabel('SNR Threshold (dB)'); ylabel('Outage Probability'); 
legend('SC','MRC','GAC'); title('Outage Probability vs SNR');

%% Assume SNR_SC, SNR_MRC, SNR_GAC are arrays of post-combining SNR (linear)
avg_SNR_SC  = mean(SNR_SC);
avg_SNR_MRC = mean(SNR_MRC);
avg_SNR_GAC = mean(SNR_GAC);

fprintf('Average post-combining SNR:\n');
fprintf('SC  : %.2f dB\n', 10*log10(avg_SNR_SC));
fprintf('MRC : %.2f dB\n', 10*log10(avg_SNR_MRC));
fprintf('GAC : %.2f dB\n', 10*log10(avg_SNR_GAC));

% Optional: bar plot
figure;
bar( [10*log10(avg_SNR_SC), 10*log10(avg_SNR_MRC), 10*log10(avg_SNR_GAC)],0.5);
ylabel('Average SNR (dB)'); xticklabels({'SC','MRC','GAC'}); grid on;
title('Average Post-Combining SNR Comparison');