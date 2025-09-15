%% ============================================================
%  MRC over Rician Fading - Complete Simulation (No Toolboxes)
%  - N branches (satellites) with independent Rician fading
%  - Per-branch K-factor (dB) configurable
%  - Outputs:
%       * Empirical CDF of post-combining SNR (dB)
%       * BER vs. average SNR for BPSK with perfect CSI & MRC
%  ============================================================
clear; clc; rng(7);

%% ------------------------ Parameters ------------------------
N_branches   = 4;                 % number of satellites/branches
K_dB         = [3, 6, 6, 10];     % per-branch Rician K-factors [dB]
K_lin        = 10.^(K_dB/10);

% Optional per-branch large-scale gains (e.g., geometry/path loss).
% Leave as ones(1,N_branches) if you don't want geometry.
G_branch     = [1.00, 0.85, 0.75, 0.60];   % relative linear gains

N_trials     = 2e5;              % Monte Carlo trials per SNR point (adjust as needed)

% For SNR-out CDF (single reference operating point)
refEbN0dB    = 10;               % average per-branch Eb/N0 (dB) at which we show CDF

% For BER curve
EbN0dB_vec   = -5:2:22;          % sweep of average per-branch Eb/N0 (dB)

%% ------------------- Helper: Rician fading -------------------
% Generate complex Rician fading samples for all branches, per trial.
% h = sqrt(K/(K+1))*e^{j*theta} + sqrt(1/(K+1))* (n_I + j n_Q)/sqrt(2)
% where theta ~ U[0,2pi), n_I,n_Q ~ N(0,1)
rician = @(Klin, M) ...
    (sqrt(Klin./(Klin+1)) .* exp(1j*2*pi*rand(1,M))) + ...
    (sqrt(1./(Klin+1))    .* (randn(1,M) + 1j*randn(1,M))/sqrt(2));

%% ----------------- One point CDF (refEbN0dB) -----------------
avgSNR_branch_lin = 10.^(refEbN0dB/10) .* G_branch;  % per-branch avg SNR (linear)

SNRout_lin_ref = zeros(N_trials,1);
for t = 1:N_trials
    gamma_sum = 0;
    for b = 1:N_branches
        h = rician(K_lin(b), 1);          % one complex fade sample
        gamma_b = avgSNR_branch_lin(b) * (abs(h)^2);
        gamma_sum = gamma_sum + gamma_b;  % MRC => sum of branch SNRs
    end
    SNRout_lin_ref(t) = gamma_sum;
end

SNRout_dB_ref = 10*log10(SNRout_lin_ref);

% Empirical CDF (no toolbox)
[SNR_sorted, idx] = sort(SNRout_dB_ref);
ecdf_y = (1:N_trials).' / N_trials;

%% ------------------- BER vs Eb/N0 (BPSK) ---------------------
% For coherent BPSK in AWGN, instantaneous BER = Q(sqrt(2*gamma_out))
% Q(x) = 0.5*erfc(x/sqrt(2)) -> so BER = 0.5*erfc(sqrt(gamma_out))
BER = zeros(size(EbN0dB_vec));
for si = 1:numel(EbN0dB_vec)
    avgSNR_branch_lin = 10.^(EbN0dB_vec(si)/10) .* G_branch;

    % Draw instantaneous post-combining SNR via Monte Carlo
    gamma_out = zeros(N_trials,1);
    for t = 1:N_trials
        gsum = 0;
        for b = 1:N_branches
            h = rician(K_lin(b), 1);
            gsum = gsum + avgSNR_branch_lin(b) * (abs(h)^2);
        end
        gamma_out(t) = gsum;
    end

    % Average BER using instantaneous SNRs
    BER(si) = mean( 0.5 * erfc( sqrt(gamma_out) ) );
end

%% ------------------------- Reporting -------------------------
meanSNRdB = mean(SNRout_dB_ref);
p10 = SNR_sorted(round(0.10*N_trials));
p50 = SNR_sorted(round(0.50*N_trials));
p90 = SNR_sorted(round(0.90*N_trials));

fprintf('MRC over Rician fading (N=%d branches)\n', N_branches);
fprintf('K-factors (dB): [%s]\n', num2str(K_dB));
fprintf('Branch gains   : [%s]\n\n', num2str(G_branch, '%.2f '));
fprintf('At Eb/N0 = %2.1f dB per branch:\n', refEbN0dB);
fprintf('  Mean SNR_out  = %5.2f dB\n', meanSNRdB);
fprintf('  10%%/50%%/90%%   = %5.2f / %5.2f / %5.2f dB (CDF quantiles)\n\n', p10, p50, p90);

%% -------------------------- Plots ----------------------------
% 1) Empirical CDF of SNR_out at refEbN0dB
figure; plot(SNR_sorted, ecdf_y, 'LineWidth', 1.5); grid on;
xlabel('Post-Combining SNR (dB)'); ylabel('CDF');
title(sprintf('MRC Output SNR CDF @ per-branch Eb/N0 = %.1f dB', refEbN0dB));

% 2) Average BER vs per-branch Eb/N0
figure; plot(EbN0dB_vec, BER, 'o-', 'LineWidth', 1.5); grid on;
xlabel('Per-branch Eb/N0 (dB)'); ylabel('Average BER (BPSK)');
title(sprintf('MRC over Rician fading (N=%d, K=[%s] dB)', N_branches, num2str(K_dB)));
ylim([1e-5 1]); set(gca, 'YScale','log');
