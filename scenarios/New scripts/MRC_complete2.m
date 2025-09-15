%% ============================================================
%  Explicit Rician MRC Simulation with Combining Vector
%  Tracks per-branch phase, magnitude, and applied MRC weights
%  ============================================================
clear; clc; rng(123);

%% PARAMETERS
N_branches   = 4;                  % number of satellites/branches
K_dB         = [6 6 6 6];          % per-branch Rician K-factor in dB
K_lin        = 10.^(K_dB/10);      % linear
G_branch     = [1 0.85 0.75 0.6];  % large-scale gain per branch
N_trials     = 5e4;                % Monte Carlo trials

%% HELPER: Generate Rician Fading Sample
rician = @(Klin) ...
    (sqrt(Klin/(Klin+1))*exp(1j*2*pi*rand)) + ...     % LOS specular
    (sqrt(1/(Klin+1))*(randn+1j*randn)/sqrt(2));      % scattered component

%% STORAGE
SNR_out_lin = zeros(N_trials,1);
H_matrix    = zeros(N_trials,N_branches); % store complex channels per trial
W_matrix    = zeros(N_trials,N_branches); % MRC combining vector per trial

%% SIMULATION LOOP
for t = 1:N_trials
    h_vec = zeros(1,N_branches);
    for b = 1:N_branches
        h_vec(b) = rician(K_lin(b));
    end
    
    H_matrix(t,:) = h_vec;                 % store complex channel
    
    % --- MRC combining vector (phase-alignment)
    w_vec = conj(h_vec);                    % conjugate = phase alignment
    W_matrix(t,:) = w_vec;
    
    % Received signals per branch (assume unit symbol x=1)
    x_vec = sqrt(G_branch) .* ones(1,N_branches); 
    y_comb = sum(w_vec .* x_vec);           % MRC combined signal
    
    % Post-combining SNR (assuming noise = 1)
    SNR_out_lin(t) = abs(y_comb)^2;
end

%% OUTPUT STATISTICS
SNR_out_dB = 10*log10(SNR_out_lin);
meanSNRdB = mean(SNR_out_dB);
fprintf('Average post-combining SNR: %.2f dB\n', meanSNRdB);

%% PLOT: CDF of SNR
figure; 
[S_sorted, idx] = sort(SNR_out_dB);
ecdf_y = (1:N_trials)/N_trials;
plot(S_sorted, ecdf_y,'LineWidth',1.5); grid on;
xlabel('Post-Combining SNR (dB)'); ylabel('CDF');
title('MRC over Rician Fading with Explicit Phase & Combining Vector');

%% PLOT: Histogram of per-branch channel phase (first branch example)
figure; 
histogram(angle(H_matrix(:,1)),50); 
xlabel('Phase (radians)'); ylabel('Counts');
title('Distribution of Phase for Branch 1');
grid on;
