%% BER Comparison: Optimal MRC vs Imperfect CSI MRC
clear; clc; rng(2025);

%% Parameters
Nrx = 4;               % number of receive antennas (SIMO)
Nbits = 1e6;           % number of transmitted bits
SNRdB = 5;            % nominal SNR per branch
sigma2 = 10^(-SNRdB/10); % noise variance

Kfactor = 5;           % Rician K-factor
errVar = 0.05;         % channel estimation error variance

%% Generate random BPSK symbols
bits = randi([0 1], Nbits, 1);
symbols = 2*bits - 1;   % BPSK mapping: 0 -> -1, 1 -> +1

%% Initialize error counters
err_opt = 0;
err_imp = 0;

%% Transmission loop
for n = 1:Nbits
    % Generate Rician fading channel (Nrx x 1)
    LOS = sqrt(Kfactor/(Kfactor+1)) * ones(Nrx,1);
    NLOS = sqrt(1/(2*(Kfactor+1))) * (randn(Nrx,1) + 1j*randn(Nrx,1));
    h = LOS + NLOS;

    % Received signal at each antenna
    noise = sqrt(sigma2/2) * (randn(Nrx,1) + 1j*randn(Nrx,1));
    r = h * symbols(n) + noise;

    % --- Optimal MRC (perfect CSI) ---
    w_opt = conj(h);
    y_opt = w_opt.' * r / norm(w_opt);
    bit_hat_opt = real(y_opt) > 0;

    % --- Imperfect MRC (with channel estimation error) ---
    h_hat = h + sqrt(errVar/2)*(randn(Nrx,1)+1j*randn(Nrx,1));
    w_imp = conj(h_hat);
    y_imp = w_imp.' * r / norm(w_imp);
    bit_hat_imp = real(y_imp) > 0;

    % Count errors
    err_opt = err_opt + (bit_hat_opt ~= bits(n));
    err_imp = err_imp + (bit_hat_imp ~= bits(n));
end

%% Compute BER
BER_opt = err_opt / Nbits;
BER_imp = err_imp / Nbits;

%% Display results
disp(['BER (Optimal MRC, perfect CSI):   ' num2str(BER_opt)]);
disp(['BER (Imperfect MRC, errVar=' num2str(errVar) '): ' num2str(BER_imp)]);

%% Sweep over SNR values for plotting
SNRdB_range = 0:30;
BER_opt_range = zeros(size(SNRdB_range));
BER_imp_range = zeros(size(SNRdB_range));

for idx = 1:length(SNRdB_range)
    SNRdB = SNRdB_range(idx);
    sigma2 = 10^(-SNRdB/10);
    Nbits_sim = 1e6;  % reduce bits for faster sweep

    err_opt = 0;
    err_imp = 0;

    for n = 1:Nbits_sim
        % Channel
        LOS = sqrt(Kfactor/(Kfactor+1)) * ones(Nrx,1);
        NLOS = sqrt(1/(2*(Kfactor+1))) * (randn(Nrx,1) + 1j*randn(Nrx,1));
        h = LOS + NLOS;

        % Transmit symbol
        bit = randi([0 1]);
        s = 2*bit - 1;

        noise = sqrt(sigma2/2) * (randn(Nrx,1) + 1j*randn(Nrx,1));
        r = h*s + noise;

        % Perfect CSI
        w_opt = conj(h);
        y_opt = w_opt.'*r / norm(w_opt);
        bit_hat_opt = real(y_opt) > 0;

        % Imperfect CSI
        h_hat = h + sqrt(errVar/2)*(randn(Nrx,1)+1j*randn(Nrx,1));
        w_imp = conj(h_hat);
        y_imp = w_imp.'*r / norm(w_imp);
        bit_hat_imp = real(y_imp) > 0;

        % Errors
        err_opt = err_opt + (bit_hat_opt ~= bit);
        err_imp = err_imp + (bit_hat_imp ~= bit);
    end

    BER_opt_range(idx) = err_opt / Nbits_sim;
    BER_imp_range(idx) = err_imp / Nbits_sim;
end

%% Plot BER curves
figure;
semilogy(SNRdB_range, BER_opt_range, 'b-o','LineWidth',1.5); hold on;
semilogy(SNRdB_range, BER_imp_range, 'r-s','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
legend('Optimal MRC (perfect CSI)','MRC with Imperfect CSI','Location','southwest');
title(['BER Performance with ' num2str(Nrx) ' Rx antennas, Rician fading']);
