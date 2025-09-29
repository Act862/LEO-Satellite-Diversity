clear; clc;

%% Simulation Parameters
Fs = 1e6;
% M = 16; % QPSK
M = 2;
k = log2(M);
bitsPerFrame = 10000;
EbN0dB = 0:2:20; % SNR range in dB
nSNR = length(EbN0dB);

ber_mrc = zeros(1, nSNR);
ber_sc = zeros(1, nSNR);

%% Channel Definition
ricianchan = comm.RicianChannel( ...
    SampleRate=Fs, ...
    PathDelays=[0.0 0.6 1.2]*1e-6, ...
    AveragePathGains=[0.1 0.5 0.2], ...
    KFactor=2.8, ...
    DirectPathDopplerShift=5e3, ...
    DirectPathInitialPhase=0.5, ...
    MaximumDopplerShift=50e3, ...
    DopplerSpectrum=doppler('Bell',8), ...
    RandomStream='mt19937ar with seed', ...
    Seed=73, ...
    PathGainsOutputPort=true);

%% Loop over SNR
for snrIdx = 1:nSNR
    % SNR and Noise Setup
    EbN0 = 10^(EbN0dB(snrIdx)/10);
    EsN0 = EbN0 * k; % symbol-level SNR
    N0 = 1/EsN0;
    sigma = sqrt(N0/2);

    % Generate bits and modulate
    bits = randi([0 1], bitsPerFrame, 1);
    tx = pskmod(bi2de(reshape(bits, [], k)), M);

    % Prepare buffers
    rxSym_mc = zeros(1, length(tx));
    rxSym_sc = zeros(1, length(tx));

    % Reset channel for repeatability
    reset(ricianchan);

    % Transmit through channel
    for i = 1:length(tx)
        [~, h] = ricianchan(tx(i));
        h_norm = sqrt(sum(abs(h).^2));
        w = h_norm .* h';

        % Noise for each branch
        n1 = sigma * (randn(1,1) + 1j*randn(1,1));
        n2 = sigma * (randn(1,1) + 1j*randn(1,1));
        n3 = sigma * (randn(1,1) + 1j*randn(1,1));

        % Received signals and co-phasing
        rx1 = h(1)*tx(i) + n1;
        rx2 = h(2)*tx(i)*exp(-1i*0.6*1e-6) + n2;
        rx3 = h(3)*tx(i)*exp(-1i*1.2*1e-6) + n3;

        % MRC combining
        w_h = conj(w)';
        rx_combined_mrc = w_h * [rx1; rx2; rx3];

        % Selection combining (SC)
        SNRs = [abs(rx1)^2, abs(rx2)^2, abs(rx3)^2];
        [~, bestIdx] = max(SNRs);
        r = [rx1, rx2, rx3];
        rx_combined_rc = r(bestIdx);

        % Demodulate
        rxSym_mc(i) = pskdemod(rx_combined_mrc, M);
        rxSym_sc(i) = pskdemod(rx_combined_rc, M);
    end

    % Reshape and count errors
    rxBits1 = reshape(de2bi(rxSym_mc, k), [], 1);
    rxBits2 = reshape(de2bi(rxSym_sc, k), [], 1);
    bits = bits(1:length(rxBits1)); % align in case of mismatch

    [errs_mrc, ber_mrc(snrIdx)] = biterr(bits, rxBits1);
    [errs_sc,  ber_sc(snrIdx)]  = biterr(bits, rxBits2);
    % 
    % fprintf("SNR = %2d dB: MRC BER = %.5g (%d errs), SC BER = %.5g (%d errs)\n", ...
    %     EbN0dB(snrIdx), ber_mrc(snrIdx), errs_mrc, ber_sc(snrIdx), errs_sc);
end

%% Plot BER curves
figure;
%   theoretic
lambda = sqrt(EbN0dB./(2+EbN0dB));
theoretical_ber = ((1-lambda)/2).^3+3*((1-lambda)/2).^4+(3/2)*((1-lambda)/2).^6;
semilogy(EbN0dB, theoretical_ber, 'LineWidth', 2); hold on;
semilogy(EbN0dB,ber_mrc,'-d','LineWidth',2);
semilogy(EbN0dB, ber_sc,'-s','LineWidth', 2);
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. E_b/N_0 for MRC and SC in Rician Fading Channel');
legend('Average BER', 'MRC', 'Selection Combining');
grid on;
