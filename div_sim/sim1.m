%   Modulation: BPSK
%   Channel: Independent Rayleigh fading for each VSAT
%   Diversity Technique: SC and MRC
%   Goal: Show BER vs SNR for each case
clc; clear; close all;
N = 1e5; % number of bits
SNR_dB = 0:2:20; % SNR range in dB
SNR = 10.^(SNR_dB/10); % linear SNR
ber_nodiv = zeros(size(SNR));
ber_sc = zeros(size(SNR));
ber_mrc = zeros(size(SNR));

%   Generate random bits
bits = randi([0 1], 1, N);
bpsk = 2*bits - 1;% BPSK mapping: 0 --> -1, 1 --> +1

for k=1:length(SNR)
    noise_var = 1/SNR(k);
    %   Rayleigh channel coefficients
    h1 = (randn(1,N) + 1i*randn(1,N))/sqrt(2);
    n1 = sqrt(noise_var/2) * (randn(1,N)+1i*randn(1,N));
    r1 = h1 .* bpsk + n1;
    r_eq = r1 ./ h1;
    bits_rx = real(r_eq) > 0;
    ber_nodiv(k) = sum(bits_rx ~= bits) / N;

    %   Selection Combining
    h2 = (randn(1,N) + 1i*randn(1,N))/sqrt(2);
    n2 = sqrt(noise_var/2) * (randn(1,N)+1i*randn(1,N));
    r2 = h2 .* bpsk + n2;

    %   Select link with higher SNR
    snr1 = abs(h1).^2;
    snr2 = abs(h2).^2;
    use_first = snr1 > snr2;

    r_sc = zeros(1,N);
    h_sc = zeros(1,N);
    r_sc(use_first) = r1(use_first);
    r_sc(~use_first) = r2(~use_first);
    h_sc(use_first) = h1(use_first);
    h_sc(~use_first) = h2(~use_first);

    r_eq_sc = r_sc ./ h_sc;
    bits_rx_sc = real(r_eq_sc) > 0;
    ber_sc(k) = sum(bits_rx_sc ~= bits)/N;

    %   Maximum Ration Combining (MRC)
    r_mrc = conj(h1).*r1 + conj(h2).*r2;
    h_eq_mrc = abs(h1).^2 + abs(h2).^2;
    r_eq_mrc = r_mrc ./ h_eq_mrc;
    bits_rx_mrc = real(r_eq_mrc) > 0;
    ber_mrc(k) = sum(bits_rx_mrc ~= bits)/N;
end

% Plot
figure;
semilogy(SNR_dB, ber_nodiv, 'r-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB, ber_sc, 'b-*', 'LineWidth', 2);
semilogy(SNR_dB, ber_mrc, 'g-s', 'LineWidth', 2);
grid on; xlabel('SNR (dB)'); ylabel('BER');
legend('No Diversity', 'Selection Combining', 'MRC');
title('Satellite Spatial Diversity - BER Comparison');

%   No diversity: Only one satellite link, susceptible to deep fades
%   SC (Selection Combining): Picks better link, improving reliability
%   MRC (Maximum Ratio Combining): Combines both links optimally for best
%   BER
