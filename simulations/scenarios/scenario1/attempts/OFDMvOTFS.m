clc; clear;

%% Parameters
N = 64;              % Number of subcarriers
M = 4;               % QPSK modulation order
numSymbols = 1000;   % Number of OFDM symbols

% Carrier frequency and sampling frequency
fc_5G = 3.5e9;       % 5G mid-band frequency
fc_ka = 30e9;        % Ka-band frequency for NTN
Fs = 15e3*N;         % Sampling frequency (assuming 15 kHz subcarrier spacing)

% Doppler
v_5G = 3;            % 3 m/s (walking speed, 5G)
v_LEO = 7500;        % 7.5 km/s typical LEO speed
c = 3e8;

fd_5G = v_5G*fc_5G/c;
fd_LEO = v_LEO*fc_ka/c;

% SNR range for simulation
SNRdB = 0:5:30;

%% QPSK Modulator/Demodulator
qpskMod = @(bits) exp(1j*pi/4) * (1 - 2*bits(1:2:end) + 1j*(1 - 2*bits(2:2:end))) / sqrt(2);
qpskDemod = @(sym) [real(sym)<0; imag(sym)<0];

%% OFDM Modulation
ofdmMod = @(x) ifft(x, N);
ofdmDemod = @(x) fft(x, N);

%% OTFS Modulation/Demodulation (simplified)
% Reference: OTFS uses ISFFT and SFFT over delay-Doppler grid

% ISFFT: Delay-Doppler to time-frequency
ISFFT = @(X) sqrt(N)*fft(ifft(X,[],1),[],2);

% SFFT: Time-frequency to delay-Doppler
SFFT = @(Y) (1/sqrt(N))*fft(ifft(Y,[],2),[],1);

%% Function to generate random bits
genBits = @(num) randi([0 1], num, 1);

%% Channel: AWGN + Doppler (simple model)
applyChannel = @(x, fd, Fs, SNR) ...
    awgn(x .* exp(1j*2*pi*fd*(0:length(x)-1)'/Fs), SNR, 'measured');

%% BER computation helper
computeBER = @(txBits, rxBits) sum(txBits~=rxBits)/length(txBits);

%% Simulation Loop

ber_ofdm_5g = zeros(length(SNRdB),1);
ber_ofdm_ntn = zeros(length(SNRdB),1);
ber_otfs_ntn = zeros(length(SNRdB),1);

for idx = 1:length(SNRdB)
    SNR = SNRdB(idx);
    
    % Generate bits for 1 OFDM symbol (N subcarriers * 2 bits/symbol)
    bits = genBits(2*N);
    
    % QPSK modulation
    qpskSyms = qpskMod(bits);
    
    % --------- 5G OFDM link ----------
    % OFDM modulation
    tx_ofdm = ofdmMod(qpskSyms);
    
    % Channel: Doppler + AWGN 5G
    rx_ofdm = applyChannel(tx_ofdm, fd_5G, Fs, SNR);
    
    % OFDM demodulation
    rx_syms = ofdmDemod(rx_ofdm);
    
    % QPSK demodulation
    rx_bits = zeros(2*N,1);
    rx_bits(1:2:end) = real(rx_syms) < 0;
    rx_bits(2:2:end) = imag(rx_syms) < 0;
    
    ber_ofdm_5g(idx) = computeBER(bits, rx_bits);
    
    % --------- NTN OFDM link ----------
    % Same bits modulated
    tx_ofdm_ntn = ofdmMod(qpskSyms);
    
    % Channel: higher Doppler LEO + AWGN
    rx_ofdm_ntn = applyChannel(tx_ofdm_ntn, fd_LEO, Fs, SNR);
    
    % OFDM demodulation
    rx_syms_ntn = ofdmDemod(rx_ofdm_ntn);
    
    % QPSK demodulation
    rx_bits_ntn = zeros(2*N,1);
    rx_bits_ntn(1:2:end) = real(rx_syms_ntn) < 0;
    rx_bits_ntn(2:2:end) = imag(rx_syms_ntn) < 0;
    
    ber_ofdm_ntn(idx) = computeBER(bits, rx_bits_ntn);
    
    % --------- NTN OTFS link ----------
    % OTFS modulation (ISFFT)
    tx_otfs = ISFFT(reshape(qpskSyms, sqrt(N), sqrt(N)));
    
    % Flatten OTFS time domain samples
    tx_otfs_vec = tx_otfs(:);
    
    % Channel: Doppler + AWGN LEO
    rx_otfs_vec = applyChannel(tx_otfs_vec, fd_LEO, Fs, SNR);
    
    % Reshape received vector
    rx_otfs = reshape(rx_otfs_vec, sqrt(N), sqrt(N));
    
    % OTFS demodulation (SFFT)
    rx_otfs_dd = SFFT(rx_otfs);
    
    % Flatten
    rx_syms_otfs = rx_otfs_dd(:);
    
    % QPSK demodulation
    rx_bits_otfs = zeros(2*N,1);
    rx_bits_otfs(1:2:end) = real(rx_syms_otfs) < 0;
    rx_bits_otfs(2:2:end) = imag(rx_syms_otfs) < 0;
    
    ber_otfs_ntn(idx) = computeBER(bits, rx_bits_otfs);
end

%% Plot results

figure; 
semilogy(SNRdB, ber_ofdm_5g, '-o', 'LineWidth', 2); hold on;
semilogy(SNRdB, ber_ofdm_ntn, '-s', 'LineWidth', 2);
semilogy(SNRdB, ber_otfs_ntn, '-d', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Comparison: 5G OFDM vs NTN OFDM vs NTN OTFS');
legend('5G OFDM (low Doppler)', 'NTN OFDM (high Doppler)', 'NTN OTFS (high Doppler)');
axis([min(SNRdB) max(SNRdB) 1e-5 1]);
