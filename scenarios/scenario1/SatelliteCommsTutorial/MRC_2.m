% Eb/N0 Vs BER and SER for PSK over Rayleigh flat fading with MRC
% Translated and extended for BER analysis

clc; clear;

%---------Input Fields------------------------
nSym = 1e6;                        % Number of symbols
EbN0dBs = -20:2:34;                % Eb/N0 range in dB
N = [1, 2, 4, 8, 10];              % Number of diversity branches
M = 4;                             % M-ary PSK (QPSK)
k = log2(M);                       % Bits per symbol
EsN0dBs = 10*log10(k) + EbN0dBs;   % Convert Eb/N0 to Es/N0

% Create figures
figure(1); hold on;
title('SER for QPSK over Rayleigh Fading with MRC');
xlabel('Eb/N0 (dB)'); ylabel('Symbol Error Rate');
grid on;

figure(2); hold on;
title('BER for QPSK over Rayleigh Fading with MRC');
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate');
grid on;

for idxN = 1:length(N)
    nRx = N(idxN);

    % Generate random input symbols
    inputSyms = randi([0 M-1], nSym, 1);
    
    % Map symbols to bits (for BER)
    inputBits = de2bi(inputSyms, k, 'left-msb'); % matrix of bits [nSym x k]
    inputBits = reshape(inputBits.', [], 1);     % serialize to [nSym*k x 1]

    % PSK modulation (unit energy)
    s = pskmod(inputSyms, M, pi/M);

    % Repeat symbols for diversity branches
    s_diversity = repmat(s.', nRx, 1);

    % Initialize error rates
    ser_sim = zeros(1, length(EbN0dBs));
    ber_sim = zeros(1, length(EbN0dBs));

    for i = 1:length(EsN0dBs)
        EsN0dB = EsN0dBs(i);

        % Rayleigh fading coefficients
        h = sqrt(1/2) * (randn(nRx, nSym) + 1j*randn(nRx, nSym));
        signal = h .* s_diversity;

        % Calculate power and noise
        P = sum(abs(signal).^2, 2) / nSym;
        gamma = 10^(EsN0dB/10);
        N0 = P / gamma;

        % Noise generation
        noise = (randn(nRx, nSym) + 1j*randn(nRx, nSym)) .* sqrt(N0/2);
        r = signal + noise;

        % MRC combining
        equalized = sum(r .* conj(h), 1);

        % Demodulation
        detectedSyms = pskdemod(equalized.', M, pi/M);
        ser_sim(i) = sum(detectedSyms ~= inputSyms) / nSym;

        % Demapped bits
        detectedBits = de2bi(detectedSyms, k, 'left-msb');
        detectedBits = reshape(detectedBits.', [], 1);

        % Bit Error Rate
        ber_sim(i) = sum(detectedBits ~= inputBits) / length(inputBits);
    end

    % Plotting
    figure(1);
    semilogy(EbN0dBs, ser_sim, 'DisplayName', ['N = ' num2str(nRx)]);

    figure(2);
    semilogy(EbN0dBs, ber_sim, 'DisplayName', ['N = ' num2str(nRx)]);
end

% Final plot styling
figure(1);
xlim([-20 35]); ylim([1e-4 1]);
legend('show');

figure(2);
xlim([-20 35]); ylim([1e-5 1]);
legend('show');
