% Eb/N0 Vs SER for PSK over Rayleigh flat fading with MRC
% Translated from Python version by Mathuranathan Viswanathan

clc; clear;

%---------Input Fields------------------------
nSym = 1e6;                        % Number of symbols
EbN0dBs = -20:2:34;                % Eb/N0 range in dB
N = [1, 2, 4, 8, 10];              % Number of diversity branches
M = 4;                             % M-ary PSK (QPSK)
k = log2(M);                       % Bits per symbol
EsN0dBs = 10*log10(k) + EbN0dBs;   % Convert Eb/N0 to Es/N0

% Create figure
figure;
hold on;

for idxN = 1:length(N)
    nRx = N(idxN);

    % Generate random input symbols
    inputSyms = randi([0 M-1], nSym, 1);
    
    % PSK modulation (unit energy)
    s = pskmod(inputSyms, M, pi/M);  % Use pi/M to rotate constellation (optional)

    % Repeat symbols for diversity branches
    s_diversity = repmat(s.', nRx, 1);

    % Initialize SER result
    ser_sim = zeros(1, length(EbN0dBs));

    for i = 1:length(EsN0dBs)
        EsN0dB = EsN0dBs(i);

        % Generate Rayleigh fading channel coefficients
        h = sqrt(1/2) * (randn(nRx, nSym) + 1j*randn(nRx, nSym));
        signal = h .* s_diversity; % Apply fading

        % Calculate signal power per branch
        P = sum(abs(signal).^2, 2) / nSym;
        gamma = 10^(EsN0dB/10); % Linear scale
        N0 = P / gamma;

        % Generate noise
        noise = (randn(nRx, nSym) + 1j*randn(nRx, nSym)) .* sqrt(N0/2);

        % Received signal
        r = signal + noise;

        % MRC equalization
        equalized = sum(r .* conj(h), 1);

        % Demodulate
        detectedSyms = pskdemod(equalized.', M, pi/M);

        % SER calculation
        ser_sim(i) = sum(detectedSyms ~= inputSyms) / nSym;
    end

    % Plot SER
    semilogy(EbN0dBs, ser_sim, 'DisplayName', ['N = ' num2str(nRx)]);
end

% Final plot formatting
grid on;
xlim([-20 35]);
ylim([1e-4 1.1]);
xlabel('Eb/N0 (dB)');
ylabel('Symbol Error Rate (P_s)');
title('SER performance for QPSK over Rayleigh fading channel with MRC');
legend('show');
