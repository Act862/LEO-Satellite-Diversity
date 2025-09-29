clc; clear;

% Parameters
numBits = 1e5;         % Number of bits to simulate
SNRdB = 0:2:20;        % SNR range in dB
numSNR = length(SNRdB);

% Pre-allocate BER vector
ber = zeros(1, numSNR);

for i = 1:numSNR
    snr = 10^(SNRdB(i)/10);
    errorCount = 0;
    
    % Loop over the bitstream in blocks of 2 (Alamouti needs 2 bits)
    for k = 1:2:numBits-1
        % Generate 2 bits
        bits = randi([0 1], 2, 1);
        % BPSK modulation: 0 -> -1, 1 -> +1
        s = 2*bits - 1;
        
        % Alamouti encoding
        % Time slot 1: transmit s(1) from Tx1, s(2) from Tx2
        % Time slot 2: transmit -conj(s(2)) from Tx1, conj(s(1)) from Tx2
        
        % Channel coefficients (Rayleigh fading)
        h = (randn(2,1) + 1j*randn(2,1))/sqrt(2); % 2x1 channel vector
        
        % Noise variance for AWGN
        N0 = 1/snr;
        noiseVar = N0/2; % per dimension
        
        % Time slot 1 received signal
        r1 = h(1)*s(1) + h(2)*s(2) + sqrt(noiseVar)*(randn + 1j*randn);
        
        % Time slot 2 received signal
        r2 = -h(1)*conj(s(2)) + h(2)*conj(s(1)) + sqrt(noiseVar)*(randn + 1j*randn);
        
        % Alamouti decoding
        s1_hat = conj(h(1))*r1 + h(2)*conj(r2);
        s2_hat = conj(h(2))*r1 - h(1)*conj(r2);
        
        % Decision (real part only since BPSK)
        bits_hat = zeros(2,1);
        bits_hat(1) = real(s1_hat) > 0;
        bits_hat(2) = real(s2_hat) > 0;
        
        % Count bit errors
        errorCount = errorCount + sum(bits ~= bits_hat);
    end
    
    ber(i) = errorCount/numBits;
    fprintf('SNR = %d dB, BER = %.5f\n', SNRdB(i), ber(i));
end

% Plot results
figure; 
semilogy(SNRdB, ber, 'b-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for 2x1 Alamouti Spatial Diversity (BPSK)');
legend('Alamouti 2 Tx, 1 Rx');
