clc; clear;

%% Part 1: Spatial Diversity WSN → Base Station

numBits = 1e5;
SNRdB = 0:5:50;      % SNR range for WSN link
numSNR = length(SNRdB);
ber = zeros(1,numSNR);

for i = 1:numSNR
    snr = 10^(SNRdB(i)/10);
    errorCount = 0;
    
    % Alamouti 2 Tx, 1 Rx
    for k = 1:2:numBits-1
        bits = randi([0 1], 2, 1);
        s = 2*bits - 1;  % BPSK
        
        % Rayleigh fading channel for 2 Tx antennas
        h = (randn(2,1) + 1j*randn(2,1))/sqrt(2);
        
        % Noise variance
        N0 = 1/snr;
        noiseVar = N0/2;
        
        % Received signals
        r1 = h(1)*s(1) + h(2)*s(2) + sqrt(noiseVar)*(randn + 1j*randn);
        r2 = -h(1)*conj(s(2)) + h(2)*conj(s(1)) + sqrt(noiseVar)*(randn + 1j*randn);
        
        % Decode
        s1_hat = conj(h(1))*r1 + h(2)*conj(r2);
        s2_hat = conj(h(2))*r1 - h(1)*conj(r2);
        
        bits_hat = zeros(2,1);
        bits_hat(1) = real(s1_hat) > 0;
        bits_hat(2) = real(s2_hat) > 0;
        
        errorCount = errorCount + sum(bits ~= bits_hat);
    end
    ber(i) = errorCount/numBits;
    fprintf('WSN link SNR = %d dB, BER = %.5f\n', SNRdB(i), ber(i));
end

%% Part 2: Base Station → LEO Satellite Uplink (Ka-band)

freq_ka = 30e9;            % 30 GHz Ka-band
c = 3e8;
lambda_ka = c/freq_ka;

d_leo = 800e3;             % 800 km distance (LEO)
Pt_bs_dBm = 40;            % BS transmit power in dBm
Gt_bs_dBi = 35;            % BS antenna gain in dBi
Gr_sat_dBi = 35;           % Satellite antenna gain in dBi

% Free Space Path Loss (FSPL) in dB
FSPL_leo_dB = 20*log10(4*pi*d_leo/lambda_ka);

% Atmospheric attenuation (rain fade)
rainAtt_dB = 5;  % typical Ka-band moderate rain attenuation

% Received power at satellite (dBm)
Pr_sat_dBm = Pt_bs_dBm + Gt_bs_dBi + Gr_sat_dBi - FSPL_leo_dB - rainAtt_dB;

% Noise power (dBm)
k = 1.38e-23;     % Boltzmann constant
T = 500;          % Noise temperature in K
B = 20e6;         % Bandwidth 20 MHz
N_dBm = 10*log10(k*T*B) + 30;

% SNR at satellite uplink (dB)
SNR_sat_dB = Pr_sat_dBm - N_dBm;

fprintf('\nLEO Satellite Uplink Parameters:\n');
fprintf('Distance: %.0f km\n', d_leo/1e3);
fprintf('FSPL: %.2f dB\n', FSPL_leo_dB);
fprintf('Rain Attenuation: %.2f dB\n', rainAtt_dB);
fprintf('Received Power at Satellite: %.2f dBm\n', Pr_sat_dBm);
fprintf('Noise Power: %.2f dBm\n', N_dBm);
fprintf('SNR at Satellite Uplink: %.2f dB\n', SNR_sat_dB);

%% Plot BER for spatial diversity link

figure;
semilogy(SNRdB, ber, 'b-o','LineWidth',2);
grid on;
xlabel('SNR at WSN Base Station Link (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR: WSN to Base Station with 2x1 Alamouti Spatial Diversity');
legend('Alamouti 2 Tx, 1 Rx');

% Show satellite link SNR on same figure as vertical line
hold on;
xline(SNR_sat_dB, 'r--', 'LineWidth', 2, 'Label', 'LEO Satellite Uplink SNR');

