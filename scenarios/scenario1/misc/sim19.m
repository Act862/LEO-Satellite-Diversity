% Filename: VSAT_LEO_BER_MRC.m
% Description: Computes and plots BER vs elevation angle and number of antennas
%              for VSAT to LEO satellite link with Rician fading and MRC.

clear; clc; close all;

%% Parameters
c = 3e8;                    % Speed of light (m/s)
fc = 20e9;                  % Carrier frequency 20 GHz
lambda = c/fc;              % Wavelength (m)
Pt_dBm = 30;                % Transmit power dBm (1 Watt)
Pt = 10^((Pt_dBm-30)/10);   % Transmit power Watts
Gt_dBi = 20;                % VSAT antenna gain dBi
Gr_dBi = 30;                % Satellite antenna element gain dBi
Gt = 10^(Gt_dBi/10);
Gr = 10^(Gr_dBi/10);
No_dBm = -174 + 10*log10(1e6); % Noise power 1 MHz bandwidth (dBm)
No = 10^((No_dBm-30)/10);       % Noise power Watts

hLEO = 600e3;               % LEO satellite altitude (m)
Re = 6371e3;                % Earth radius (m)

elev_angles = 10:5:90;      % Elevation angles (degrees)
index = 0:4;
M_antennas = 2.^index;          % Number of antennas

K = 10;                     % Rician K-factor (LOS dominance)

numTrials = 1e5;            % Monte Carlo trials per scenario

%% Preallocate BER matrix
BER = zeros(length(M_antennas), length(elev_angles));

%% Rician power sample generator function
% Returns vector of combined instantaneous SNR samples with M antennas MRC
rician_mrc_snr = @(K, avg_snr, M, N) ...
    sum(abs(sqrt(avg_snr/(K+1))*(randn(M,N) + 1i*randn(M,N)) + sqrt(K*avg_snr/(K+1))).^2,1);

%% Simulation loop
for idxM = 1:length(M_antennas)
    M = M_antennas(idxM);
    
    for idxEl = 1:length(elev_angles)
        theta = elev_angles(idxEl);
        elev_rad = deg2rad(theta);
        
        % Calculate slant range distance
        d = sqrt((Re+hLEO)^2 - (Re*cos(elev_rad))^2) - Re*sin(elev_rad);
        
        % Free space path loss
        FSPL = (4*pi*d/lambda)^2;
        
        % Average received power per antenna
        Pr = Pt * Gt * Gr / FSPL;
        
        % Average SNR per antenna
        avg_snr = Pr / No;
        
        % Generate instantaneous SNR samples after MRC
        gamma_samples = rician_mrc_snr(K, avg_snr, M, numTrials);
        
        % Theoretical BPSK BER for AWGN with instantaneous SNR gamma
        % BER per sample: Q(sqrt(2*gamma))
        ber_samples = qfunc(sqrt(2*gamma_samples));
        
        % Average BER over trials
        BER(idxM, idxEl) = mean(ber_samples);
        
        fprintf('M=%d, Elev=%d°, BER=%.4e\n', M, theta, BER(idxM, idxEl));
    end
end

%% 3D Surface plot of BER
[X, Y] = meshgrid(elev_angles, M_antennas);

figure('Color','w');
surf(X, Y, log10(BER), 'EdgeColor', 'none');
colormap(jet);
colorbar;
caxis([-6 0]) % Show BER from 10^-6 to 1
colorbar('Ticks', -6:0, 'TickLabels', {'1e-6','1e-5','1e-4','1e-3','1e-2','1e-1','1'});
xlabel('Elevation Angle (degrees)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Number of Antennas (M)', 'FontWeight', 'bold', 'FontSize', 12);
zlabel('Log_{10}(BER)', 'FontWeight', 'bold', 'FontSize', 12);
title('BER vs Elevation Angle and Number of Antennas (Rician fading + MRC)', ...
      'FontWeight', 'bold', 'FontSize', 14);
view(45,30);
grid on;
set(gca,'FontSize',12);
