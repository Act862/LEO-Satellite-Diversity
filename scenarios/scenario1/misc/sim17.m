% VSAT <-> LEO Satellite communication with Doppler + MRC

clear; clc; close all;

% Constants
c = 3e8;                 % Speed of light (m/s)
fc = 20e9;               % Carrier frequency 20 GHz
lambda = c/fc;           % Wavelength (m)
Pt_dBm = 30;             % Transmit power in dBm (1 Watt = 30 dBm)
Pt = 10^((Pt_dBm-30)/10);% Transmit power in Watts
Gt_dBi = 20;             % VSAT antenna gain in dBi
Gr_dBi = 30;             % Satellite antenna element gain in dBi
Gt = 10^(Gt_dBi/10);
Gr = 10^(Gr_dBi/10);
No_dBm = -174 + 10*log10(1e6); % Noise power for 1 MHz BW in dBm
No = 10^((No_dBm-30)/10);

% LEO satellite parameters
hLEO = 600e3;            % Altitude of LEO satellite in meters (600 km)
vLEO = 7.5e3;            % Satellite velocity (m/s), typical LEO speed

% Simulation parameters
elev_angles = 10:5:90;   % Elevation angles in degrees
exponents = 0:4;
M_antennas = 2.^exponents;       % Number of receive antennas

% Preallocate SNR matrix
SNR = zeros(length(M_antennas), length(elev_angles));

for idxM = 1:length(M_antennas)
    M = M_antennas(idxM);
    
    for idxEl = 1:length(elev_angles)
        theta = elev_angles(idxEl);
        
        % Calculate distance from VSAT to LEO satellite at elevation theta
        % Using simplified Earth geometry:
        Re = 6371e3;  % Earth radius (m)
        elev_rad = deg2rad(theta);
        
        % Range calculation (slant range)
        d = sqrt( (Re+hLEO)^2 - (Re*cos(elev_rad))^2 ) - Re*sin(elev_rad);
        
        % Free space path loss (FSPL)
        FSPL = (4*pi*d/lambda)^2;
        
        % Doppler shift
        % Assuming relative velocity component along line of sight:
        v_rel = vLEO * cos(elev_rad);
        fd = (v_rel/lambda);
        
        % Received power (single antenna)
        Pr = Pt * Gt * Gr / FSPL;
        
        % Noise power (Watts)
        % Already defined as No
        
        % MRC: SNR improves linearly with number of antennas M
        SNR_single = Pr / No;
        SNR_MRC = M * SNR_single;
        
        % Store in matrix (in dB)
        SNR(idxM, idxEl) = 10*log10(SNR_MRC);
    end
end

% Plot SNR vs Elevation Angle for different antenna numbers
figure;
hold on; grid on; box on;
colors = jet(length(M_antennas));
for idxM = 1:length(M_antennas)
    plot(elev_angles, SNR(idxM,:), 'LineWidth', 1.5, 'Color', colors(idxM,:));
end
xlabel('Elevation Angle (degrees)');
ylabel('Received SNR (dB)');
title('Received SNR vs Elevation Angle for Different Number of Antennas (MRC)');
legendStrings = arrayfun(@(x) sprintf('M = %d', x), M_antennas, 'UniformOutput', false);
legend(legendStrings, 'Location', 'bestoutside');
set(gca,'FontSize',12);
