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

% Simulation parameters
elev_angles = 10:5:90;   % Elevation angles in degrees
exponents = 0:4;
M_antennas = 2.^exponents;       % Number of receive antennas

% Rician fading parameters
K = 0;                  % Rician K factor (linear)
gamma_th_dB = 10;        % Outage threshold SNR in dB
gamma_th = 10^(gamma_th_dB/10);

% Preallocate outage probability matrix
P_out = zeros(length(M_antennas), length(elev_angles));

% Helper function: Marcum Q-function Q1(a,b) for outage calc
% Use MATLAB's marcumq function if available (R2021a+),
% otherwise approximate below
marcumQ = @(a,b) marcumq(a,b);

for idxM = 1:length(M_antennas)
    M = M_antennas(idxM);
    
    for idxEl = 1:length(elev_angles)
        theta = elev_angles(idxEl);
        
        % Distance calculation (slant range)
        Re = 6371e3;  % Earth radius (m)
        elev_rad = deg2rad(theta);
        d = sqrt( (Re+hLEO)^2 - (Re*cos(elev_rad))^2 ) - Re*sin(elev_rad);
        
        % Free space path loss (FSPL)
        FSPL = (4*pi*d/lambda)^2;
        
        % Average received SNR per antenna (linear)
        Pr = Pt * Gt * Gr / FSPL;
        avg_snr = Pr / No;
        
        % Parameters for Rician fading combined SNR
        % Mean SNR per branch (linear scale)
        Omega = avg_snr; % average received power (per branch)
        % Scale factor for the noncentral chi-square distribution in MRC
        
        % Approximate outage probability with MRC and Rician fading:
        % According to Simon & Alouini (Digital Comm. over Fading Channels):
        % P_out = 1 - Q_M(sqrt(2*K*M), sqrt(2*gamma_th/Omega * (1+K)))
        % where Q_M is the generalized Marcum Q-function of order M.
        % We'll use MATLAB's marcumq in a loop for order M.
        
        % Calculate the argument of Marcum Q-function
        a = sqrt(2*K*M);
        b = sqrt(2*gamma_th/Omega * (1+K));
        
        % Generalized Marcum Q for order M can be computed by
        % repeated calls or approximations.
        % Here, we'll approximate for integer M by:
        % P_out = 1 - marcumq(a,b,M)
        % MATLAB marcumq supports the order argument from R2021a onwards.
        
        try
            % MATLAB built-in function with order M
            outage_prob = 1 - marcumq(a,b,M);
        catch
            % If marcumq with order not available, approximate M=1 case only
            % or fallback to single antenna outage (conservative)
            if M == 1
                outage_prob = 1 - marcumq(a,b);
            else
                % Approximate using single antenna outage to avoid error
                outage_prob = 1 - marcumq(a,b);
                % Warn user
                warning('marcumq with order not supported. Approximating MRC with single antenna outage.');
            end
        end
        
        % Store outage probability
        P_out(idxM, idxEl) = outage_prob;
    end
end

% Plot outage probability vs elevation angle for different antennas
figure;
hold on; grid on; box on;
colors = jet(length(M_antennas));
for idxM = 1:length(M_antennas)
    semilogy(elev_angles, P_out(idxM,:), 'LineWidth', 1.5, 'Color', colors(idxM,:));
end
xlabel('Elevation Angle (degrees)');
ylabel('Outage Probability');
title(sprintf('Outage Probability vs Elevation Angle for Different Number of Antennas (Rician K=%d)', K));
legendStrings = arrayfun(@(x) sprintf('M = %d', x), M_antennas, 'UniformOutput', false);
legend(legendStrings, 'Location', 'bestoutside');
set(gca,'FontSize',12);
ylim([1e-5 1]);
