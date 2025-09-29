% Filename: VSAT_LEO_Outage_Analytical.m
% Analytical outage probability for VSAT-LEO MRC Rician fading link

clear; clc; close all;

%% Parameters (same as before)
c = 3e8; fc = 20e9; lambda = c/fc;
Pt_dBm = 30; Pt = 10^((Pt_dBm-30)/10);
Gt_dBi = 20; Gr_dBi = 30;
Gt = 10^(Gt_dBi/10);
Gr = 10^(Gr_dBi/10);
No_dBm = -174 + 10*log10(1e6);
No = 10^((No_dBm-30)/10);

hLEO = 600e3; Re = 6371e3;

elev_angles = 0:5:180;
M_antennas = 1:2:8;

K = 10; % Rician K factor
gamma_th_dB = 10;
gamma_th = 10^(gamma_th_dB/10);

%% Preallocate outage probability matrix
P_out_analytical = zeros(length(M_antennas), length(elev_angles));

for idxM = 1:length(M_antennas)
    M = M_antennas(idxM);
    
    for idxEl = 1:length(elev_angles)
        theta = elev_angles(idxEl);
        elev_rad = deg2rad(theta);
        
        % Distance and FSPL
        d = sqrt((Re+hLEO)^2 - (Re*cos(elev_rad))^2) - Re*sin(elev_rad);
        FSPL = (4*pi*d/lambda)^2;
        
        % Average SNR per branch
        Pr = Pt * Gt * Gr / FSPL;
        avg_snr = Pr / No;
        
        % Arguments for Marcum Q-function
        a = sqrt(2*K*M);
        b = sqrt(2*(M+K*M)*gamma_th/avg_snr);
        
        % Calculate outage probability using generalized Marcum Q-function
        % MATLAB's marcumq supports order M (R2021a+)
        try
            Qval = marcumq(a,b,M);
        catch
            % Fallback for older versions: approximate with order=1
            warning('marcumq with order not supported, approximating with order 1');
            Qval = marcumq(a,b);
        end
        
        P_out_analytical(idxM, idxEl) = 1 - Qval;
    end
end

%% Plot Analytical Outage Probability
[X, Y] = meshgrid(elev_angles, M_antennas);

figure('Color','w');
surf(X, Y, log10(P_out_analytical), 'EdgeColor', 'none');
colormap(jet);
colorbar;
caxis([-5 0])
colorbar('Ticks', -5:0, 'TickLabels', {'1e-5','1e-4','1e-3','1e-2','1e-1','1'});
xlabel('Elevation Angle (degrees)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Number of Antennas (M)', 'FontWeight', 'bold', 'FontSize', 12);
zlabel('Log_{10}(Outage Probability)', 'FontWeight', 'bold', 'FontSize', 12);
title(sprintf('Analytical Outage Probability vs Elevation Angle and Antennas\nRician K=%d, Threshold=%d dB', K, gamma_th_dB), ...
      'FontWeight', 'bold', 'FontSize', 14);
view(45,30);
grid on;
set(gca,'FontSize',12);
