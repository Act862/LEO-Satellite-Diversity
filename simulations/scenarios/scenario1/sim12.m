% Constants
fc = 30e9;                  % Carrier frequency (Hz)
c = 3e8;                    % Speed of light (m/s)
lambda = c / fc;            % Wavelength (m)
Re = 6371e3;                % Earth radius (m)
h_sat = 600e3;              % Satellite altitude (m)
vsat_speed = 0;             % VSAT is stationary
Gt_dB = 35;
Gr_dB = 35;
Gt = 10^(Gt_dB/10);
Gr = 10^(Gr_dB/10);
Pt_dBm = 30;                % Transmit power in dBm
Pt = 10^((Pt_dBm - 30)/10); % Convert to Watts

% Log-logistic fading parameters (foliage environment)
mu = 0;                     % Median of log-scale
sigma = 1.2;                % Spread (scale parameter)

% Elevation angles
theta_deg = linspace(5, 90, 100);
theta_rad = deg2rad(theta_deg);

% Slant range
slant_range = @(theta) sqrt((Re + h_sat).^2 - (Re * cos(theta)).^2) - Re * sin(theta);
d = slant_range(theta_rad);

% Satellite angular velocity
G = 6.67430e-11;
M = 5.972e24;
w_sat = sqrt(G * M / (Re + h_sat)^3); % rad/s
v_rel = w_sat * Re * cos(theta_rad);  % relative velocity at surface projection

% Doppler shift
fd = (fc / c) * v_rel;

% Doppler spread (foliage: assume local scatterers ~ λ/2)
doppler_spread = abs(fd) * 0.5;  % Assume range ±fd/2

% Free-space path loss
Lfs = (4 * pi * d * fc / c).^2;
PL_dB = 10 * log10(Lfs);

% Log-logistic fading gain (in dB domain)
% Generate 100 samples and average for stable estimate
samples = 1000;
fading_linear = (1 ./ (1 + (randn(samples,1)/sigma).^2)).^(1/sigma); % log-logistic samples
fading_dB = 10 * log10(mean(fading_linear));

% Received power in dBm
Pr_dBm = Pt_dBm + Gt_dB + Gr_dB - PL_dB + fading_dB;

% Plot Doppler Shift
figure;
plot(theta_deg, fd/1e3, 'b', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Doppler Shift (kHz)');
title('Doppler Shift vs Elevation Angle');
grid on;

% Plot Doppler Spread
figure;
plot(theta_deg, doppler_spread/1e3, 'm', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Doppler Spread (kHz)');
title('Doppler Spread vs Elevation Angle (Log-Logistic)');
grid on;

% Plot Path Loss
figure;
plot(theta_deg, PL_dB, 'k', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Path Loss (dB)');
title('Free-Space Path Loss vs Elevation Angle');
grid on;

% Plot Received Power
figure;
plot(theta_deg, Pr_dBm, 'g', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Received Power (dBm)');
title('Received Power vs Elevation Angle (Log-Logistic Fading)');
grid on;
