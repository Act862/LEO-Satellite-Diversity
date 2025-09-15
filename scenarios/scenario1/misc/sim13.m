% Constants
fc = 30e9;                  % Carrier frequency (Hz)
c = 3e8;                    % Speed of light (m/s)
lambda = c / fc;            % Wavelength (m)
Re = 6371e3;                % Earth radius (m)
h_sat = 600e3;              % Satellite altitude (m)
Gt_dB = 35; Gr_dB = 35;
Gt = 10^(Gt_dB/10); Gr = 10^(Gr_dB/10);
Pt_dBm = 30;                % Transmit power (dBm)
Pt = 10^((Pt_dBm - 30)/10); % Transmit power in Watts

% Log-logistic fading parameters
sigma = 1.2;                % Scale parameter
samples = 1000;             % Fading samples per point
M = 4;                      % Number of receiver antennas (satellite)

% Elevation angle sweep
theta_deg = linspace(5, 90, 100);
theta_rad = deg2rad(theta_deg);

% Slant range function
slant_range = @(theta) sqrt((Re + h_sat).^2 - (Re * cos(theta)).^2) - Re * sin(theta);
d = slant_range(theta_rad);

% Path loss
Lfs = (4 * pi * d * fc / c).^2;
PL_dB = 10 * log10(Lfs);

% Doppler shift and spread
G = 6.67430e-11; M_E = 5.972e24;
w_sat = sqrt(G * M_E / (Re + h_sat)^3); % rad/s
v_rel = w_sat * Re * cos(theta_rad);    % relative velocity at surface
fd = (fc / c) * v_rel;                  % Doppler shift (Hz)
doppler_spread = abs(fd) * 0.5;         % Spread estimate

% Preallocate received power
Pr_MRC_dBm = zeros(size(theta_deg));

for i = 1:length(theta_deg)
    d_i = d(i);
    pr = Pt * Gt * Gr / ( (4 * pi * d_i * fc / c)^2 );

    % M fading branches (log-logistic)
    fading = zeros(M, samples);
    for m = 1:M
        % Log-logistic fading samples per antenna
        X = randn(samples,1);   % Standard normal
        fading(m,:) = (1 ./ (1 + (X / sigma).^2)).^(1/sigma);
    end

    % MRC combining: sum of powers
    total_gain = sum(fading, 1);
    Pr_samples = pr * total_gain;      % Received power samples (Watts)
    Pr_mean = mean(Pr_samples);        % Average over samples
    Pr_MRC_dBm(i) = 10 * log10(Pr_mean) + 30; % Convert to dBm
end

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
title('Doppler Spread vs Elevation Angle');
grid on;

% Plot Path Loss
figure;
plot(theta_deg, PL_dB, 'k', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Path Loss (dB)');
title('Free-Space Path Loss vs Elevation Angle');
grid on;

% Plot Received Power with Spatial Diversity
figure;
plot(theta_deg, Pr_MRC_dBm, 'g', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Received Power (dBm)');
title(sprintf('Received Power vs Elevation Angle (Log-Logistic, M = %d)', M));
grid on;
