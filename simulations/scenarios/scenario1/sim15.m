% Constants
c = 3e8; % Speed of light (m/s)
R = 6371e3; % Earth radius in meters
h = 1200e3; % Satellite altitude in meters

% Elevation angles (in degrees)
elevation_deg = 10:5:90;
elevation_rad = deg2rad(elevation_deg);
num_angles = length(elevation_deg);

% Band parameters
bands = {'C', 'S', 'Ka'};
data_rates = [5e6, 2e6, 50e6]; % bits/sec
extra_band_latency = [100e-6, 80e-6, 150e-6]; % seconds

% Other delays (fixed)
wsn_to_vsat_delay = 5e-3; % seconds
vsat_processing_delay = 2e-3; % seconds
cloud_processing_delay = 10e-3; % seconds
packet_size_bits = 1024 * 8;

% Initialize delay matrix
total_delay = zeros(length(bands), num_angles);

for i = 1:num_angles
    theta = elevation_rad(i);
    
    % Slant range (approximation using spherical earth model)
    d = sqrt((R + h)^2 - (R * cos(theta))^2) - R * sin(theta); % one-way
    prop_delay = 2 * d / c; % up and down link
    
    for b = 1:length(bands)
        tx_delay = packet_size_bits / data_rates(b);
        total_delay(b, i) = wsn_to_vsat_delay + vsat_processing_delay + ...
                            tx_delay + prop_delay + cloud_processing_delay + ...
                            extra_band_latency(b);
    end
end

% Plotting
figure;
plot(elevation_deg, total_delay(1,:)*1e3, '-o', ...
     elevation_deg, total_delay(2,:)*1e3, '-s', ...
     elevation_deg, total_delay(3,:)*1e3, '-^', 'LineWidth', 1.5);
xlabel('Elevation Angle (degrees)');
ylabel('Total Delay (ms)');
legend(bands, 'Location', 'northeast');
title('End-to-End Delay vs. Elevation Angle (WSN–LEO–Cloud)');
grid on;
