% Constants
c = 3e8; % Speed of light (m/s)
h = 1200e3; % LEO altitude (m)
wsn_to_vsat_delay = 5e-3; % Assume 5ms from WSN to VSAT_1
vsat_processing_delay = 2e-3; % VSAT processing delay
cloud_processing_delay = 10e-3; % Delay in cloud server
packet_size_bytes = 1024;
packet_size_bits = packet_size_bytes * 8;

% Frequency bands
bands = {'C', 'S', 'Ka'};
frequencies = [4e9, 2.5e9, 20e9]; % Hz
bandwidths = [45e6, 20e6, 250e6]; % Hz
data_rates = [5e6, 2e6, 50e6]; % bits/sec
extra_band_latency = [100e-6, 80e-6, 150e-6]; % seconds

% Initialize delay storage
total_delay = zeros(1, length(bands));

for i = 1:length(bands)
    % Transmission delay (bits / rate)
    tx_delay = packet_size_bits / data_rates(i);
    
    % Propagation delay: VSAT_1 -> Sat + Sat -> VSAT_2
    prop_delay = 2 * (h / c); % up + downlink
    
    % Add up all delays
    total_delay(i) = wsn_to_vsat_delay + vsat_processing_delay + ...
                     tx_delay + prop_delay + cloud_processing_delay + ...
                     extra_band_latency(i);
    
    fprintf('%s Band Total Delay: %.2f ms\n', bands{i}, total_delay(i)*1e3);
end

% Plotting
figure;
bar(total_delay*1e3);
set(gca, 'XTickLabel', bands);
ylabel('Total Delay (ms)');
title('End-to-End Delay for Different Frequency Bands (LEO Relay)');
grid on;
