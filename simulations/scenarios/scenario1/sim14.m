%% --- System Parameters (Define these at the beginning of your script) ---
clear; clc; close all;

% General Parameters
dataRate = 1e6; % bps
symbolRate = dataRate; % For BPSK, symbol rate = bit rate
numBits = 1e5; % Number of bits to simulate for BER calculation
fc = 2.5e9; % Carrier frequency (2.5 GHz for S-band LEO, adjust as needed)
c = physconst('lightspeed'); % Speed of light

% VSAT 1 (Transmitter - Uplink)
TxPower_VSAT1 = 50; % dBm (convert to Watts for calculations)
Gt_VSAT1 = 30; % dBi (Antenna gain of VSAT 1)
TxPower_VSAT1_W = 10^((TxPower_VSAT1 - 30)/10); % Watts
numTxAntennas_VSAT1 = 1;

% LEO Satellite (Relay)
Gt_LEO_Uplink = 15; % dBi (LEO satellite receive antenna gain for uplink)
Gr_LEO_Downlink = 15; % dBi (LEO satellite transmit antenna gain for downlink)
TxPower_LEO = 20; % dBm (LEO satellite transmit power for downlink)
TxPower_LEO_W = 10^((TxPower_LEO - 30)/10); % Watts
NF_LEO = 3; % dB (Noise Figure of LEO receiver)
T_LEO = 290; % Kelvin (System Noise Temperature of LEO receiver)

% VSAT 2 (Receiver - Downlink)
Gr_VSAT2 = 30; % dBi (Antenna gain of VSAT 2)
numRxAntennas_VSAT2 = 4;
NF_VSAT2 = 5; % dB (Noise Figure of VSAT 2 receiver)
T_VSAT2 = 290; % Kelvin (System Noise Temperature of VSAT 2 receiver)

% Channel Parameters
K_rician = 10; % Rician K-factor (dB) for satellite links (10 dB is strong LoS)
K_rician_linear = 10^(K_rician/10);

% Foliage Attenuation Model (ITU-R P.833-9 simplified model for basic simulation)
% This is a highly simplified model. For more accuracy, consider specific
% tree types, density, and path length through foliage.
attenuation_per_meter_foliage = 0.5; % dB/meter (example value, varies significantly)
path_length_through_foliage = 100; % meters (example, adjust based on your scenario)
A_foliage = attenuation_per_meter_foliage * path_length_through_foliage; % Total foliage attenuation in dB

% Noise Calculations
k = 1.38e-23; % Boltzmann's constant
B = dataRate; % Bandwidth (assuming equal to data rate for simplicity)

% Calculate thermal noise power for LEO and VSAT2
N0_LEO = k * T_LEO * B; % Watts
N0_VSAT2 = k * T_VSAT2 * B; % Watts

%% --- LEO Satellite Orbit Simulation (Simplified) ---
% For a full simulation, you'd use SGP4 or similar models for orbital propagation.
% Here, we'll just use a simplified distance model.

altitude_LEO = 500e3; % meters (LEO altitude)
earth_radius = 6371e3; % meters

% Assume VSATs are at a fixed location relative to the ground track.
% Distance will vary based on elevation angle. For simplicity, let's assume
% a worst-case or average elevation angle.
elevation_angle_deg = 30; % degrees
elevation_angle_rad = deg2rad(elevation_angle_deg);

% Calculate slant range (simplified geometry)
% R = sqrt(Re^2 + (Re+h)^2 - 2*Re*(Re+h)*cos(90+El)) - using law of cosines
% A simpler approximation for reasonable elevation angles:
slant_range_uplink = altitude_LEO / sin(elevation_angle_rad); % meters
slant_range_downlink = slant_range_uplink; % Assuming similar geometry

%% --- Uplink Evaluation (VSAT 1 to LEO Satellite) ---

fprintf('--- Uplink (VSAT 1 to LEO Satellite) Evaluation ---\n');

% 1. Path Loss (Free Space Path Loss)
Lp_uplink_dB = 20*log10(4*pi*slant_range_uplink*fc/c); % dB

% 2. Total Channel Attenuation (includes foliage and fading)
% For simulation, we'll simulate many iterations to get average BER/SNR
num_iterations = 1000; % Number of channel realizations

SNR_uplink_linear_iterations = zeros(1, num_iterations);
BER_uplink_iterations = zeros(1, num_iterations);

for i = 1:num_iterations
    % Rician Fading Channel for Uplink
    % Generate Rician fading coefficients for 1 Tx, 1 Rx antenna
    % h = sqrt(K_rician_linear/(K_rician_linear+1)) * LoS_component + sqrt(1/(K_rician_linear+1)) * NLoS_component
    LoS_component = 1; % Assume perfect LoS for simplicity of direction
    NLoS_component = (randn + 1i*randn)/sqrt(2); % Rayleigh fading component
    h_uplink = sqrt(K_rician_linear/(K_rician_linear+1)) * LoS_component + ...
               sqrt(1/(K_rician_linear+1)) * NLoS_component;
    
    channel_gain_uplink_linear = abs(h_uplink)^2; % Power gain of fading channel
    channel_gain_uplink_dB = 10*log10(channel_gain_uplink_linear);

    % Total effective path loss and attenuation
    Total_Att_uplink_dB = Lp_uplink_dB + A_foliage - channel_gain_uplink_dB; % Note: channel_gain_dB is positive for gain, so subtract it
    Total_Att_uplink_linear = 10^(Total_Att_uplink_dB/10);

    % 3. Received Power at LEO Satellite
    Pr_LEO_W = TxPower_VSAT1_W * (10^(Gt_VSAT1/10)) * (10^(Gt_LEO_Uplink/10)) / Total_Att_uplink_linear; % Watts

    % 4. Calculate SNR at LEO Satellite
    SNR_LEO_linear = Pr_LEO_W / (N0_LEO * 10^(NF_LEO/10)); % N0 is thermal noise, NF is noise figure
    SNR_LEO_dB = 10*log10(SNR_LEO_linear);
    SNR_uplink_linear_iterations(i) = SNR_LEO_linear;

    % 5. BER Calculation (example for BPSK)
    % For BPSK, BER = 0.5 * erfc(sqrt(SNR))
    BER_uplink = 0.5 * erfc(sqrt(SNR_LEO_linear));
    BER_uplink_iterations(i) = BER_uplink;
end

Avg_SNR_uplink_dB = 10*log10(mean(SNR_uplink_linear_iterations));
Avg_BER_uplink = mean(BER_uplink_iterations);

fprintf('Average Uplink SNR (at LEO): %.2f dB\n', Avg_SNR_uplink_dB);
fprintf('Average Uplink BER (at LEO): %.2e\n', Avg_BER_uplink);

% Outage Probability Uplink
SNR_threshold_uplink_dB = 5; % dB (example threshold for acceptable performance)
outage_count_uplink = sum(SNR_uplink_linear_iterations < 10^(SNR_threshold_uplink_dB/10));
outage_probability_uplink = outage_count_uplink / num_iterations;
fprintf('Uplink Outage Probability (SNR < %.1f dB): %.4f\n', SNR_threshold_uplink_dB, outage_probability_uplink);


%% --- Downlink Evaluation (LEO Satellite to VSAT 2) ---

fprintf('\n--- Downlink (LEO Satellite to VSAT 2) Evaluation ---\n');

% 1. Path Loss (Free Space Path Loss)
Lp_downlink_dB = 20*log10(4*pi*slant_range_downlink*fc/c); % dB

SNR_downlink_linear_iterations = zeros(1, num_iterations);
BER_downlink_iterations = zeros(1, num_iterations);

for i = 1:num_iterations
    % Rician Fading Channel for Downlink (with 4 Rx antennas at VSAT 2)
    h_downlink = zeros(1, numRxAntennas_VSAT2); % Fading coefficients for each Rx antenna
    for ant = 1:numRxAntennas_VSAT2
        LoS_component_dl = 1;
        NLoS_component_dl = (randn + 1i*randn)/sqrt(2);
        h_downlink(ant) = sqrt(K_rician_linear/(K_rician_linear+1)) * LoS_component_dl + ...
                          sqrt(1/(K_rician_linear+1)) * NLoS_component_dl;
    end

    % 2. Maximal Ratio Combining (MRC) for 4 Rx antennas
    % The effective channel gain for MRC is the sum of squared magnitudes of individual channel gains.
    channel_gain_downlink_MRC_linear = sum(abs(h_downlink).^2);
    channel_gain_downlink_MRC_dB = 10*log10(channel_gain_downlink_MRC_linear);

    % Total effective path loss and attenuation
    Total_Att_downlink_dB = Lp_downlink_dB + A_foliage - channel_gain_downlink_MRC_dB;
    Total_Att_downlink_linear = 10^(Total_Att_downlink_dB/10);

    % 3. Received Power at VSAT 2
    Pr_VSAT2_W = TxPower_LEO_W * (10^(Gr_LEO_Downlink/10)) * (10^(Gr_VSAT2/10)) / Total_Att_downlink_linear;

    % 4. Calculate SNR at VSAT 2 (after MRC)
    SNR_VSAT2_linear = Pr_VSAT2_W / (N0_VSAT2 * 10^(NF_VSAT2/10));
    SNR_downlink_linear_iterations(i) = SNR_VSAT2_linear;

    % 5. BER Calculation (example for BPSK with MRC)
    % For BPSK with MRC, BER = (0.5 * erfc(sqrt(SNR))) with effective SNR
    % Effective SNR with MRC for BPSK is just the sum of individual SNRs if noise is i.i.d.
    % Here, we directly use the combined SNR.
    BER_downlink = 0.5 * erfc(sqrt(SNR_VSAT2_linear));
    BER_downlink_iterations(i) = BER_downlink;
end

Avg_SNR_downlink_dB = 10*log10(mean(SNR_downlink_linear_iterations));
Avg_BER_downlink = mean(BER_downlink_iterations);

fprintf('Average Downlink SNR (at VSAT 2): %.2f dB\n', Avg_SNR_downlink_dB);
fprintf('Average Downlink BER (at VSAT 2): %.2e\n', Avg_BER_downlink);

% Outage Probability Downlink
SNR_threshold_downlink_dB = 5; % dB
outage_count_downlink = sum(SNR_downlink_linear_iterations < 10^(SNR_threshold_downlink_dB/10));
outage_probability_downlink = outage_count_downlink / num_iterations;
fprintf('Downlink Outage Probability (SNR < %.1f dB): %.4f\n', SNR_threshold_downlink_dB, outage_probability_downlink);

%% --- Display Procedures (Conceptual) ---
fprintf('\n--- Uplink and Downlink Procedures (Conceptual Display) ---\n');

fprintf('Uplink Procedure (VSAT 1 to LEO Satellite):\n');
fprintf('  1. WSN data collected by VSAT 1.\n');
fprintf('  2. VSAT 1 modulates data (e.g., BPSK).\n');
fprintf('  3. Signal transmitted from single antenna, passes through wooded foliage and free space.\n');
fprintf('  4. LEO Satellite receives signal, affected by Rician fading and noise.\n');
fprintf('  5. SNR and BER evaluated at LEO receiver.\n');

fprintf('\nDownlink Procedure (LEO Satellite to VSAT 2):\n');
fprintf('  1. LEO Satellite processes/relays data.\n');
fprintf('  2. LEO Satellite modulates data and transmits.\n');
fprintf('  3. Signal passes through free space and wooded foliage near VSAT 2.\n');
fprintf('  4. VSAT 2 receives signal on 4 antennas.\n');
fprintf('  5. Maximal Ratio Combining (MRC) is applied to combine signals from 4 antennas, maximizing SNR.\n');
fprintf('  6. Demodulation and BER evaluation at VSAT 2.\n');

%% --- Visualization (Basic Histograms) ---
figure;
subplot(2,1,1);
histogram(10*log10(SNR_uplink_linear_iterations), 50);
title('Uplink SNR Distribution (VSAT 1 to LEO)');
xlabel('SNR (dB)');
ylabel('Frequency');
grid on;

subplot(2,1,2);
histogram(10*log10(SNR_downlink_linear_iterations), 50);
title('Downlink SNR Distribution (LEO to VSAT 2 - with MRC)');
xlabel('SNR (dB)');
ylabel('Frequency');
grid on;