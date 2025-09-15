% --- Parameters ---
numSymbols = 200;
sps = 8;                 % Samples per symbol
rolloff = 0.35;          % RRC rolloff
span = 6;                % RRC filter span
Fs = 1e6;                % Sampling frequency
Ts = 1 / Fs;

% --- Generate QPSK symbols ---
data = randi([0 3], numSymbols, 1);
symbols = pskmod(data, 4, pi/4);

% --- Upsample ---
tx_raw = upsample(symbols, sps);  % No pulse shaping

% --- Pulse shaping filter (Root Raised Cosine) ---
rrc = rcosdesign(rolloff, span, sps, 'sqrt');
tx_shaped = conv(upsample(symbols, sps), rrc, 'same');

% --- Simulate channel: simple low-pass filter (band-limited channel) ---
channel_filter = fir1(64, 0.5);  % Simple LPF to simulate ISI
rx_raw = conv(tx_raw, channel_filter, 'same');
rx_shaped = conv(tx_shaped, channel_filter, 'same');

% --- Downsample and normalize ---
rx_raw_down = downsample(rx_raw, sps);
rx_shaped_down = downsample(rx_shaped, sps);
rx_raw_down = rx_raw_down(10:end-10);       % Trim edges
rx_shaped_down = rx_shaped_down(10:end-10); % Trim edges

% --- Plot Eye Diagrams ---
eyediagram(rx_raw(100:end-100), 2*sps);    % Before pulse shaping
figure;
eyediagram(rx_shaped(100:end-100), 2*sps); % After pulse shaping

% --- Plot Constellation (Decision Quality) ---
figure;
subplot(1,2,1);
plot(real(rx_raw_down), imag(rx_raw_down), 'bo');
title('Constellation: Without Pulse Shaping');
axis equal; grid on;

subplot(1,2,2);
plot(real(rx_shaped_down), imag(rx_shaped_down), 'r+');
title('Constellation: With Pulse Shaping');
axis equal; grid on;
