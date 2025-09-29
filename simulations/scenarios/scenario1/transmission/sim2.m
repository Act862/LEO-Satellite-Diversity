clc; clear; close all;

%% === 1. Load Audio Sample ===
% You can replace this with 'load speech_dft.mat' or similar if using MAT
[audio, Fs] = audioread('handel.wav');  % Mono WAV file
audio = audio(:,1);                     % Ensure mono
soundsc(audio, Fs);                    % Play original audio
disp('Original audio playing...');

pause(length(audio)/Fs + 1);           % Wait for playback to finish

%% === 2. MAC Layer: Packetize ===
bits_per_sample = 8;
audio_q = uint8((audio + 1) * 127.5);         % Quantize to 8-bit unsigned
bitstream = reshape(de2bi(audio_q, bits_per_sample, 'left-msb').', [], 1);  % Serialize

packet_size = 256;                            % bits per MAC packet
num_packets = ceil(length(bitstream) / packet_size);
bitstream = [bitstream; zeros(num_packets * packet_size - length(bitstream), 1)];  % Padding

%% === 3. PHY Layer: BPSK Modulation with Header + Cyclic Prefix ===
header = [1 0 1 0];  % Simple known pattern for sync (e.g., preamble)
cp_len = 16;
mod_signal = [];

for i = 1:num_packets
    pkt = bitstream((i-1)*packet_size + (1:packet_size));
    pkt_full = [header.'; pkt];  % Add header
    symbols = 2*pkt_full - 1;    % BPSK: 0→-1, 1→+1

    % Add cyclic prefix
    cp = symbols(end-cp_len+1:end);
    symbols_cp = [cp; symbols];

    mod_signal = [mod_signal; symbols_cp];  % Concatenate to overall signal
end

%% === 4. Channel Model: Multipath + Doppler ===
fc = 2e9;              % Carrier frequency (2 GHz)
v = 30;                % Relative velocity (m/s)
c = 3e8;               % Speed of light
Ts = 1/Fs;
fD = (v/c) * fc;

% Multipath channel parameters
path_delays = [0 4 8]*Ts;               % in seconds
path_gains_dB = [0 -4 -8];
path_gains = 10.^(path_gains_dB/20);
num_paths = length(path_delays);
delay_samples = round(path_delays / Ts);

% Doppler: generate time-varying phase shifts
n = (0:length(mod_signal)-1).';
doppler = exp(1j * 2 * pi * fD * Ts * n);

% Apply multipath + Doppler
tx_signal = zeros(size(mod_signal));
for i = 1:num_paths
    delayed = [zeros(delay_samples(i),1); mod_signal(1:end-delay_samples(i))];
    tx_signal = tx_signal + double(path_gains(i)) * double(delayed) .* double(doppler);

end

%% === 5. Receiver: Remove CP, Demodulate, Reconstruct ===
rx_bits = [];

offset = length(header) + cp_len;
pkt_len = packet_size + length(header) + cp_len;

for i = 0:num_packets-1
    pkt_rx = tx_signal(i*pkt_len + (1:pkt_len));

    % Remove cyclic prefix
    symbols_rx = pkt_rx(cp_len+1:end);

    % Simple symbol decision (no channel estimation)
    bits_rx = real(symbols_rx) > 0;
    rx_bits = [rx_bits; bits_rx(length(header)+1:end)];
end

% BER
bitstream_trimmed = bitstream(1:length(rx_bits));
BER = sum(rx_bits ~= bitstream_trimmed) / length(rx_bits);
disp(['Bit Error Rate (BER): ' num2str(BER)]);

%% === 6. Reconstruct Audio ===
rx_bytes = reshape(rx_bits, bits_per_sample, []).';
audio_rx_q = bi2de(rx_bytes, 'left-msb');
audio_rx = double(audio_rx_q) / 127.5 - 1;
audio_rx = audio_rx(1:length(audio));  % Trim to original length

soundsc(audio_rx, Fs);
disp('Received audio playing...');
pause(length(audio_rx)/Fs + 1);

%% === 7. Plot ===
figure;
subplot(2,1,1);
plot(real(mod_signal)); title('Transmitted BPSK Waveform'); xlabel('Sample'); ylabel('Amplitude');
subplot(2,1,2);
spectrogram(tx_signal, 128, 120, 128, Fs, 'yaxis');
title('Spectrogram of Transmitted Signal with Doppler + Multipath');
