% --- Parameters ---
numSymbols = 100;        % Number of QPSK symbols
sps = 8;                 % Samples per symbol (upsampling factor)
rolloff = 0.35;          % RRC roll-off factor
span = 6;                % RRC filter span (in symbols)
Fs = 1e6;                % Sampling frequency for plotting

% --- Generate QPSK Symbols ---
data = randi([0 3], numSymbols, 1);
symbols = pskmod(data, 4, pi/4);  % QPSK with π/4 offset

% --- Upsample ---
upsampled = upsample(symbols, sps);  % Insert zeros

% --- RRC Filter Design ---
rrcFilter = rcosdesign(rolloff, span, sps, 'sqrt');

% --- Apply Pulse Shaping ---
shapedSignal = conv(upsampled, rrcFilter, 'same');

% --- Time Axes ---
t_original = (0:length(upsampled)-1)/Fs;
t_shaped = (0:length(shapedSignal)-1)/Fs;

% --- Frequency Axes ---
Nfft = 2048;
f = linspace(-Fs/2, Fs/2, Nfft);

spectrum_raw = fftshift(abs(fft(upsampled, Nfft)));
spectrum_shaped = fftshift(abs(fft(shapedSignal, Nfft)));

% --- Plot Time Domain ---
figure;
subplot(2,1,1);
plot(t_original*1e3, real(upsampled), 'b--', 'DisplayName', 'Unshaped');
hold on;
plot(t_shaped*1e3, real(shapedSignal), 'r', 'LineWidth', 1.5, 'DisplayName', 'Pulse Shaped');
xlabel('Time (ms)');
ylabel('Amplitude');
title('Time Domain: Unshaped vs Pulse Shaped Signal');
legend;
grid on;

% --- Plot Frequency Domain ---
subplot(2,1,2);
plot(f/1e3, 20*log10(spectrum_raw/max(spectrum_raw)), 'b--', 'DisplayName', 'Unshaped');
hold on;
plot(f/1e3, 20*log10(spectrum_shaped/max(spectrum_shaped)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Pulse Shaped');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Frequency Domain: Spectrum of Unshaped vs Pulse Shaped');
legend;
grid on;
ylim([-80 5]);
