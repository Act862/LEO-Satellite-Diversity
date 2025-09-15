Fc = 2e9;
Fs = 1e6;
Rs = 100e3;
M = 4;
k = log2(M);
N = 10000;
% LEO Speed of satellite
v = 7.5e3;
maxDoppler = (v/3e8)*Fc;
passDuration = 0.6; % (600 ms)

samplesPerSymbol = Fs/Rs;
totalSamples = ceil(N/k*samplesPerSymbol);
t = (0:totalSamples-1)/Fs;
t = t(:);

dopplerShift = maxDoppler*cos(2*pi*t/passDuration);

dataIn = randi([0 1], N, 1);
dataSym = bi2de(reshape(dataIn, [], k));
modSignal = pskmod(dataSym, M, pi/4);

txSignal = upfirdn(modSignal, ones(samplesPerSymbol,1), samplesPerSymbol);
t = t(1:length(txSignal));

phaseOffset = 2*pi*cumsum(dopplerShift)/Fs;
dopplerSignal = txSignal.*exp(1j*phaseOffset);
EbNo = 10;
rxSignal = awgn(dopplerSignal, EbNo + 10*log10(k*samplesPerSymbol), 'measured');
rxCorrected = rxSignal .* exp(-1j*phaseOffset);
% rxCorrected = rxSignal;
rxDown = downsample(rxCorrected, samplesPerSymbol);
rxDataSym = pskdemod(rxDown, M, pi/4);
rxData = de2bi(rxDataSym, k);
rxData = rxData(:);

[numErrors, ber] = biterr(dataIn, rxData);
fprintf('BER with time-varying Doppler: %f (Errors: %d / %d\n', ber, numErrors, N);

figure;
plot(t, dopplerShift/1e3);
xlabel('Time (s)');
ylabel('Doppler Shift (kHz)');
title('Time-Varying Doppler Shift Profile');
grid on;

% Take samples before and after compensation
rxDown_raw = downsample(rxSignal, samplesPerSymbol);

% Plot constellation animation
figure;
for k = 1:10:length(rxDown)
    clf;
    subplot(1,2,1);
    plot(rxDown_raw(1:k), 'ro');
    title('Before Doppler Compensation');
    axis([-2 2 -2 2]);
    grid on;
    xlabel('In-phase'); ylabel('Quadrature');

    subplot(1,2,2);
    plot(rxDown(1:k), 'go');
    title('After Doppler Compensation');
    axis([-2 2 -2 2]);
    grid on;
    xlabel('In-phase'); ylabel('Quadrature');

    drawnow;
end