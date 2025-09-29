Fc = 2e9;   %   Carrier Frequency (2 GHz)
Fd = 50e3;  %   Doppler Shift (50 kHz)
Fs = 1e6;   %   Sampling Frequency
Rs = 100e3; %   Symbol rate from the source
M = 4;      %   Modulation Order (QPSK)
k = log2(M);%   bits per symbol
N = 1024;  %   Number of bits

%   Generate random bits and modulate
dataIn = randi([0 1], N, 1);
dataSym = bi2de(reshape(dataIn, [], k));
%   offset by 45 degrees
modSignal = pskmod(dataSym, M, pi/4);

%   upsample
samplesPerSymbol = Fs/Rs;
txSignal = upfirdn(modSignal, ones(samplesPerSymbol, 1), samplesPerSymbol);

%   time vector
t = (0:length(txSignal)-1)/Fs;
subplot(1,2,1);
title('txSignal');
plot(t,txSignal);
xlabel('t');
ylabel('Px');
grid on;
plot(t,txSignal);
%   apply doppler shift
dopplerSignal = txSignal .* exp(1j*2*pi*Fd*t).';
subplot(1,2,2);
title('Doppler Signal');
plot(t,dopplerSignal);
xlabel('t');
ylabel('Px');
grid on;

%   Pass through the channel, add noise
EbNo = 10;
rxSignal = awgn(dopplerSignal, EbNo + 10*log10(k*samplesPerSymbol), 'measured');

rxCorrected = rxSignal .* exp(-1j*2*pi*Fd*t).';

%   Downsample and Demodulate
rxDown = downsample(rxCorrected, samplesPerSymbol);
rxDataSym = pskdemod(rxDown, M, pi/4);
rxData = de2bi(rxDataSym, k);
rxData = rxData(:);

%   BER calculation
[numErrors, ber] = biterr(dataIn, rxData);
fprintf('BER: %f, Errors: %d out of %d\n', ber, numErrors, N);
