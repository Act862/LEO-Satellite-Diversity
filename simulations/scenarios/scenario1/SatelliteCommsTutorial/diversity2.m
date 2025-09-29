%   Spatial Diversity
clear;clc;
Fs = 1e6;
bits = randi([0 1],10000,1);
tx = pskmod(bi2de(reshape(bits,[],2)),4,pi/4,'gray');
% tx = upsample(symbols,Fs/1e5); % 100 kSymbols/sec
ricianchan = comm.RicianChannel( ...
    SampleRate=1e6, ...
    PathDelays=[0.0 0.6 1.2]*1e-6, ...
    AveragePathGains=[0.1 0.5 0.2], ...
    KFactor=2.8, ...
    DirectPathDopplerShift=5e3, ...
    DirectPathInitialPhase=0.5, ...
    MaximumDopplerShift=50e3, ...
    DopplerSpectrum=doppler('Bell',0.5), ...
    RandomStream='mt19937ar with seed', ...
    Seed=73, ...
    PathGainsOutputPort=true);

% Info=info(ricianchan)
% ChanCoeff=Info.ChannelFilterCoefficients
rxSym_mc = zeros(1,length(tx));
rxSym_sc = zeros(1,length(tx));
n0 = 1/2;
sigma = sqrt(n0/2);
for i=1:length(tx)
    % send a symbol from the transmit signal
    [~,h] = ricianchan(tx(i));
    h_norm = sqrt(sum(abs(h).^2));
    w = h_norm.*h';
    % rx1 = h(1)*tx(i)*exp(-2*1i*pi*300)+n0;
    % rx2 = h(2)*tx(i)*exp(-2*1i*pi*450)+n0;
    % rx3 = h(3)*tx(i)*exp(-2*1i*pi*20)+n0;
    n1 = sigma * (randn(1,1) + 1j*randn(1,1));
    n2 = sigma * (randn(1,1) + 1j*randn(1,1));
    n3 = sigma * (randn(1,1) + 1j*randn(1,1));
    rx1 = h(1)*tx(i)+n1;
    rx2 = h(2)*tx(i)+n2;
    rx3 = h(3)*tx(i)+n3;
    w_h = conj(w)';
    rx_combined_mrc = w_h*[rx1;rx2;rx3];
    SNR1 = abs(rx1).^2;
    SNR2 = abs(rx2).^2;
    SNR3 = abs(rx3).^2;
    [~,bestIndex] = max([SNR1 SNR2 SNR3],[],2);
    switch bestIndex
        case 1, rx_combined_rc = rx1;
        case 2, rx_combined_rc = rx2;
        case 3, rx_combined_rc = rx3;
    end
    % rx = downsample(rx_combined,Fs/1e5);
    rxSym_mc(i) = pskdemod(rx_combined_mrc,4,pi/4,'gray');
    rxSym_sc(i) = pskdemod(rx_combined_rc,4,pi/4,'gray');
end
rxSym_mc = rxSym_mc';
rxSym_sc = rxSym_sc';
rxBits1 = reshape(de2bi(rxSym_mc,2),[],1);
rxBits2 = reshape(de2bi(rxSym_sc,2),[],1);
bits = bits(1:length(rxBits1));
[errs,ber] = biterr(bits,rxBits1);
fprintf("BER with Spatial Diversity (MRC): %.5f (%d errors)\n", ber,errs);
[errs,ber] = biterr(bits,rxBits2);
fprintf("BER with Spatial Diversity (SC): %.5f (%d errors)\n", ber,errs);
