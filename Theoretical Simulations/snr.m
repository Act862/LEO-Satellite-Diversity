clear;clc;
symbolsPerFrame = 400;
M = 4;
k = log2(M);
%   Generate bits for transmission
bits = randi([0 1], symbolsPerFrame*k,1);
txSyms = bit2int(bits,k);
txSig = pskmod(txSyms,M);
%   Transmit the symbols
[rxSig,No] = awgn(txSig,8,"measured");
%   demodulate
rxSyms = pskdemod(rxSig,M);
rxBits = int2bit(rxSyms,k);
errors = biterr(rxBits,bits);
ber = errors/(symbolsPerFrame*k)

%   check the SNR
Ptx = mean(abs(txSig).^2);
snr_watts = Ptx/No;
snr_dB = pow2db(snr_watts)

