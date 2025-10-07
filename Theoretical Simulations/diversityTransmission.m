clear;clc;
N = 256;
img = imread('Transmissions/cameraman.tif');
img = imresize(img,[N N]);
%   serialize the image
img_Linear = reshape(img,[1 N*N]);
%   turn the image data into bits
%   Assume 8bits per pixel value (image entry)
dataIn = reshape(int2bit(img_Linear,8),...
    [1 8*N*N]);
%   modulate data
M = 512;
modtype = "QAM";
k = log2(M);
K = 3;
L = 4;  % diversity order
KdB = 10*log10(K);
ebnrdB = -10:15;
nof_iterations = 1e7;
[dataIn,paddedBits] = padSequence(dataIn,k);
snrdB = convertSNR(ebnrdB,"ebno","snr","BitsPerSymbol",k);
dataSym = bit2int(dataIn', k);

if modtype == "PSK"
    txSig = pskmod(dataSym,M);
elseif modtype == "QAM"
    txSig = qammod(dataSym,M,'gray','UnitAveragePower',false);
end

bitErrorRate_sc = zeros(1,length(snrdB));
bitErrorRate_mrc = zeros(1,length(snrdB));
bitErrorRate_egc = zeros(1,length(snrdB));
for n=1:length(snrdB)
    nBits = 0;
    nErrs_sc = 0;
    nErrs_mrc = 0;
    nErrs_egc = 0;
    while nBits < nof_iterations
        h1 = rice_fading(KdB,length(txSig),1)';
        h2 = rice_fading(KdB,length(txSig),1)';
        h3 = rice_fading(KdB,length(txSig),1)';
        h4 = rice_fading(KdB,length(txSig),1)';

        channelOut1 = h1.*txSig;
        channelOut2 = h2.*txSig;
        channelOut3 = h3.*txSig;
        channelOut4 = h4.*txSig;

        rxSig1 = awgn(channelOut1,snrdB(n),"measured");
        rxSig2 = awgn(channelOut2,snrdB(n),"measured");
        rxSig3 = awgn(channelOut3,snrdB(n),"measured");
        rxSig4 = awgn(channelOut4,snrdB(n),"measured");

        rxSig = [rxSig1 rxSig2 rxSig3 rxSig4];
        h = [h1 h2 h3 h4];
        r = selectionCombining(rxSig, h);
        r2 = maximalRatioCombining(rxSig,h);
        r3 = equalGainCombining(rxSig, h);
        if modtype == "PSK"
            rxSym_sc = pskdemod(r,M);
            rxSym_mrc = pskdemod(r2,M);
            rxSym_egc = pskdemod(r3,M);
        elseif modtype == "QAM"
            rxSym_sc = qamdemod(r,M,"gray",'OutputType','integer');
            rxSym_mrc = qamdemod(r2,M,"gray",'OutputType','integer');
            rxSym_egc = qamdemod(r3,M,"gray",'OutputType','integer');
        end
        %   back to bits
        dataOut_sc = int2bit(rxSym_sc,k);
        dataOut_mrc = int2bit(rxSym_mrc,k);
        dataOut_egc = int2bit(rxSym_egc,k);
        errors_sc = biterr(dataIn,dataOut_sc');
        errors_mrc = biterr(dataIn,dataOut_mrc');
        errors_egc = biterr(dataIn,dataOut_egc');
        nBits = nBits + length(dataIn);
        nErrs_sc = nErrs_sc + errors_sc;
        nErrs_mrc = nErrs_mrc + errors_mrc;
        nErrs_egc = nErrs_egc + errors_egc;
    end
    bitErrorRate_sc(n) = nErrs_sc/nBits;
    bitErrorRate_mrc(n) = nErrs_mrc/nBits;
    bitErrorRate_egc(n) = nErrs_egc/nBits;
end
img_Linear = bit2int(dataOut_sc(1:length(dataOut_sc)-paddedBits),8);
img_Linear_mrc = bit2int(dataOut_mrc(1:length(dataOut_mrc)-paddedBits),8);
img_Linear_egc = bit2int(dataOut_egc(1:length(dataOut_egc)-paddedBits),8);
img_reconstructed = reshape(img_Linear,[N N]);
img_reconstructed_mrc = reshape(img_Linear_mrc,[N N]);
img_reconstructed_egc = reshape(img_Linear_egc,[N N]);

%%  Theoretical BER Calculation
berTheory = berfading(ebnrdB,lower(modtype),M,L,K);

%%  Results Demonstration
figure;
subplot(1,3,1);
imshow(uint8(img_reconstructed));
title('Received Image with SC');
subplot(1,3,2);
imshow(uint8(img_reconstructed_mrc));
title('Received Image with MRC');
subplot(1,3,3);
imshow(uint8(img_reconstructed_egc));
title('Received Image with EGC');

figure;
semilogy(ebnrdB,bitErrorRate_sc,'*','LineWidth',1.5);
hold on;
semilogy(ebnrdB,bitErrorRate_mrc,'*','LineWidth',1.5);
semilogy(ebnrdB,bitErrorRate_egc,'*','LineWidth',1.5);
semilogy(ebnrdB,berTheory);
hold off;
grid on;
legend('SC','MRC','EGC');
if M==2
    text = sprintf("BER for B%s",modtype);
elseif M==4
    text = sprintf("BER for Q%s",modtype);
else
    text = sprintf("BER for %d-%s",M,modtype);
end
title(text);

%%  Helper Functions

function [paddedSequence,paddedBits] = padSequence(seq,k)
%   pad the sequence to me of length multiple of k
N = length(seq);
paddedSequence = seq;
while mod(N,k)~=0
    paddedSequence = [paddedSequence 0];
    N = N + 1;
end
paddedBits = length(paddedSequence) - length(seq);
end

function r = rice_fading(Kdb,N,Mi)
K = 10^(Kdb/10);
mu = sqrt(K/(2*(K+1)));
sigma = sqrt(1/(2*(K+1)));
% const = 1/(2*(K+1));
x = randn(1,N);
y = randn(1,N);
% r = sqrt(const*((x+sqrt(2*K)).^2 + y.^2));
r = (sigma*x+mu) + 1i*(sigma*y); % complex mu = 0
rt = zeros(1,Mi*length(r));
ki = 1;
for i=1:length(r)
    %   Rician Fading channel samples
    rt(ki:i*Mi) = r(i);
    ki = ki+Mi;
end
r = rt;
r = r/sqrt(mean(abs(r).^2));
end

function r = selectionCombining(rxSig, h)
% rxSig: [N x L]
% h:     [N x L]
[N, L] = size(h);
r = zeros(N,1);
for n = 1:N
    [~, idx] = max(abs(h(n,:)).^2);  % choose best branch
    r(n) = rxSig(n, idx) ./ h(n, idx); % equalize the chosen branch
end
end

function r = maximalRatioCombining(rxSig, h)
% rxSig and h are [N x L] matrices
[N, L] = size(rxSig);
r = zeros(N,1);
for n = 1:N
    num = 0;
    den = 0;
    for l = 1:L
        num = num + conj(h(n,l)) * rxSig(n,l);
        den = den + abs(h(n,l))^2;
    end
    r(n) = num / den;   % proper normalization
end
end

function r = equalGainCombining(rxSig, h)
[N, L] = size(rxSig);
r = zeros(N,1);
for n = 1:N
    num = 0;
    den = 0;
    for l = 1:L
        num = num + exp(-1i*angle(h(n,l)))*rxSig(n,l);
        den = den + abs(h(n,l));
    end
    r(n) = num / den;   % proper normalization
end
end