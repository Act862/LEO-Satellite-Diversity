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
M = 256;
k = log2(M);
K = 10;
KdB = 10*log10(K);
dataSym = bit2int(dataIn', k);
txSig = qammod(dataSym,M);
h1 = rice_fading(KdB,length(txSig),1)';
h2 = rice_fading(KdB,length(txSig),1)';
h3 = rice_fading(KdB,length(txSig),1)';
h4 = rice_fading(KdB,length(txSig),1)';

channelOut1 = h1.*txSig;
channelOut2 = h2.*txSig;
channelOut3 = h3.*txSig;
channelOut4 = h4.*txSig;

rxSig1 = awgn(channelOut1,20,"measured");
rxSig2 = awgn(channelOut2,20,"measured");
rxSig3 = awgn(channelOut3,20,"measured");
rxSig4 = awgn(channelOut4,20,"measured");

rxSig = [rxSig1 rxSig2 rxSig3 rxSig4];
h = [h1 h2 h3 h4];
r = selectionCombining(rxSig, h);
r2 = maximalRatioCombining(rxSig,h);
rxSym_sc = qamdemod(r,M);
rxSym_mrc = qamdemod(r2,M);
%   back to bits
dataOut_sc = int2bit(rxSym_sc,k);
dataOut_mrc = int2bit(rxSym_mrc,k);
biterr(dataIn,dataOut_sc')/length(dataIn)
biterr(dataIn,dataOut_mrc')/length(dataIn)

img_Linear = bit2int(dataOut_sc,8);
img_Linear_mrc = bit2int(dataOut_mrc,8);
img_reconstructed = reshape(img_Linear,[N N]);
img_reconstructed_mrc = reshape(img_Linear_mrc,[N N]);
subplot(1,2,1);
imshow(uint8(img_reconstructed));
subplot(1,2,2);
imshow(uint8(img_reconstructed_mrc));


%%  Helper Functions
function r = rice_fading(Kdb,N,Mi)
    K = 10^(Kdb/10);
    const = 1/(2*(K+1));
    x = randn(1,N);
    y = randn(1,N);
    r = sqrt(const*((x+sqrt(2*K)).^2 + y.^2));
    rt = zeros(1,Mi*length(r));
    ki = 1;
    for i=1:length(r)
        rt(ki:i*Mi) = r(i);
        ki = ki+Mi;
    end
    r = rt;
end

function r = selectionCombining(rxSig, h)
    % rxSig: [N x L]
    % h:     [N x L]
    [N, L] = size(h)
    r = zeros(N,1);
    for n = 1:N
        [~, idx] = max(abs(h(n,:)).^2);  % choose best branch
        r(n) = rxSig(n, idx) ./ h(n, idx); % equalize the chosen branch
    end
end

function r = maximalRatioCombining(rxSig,h)
    %   Generate precoding vector
    h_norm = zeros(1,size(rxSig,1));
    for n=1:size(h,1)
        s = 0;
        for m=1:size(h,2)
            s = s + abs(h(n,m)).^2;
        end
        h_norm(n) = 1/sqrt(s);
    end
    size(h_norm)
    w = h_norm'.*h;
    conj(w)
    %   Apply the combining
    r = sum(conj(w).*rxSig,2);
end