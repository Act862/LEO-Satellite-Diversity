%%  Create the data
%   create an image
% img = 255*ones(128,256);
img = imread("cameraman.tif");
% img(3,:) = 0;

%   Serialize the image
img_serialized = reshape(img,[1 size(img,1)*size(img,2)]);

%   Turn the data into bits
l = 8;
dataIn = reshape(int2bit(img_serialized,l),[1 size(img,1)*size(img,2)*l]);
M = 4;
k = log2(M);
ricianK = 10;
ricianK_dB = 10*log10(ricianK);
ebnr = -10:0.1:0;
snr = convertSNR(ebnr,"ebno","snr","BitsPerSymbol",k);
ber = zeros(1,length(snr));

for n=1:length(snr)
    numErrors = 0;
    numBits = 0;

    while numErrors < 200 && numBits < 1e6
        %   Enter the modulation system - QPSK
        %   pad to reach a power of two
        txSym = bit2int(dataIn',k);
        txSig = pskmod(txSym,M);
        h = rice_fading(ricianK_dB,length(txSig),1)';
        %   Transmit the signal
        %   Coherent detection with removing the channel coefficient
        rxSig = awgn(h.*txSig,snr(n),"measured")./h;
        rxSym = pskdemod(rxSig,M);
        dataOut = int2bit(rxSym,k);

        % BER
        err = biterr(dataIn',dataOut);
        numErrors = numErrors  +err;
        numBits = numBits + length(dataIn);
    end
    ber(n) = numErrors/numBits;
end
imgData = bit2int(dataOut,8);
receivedImage = reshape(imgData',[size(img,1) size(img,2)]);

ber_theory = berawgn(ebnr,"psk",M,'nondiff');
ber_theory_fading = berfading(ebnr,"psk",M,1,ricianK);
figure;
semilogy(ebnr,ber_theory);hold on;
semilogy(ebnr,ber_theory_fading);
semilogy(ebnr,ber,'*');hold off;
axis tight;
grid on;
text = sprintf('rician K = %d',ricianK);
legend('awgn',text,'simulation');

figure;
subplot(2,1,1);
imshow(img);
subplot(2,1,2);
imshow(uint8(receivedImage));