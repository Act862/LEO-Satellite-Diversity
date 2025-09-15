%   QPSK Modulation
%   Assume a 2-branch system Rx, and one transmitted signal
M = 4;
k = log2(M);
EbNoVec = -10:30;
EbNoVec2 = -20:20;
K = 10;
numSymPerFrame = 100;
snrdB = convertSNR(EbNoVec,"ebno","snr",BitsPerSymbol=k);
snrdB2 = convertSNR(EbNoVec2, "ebno","snr",BitsPerSymbol=k);
berEst = zeros(size(EbNoVec));

for n=1:length(snrdB)
    numErrs = 0;
    numBits = 0;

    while numErrs < 200 && numBits < 1e7
        dataIn = randi([0 1], numSymPerFrame*k,1);
        dataSym = bit2int(dataIn,k);
        txSig = pskmod(dataSym,M);
        % h1 = getRicianCCoef(K, length(txSig));
        fadedTxSig1 = txSig;
        % h2 = getRicianCCoef(K, length(txSig));
        fadedTxSig2 = txSig;
        rxSig1 = awgn(fadedTxSig1,snrdB(n),'measured');
        rxSig2 = awgn(fadedTxSig2,snrdB2(n),'measured');

        rxSig = max(rxSig1,rxSig2);

        rxSym = pskdemod(rxSig,M);
        dataOut = int2bit(rxSym,k);
        nErrors = biterr(dataIn,dataOut);
        numErrs = numErrs + nErrors;
        numBits = numBits + numSymPerFrame*k;
    end
    berEst(n) = numErrs/numBits;
end

berTheory = berfading(EbNoVec,'psk',M,1,200);

semilogy(EbNoVec2,berEst,'*');
hold on;
semilogy(EbNoVec,berTheory);
grid on;
legend('Estimated BER','Theoretical BER')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')

%%  Function Definitions
function h = getRicianCCoef(K,N)
    hLOS = exp(1i*2*pi*rand(1,N));
    hNLOS = sqrt(1/2)*(randn(1,N) + 1i*randn(1,N));
    h = sqrt(K/(K + 1))*hLOS + sqrt(1/(K + 1))*hNLOS;
end