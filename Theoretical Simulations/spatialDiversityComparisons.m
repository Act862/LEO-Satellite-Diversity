%   Selection Combining Example
%   Calculating BER under Rician Fading
%   Calculating Outage under Rician Fading

clear; clc;

%% Modulation parameters
M = 16;                          % M-PSK
K = 3;
k = log2(M);                    % bits per symbol
ebnr = -20 + (10 - (-20)) * rand(1,100);   % random Eb/No samples
numSymsPerFrame = 100;
snrdB = convertSNR(ebnr,"ebno","snr",BitsPerSymbol=k);

berEst = zeros(1,length(snrdB));
berEst_mrc = zeros(1,length(snrdB));
%% Channel definitions: two independent Rician fading branches
ricianchan1 = comm.RicianChannel( ...
    SampleRate=1, ...
    KFactor=K, ...
    AveragePathGains=0,...
    PathDelays=0,...
    DirectPathDopplerShift=0,...
    DirectPathInitialPhase=0,...
    NormalizePathGains=true,...
    MaximumDopplerShift=0, ...
    RandomStream='mt19937ar with seed', ...
    Seed=randi(1e6), ...
    PathGainsOutputPort=true);

ricianchan2 = comm.RicianChannel( ...
    SampleRate=1, ...
    KFactor=K, ...
    AveragePathGains=0,...
    PathDelays=0,...
    DirectPathDopplerShift=0,...
    DirectPathInitialPhase=0,...
    NormalizePathGains=true,...
    MaximumDopplerShift=0, ...
    RandomStream='mt19937ar with seed', ...
    Seed=randi(1e6), ...
    PathGainsOutputPort=true);


%% Simulation loop
for n = 1:length(snrdB)
    numErrs_sc = 0;
    numBits = 0;
    numErrs_mrc = 0;
    while numBits < 1e6
        % Generate random bits
        dataIn = randi([0 1], numSymsPerFrame*k, 1);

        % Map to symbols
        dataSym = bit2int(dataIn,k);
        txSig = qammod(dataSym,M);

        % Pass through two independent Rician channels
        [rx1, h1] = ricianchan1(txSig);
        [rx2, h2] = ricianchan2(txSig);
        
        %   selection signal combining
        selSig = zeros(1,length(rx1))';
        for s=1:length(rx1)
            if abs(h1(s))^2 > abs(h2(s))^2
                selSig(s) = rx1(s);
            else
                selSig(s) = rx2(s);
            end
        end
        selH = max(abs(h1).^2,abs(h2).^2);
        
        %   maximal ratio combining
        w1 = (1./sqrt(abs(h1).^2 + abs(h2).^2)).*h1;
        w2 = (1./sqrt(abs(h1).^2 + abs(h2).^2)).*h2;
        mrcSig = conj(w1).*rx1 + conj(w2).*rx2;
        
        % Add AWGN after selection
        rxSig = awgn(selSig, snrdB(n), 'measured');
        rxSig_mrc = awgn(mrcSig,snrdB(n),'measured');
        % Demodulate
        rxSym = qamdemod(rxSig, M);
        rxSym_mrc = qamdemod(rxSig_mrc,M);
        % Turn to bits
        dataOut = int2bit(rxSym, k);
        dataOut_mrc = int2bit(rxSym_mrc,k);
        
        % Count errors
        numErrors = biterr(dataIn, dataOut);
        numErrors_mrc = biterr(dataIn,dataOut_mrc);
        numErrs_sc = numErrs_sc + numErrors;
        numErrs_mrc = numErrs_mrc + numErrors_mrc;
        numBits = numBits + numSymsPerFrame*k;
    end
    berEst(n) = numErrs_sc/numBits;
    berEst_mrc(n) = numErrs_mrc/numBits;
end

%% Theoretical curve
ebnr2 = -20:10;
berTheory_noDiv = berfading(ebnr2,'qam',M,1,K);
berTheory = berfading(ebnr2,'qam',M,2,K);

%% Plot results
semilogy(ebnr, berEst, '*');
hold on;
semilogy(ebnr,berEst_mrc,'*');
semilogy(ebnr2, berTheory);
semilogy(ebnr2, berTheory_noDiv);
hold off;
grid on;
legend('Estimated BER (SC)','Estimated BER (MRC)','Theoretical BER', 'Theoretical BER (no diversity)');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate');
title('Selection and Maximal Ratio Combining over Rician Fading (16-QAM)');
