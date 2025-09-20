clear;clc;

%   Rician Channel System Object
ricianchan = comm.RicianChannel(...
    "SampleRate",10e6,...                   
    "PathDelays",[6 8 12]*1e-3,...       % seconds
    "AveragePathGains",[0.1 0.2 0.5],... % dB
    KFactor=3,...
    DirectPathDopplerShift=5.0,...
    DirectPathInitialPhase=0,...
    MaximumDopplerShift=720e3,...
    DopplerSpectrum=doppler('Bell',8),...
    RandomStream='mt19937ar with seed',...
    Seed=73,...
    PathGainsOutputPort=true);

%   Modulation and transmitted signal parameters
M = 2;
k = log2(M);
nof_syms = 1000;
snr = -10:30;
berEst = zeros(1,length(snr));
nof_bits = nof_syms*k;
%   Simulation for various SNRs
for n=1:length(snr)
    total_bits = 0;
    total_errors = 0;
    while total_bits < 1e6
        dataIn = randi([0 1], nof_syms*k, 1);
        dataSym = bit2int(dataIn,k);
        txSig = pskmod(dataSym,M,0,"gray");
        [~, h] = ricianchan(txSig);
        path1 = txSig.*h(:,1);
        path2 = txSig.*h(:,2);
        path3 = txSig.*h(:,3);
        %   choose the index of the maximum h
        [~,idx] = max(abs(h).^2,[],2);
        chosen = (idx==1).*path1 + (idx==2).*path2 + (idx==3).*path3;
        channelOut = chosen;
        %   Add the thermal noise at the receiver
        %   Assume noise varies based on snr vector
        sigpower = mean(abs(txSig).^2);
        sigpowerdB = pow2db(sigpower);
        rxSig = awgn(channelOut,snr(n),sigpower);
        noise = rxSig - channelOut;
        noisePower = mean(abs(noise).^2);
        SNR_linear = sigpower/noisePower;
        SNR_db = 10*log10(SNR_linear);
        %   Coherent detection : full channel knowdedge
        %   Choose the maximum of the three channels
        h_detection = max(h,[],2);
        rxSig_coherent = rxSig./h_detection;
        rxSym = pskdemod(rxSig_coherent,M,0,"gray");
        dataOut = int2bit(rxSym,k);
        
        nof_errors = biterr(dataIn,dataOut);
        total_errors = total_errors + nof_errors;
        total_bits = nof_bits + total_bits;
    end
    
    berEst(n) = total_errors/total_bits;
end

berTheory = berfading(snr,'psk',M,3,3);

figure;
semilogy(snr,berEst,'r*');
hold on;
semilogy(snr,berTheory);
hold off;

hp = scatterplot(rxSig_coherent);
hold on;
scatterplot(txSig,[],[],'r*',hp);
grid;
hold off;