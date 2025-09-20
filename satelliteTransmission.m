clear; clc;

%% Channel definition: Ka-band LEO satellite
ricianchan = comm.RicianChannel(...
    "SampleRate",10e6,...                   % Symbol/sample rate (10 MHz)
    "PathDelays",[6 8 12]*1e-3,...          % Multipath delays [s]
    "AveragePathGains",[0.1 0.2 0.5],...    % Path average gains [dB]
    KFactor=3,...                           % Rician K-factor (LoS vs scatter)
    DirectPathDopplerShift=5.0,...          % Doppler on LoS
    DirectPathInitialPhase=0,...
    MaximumDopplerShift=720e3,...           % LEO Doppler (≈ 500–700 kHz @ 20 GHz, 7.5 km/s)
    DopplerSpectrum=doppler('Bell',8),...
    RandomStream='mt19937ar with seed',...
    Seed=73,...
    PathGainsOutputPort=true);

%% Modulation parameters
M = 4;                         % QPSK
k = log2(M);                   % Bits per symbol
nof_syms = 2000;                % Symbols per frame (longer for fading)
snrVec = -10:5:30;              % SNR range [dB]
berEst_SC = zeros(1,length(snrVec));
berEst_MRC = zeros(1,length(snrVec));
nof_bits = nof_syms*k;

%% Simulation
for n = 1:length(snrVec)
    total_bits = 0;
    total_errors_SC = 0;
    total_errors_MRC = 0;
    
    while total_bits < 1e6   % ensure enough errors for statistics
        % --- Bit generation
        dataIn = randi([0 1], nof_bits, 1);
        dataSym = bit2int(dataIn,k);       % group bits into integers
        txSig = pskmod(dataSym,M,0,"gray"); % QPSK modulation

        % --- Channel
        [fadedSig,h] = ricianchan(txSig);
        path1 = txSig .* h(:,1);
        path2 = txSig .* h(:,2);
        path3 = txSig .* h(:,3);

        % --- Add AWGN (measured wrt channel output power)
        rxSig = awgn(fadedSig, snrVec(n), 'measured');

        % =============================
        % Selection Combining (SC)
        % =============================
        [~,idx] = max(abs(h).^2,[],2);   % pick strongest path index
        chosen = (idx==1).*path1 + (idx==2).*path2 + (idx==3).*path3;
        chosen_h = (idx==1).*h(:,1) + (idx==2).*h(:,2) + (idx==3).*h(:,3);
        rxSig_SC = awgn(chosen, snrVec(n), 'measured');
        rxSig_SC_eq = rxSig_SC ./ chosen_h; % equalize
        rxSym_SC = pskdemod(rxSig_SC_eq,M,0,"gray");
        dataOut_SC = int2bit(rxSym_SC,k);
        total_errors_SC = total_errors_SC + biterr(dataIn,dataOut_SC);

        % =============================
        % Maximum Ratio Combining (MRC)
        % =============================
        % MRC combining (normalize by sum of channel powers)
        mrc_num = conj(h(:,1)).*path1 + conj(h(:,2)).*path2 + conj(h(:,3)).*path3;
        mrc_den = sum(abs(h).^2,2);
        rxSig_MRC = mrc_num ./ mrc_den;
        rxSig_MRC_noisy = awgn(rxSig_MRC, snrVec(n), 'measured');
        rxSym_MRC = pskdemod(rxSig_MRC_noisy,M,0,"gray");
        dataOut_MRC = int2bit(rxSym_MRC,k);
        total_errors_MRC = total_errors_MRC + biterr(dataIn,dataOut_MRC);

        % --- Update counter
        total_bits = total_bits + nof_bits;
    end
    
    berEst_SC(n) = total_errors_SC / total_bits;
    berEst_MRC(n) = total_errors_MRC / total_bits;
end

%% Theoretical reference (QPSK in Rician fading, no diversity)
berTheory = berfading(snrVec,'psk',M,1,3);  % single-branch Rician
%% Plot results
figure;
semilogy(snrVec,berEst_SC,'r*-','LineWidth',1.5); hold on;
semilogy(snrVec,berEst_MRC,'bo-','LineWidth',1.5);
semilogy(snrVec,berTheory,'k--','LineWidth',1.5);
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
legend('Selection Combining (SC)','Max Ratio Combining (MRC)','Theory (Rician, 3 branch)');
title('QPSK over Ka-band LEO Rician Channel with Diversity');
grid on;
