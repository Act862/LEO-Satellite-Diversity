%% ============================================================
% SC / MRC / GAC with Time-Correlated Rician Fading
clear; clc; rng(42);

%% ================= PARAMETERS ===============================
N_branches   = 4;
N_symbols    = 1e4;
modSchemes   = {'BPSK','QPSK','16QAM'};
avgSNRdB    = [10 8 6 4];
avgSNR_lin  = 10.^(avgSNRdB/10);

% GAC geometry weighting
slantRanges = [800 1000 1200 1500]; % km
Rmax = max(slantRanges);
geoFactor = Rmax ./ slantRanges;

% Rician parameters
K_factor = 5; 
Ts = 1e-3;          % symbol duration [s]
Tc = 0.01;          % coherence time [s]
rho = exp(-Ts/Tc);  % correlation coefficient

%% ================= MAIN LOOP OVER MODULATIONS =================
for mIdx = 1:length(modSchemes)
    modScheme = modSchemes{mIdx};
    
    % --- Modulation
    switch modScheme
        case 'BPSK'
            M = 2; x = randi([0 M-1],1,N_symbols); s = 2*x-1;
        case 'QPSK'
            M = 4; x = randi([0 M-1],1,N_symbols); s = exp(1j*pi/2*x)/sqrt(2);
        case '16QAM'
            M = 16; x = randi([0 M-1],1, N_symbols);
            a = [-3 -1 1 3]; s = (a(mod(x,4)+1) + 1j*a(floor(x/4)+1))/sqrt(10);
    end
    
    % --- Initialize fading and storage
    h_prev = sqrt(K_factor/(K_factor+1))*exp(1j*2*pi*rand(1,N_branches));
    SNR_SC_time  = zeros(1,N_symbols); SNR_MRC_time = zeros(1,N_symbols); SNR_GAC_time = zeros(1,N_symbols);
    BER_SC = zeros(1,N_symbols); BER_MRC = zeros(1,N_symbols); BER_GAC = zeros(1,N_symbols);
    
    for n = 1:N_symbols
        % --- Time-correlated Rician fading per branch
        w = sqrt(1/(K_factor+1))*(randn(1,N_branches)+1j*randn(1,N_branches))/sqrt(2);
        h = rho*h_prev + sqrt(1-rho^2)*(sqrt(K_factor/(K_factor+1))*exp(1j*2*pi*rand(1,N_branches)) + w);
        h_prev = h;
        
        % --- AWGN
        n_vec = (randn(1,N_branches)+1j*randn(1,N_branches))/sqrt(2);
        
        % --- Received symbols
        y = sqrt(avgSNR_lin).*h.*s(n) + n_vec;
        
        % --- Instantaneous SNR per branch
        gamma = avgSNR_lin .* abs(h).^2;
        
        % --- SC
        [gamma_SC, idx_SC] = max(gamma); SNR_SC_time(n)=10*log10(gamma_SC); y_SC=y(idx_SC);
        % --- MRC
        w_MRC = conj(h); y_MRC=sum(w_MRC.*y); gamma_MRC=sum(gamma); SNR_MRC_time(n)=10*log10(gamma_MRC);
        % --- GAC
        w_GAC = conj(h).*geoFactor; y_GAC=sum(w_GAC.*y); gamma_GAC=sum(gamma.*geoFactor); SNR_GAC_time(n)=10*log10(gamma_GAC);
        
        % --- Demodulation & BER
        switch modScheme
            case 'BPSK'
                BER_SC(n)=real(y_SC)>0 ~= x(n);
                BER_MRC(n)=real(y_MRC)>0 ~= x(n);
                BER_GAC(n)=real(y_GAC)>0 ~= x(n);
            case 'QPSK'
                BER_SC(n)=mod(round(2*angle(y_SC)/pi),4)~=x(n);
                BER_MRC(n)=mod(round(2*angle(y_MRC)/pi),4)~=x(n);
                BER_GAC(n)=mod(round(2*angle(y_GAC)/pi),4)~=x(n);
            case '16QAM'
                I=[-3 -1 1 3]; Q=[-3 -1 1 3];
                [~,idxI_SC]=min(abs(real(y_SC)-I)); [~,idxQ_SC]=min(abs(imag(y_SC)-Q)); x_hat_SC=idxQ_SC*4+idxI_SC-4;
                [~,idxI_MRC]=min(abs(real(y_MRC)-I)); [~,idxQ_MRC]=min(abs(imag(y_MRC)-Q)); x_hat_MRC=idxQ_MRC*4+idxI_MRC-4;
                [~,idxI_GAC]=min(abs(real(y_GAC)-I)); [~,idxQ_GAC]=min(abs(imag(y_GAC)-Q)); x_hat_GAC=idxQ_GAC*4+idxI_GAC-4;
                BER_SC(n)=x_hat_SC~=x(n); BER_MRC(n)=x_hat_MRC~=x(n); BER_GAC(n)=x_hat_GAC~=x(n);
        end
    end
    
    %% =============== PLOTS =================================
    figure('Name',modScheme);
    timeVec=1:N_symbols;
    
    subplot(3,1,1); hold on; grid on;
    plot(timeVec,SNR_SC_time,'r'); plot(timeVec,SNR_MRC_time,'b'); plot(timeVec,SNR_GAC_time,'g');
    ylabel('SNR (dB)'); title([modScheme ': Post-Combining SNR vs Time']); legend('SC','MRC','GAC');
    
    subplot(3,1,2); hold on; grid on;
    plot(timeVec,BER_SC,'r'); plot(timeVec,BER_MRC,'b'); plot(timeVec,BER_GAC,'g');
    ylabel('Instantaneous BER'); set(gca,'YScale','log'); legend('SC','MRC','GAC'); title([modScheme ': Instantaneous BER vs Time']);
    
    subplot(3,1,3); hold on; grid on;
    % --- Outage probability for different SNR thresholds
    gamma_th_dB = 0:2:20; % threshold in dB
    outage_SC  = arrayfun(@(th) mean(SNR_SC_time<th),gamma_th_dB);
    outage_MRC = arrayfun(@(th) mean(SNR_MRC_time<th),gamma_th_dB);
    outage_GAC = arrayfun(@(th) mean(SNR_GAC_time<th),gamma_th_dB);
    plot(gamma_th_dB,outage_SC,'r-o'); plot(gamma_th_dB,outage_MRC,'b-s'); plot(gamma_th_dB,outage_GAC,'g-^');
    xlabel('SNR Threshold (dB)'); ylabel('Outage Probability'); title([modScheme ': Outage vs Threshold']);
    legend('SC','MRC','GAC'); grid on;
end
