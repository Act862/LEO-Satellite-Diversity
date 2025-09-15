% Script with reinforcement learning adaptive combining (Q-learning)
clear; clc;

%% Simulation time settings
startTime = datetime(2025,09,11,12,01,0);
endTime = startTime + minutes(5);
sampleTime = 10;                    % seconds
t = startTime:seconds(sampleTime):endTime;

%% Scenario & satellites
sc = satelliteScenario(startTime,endTime,sampleTime);
sat = satellite(sc,"leoSatelliteConstellation.tle");

% Transmitters on satellites
satTxs = [];
gTxs = [];
systemLoss = 3;
dishDiameter = 0.5;
apertureEfficiency = 0.5;
for s = sc.Satellites
    gimbalTx = gimbal(s);
    gTxs = [gTxs gimbalTx];
    tx = transmitter(gimbalTx,Name="TxSat"+s.ID,Frequency=30e9,Power=20,BitRate=20,SystemLoss=systemLoss);
    gaussianAntenna(tx,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);
    satTxs = [satTxs tx];
end

% Ground station
gsLat = 37.9838; gsLon = 23.7275; gsAlt = 0.1; % km
gs = groundStation(sc,gsLat,gsLon,"Altitude",gsAlt,"MinElevationAngle",5,"Name","Base Station");
gainToNoiseTemperatureRatio = 5;

% Point sats at GS
for s = 1:length(sc.Satellites)
    pointAt(gTxs(s), gs);
end

% Build links
links = [];
for s = 1:length(sc.Satellites)
    gimbalRx_s = gimbal(gs);
    pointAt(gimbalRx_s, sat(s));
    rx_s = receiver(gimbalRx_s,Name="Rx_GS_Sat"+sat(s).ID, GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio, SystemLoss=systemLoss);
    lnk = link(satTxs(s), rx_s);
    links = [links lnk];
end

%% Output arrays
average_ebnr = zeros(1,length(t));
average_rp   = zeros(1,length(t));
snr_mrc      = zeros(1,length(t));
snr_sc       = zeros(1,length(t));
snr_adaptive = zeros(1,length(t));   % RL-driven scheme

%% Channel parameters
K = 5;                        
pilotSNR_dB = 10;              
Npilots = 4;                   
pilotSNR_lin = 10^(pilotSNR_dB/10);
sigma_e2 = 1/(pilotSNR_lin * Npilots);

%% --- Reinforcement Learning Setup ---
% State space = discretized uncertainty (U)
U_bins = linspace(0,1,10);   % 10 bins for uncertainty
nStates = length(U_bins);

% Actions: 1=MRC, 0=SC
nActions = 2;

% Initialize Q-table
Q = zeros(nStates,nActions);

% Hyperparameters
alpha = 0.4;       % learning rate
gammaRL = 0.9;     % discount factor
epsilon = 0.75;     % exploration rate

numEpisodes = 1000;   % training episodes

%% --- Training Phase (offline learning) ---
for ep = 1:numEpisodes
    for k = 1:length(t)
        % ----- Generate one state (uncertainty U) -----
        % For training, simulate one satellite channel only
        phi = 2*pi*rand;
        LOS = sqrt(K/(K+1));
        sigma = sqrt(1/(2*(K+1)));
        X = randn; Y = randn;
        h_complex = LOS*exp(1j*phi) + sigma*(X + 1j*Y);
        h_abs2 = abs(h_complex)^2;

        % Channel estimate
        e_real = sqrt(sigma_e2)*randn;
        e_imag = sqrt(sigma_e2)*randn;
        e_complex = (e_real + 1j*e_imag)/sqrt(2);
        h_hat = h_complex + e_complex;
        U = sigma_e2 / (abs(h_hat)^2 + eps);

        % Discretize state
        [~,iU] = min(abs(U_bins - U));

        % Choose action (epsilon-greedy)
        if rand < epsilon
            action = randi([0 1]); % explore
        else
            [~,action] = max(Q(iU,:)); % exploit
            action = action-1;
        end

        % --- Reward: which method gives higher SNR? ---
        ebno_nom_lin = 10^(10*rand/10); % random nominal Eb/N0
        ebno_faded_lin = ebno_nom_lin * h_abs2;
        snr_sc_sample = 10*log10(ebno_faded_lin);

        ebno_est_faded_lin = ebno_nom_lin * abs(h_hat)^2;
        snr_mrc_sample = 10*log10(ebno_est_faded_lin);

        if action==1
            reward = snr_mrc_sample; % reward = achieved SNR (dB)
        else
            reward = snr_sc_sample;
        end

        % Next state
        phi2 = 2*pi*rand;
        h2 = LOS*exp(1j*phi2) + sigma*(randn + 1j*randn);
        hhat2 = h2 + (sqrt(sigma_e2)/sqrt(2))*(randn+1j*randn);
        U2 = sigma_e2 / (abs(hhat2)^2 + eps);
        [~,iU2] = min(abs(U_bins - U2));

        % Update Q
        Q(iU,action+1) = Q(iU,action+1) + ...
            alpha * (reward + gammaRL * max(Q(iU2,:)) - Q(iU,action+1));
    end
end

%% --- Deployment Phase: Use learned Q for decisions ---
for k = 1:length(t)
    ebnr_sum_lin = 0;
    pr_sum_lin = 0;
    count_non_zeros_ebnr = 0;
    count_non_zeros_pr = 0;
    numSat_underService = 0;
    max_ebnr_lin = 0;

    ebno_faded_lin_vec = zeros(1,length(sc.Satellites));
    h_vec = zeros(1,length(sc.Satellites));
    hhat_vec = zeros(1,length(sc.Satellites));
    U_vec = zeros(1,length(sc.Satellites));
    pr_lin_vec = zeros(1,length(sc.Satellites));

    for s = 1:length(sc.Satellites)
        % Fading
        phi = 2*pi*rand;
        LOS = sqrt(K/(K+1));
        sigma = sqrt(1/(2*(K+1)));
        X = randn; Y = randn;
        h_complex = LOS*exp(1j*phi) + sigma*(X + 1j*Y);
        h_abs2 = abs(h_complex)^2;

        ebno_nom_dB = ebno(links(s), t(k));
        if ebno_nom_dB == -Inf || isnan(ebno_nom_dB)
            ebno_nom_lin = 0;
        else
            ebno_nom_lin = 10^(ebno_nom_dB/10);
        end

        ebno_faded_lin = ebno_nom_lin * h_abs2;

        % Estimation
        e_complex = (sqrt(sigma_e2)/sqrt(2))*(randn+1j*randn);
        h_hat = h_complex + e_complex;
        U = sigma_e2 / (abs(h_hat)^2 + eps);

        [~, pr_dB] = sigstrength(links(s), t(k));
        if pr_dB == -Inf || isnan(pr_dB)
            pr_lin = 0;
        else
            pr_lin = 10^(pr_dB/10);
        end

        ebno_faded_lin_vec(s) = ebno_faded_lin;
        h_vec(s) = h_abs2;
        hhat_vec(s) = h_hat;
        U_vec(s) = U;
        pr_lin_vec(s) = pr_lin;

        if ebno_faded_lin > 0, count_non_zeros_ebnr = count_non_zeros_ebnr + 1; numSat_underService=numSat_underService+1; end
        if pr_lin > 0, count_non_zeros_pr = count_non_zeros_pr + 1; end

        ebnr_sum_lin = ebnr_sum_lin + ebno_faded_lin;
        pr_sum_lin = pr_sum_lin + pr_lin;
        max_ebnr_lin = max(max_ebnr_lin, ebno_faded_lin);
    end

    % True MRC
    if ebnr_sum_lin > 0, snr_mrc(k) = 10*log10(ebnr_sum_lin); else, snr_mrc(k)=-Inf; end
    if max_ebnr_lin > 0, snr_sc(k) = 10*log10(max_ebnr_lin); else, snr_sc(k)=-Inf; end

    % Imperfect MRC
    ebno_est_faded_lin = ebno_faded_lin_vec ./ (h_vec+eps) .* abs(hhat_vec).^2;
    imperfect_mrc_lin = sum(ebno_est_faded_lin);
    if imperfect_mrc_lin > 0, imperfect_mrc_dB = 10*log10(imperfect_mrc_lin); else, imperfect_mrc_dB=-Inf; end

    % --- RL Decision ---
    U_global = max(U_vec);
    [~,iU] = min(abs(U_bins - U_global));
    [~,action] = max(Q(iU,:));
    action = action-1;

    if action==1
        snr_adaptive(k) = imperfect_mrc_dB;
    else
        snr_adaptive(k) = snr_sc(k);
    end

    % Averages
    if count_non_zeros_ebnr > 0
        avg_ebnr_lin = sum(ebno_faded_lin_vec) / count_non_zeros_ebnr;
        average_ebnr(k) = 10*log10(max(avg_ebnr_lin, eps));
    else
        average_ebnr(k) = -Inf;
    end
    if count_non_zeros_pr > 0
        average_rp_lin = sum(pr_lin_vec) / count_non_zeros_pr;
        average_rp(k) = 10*log10(max(average_rp_lin, eps));
    else
        average_rp(k) = -Inf;
    end

    disp("Satellites in service at "+string(t(k))+" : "+numSat_underService);
end

%% Plots
figure;
subplot(1,2,1);
plot(t,10.^(average_ebnr./10),'LineWidth',2);
title('average Eb/N_0 vs time'); ylabel('Eb/N_0 (linear)'); xlabel('time'); grid on;
subplot(1,2,2);
plot(t,10.^(average_rp./10),'LineWidth',2);
title('average receiver power vs time'); ylabel('P (linear)'); xlabel('time'); grid on;

figure;
semilogy(t, 10.^(snr_mrc./10),'LineWidth',1.5); hold on;
semilogy(t, 10.^(snr_sc./10),'LineWidth',1.5);
semilogy(t, 10.^(snr_adaptive./10),'LineWidth',1.5); hold off;
legend('MRC (true)','SC','RL-Adaptive');
title('MRC vs SC vs RL Adaptive (linear)'); xlabel('time'); ylabel('Eb/N_0 (linear)'); grid on;

figure;
plot(t, snr_mrc,'-','LineWidth',1.2); hold on;
plot(t, snr_sc,'--','LineWidth',1.2);
plot(t, snr_adaptive,':','LineWidth',1.4); hold off;
legend('MRC (true)','SC','RL-Adaptive');
title('MRC vs SC vs RL Adaptive (dB)'); xlabel('time'); ylabel('Eb/N_0 (dB)'); grid on;
