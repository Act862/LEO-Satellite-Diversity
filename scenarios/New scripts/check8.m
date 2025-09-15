% Script with reinforcement learning adaptive combining (Threshold Q-learning)
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
uncertainty  = zeros(1,length(t));   % global uncertainty
rl_action    = zeros(1,length(t));   % 1=MRC, 0=SC

%% Channel parameters
K = 5;                        
pilotSNR_dB = 10;              
Npilots = 4;                   
pilotSNR_lin = 10^(pilotSNR_dB/10);
sigma_e2 = 1/(pilotSNR_lin * Npilots);

%% --- Reinforcement Learning Setup (Threshold Tuning) ---
nBins = 21;                           % discretize threshold [0,1]
thresholds = linspace(0,1,nBins);
actions = [-1 0 1];                   % shift threshold left, stay, right
nActions = numel(actions);

alpha = 0.1; gamma = 0.95;
epsilon = 1.0; epsilonMin = 0.05; epsDecay = 0.999;

Q = zeros(nBins,nActions);
numEpisodes = 500; stepsPerEp = 50;

% --- Training phase (offline) ---
for ep = 1:numEpisodes
    state = randi(nBins);             % start threshold index
    prevScheme = '';
    for step = 1:stepsPerEp
        % Epsilon-greedy
        if rand < epsilon
            aIdx = randi(nActions);
        else
            [~,aIdx] = max(Q(state,:));
        end
        action = actions(aIdx);
        nextState = max(1,min(nBins,state+action));
        th = thresholds(nextState);

        % --- Simulate channel
        phi = 2*pi*rand;
        LOS = sqrt(K/(K+1));
        sigma = sqrt(1/(2*(K+1)));
        h = LOS*exp(1j*phi) + sigma*(randn+1j*randn);
        h_hat = h + (sqrt(sigma_e2)/sqrt(2))*(randn+1j*randn);
        U = sigma_e2/(abs(h_hat)^2 + eps);

        % Decide SC/MRC
        if U < th
            scheme = 'MRC';
            reward = 10*log10(abs(h_hat)^2) - U;  % reward ~ SNR - penalty
        else
            scheme = 'SC';
            reward = 10*log10(abs(h)^2) - 0.5*U;
        end

        % Switching penalty
        if ~isempty(prevScheme) && ~strcmp(scheme,prevScheme)
            reward = reward - 0.5;
        end

        % Q-update
        Q(state,aIdx) = Q(state,aIdx) + ...
            alpha*(reward + gamma*max(Q(nextState,:)) - Q(state,aIdx));

        state = nextState; prevScheme = scheme;
    end
    epsilon = max(epsilon*epsDecay,epsilonMin);
end

%% --- Deployment phase ---
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
        phi = 2*pi*rand;
        LOS = sqrt(K/(K+1));
        sigma = sqrt(1/(2*(K+1)));
        h = LOS*exp(1j*phi) + sigma*(randn+1j*randn);
        h_abs2 = abs(h)^2;

        ebno_nom_dB = ebno(links(s), t(k));
        if ebno_nom_dB==-Inf || isnan(ebno_nom_dB)
            ebno_nom_lin=0;
        else
            ebno_nom_lin=10^(ebno_nom_dB/10);
        end
        ebno_faded_lin = ebno_nom_lin * h_abs2;

        % Estimation
        h_hat = h + (sqrt(sigma_e2)/sqrt(2))*(randn+1j*randn);
        U = sigma_e2/(abs(h_hat)^2+eps);

        [~, pr_dB] = sigstrength(links(s), t(k));
        if pr_dB==-Inf || isnan(pr_dB)
            pr_lin=0;
        else
            pr_lin=10^(pr_dB/10);
        end

        ebno_faded_lin_vec(s)=ebno_faded_lin;
        h_vec(s)=h_abs2;
        hhat_vec(s)=h_hat;
        U_vec(s)=U;
        pr_lin_vec(s)=pr_lin;

        if ebno_faded_lin>0, count_non_zeros_ebnr=count_non_zeros_ebnr+1; numSat_underService=numSat_underService+1; end
        if pr_lin>0, count_non_zeros_pr=count_non_zeros_pr+1; end
        ebnr_sum_lin=ebnr_sum_lin+ebno_faded_lin;
        pr_sum_lin=pr_sum_lin+pr_lin;
        max_ebnr_lin=max(max_ebnr_lin,ebno_faded_lin);
    end

    % True MRC
    if ebnr_sum_lin>0, snr_mrc(k)=10*log10(ebnr_sum_lin); else, snr_mrc(k)=-Inf; end
    if max_ebnr_lin>0, snr_sc(k)=10*log10(max_ebnr_lin); else, snr_sc(k)=-Inf; end

    % Imperfect MRC
    ebno_est_faded_lin = ebno_faded_lin_vec ./ (h_vec+eps) .* abs(hhat_vec).^2;
    imperfect_mrc_lin = sum(ebno_est_faded_lin);
    if imperfect_mrc_lin>0, imperfect_mrc_dB=10*log10(imperfect_mrc_lin); else, imperfect_mrc_dB=-Inf; end

    % --- RL Decision ---
    U_global = max(U_vec);
    [~,state] = min(abs(thresholds-U_global));
    [~,aIdx] = max(Q(state,:));
    chosenTh = thresholds(state);
    if U_global < chosenTh
        snr_adaptive(k)=imperfect_mrc_dB;
        rl_action(k)=1;
    else
        snr_adaptive(k)=snr_sc(k);
        rl_action(k)=0;
    end
    uncertainty(k)=U_global;

    % Averages
    if count_non_zeros_ebnr>0
        avg_ebnr_lin=sum(ebno_faded_lin_vec)/count_non_zeros_ebnr;
        average_ebnr(k)=10*log10(max(avg_ebnr_lin,eps));
    else
        average_ebnr(k)=-Inf;
    end
    if count_non_zeros_pr>0
        average_rp_lin=sum(pr_lin_vec)/count_non_zeros_pr;
        average_rp(k)=10*log10(max(average_rp_lin,eps));
    else
        average_rp(k)=-Inf;
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
semilogy(t,10.^(snr_mrc./10),'LineWidth',1.5); hold on;
semilogy(t,10.^(snr_sc./10),'LineWidth',1.5);
semilogy(t,10.^(snr_adaptive./10),'LineWidth',1.5); hold off;
legend('MRC (true)','SC','RL-Adaptive');
title('MRC vs SC vs RL Adaptive (linear)'); xlabel('time'); ylabel('Eb/N_0 (linear)'); grid on;

figure;
plot(t,snr_mrc,'-','LineWidth',1.2); hold on;
plot(t,snr_sc,'--','LineWidth',1.2);
plot(t,snr_adaptive,':','LineWidth',1.4); hold off;
legend('MRC (true)','SC','RL-Adaptive');
title('MRC vs SC vs RL Adaptive (dB)'); xlabel('time'); ylabel('Eb/N_0 (dB)'); grid on;

figure;
yyaxis left
plot(t,uncertainty,'LineWidth',1.5);
ylabel('Uncertainty U');
yyaxis right
stairs(t,rl_action,'r','LineWidth',1.2);
ylabel('RL Action (0=SC,1=MRC)');
xlabel('Time');
title('Uncertainty and RL Action over Time');
grid on;
legend('Uncertainty','RL Action');

