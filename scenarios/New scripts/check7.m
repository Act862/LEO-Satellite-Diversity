% RL threshold-based adaptive combining (Q-learning with threshold states)
clear; clc; rng(0);

%% Simulation time settings
startTime = datetime(2025,09,11,12,01,0);
endTime = startTime + minutes(5);
sampleTime = 10;                    % seconds
t = startTime:seconds(sampleTime):endTime;

%% Scenario & satellites (unchanged)
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

%% Channel parameters (same as your script)
K = 5;                        
pilotSNR_dB = 10;              
Npilots = 4;                   
pilotSNR_lin = 10^(pilotSNR_dB/10);
sigma_e2 = 1/(pilotSNR_lin * Npilots);

%% --- NEW Reinforcement Learning Setup (threshold-as-state Q-learning) ---
nThresholdBins = 21;                 % bins for threshold in [0,1]
thresholds = linspace(0,1,nThresholdBins);

actions = [-1, 0, 1];                % move left, stay, move right
nActions = numel(actions);

% Q-table: rows = threshold bins (states), cols = actions
Q = zeros(nThresholdBins, nActions);

% Q-learning hyperparameters
alpha = 0.1;          % learning rate
gammaRL = 0.99;       % discount factor
epsilonStart = 1.0;   % initial exploration
epsilonEnd   = 0.05;  % final exploration
epsilonDecaySteps = 10000;

numEpisodes = 2000;   % training episodes
stepsPerEpisode = min(50, length(t));  % steps per episode for training

% Reward shaping
lambda_uncertainty = 1.0;   % penalty factor for U in reward
switch_cost = 0.5;          % penalty if chosen combining changes from previous step

% Channel model params used inside train (same idea as your generation)
LOS_factor = sqrt(K/(K+1));
sigma_factor = sqrt(1/(2*(K+1)));

% helper for epsilon scheduling
epsilonByStep = @(step) (epsilonStart*(1 - min(1,step/epsilonDecaySteps)) + epsilonEnd*min(1,step/epsilonDecaySteps));

%% --- Offline Training (simulate channel & estimation) ---
globalStep = 0;
for ep = 1:numEpisodes
    % randomize starting threshold (state)
    stateIdx = randi(nThresholdBins);
    prev_scheme = NaN;  % used to penalize switching between scheme choices
    for step = 1:stepsPerEpisode
        eps = epsilonByStep(globalStep);
        % epsilon-greedy: choose action index
        if rand < eps
            aIdx = randi(nActions);
        else
            [~, aIdx] = max(Q(stateIdx,:));
        end
        action = actions(aIdx);
        % next state index after applying action (clip)
        nextStateIdx = max(1, min(nThresholdBins, stateIdx + action));
        threshold = thresholds(nextStateIdx);
        
        % --- generate channel sample & estimate (same model you used) ---
        phi = 2*pi*rand;
        X = randn; Y = randn;
        h_complex = LOS_factor*exp(1j*phi) + sigma_factor*(X + 1j*Y);
        h_abs2 = abs(h_complex)^2;
        
        % estimation noise
        e_real = sqrt(sigma_e2)*randn;
        e_imag = sqrt(sigma_e2)*randn;
        e_complex = (e_real + 1j*e_imag)/sqrt(2);
        h_hat = h_complex + e_complex;
        U = sigma_e2 / (abs(h_hat)^2 + eps);  % uncertainty measure
        
        % Pick a nominal Eb/N0 (use a random sample of link conditions, as before)
        ebno_nom_dB = 10*rand; % small randomization (0..10 dB); adapt to your system
        ebno_nom_lin = 10^(ebno_nom_dB/10);
        ebno_faded_lin = ebno_nom_lin * h_abs2;
        
        % SNR for SC (best branch in single-sat toy; here we use faded link as SC)
        snr_sc_sample_dB = 10*log10(max(ebno_faded_lin, eps));
        
        % For MRC the agent uses estimated channel for weighting; emulate imperfect MRC:
        ebno_est_faded_lin = ebno_nom_lin * abs(h_hat)^2;
        snr_mrc_sample_dB = 10*log10(max(ebno_est_faded_lin, eps));
        
        % Decide which combining would be used based on threshold (this sim is to compute reward)
        chosen_scheme = (U < threshold); % 1 -> MRC, 0 -> SC
        if chosen_scheme == 1
            reward_snr = snr_mrc_sample_dB;
        else
            reward_snr = snr_sc_sample_dB;
        end
        
        switchingPenalty = 0;
        if ~isnan(prev_scheme)
            if prev_scheme ~= chosen_scheme
                switchingPenalty = switch_cost;
            end
        end
        
        reward = reward_snr - lambda_uncertainty * U - switchingPenalty;
        
        % Simulate a next-state measurement (one-step ahead) for Q update:
        phi2 = 2*pi*rand;
        X2 = randn; Y2 = randn;
        h2 = LOS_factor*exp(1j*phi2) + sigma_factor*(X2 + 1j*Y2);
        e2 = (sqrt(sigma_e2)/sqrt(2))*(randn + 1j*randn);
        hhat2 = h2 + e2;
        U2 = sigma_e2 / (abs(hhat2)^2 + eps);
        % discretize next state (threshold remains the agent's state variable; next state is threshold after action)
        % but the discretization of U2 is not needed because our states are thresholds.
        % We use nextStateIdx for the TD update
        bestNext = max(Q(nextStateIdx,:));
        Q(stateIdx,aIdx) = Q(stateIdx,aIdx) + alpha * (reward + gammaRL * bestNext - Q(stateIdx,aIdx));
        
        % update
        prev_scheme = chosen_scheme;
        stateIdx = nextStateIdx;
        globalStep = globalStep + 1;
    end
end

%% --- Deployment Phase: Use learned Q for decisions on the actual satellite links ---
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
        X = randn; Y = randn;
        h_complex = LOS_factor*exp(1j*phi) + sigma_factor*(X + 1j*Y);
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

    % True MRC (ideal): sum of faded linear Eb/N0
    if ebnr_sum_lin > 0, snr_mrc(k) = 10*log10(ebnr_sum_lin); else, snr_mrc(k)=-Inf; end
    if max_ebnr_lin > 0, snr_sc(k) = 10*log10(max_ebnr_lin); else, snr_sc(k)=-Inf; end

    % Imperfect MRC using estimated channel coefficients (same as before)
    ebno_est_faded_lin = ebno_faded_lin_vec ./ (h_vec+eps) .* abs(hhat_vec).^2;
    imperfect_mrc_lin = sum(ebno_est_faded_lin);
    if imperfect_mrc_lin > 0, imperfect_mrc_dB = 10*log10(imperfect_mrc_lin); else, imperfect_mrc_dB=-Inf; end

    % --- RL Decision: choose threshold bin from Q-table greedy policy ---
    % For the deployment we use the *global* uncertainty measure (as your script used)
    U_global = max(U_vec);
    % find nearest threshold index to current threshold state (we will act greedily on that state)
    % But RL state is threshold — we must pick an index representing the current threshold.
    % To do that we pick the best action chain fixed-point: simplest is to pick the threshold state
    % whose greedy action points to itself (fixed-points), otherwise pick argmax of state-value.
    % For simplicity here we choose the greedy action from the state with highest state-value.
    stateValues = max(Q,[],2);   % value per threshold state
    [~, bestStateIdx] = max(stateValues); 
    % Alternatively we could map the current U_global to a threshold candidate; here we use greedy bestStateIdx.
    greedyActionIdx = find(Q(bestStateIdx,:) == max(Q(bestStateIdx,:)),1,'first');
    greedyAction = actions(greedyActionIdx);
    % next threshold index if the greedy action is taken from bestStateIdx
    chosenStateIdx = max(1, min(nThresholdBins, bestStateIdx + greedyAction));
    chosenThreshold = thresholds(chosenStateIdx);

    % Another valid approach: choose the threshold whose index is argmax over state-value and use it directly.
    % The chosenThreshold is now used for decision:
    if U_global < chosenThreshold
        snr_adaptive(k) = imperfect_mrc_dB;
    else
        snr_adaptive(k) = snr_sc(k);
    end

    % Averages (unchanged)
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

%% Plots (unchanged)
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
