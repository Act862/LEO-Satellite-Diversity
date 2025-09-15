%% LEO Constellation + RL-aided GAC vs SC/MRC
% Satellite Toolbox + Reinforcement Learning Toolbox
% Author: (you)
clear; clc; close all;

%% ===================== USER INPUTS =====================
tlePath   = "leoSatelliteConstellation.tle";         % <-- your LEO constellation TLE file
gsName    = "Lamia Base Station";
gsLat     = 38.875941904578895;
gsLon     = 22.437526584061686;
gsAlt     = 0;                           % meters
startTime = datetime(2020,1,11,14,0,0);  % simulation start
stopTime  = startTime + hours(10);         % simulation span
sampleSec = 60;                          % scenario sample time (s)

% RF/link setup (simple/consistent across all links)
txFrequencyHz = 20e9;          % downlink frequency (Hz)
txPowerW      = 15;            % per-satellite TX power (W)
satDish_m     = 0.5;           % sat antenna (gaussian dish)
gsDish_m      = 2.0;           % GS antenna (gaussian dish)
gsG_T_dB      = 20;            % GS G/T (dB/K), a single RX at GS
reqEbNo_dB    = 1;             % Required Eb/No (dB) if you use margin calcs
bitRate       = 20e6;          % just used in Eb/No calc if needed
c             = physconst('LightSpeed');

% RL settings
maxEpisodes   = 80;            % keep modest so it runs fast; increase for better fit
maxStepsPerEp = 300;           % steps per episode (time steps sampled)
useGPU        = false;         % set true if you want RL training on GPU

rng(1);

%% ===================== SCENARIO =====================
sc = satelliteScenario(startTime, stopTime, sampleSec);

% Create the constellation from the TLE file
sats = satellite(sc, tlePath);

% One ground station + a single tracking gimbal + receiver
gs = groundStation(sc, Name=gsName, Latitude=gsLat, Longitude=gsLon, Altitude=gsAlt);

% A single RX chain at the GS
gGs = gimbal(gs, "MountingAngles", [0;180;0], "MountingLocation", [0; 0; -5]); % simple mount
gsRx = receiver(gGs, Name="GS_RX", MountingLocation=[0;0;1], ...
    GainToNoiseTemperatureRatio=gsG_T_dB, RequiredEbNo=reqEbNo_dB);
gaussianAntenna(gsRx, "DishDiameter", gsDish_m);

% Add a TX on every satellite, each with its own pointing gimbal
nSats = numel(sats);
satTx  = cell(1, nSats);
gSats  = cell(1, nSats);
links  = cell(1, nSats);
for k = 1:nSats
    gSats{k} = gimbal(sats(k), "MountingLocation", [0; 0; 2]);
    satTx{k} = transmitter(gSats{k}, "Name", "SatTx"+k, ...
        Frequency=txFrequencyHz, Power=txPowerW, BitRate=bitRate);
    gaussianAntenna(satTx{k}, "DishDiameter", satDish_m);
    % Point both ways
    pointAt(gSats{k}, gs);
    pointAt(gGs, sats(k)); % GS tracks whichever is set active (we’ll compute all anyway)

    % Link: sat k downlink to GS
    links{k} = link(satTx{k}, gsRx);
end

% Pre-collect all scenario sample times
T = sc.StartTime : seconds(sc.SampleTime) : sc.StopTime;
NT = numel(T);

%% ===================== TRUE SLANT RANGES & VISIBILITY =====================
% We'll compute truth slant ranges (ECEF) at each time to the GS
% Convert GS to ECEF
wgs = wgs84Ellipsoid("meter");
[gsX, gsY, gsZ] = geodetic2ecef(wgs, gs.Latitude, gs.Longitude, gs.Altitude);

trueMaxSlant_m   = zeros(NT,1);
visibleMask      = false(NT, nSats);   % visibility flag (access)
slantAll_m       = NaN(NT, nSats);     % per-sat slant range (truth)
ebnoAll_dB       = -Inf(NT, nSats);    % per-sat Eb/No from link budget (dB)

% We’ll also grab access flags and Eb/No via Toolbox for fairness
for k = 1:nSats
    % Evaluate Eb/No and range using link object
    % (ebno() and range() return time-synchronized vectors over scenario sample times)
    [ebno_k, t_k] = ebno(links{k});
    % Ensure alignment to T
    % Most Satellite Toolbox functions are aligned to scenario SampleTime already,
    % but we’ll map indices robustly:
    [~, ia, ib] = intersect(T, t_k);
    ebnoAll_dB(ia, k) = ebno_k(ib);

    % Access intervals -> boolean visibility
    acc = access(gs, sats(k));
    ints = accessIntervals(acc); %#ok<NASGU> (optional: inspect/plot)

    % Compute ECEF state and slant range
    % (states returns position in ECEF if requested)
    pos_k = zeros(3, NT);
    for ii = 1:NT
        [pECEF, ~] = states(sats(k), T(ii), "CoordinateFrame","ecef");
        pos_k(:,ii) = pECEF; % meters
        % line-of-sight visibility if Eb/No finite and not -Inf (simple proxy)
        if isfinite(ebnoAll_dB(ii,k)) && (ebnoAll_dB(ii,k) > -300)
            visibleMask(ii,k) = true;
        end
    end
    % slant range
    dx = pos_k(1,:)' - gsX;
    dy = pos_k(2,:)' - gsY;
    dz = pos_k(3,:)' - gsZ;
    slantAll_m(:,k) = sqrt(dx.^2 + dy.^2 + dz.^2);
end

% Truth max slant range among visible satellites; if none visible, set NaN
for ii = 1:NT
    visIdx = find(visibleMask(ii,:));
    if ~isempty(visIdx)
        trueMaxSlant_m(ii) = max(slantAll_m(ii, visIdx));
    else
        trueMaxSlant_m(ii) = NaN;
    end
end

%% ===================== FEATURE ENGINEERING FOR RL =====================
% Observation at time t:
%   [ sin(2πtod), cos(2πtod), sin(2πdoy), cos(2πdoy), #visible/nSats ]
% Minimal, time/context-only (no future/cheating).
tod = seconds(timeofday(T)) / (24*3600); % fractional time-of-day
doy = day(T, "dayofyear") / 366;        % fractional day-of-year (approx)
fracVis = sum(visibleMask,2) / max(1,nSats);

obsAll = [sin(2*pi*tod(:)), cos(2*pi*tod(:)), ...
          sin(2*pi*doy(:)), cos(2*pi*doy(:)), ...
          fracVis(:)];

% Targets (truth): normalized max slant range in [0,1] for RL
% normalize by 5000 km (robust upper bound for LEO slant to GS) — tweak as you like
normDen_m = 5e6; % 5000 km
yTrue = trueMaxSlant_m / normDen_m;

% Mask out timesteps with no visibility (no targets)
valid = ~isnan(yTrue);
obsValid = obsAll(valid,:);
yValid  = yTrue(valid);

%% ===================== RL ENVIRONMENT (Continuous Action: predict scalar) =====================
% Create an index mapping so env steps across valid rows
validIdx = find(valid);
Nvalid   = numel(validIdx);

% Guard for tiny runs (e.g., 8–10 minutes)
if Nvalid < 2
    error('RL:NotEnoughData','Not enough valid samples for RL (Nvalid=%d). Increase sim time or relax visibility.',Nvalid);
end

% Specs
obsInfo = rlNumericSpec([5 1], LowerLimit=-inf(5,1), UpperLimit=inf(5,1));
obsInfo.Name = "Observation";
actInfo = rlNumericSpec([1 1], LowerLimit=0, UpperLimit=1); % predicted normalized max slant range
actInfo.Name = "PredictedNormMaxSlant";

% Wrap local functions so signatures match:
%   resetFcn: 0 inputs
%   stepFcn:  (Action, LoggedSignals)
resetFcn = @() localReset(obsValid, Nvalid, maxStepsPerEp);
stepFcn  = @(Action,LoggedSignals) localStep(Action, LoggedSignals, obsValid, yValid, Nvalid, maxStepsPerEp);

% Create env (NOTE: only once)
env = rlFunctionEnv(obsInfo, actInfo, stepFcn, resetFcn);


%% ===================== RL AGENT (DDPG) =====================
% Simple actor-critic with tiny networks
statePath = [
    featureInputLayer(5,"Normalization","none","Name","obs")
    fullyConnectedLayer(32,"Name","fc1")
    reluLayer("Name","relu1")
    fullyConnectedLayer(32,"Name","fc2")
    reluLayer("Name","relu2")];

actorPath = [
    statePath
    fullyConnectedLayer(1,"Name","actFC")
    sigmoidLayer("Name","sig")]; % outputs [0,1]

criticPath = [
    featureInputLayer(5,"Normalization","none","Name","obs")
    fullyConnectedLayer(32,"Name","cfc1")
    reluLayer("Name","crelu1")
    concatenationLayer(1,2,"Name","concat") % concat obs features with action
    fullyConnectedLayer(32,"Name","cfc2")
    reluLayer("Name","crelu2")
    fullyConnectedLayer(1,"Name","qout")];

criticObsPath = layerGraph(criticPath);
criticActPath = featureInputLayer(1,"Normalization","none","Name","act");
criticLG = addLayers(criticObsPath, criticActPath);

% Only connect the action to concat/in2
criticLG = connectLayers(criticLG,"act","concat/in2");


actorNet  = dlnetwork(layerGraph(actorPath));
criticNet = dlnetwork(criticLG);

actorOpts = rlOptimizerOptions(LearnRate=1e-3, GradientThreshold=1);
criticOpts= rlOptimizerOptions(LearnRate=5e-3, GradientThreshold=1);

agentOpts = rlDDPGAgentOptions( ...
    SampleTime = sc.SampleTime, ...
    DiscountFactor = 0.999, ...
    ExperienceBufferLength = 1e5, ...
    MiniBatchSize = 256, ...
    TargetSmoothFactor = 1e-3);


actor = rlContinuousDeterministicActor(actorNet, obsInfo, actInfo);
critic = rlQValueFunction(criticNet, obsInfo, actInfo);

agentOpts.ActorOptimizerOptions  = actorOpts;
agentOpts.CriticOptimizerOptions = criticOpts;

agent = rlDDPGAgent( ...
    rlContinuousDeterministicActor(actorNet, obsInfo, actInfo), ...
    rlQValueFunction(criticNet, obsInfo, actInfo), ...
    agentOpts);

trainOpts = rlTrainingOptions( ...
    MaxEpisodes=maxEpisodes, ...
    MaxStepsPerEpisode=maxStepsPerEp, ...
    StopTrainingCriteria="AverageReward", ...
    StopTrainingValue=-5e-4, ... % loose stop
    ScoreAveragingWindowLength=10, ...
    Plots="training-progress", ...
    UseParallel=false);

if useGPU, trainOpts.UseDevice = "gpu"; end %#ok<UNRCH>

trainingStats = train(agent, env, trainOpts); %#ok<NASGU>

%% ===================== INFERENCE: RL PREDICTIONS OVER ALL TIMES =====================
predNorm = NaN(NT,1);
for ii = 1:NT
    o = obsAll(ii,:)';
    if ~valid(ii)
        predNorm(ii) = NaN;
    else
        a = getAction(agent, o);
        if iscell(a), a = a{1}; end
        predNorm(ii) = double(a);
    end
end
predMaxSlant_m = predNorm * normDen_m;

%% ===================== COMBINING: GAC vs SC vs MRC =====================
% We use Toolbox Eb/No per-satellite per-time (ebnoAll_dB).
% Convert to linear SNR per link (assuming Eb/No ~ SNR per symbol for simplicity here).
SNRlin = 10.^(ebnoAll_dB/10);

% --- Selection Combining (SC): max SNR among visible
SNR_SC = zeros(NT,1);
for ii = 1:NT
    s = SNRlin(ii,:);
    s(~visibleMask(ii,:)) = -Inf;
    SNR_SC(ii) = max(s);
    if ~isfinite(SNR_SC(ii)), SNR_SC(ii) = NaN; end
end

% --- Maximal Ratio Combining (MRC): sum of SNRs (coherent)
SNR_MRC = zeros(NT,1);
for ii = 1:NT
    s = SNRlin(ii,:);
    s(~visibleMask(ii,:)) = 0;
    SNR_MRC(ii) = sum(s);
    if SNR_MRC(ii)==0, SNR_MRC(ii)=NaN; end
end

% --- Geometry-Assisted Combining (GAC):
% Build geometry weights from PREDICTED slant range.
% Free-space path loss ~ R^2, so an intuitive geometry weight ∝ 1/R^2.
% Normalize weights to sum to 1 across *visible* sats at each time.
SNR_GAC = zeros(NT,1);
lambda = c / txFrequencyHz;

for ii = 1:NT
    vis = visibleMask(ii,:);
    if ~any(vis) || isnan(predMaxSlant_m(ii))
        SNR_GAC(ii) = NaN; continue;
    end

    % For GAC weights, we need a per-satellite geometry surrogate.
    % We’ll map the predicted MAX slant to a per-satellite factor by
    % using each sat’s ACTUAL slant but clipping it to the predicted max
    % (i.e., avoid over-trusting geometry when predict says “long”). Alternatives exist.
    Ri = slantAll_m(ii, :);
    Ri(~vis) = NaN;
    Ri_eff = min(Ri, predMaxSlant_m(ii)); % geometry-assist using the RL prediction cap
    % Weight ∝ 1/Ri_eff^2 (Free-space). Avoid divide-by-zero:
    wi = 1 ./ max(Ri_eff.^2, 1);
    wi(~vis | isnan(wi)) = 0;
    if sum(wi)==0, SNR_GAC(ii)=NaN; continue; end
    wi = wi / sum(wi);

    % Combine linearly in SNR domain with geometry weights.
    % (If you prefer amplitude-domain combining, replace with
    %   (sum(wi .* sqrt(SNRi)))^2  which mimics coherent weighting.)
    s = SNRlin(ii,:);
    s(~vis) = 0;
    SNR_GAC(ii) = sum(wi .* s);
    if SNR_GAC(ii)==0, SNR_GAC(ii)=NaN; end
end

%% ===================== PLOTS & METRICS =====================
figure; 
plot(T, trueMaxSlant_m/1e3, '.', 'DisplayName','Truth max slant (km)'); hold on;
plot(T, predMaxSlant_m/1e3, '-', 'DisplayName','RL predicted max slant (km)');
grid on; xlabel('Time'); ylabel('km'); title('Max Slant Range: Truth vs RL Prediction');
legend('Location','best');

% Compare SNRs
figure;
plot(T, 10*log10(SNR_SC), '.', 'DisplayName','SC'); hold on;
plot(T, 10*log10(SNR_MRC), '.', 'DisplayName','MRC');
plot(T, 10*log10(SNR_GAC), '.', 'DisplayName','GAC (RL-assisted)');
grid on; xlabel('Time'); ylabel('Combined SNR (dB)');
title('Combining Gain: SC vs MRC vs GAC (RL-assisted)');
legend('Location','best');

% (Optional) BER for QPSK using simple approximations
% (Using SNR as Eb/No here for illustration.)
Q = @(x) 0.5*erfc(x/sqrt(2));
berQPSK = @(snrLin) 0.5*(2*Q(sqrt(2*snrLin)) - Q(sqrt(2*snrLin)).^2);

figure;
plot(T, berQPSK(SNR_SC), '.', 'DisplayName','SC'); hold on;
plot(T, berQPSK(SNR_MRC), '.', 'DisplayName','MRC');
plot(T, berQPSK(SNR_GAC), '.', 'DisplayName','GAC (RL-assisted)');
set(gca,'YScale','log'); grid on;
xlabel('Time'); ylabel('P_b (log)'); title('QPSK BER (Illustrative)');
legend('Location','best');

%% ===================== OPTIONAL: Extend time & 3 closest sats =====================
% Example: for each time, list the 3 nearest satellites (smallest slant)
Nclose = 3;
closestIdx = NaN(NT, Nclose);
for ii = 1:NT
    vis = visibleMask(ii,:);
    if ~any(vis), continue; end
    [~, order] = sort(slantAll_m(ii, vis), 'ascend'); % nearest
    visList = find(vis);
    kk = min(Nclose, numel(order));
    closestIdx(ii,1:kk) = visList(order(1:kk));
end
% closestIdx(ii,:) has indices of nearest sats at time i (NaN where not enough)

disp('Done.')
