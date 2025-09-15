% Script with adaptive combining (switch MRC/SC) keeping original structure
% Assuming downlink transmission.
% For uplink transmission we need to add transmitters on the base station
% and one receiver on each satellite.
clear; clc;

%% Simulation time settings
startTime = datetime(2025,09,11,12,01,0);
endTime = startTime + minutes(5);
sampleTime = 10;                    % seconds
t = startTime:seconds(sampleTime):endTime;

%% Scenario & sats
sc = satelliteScenario(startTime,endTime,sampleTime);
sat = satellite(sc,"leoSatelliteConstellation.tle");

% Transmitters on sats (downlink)
satTxs = [];
gTxs = [];
systemLoss = 3;
dishDiameter = 0.6;
apertureEfficiency = 0.5881;
for s = sc.Satellites
    gimbalTx = gimbal(s);
    gTxs = [gTxs gimbalTx];
    tx = transmitter(gimbalTx,Name="TxSat"+s.ID,Frequency=20e9,Power=20,BitRate=10,SystemLoss=systemLoss);
    gaussianAntenna(tx,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);
    satTxs = [satTxs tx];
end

% Ground station
gsLat = 37.9838; gsLon = 23.7275; gsAlt = 0.1; % km
gs = groundStation(sc,gsLat,gsLon,"Altitude",gsAlt,"MinElevationAngle",10,"Name","Base Station");
gainToNoiseTemperatureRatio = 18.5;

% Point sats at GS
for s = 1:length(sc.Satellites)
    pointAt(gTxs(s), gs);
end

% Build link objects (one RX gimbal per satellite)
links = [];
for s = 1:length(sc.Satellites)
    gimbalRx_s = gimbal(gs);
    pointAt(gimbalRx_s, sat(s));
    rx_s = receiver(gimbalRx_s,Name="Rx_GS_Sat"+sat(s).ID, GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio, SystemLoss=systemLoss);
    lnk = link(satTxs(s), rx_s);
    links = [links lnk];
end

%% Simulation
% --- Parameters for fading & estimator ---
K = 1;                         % Rician K-factor (linear)
pilotSNR_dB = 10;              % pilot SNR (dB) used for estimator quality
Npilots = 4;                   % number of pilot symbols
pilotSNR_lin = 10^(pilotSNR_dB/10);
sigma_e2 = 1/(pilotSNR_lin * Npilots);   % simple estimator variance model
gammaVec = [0 0.05 0.1 0.5 1];                   % uncertainty threshold for switching (tune this)

for gamma=gammaVec
    % output arrays
    average_ebnr = zeros(1,length(t));
    average_rp   = zeros(1,length(t));
    snr_mrc      = zeros(1,length(t));   % will store dB
    snr_sc       = zeros(1,length(t));   % dB
    snr_adaptive = zeros(1,length(t));   % dB (our switching scheme)
    for k = 1:length(t)
        ebnr_sum_lin = 0;            % linear sum for MRC (sum of Eb/N0 linear)
        pr_sum_lin = 0;              % linear sum of received power
        count_non_zeros_ebnr = 0;
        count_non_zeros_pr = 0;
        numSat_underService = 0;
        max_ebnr_lin = 0;            % linear max for SC
        
        % for adaptive: collect estimate variables locally
        ebno_nom_lin_vec = zeros(1,length(sc.Satellites)); % store nominal Eb/N0 linear
        ebno_faded_lin_vec = zeros(1,length(sc.Satellites)); % true faded Eb/N0 linear
        h_vec = zeros(1,length(sc.Satellites));   % true complex h (we store mag^2)
        hhat_vec = zeros(1,length(sc.Satellites));% estimated complex h_hat (store complex)
        U_vec = zeros(1,length(sc.Satellites));   % uncertainty per link
        pr_lin_vec = zeros(1,length(sc.Satellites)); % pr linear
        
        for s = 1:length(sc.Satellites)
            % --- Rician fading sample (complex) ---
            phi = 2*pi*rand;
            LOS = sqrt(K/(K+1));
            sigma = sqrt(1/(2*(K+1)));
            X = randn; Y = randn;
            h_complex = LOS*exp(1j*phi) + sigma*(X + 1j*Y);
            h_abs2 = abs(h_complex)^2;
            
            % --- Nominal Eb/N0 (dB -> linear) ---
            ebno_nom_dB = ebno(links(s), t(k)); % could be -Inf
            if ebno_nom_dB == -Inf || isnan(ebno_nom_dB)
                ebno_nom_lin = 0;
            else
                ebno_nom_lin = 10^(ebno_nom_dB/10);
            end
            
            % --- Faded Eb/N0 (linear) ---
            ebno_faded_lin = ebno_nom_lin * h_abs2; % linear multiplication by |h|^2
            
            % --- Channel estimation (h_hat) and uncertainty metric U ---
            % e ~ CN(0, sigma_e2)
            e_real = sqrt(sigma_e2)*randn;
            e_imag = sqrt(sigma_e2)*randn;
            e_complex = (e_real + 1j*e_imag)/sqrt(2); % complex normal scaling
            h_hat = h_complex + e_complex;
            U = sigma_e2 / (abs(h_hat)^2 + eps);  % relative uncertainty metric
            
            % --- Received power (sigstrength). Convert to linear (W) if available ---
            [~, pr_dB] = sigstrength(links(s), t(k)); % may be -Inf
            if pr_dB == -Inf || isnan(pr_dB)
                pr_lin = 0;
            else
                pr_lin = 10^(pr_dB/10);
            end
            
            % --- Collect per-satellite for later combining and statistics ---
            ebno_nom_lin_vec(s) = ebno_nom_lin;
            ebno_faded_lin_vec(s) = ebno_faded_lin;
            h_vec(s) = h_abs2;         % store magnitude^2 (real)
            hhat_vec(s) = h_hat;       % complex estimated
            U_vec(s) = U;
            pr_lin_vec(s) = pr_lin;
            
            % --- Local bookkeeping (counts / sums) ---
            if ebno_faded_lin > 0
                count_non_zeros_ebnr = count_non_zeros_ebnr + 1;
                numSat_underService = numSat_underService + 1;
            end
            if pr_lin > 0
                count_non_zeros_pr = count_non_zeros_pr + 1;
            end
            
            % For MRC true (we will sum later); still accumulate sums for quick checks:
            ebnr_sum_lin = ebnr_sum_lin + ebno_faded_lin;   % linear sum for true MRC
            pr_sum_lin = pr_sum_lin + pr_lin;
            
            % track SC max (linear)
            if ebno_faded_lin > max_ebnr_lin
                max_ebnr_lin = ebno_faded_lin;
            end
        end % end per-satellite loop
        
        % --- Compute MRC (true) and SC using linear domain, then convert to dB ---
        % MRC true (uses true h)
        if ebnr_sum_lin > 0
            snr_mrc(k) = 10*log10(ebnr_sum_lin);
        else
            snr_mrc(k) = -Inf;
        end
        
        % SC (best single link by true faded Eb/N0)
        if max_ebnr_lin > 0
            snr_sc(k) = 10*log10(max_ebnr_lin);
        else
            snr_sc(k) = -Inf;
        end
        
        % --- Imperfect MRC (using estimates h_hat) approximation ---
        % Reconstruct nominal Eb/N0 per link (avoid dividing by zero)
        ebno_nom_lin_vec2 = zeros(size(ebno_nom_lin_vec)); % same length
        for s = 1:length(ebno_nom_lin_vec)
            if h_vec(s) > 0
                ebno_nom_lin_vec2(s) = ebno_faded_lin_vec(s) / h_vec(s);
            else
                ebno_nom_lin_vec2(s) = 0;
            end
        end
        % estimated faded Eb/N0 using |h_hat|^2
        ebno_est_faded_lin = ebno_nom_lin_vec2 .* (abs(hhat_vec).^2);
        imperfect_mrc_lin = sum(ebno_est_faded_lin);
        if imperfect_mrc_lin > 0
            imperfect_mrc_dB = 10*log10(imperfect_mrc_lin);
        else
            imperfect_mrc_dB = -Inf;
        end
        
        % --- Adaptive switching: use imperfect MRC if uncertainty small, else SC ---
        U_global = median(U_vec); % decision metric (you may choose mean or percentile)
        if U_global <= gamma
            % choose imperfect MRC
            snr_adaptive(k) = imperfect_mrc_dB;
        else
            % fallback to SC (true SC uses true h; practically one would choose best estimate)
            snr_adaptive(k) = snr_sc(k);
        end
        
        % --- Averages for plotting (match your original outputs) ---
        if count_non_zeros_ebnr > 0
            % average_ebnr in your original code was in dB; keep same meaning but compute correctly
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
        
        disp("The number of satellites servicing the base station at time "+ string(t(k)) +" is: "+numSat_underService);
    end % end time loop
    
    %% Plots (kept in same style as your original script)
    % Convert average_ebnr/average_rp back to linear for the left plots (you used 10.^(.../10))
    figure;
    subplot(1,2,1);
    plot(t,10.^(average_ebnr./10),'LineWidth',2);
    title('average Eb/N_0 against time');
    ylabel('Eb/N_0 (linear scale) Watts');
    xlabel('simulation time (datetime)');
    grid on;
    
    subplot(1,2,2);
    plot(t,10.^(average_rp./10),'LineWidth',2);
    title('average receiver power against time');
    ylabel('P (linear scale) Watts');
    xlabel('simulation time (datetime)');
    grid on;
    
    % Comparison SNR traces
    figure;
    semilogy(t, 10.^(snr_mrc./10),'LineWidth',1.5); % plotting in linear on semilogy matches your earlier intent
    hold on;
    semilogy(t, 10.^(snr_sc./10),'LineWidth',1.5);
    semilogy(t, 10.^(snr_adaptive./10),'LineWidth',1.5);
    hold off;
    grid on;
    legend('MRC (true)','SC','Adaptive (switch MRC/SC)');
    text = sprintf("MRC vs SC vs Adaptive (linear scale),\\gamma=%.2f",gamma);
    title(text);
    xlabel('simulation time (datetime)');
    ylabel('Eb/N_0 (linear Watts)');
    
    % Optional: also plot dB traces for clearer reading
    figure;
    plot(t, snr_mrc, '-','LineWidth',1.2); hold on;
    plot(t, snr_sc, '--','LineWidth',1.2);
    plot(t, snr_adaptive, ':','LineWidth',1.4);
    hold off;
    grid on;
    legend('MRC (true)','SC','Adaptive (dB)');
    text = sprintf('MRC vs SC vs Adaptive (dB), \\gamma=%.2f',gamma);
    title(text);
    xlabel('simulation time (datetime)');
    ylabel('Eb/N_0 (dB)');
end