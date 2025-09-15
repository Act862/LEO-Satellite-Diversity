%   Simple channel simulation
%   Satellite is the transmitter and base station is the receiver
clear;clc;

startTime = datetime(2025,8,30,17,31,0);
endTime = startTime + minutes(10);
sampleTime = 30;                    % seconds

t = startTime:seconds(sampleTime):endTime;

%   Initialize the satellite scenario
sc = satelliteScenario(startTime,endTime,sampleTime);

%   Add satellites
sat = satellite(sc,"leoSatelliteConstellation.tle");

%   add tx on each satellite
satTxs = [];
gTxs = [];
systemLoss = 3;
dishDiameter = 0.5;
apertureEfficiency = 0.5;

for s=sc.Satellites
    gimbalTx = gimbal(s);
    gTxs = [gTxs gimbalTx];
    tx = transmitter(gimbalTx,Name="TxSat"+s.ID,Frequency=30e9,...
        Power=20,BitRate=20,SystemLoss=systemLoss);
    gaussianAntenna(tx,DishDiameter=dishDiameter,...
        ApertureEfficiency=apertureEfficiency);
    satTxs = [satTxs tx];
end

%   Add ground station
gsLat = 37.9838;     % deg
gsLon = 23.7275;     % deg
gsAlt = 0.1;         % km
gs = groundStation(sc,gsLat,gsLon,"Altitude",gsAlt,"MinElevationAngle",...
    5,"Name","Base Station");

gainToNoiseTemperatureRatio = 5;

%   make all the satellites point at the base station
for s=1:length(sc.Satellites)
    pointAt(gTxs(s), gs);
end

%   generate all the links (multi-gimbal GS side)
links = [];
for s=1:length(sc.Satellites)
    % independent gimbal for this satellite
    gimbalRx_s = gimbal(gs);
    pointAt(gimbalRx_s, sat(s));

    % independent receiver for this satellite
    rx_s = receiver(gimbalRx_s,Name="Rx_GS_Sat"+sat(s).ID, ...
        GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio,...
        SystemLoss = systemLoss);

    % build link
    lnk = link(satTxs(s),rx_s);
    links = [links lnk];
end

% geometric access check (for comparison)
% ac = access(sat,gs);
% play(sc);

% Preallocate arrays
average_ebnr = zeros(1,length(t));
average_rp = zeros(1,length(t));
snr_mrc = zeros(1,length(t));
snr_sc = zeros(1,length(t));
snr_hybrid_std = zeros(1,length(t));
snr_hybrid_minr = zeros(1,length(t));

for k = 1:length(t)
    ebnr_lin_vec = zeros(1,length(sc.Satellites));
    pr_vec = zeros(1,length(sc.Satellites));
    r_vec = zeros(1,length(sc.Satellites));
    h_vec = zeros(1,length(sc.Satellites));
    numSat_underService = 0;
    
    % --- Step 1: collect Eb/N0 (linear), received power, slant range, fading ---
    for s = 1:length(sc.Satellites)
        % Rician fading
        K = 5;
        phi = 2*pi*rand;
        LOS = sqrt(K/(K+1));
        sigma = sqrt(1/(2*(K+1)));
        X = randn;
        Y = randn;
        h = LOS*exp(1j*phi) + sigma*(X + 1j*Y);
        h_vec(s) = abs(h)^2;
        
        % Nominal Eb/N0 from satellite toolbox
        ebno_nom_dB = ebno(links(s), t(k));
        if ebno_nom_dB == -Inf
            ebnr_lin_vec(s) = 0;
        else
            ebnr_lin_vec(s) = 10^(ebno_nom_dB/10) * abs(h)^2; % linear with fading
            numSat_underService = numSat_underService + 1;
        end
        
        % Received power and slant range
        [~,pr] = sigstrength(links(s), t(k));
        pr_vec(s) = pr;
        [~,~,r] = aer(sc.Satellites(s), gs, t(k));
        r_vec(s) = r;
    end
    
    % --- Step 2: compute combining methods ---
    if numSat_underService > 0
        % MRC (linear sum)
        snr_mrc_lin = sum(ebnr_lin_vec);
        snr_mrc(k) = 10*log10(snr_mrc_lin);
        
        % SC (select max)
        snr_sc_lin = max(ebnr_lin_vec);
        snr_sc(k) = 10*log10(snr_sc_lin);
        
        % Hybrid MRC (standard normalization)
        weights_std = (h_vec ./ r_vec.^2);
        weights_std = weights_std / sum(weights_std);
        snr_hybrid_std_lin = sum(weights_std .* ebnr_lin_vec);
        snr_hybrid_std(k) = 10*log10(snr_hybrid_std_lin);
        
        % Hybrid MRC with min(r) normalization
        snr_hybrid_minr_lin = sum((h_vec ./ r_vec.^2) .* ebnr_lin_vec) / min(r_vec)^2;
        snr_hybrid_minr(k) = 10*log10(snr_hybrid_minr_lin);
        
        % Average Eb/N0 and received power
        average_ebnr(k) = 10*log10(sum(ebnr_lin_vec)/numSat_underService);
        average_rp(k) = sum(pr_vec)/numSat_underService;
    end
    
    disp("Satellites servicing base station at time " + string(t(k)) + ": " + numSat_underService);
end

% --- Step 3: plotting ---
figure;
subplot(1,2,1);
plot(t,10.^(average_ebnr./10),'LineWidth',2);
title('Average Eb/N_0 (linear) vs time');
ylabel('Eb/N_0 (linear Watts)');
xlabel('Time');
grid on;

subplot(1,2,2);
plot(t,average_rp,'LineWidth',2);
title('Average received power vs time');
ylabel('Power (Watts)');
xlabel('Time');
grid on;

figure;
semilogy(t,snr_mrc,'LineWidth',1.5); hold on;
semilogy(t,snr_sc,'LineWidth',1.5);
semilogy(t,snr_hybrid_std,'LineWidth',1.5);
semilogy(t,snr_hybrid_minr,'--','LineWidth',1.5);
hold off;
grid on;
legend('MRC','SC','Hybrid std','Hybrid min(r)');
title('Comparison of combining methods');
xlabel('Time');
ylabel('Eb/N_0 [dB]');
