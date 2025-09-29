%   Simple channel simulation
%   Satellite is the transmitter and base station is the receiver
clear;clc;

startTime = datetime(2025,8,30,17,31,0);
endTime = startTime + minutes(30);
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

%   check during each time slot k, the satellites that are visible.
average_ebnr = zeros(1,length(t));
average_rp = zeros(1,length(t));
snr_mrc = zeros(1,length(t));
snr_sc = zeros(1,length(t));
snr_gac = zeros(1,length(t));
for k=1:length(t)
    ebnr_sum = 0;
    pr_sum = 0;
    ebnr_weighted_sum = 0;
    count_non_zeros_ebnr = 0;
    count_non_zeros_pr = 0;
    numSat_underService = 0;
    max_ebnr = -Inf;
    max_r = -Inf;
    for s=1:length(sc.Satellites)
        ebnr = ebno(links(s),t(k));
        [~,pr] = sigstrength(links(s),t(k));
        [~,~,r] = aer(sc.Satellites(s),gs,t(k));
        if ebnr == -Inf
            ebnr = 0;
        else
            count_non_zeros_ebnr = count_non_zeros_ebnr + 1;
            numSat_underService = numSat_underService + 1;
        end

        if pr == -Inf
            pr = 0;
        else
            count_non_zeros_pr = count_non_zeros_pr + 1;
        end
        if ebnr > max_ebnr
            max_ebnr = ebnr;
        end
        if r > max_r
            max_r = r;
        end
        pr_sum = pr_sum + pr;
        ebnr_sum = ebnr_sum + ebnr;
        ebnr_weighted_sum = ebnr_weighted_sum + r*ebnr;
    end
    snr_mrc(k) = ebnr_sum;
    snr_sc(k) = max_ebnr;
    snr_gac(k) = (1/max_r)*ebnr_weighted_sum;
    disp("The number of satellites servicing the base station at time "+ string(t(k)) +" is: "+numSat_underService);
    if count_non_zeros_ebnr > 0
        average_ebnr(k) = ebnr_sum/count_non_zeros_ebnr;
    end
    if count_non_zeros_pr > 0
        average_rp(k) = pr_sum/count_non_zeros_pr;
    end
end

%   display data
%   on linear scale
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

figure;
semilogy(t,snr_mrc,'LineWidth',1.5);
hold on;
semilogy(t,snr_sc,'LineWidth',1.5);
semilogy(t,snr_gac,'LineWidth',1.5);
hold off;
grid on;
legend('MRC','SC', 'GAC');
title('MRC vs SC vs GAC');
xlabel('simulation time (datetime)');
ylabel('Eb/N_0 [dBW]');