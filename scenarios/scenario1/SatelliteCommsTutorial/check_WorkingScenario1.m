clear; clc;
%   Two satellites, two base stations.
%   GS1 transmits to both Sats (needs 2 gimbals).
%   Both Sats forward to GS2 (1 gimbal each).
%   We analyze link Eb/N0 and BER, pass by pass.

%% Simulation Setup
startTime = datetime(2020,1,11,14,00,0);
stopTime = startTime + days(1);
sampleTime = 10;  % Coarse sim time (fast). We'll refine with interpolation later.

sc = satelliteScenario(startTime, stopTime, sampleTime);

%% Add Satellites
earthRadius = 6371e3;
sat1 = satellite(sc, earthRadius+600e3,0,60,0,0,0, ...
    "Name","Sat1","OrbitPropagator","two-body-keplerian");
sat2 = satellite(sc, earthRadius+600e3,0,50,0,0,0, ...
    "Name","Sat2","OrbitPropagator","two-body-keplerian");

%% Add Gimbals (Satellites)
g1Sat1 = gimbal(sat1,"MountingLocation",[0;1;2]);
g2Sat1 = gimbal(sat1,"MountingLocation",[0;-1;2]);
g1Sat2 = gimbal(sat2,"MountingLocation",[0;1;2]);
g2Sat2 = gimbal(sat2,"MountingLocation",[0;-1;2]);

%% Tx/Rx on Satellites
sat1Tx = transmitter(g1Sat1,"MountingLocation",[0;0;1], ...
    Frequency=20e9,Power=15);
sat2Tx = transmitter(g1Sat2,"MountingLocation",[0;0;1], ...
    Frequency=20e9,Power=15);
gaussianAntenna(sat1Tx,"DishDiameter",0.5);
gaussianAntenna(sat2Tx,"DishDiameter",0.5);

sat1Rx = receiver(g2Sat1,"MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=3,RequiredEbNo=4);
sat2Rx = receiver(g2Sat2,"MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=3,RequiredEbNo=4);
gaussianAntenna(sat1Rx,"DishDiameter",0.5);
gaussianAntenna(sat2Rx,"DishDiameter",0.5);

%% Ground Stations
gs1 = groundStation(sc,Name="Lamia Base Station", ...
    Latitude=38.8759419,Longitude=22.4375266);

gs2 = groundStation(sc,Latitude=52.2294963,Longitude=0.1487094, ...
    Name="Ground Station 2");

%% Gimbals (Ground Stations)
% GS1: transmitter
g1Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);

gs1Tx1 = transmitter(g1Gs1,"Name","GS1Tx1","MountingLocation",[0;0;1], ...
    "Frequency",30e9,"Power",40,"BitRate",20);
gs1Tx2 = transmitter(g2Gs1,"Name","GS1Tx2","MountingLocation",[0;0;1], ...
    "Frequency",30e9,"Power",40,"BitRate",20);
gaussianAntenna(gs1Tx1,"DishDiameter",2);
gaussianAntenna(gs1Tx2,"DishDiameter",2);

% GS2: receiver
g1Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);

gs2Rx1 = receiver(g1Gs2,"Name","GS2Rx1","MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",3,"RequiredEbNo",1);
gs2Rx2 = receiver(g2Gs2,"Name","GS2Rx2","MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",3,"RequiredEbNo",1);
gaussianAntenna(gs2Rx1,"DishDiameter",2);
gaussianAntenna(gs2Rx2,"DishDiameter",2);

%% Pointing
pointAt(g1Gs1,sat1); pointAt(g2Gs1,sat2);
pointAt(g2Sat1,gs1); pointAt(g2Sat2,gs1);

pointAt(g1Gs2,sat1); pointAt(g2Gs2,sat2);
pointAt(g1Sat1,gs2); pointAt(g1Sat2,gs2);

%% Links
lnk1 = link(gs1Tx1,sat1Rx,sat1Tx,gs2Rx1);
lnk2 = link(gs1Tx2,sat2Rx,sat2Tx,gs2Rx2);

%% Loop over passes
lint1 = linkIntervals(lnk1);
lint2 = linkIntervals(lnk2);

numPasses = min(height(lint1),height(lint2));
avgBER_SC  = zeros(numPasses,1);
avgBER_MRC = zeros(numPasses,1);

% Create figures once
fig1 = figure; hold on; grid on; % Eb/No plot
fig2 = figure; hold on; grid on; set(gca,'YScale','log'); % BER plot

for i = 1:numPasses
    % Find overlap of i-th pass
    passStart = max(lint1.StartTime(i), lint2.StartTime(i));
    passEnd   = min(lint1.EndTime(i),   lint2.EndTime(i));
    if passStart >= passEnd
        continue; % no overlap
    end
    
    % Fine time grid
    common_time = (passStart:seconds(1):passEnd)';
    
    % Evaluate Eb/No for both links in this interval
    [e1, t1] = ebno(lnk1, passStart, passEnd);
    [e2, t2] = ebno(lnk2, passStart, passEnd);
    
    % Interpolate to common grid
    t1n = datenum(t1); t2n = datenum(t2); tcn = datenum(common_time);
    e1i_dB = interp1(t1n, e1, tcn, 'pchip', -Inf);
    e2i_dB = interp1(t2n, e2, tcn, 'pchip', -Inf);
    
    % Convert to linear
    E1_lin = 10.^(e1i_dB/10);
    E2_lin = 10.^(e2i_dB/10);
    
    % Combining
    SC_lin  = max(E1_lin,E2_lin);
    MRC_lin = E1_lin + E2_lin;
    SC_dB  = 10*log10(SC_lin);  SC_dB(SC_lin==0)=-Inf;
    MRC_dB = 10*log10(MRC_lin); MRC_dB(MRC_lin==0)=-Inf;
    
    % BER (BPSK)
    berSC  = qfunc(sqrt(2*SC_lin));
    berMRC = qfunc(sqrt(2*MRC_lin));
    
    % Store averages
    avgBER_SC(i)  = mean(berSC(berSC>0));
    avgBER_MRC(i) = mean(berMRC(berMRC>0));
    
    % --- Plot Eb/No for this pass
    figure(fig1);
    plot(common_time, SC_dB,'-','DisplayName',sprintf('SC pass %d',i));
    plot(common_time, MRC_dB,'--','DisplayName',sprintf('MRC pass %d',i));
    
    % --- Plot BER for this pass
    figure(fig2);
    semilogy(common_time, berSC,'-','DisplayName',sprintf('SC pass %d',i));
    semilogy(common_time, berMRC,'--','DisplayName',sprintf('MRC pass %d',i));
end

% Finalize Eb/No figure
figure(fig1);
ylabel('Eb/N_0 (dB)'); xlabel('Time');
title('Eb/N0 per pass');
legend('Location','best'); datetick('x','keeplimits');

% Finalize BER figure
figure(fig2);
ylabel('BER'); xlabel('Time');
title('BER per pass (BPSK)');
legend('Location','best'); datetick('x','keeplimits');

% Plot average BER per pass
figure;
bar(1:numPasses, [avgBER_SC avgBER_MRC]);
legend('SC','MRC'); xlabel('Pass #'); ylabel('Avg BER');
title('Average BER per pass');
grid on;