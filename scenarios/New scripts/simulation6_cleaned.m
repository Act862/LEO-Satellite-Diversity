clear; clc;
%   Reinvention of the Scenario6: Two Satellites and Two Ground Stations
%   The two satellites work as spatial diversity redundant resources.
%   The one base station sends data to both satellites, meaning it needs two gimbals.
%   Then both satellites send the data to the one gimbal of the second base station.
%   Setting up the simulation scenario.

startTime = datetime(2020,1,11,14,00,0);
stopTime = startTime + days(1);
sampleTime = 10;

sc = satelliteScenario(startTime,stopTime,sampleTime);

% Add the satellites to the scenarios
sat1 = satellite(sc,earthRadius+600e3,0,60,0,0,0, ...
    "Name","Sat1","OrbitPropagator","two-body-keplerian");
sat2 = satellite(sc,earthRadius+600e3,0,50,0,0,0, ...
    "Name","Sat2","OrbitPropagator","two-body-keplerian");

% Add the gimbals on the satellites
g1Sat1 = gimbal(sat1,"MountingLocation",[0;1;2]);
g2Sat1 = gimbal(sat1,"MountingLocation",[0;-1;2]);
g1Sat2 = gimbal(sat2,"MountingLocation",[0;1;2]);
g2Sat2 = gimbal(sat2,"MountingLocation",[0;-1;2]);

% Add the transmitters and receivers on each gimbal
% Gimbal 1 on each satellite is the transmitter
% Gimbal 2 on each satellite is the receiver
% We choose to have the same frequency on both satellites
sat1Tx = transmitter(g1Sat1,"MountingLocation",[0;0;1], ...
    Frequency=20e9,Power=20);
sat2Tx = transmitter(g1Sat2,"MountingLocation",[0;0;1], ...
    Frequency=20e9,Power=20);
gaussianAntenna(sat1Tx,"DishDiameter",0.5,"ApertureEfficiency",0.6447);
gaussianAntenna(sat2Tx,"DishDiameter",0.5,"ApertureEfficiency",0.6447);

sat1Rx = receiver(g2Sat1,"MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=13,RequiredEbNo=0);
sat2Rx = receiver(g2Sat2,"MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=13,RequiredEbNo=0);
gaussianAntenna(sat1Rx,"DishDiameter",0.33,"ApertureEfficiency",0.6578);
gaussianAntenna(sat2Rx,"DishDiameter",0.33,"ApertureEfficiency",0.6578);

% Adding the ground stations
name = "Lamia Base Station";
lat = 38.875941904578895;
lon = 22.437526584061686;
gs1 = groundStation(sc,Name=name,Latitude=lat,Longitude=lon);

latitude = 52.2294963;                                             
longitude = 0.1487094;                                             
gs2 = groundStation(sc,latitude,longitude,Name="Ground Station 2");

% Add two gimbals on each ground station
g1Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);

gs1Tx1 = transmitter(g1Gs1,"Name","GS1Tx1","MountingLocation",[0;0;1], ...
    "Frequency",30e9,"Power",3,"BitRate",10);
gs1Tx2 = transmitter(g2Gs1,"Name","GS1Tx2","MountingLocation",[0;0;1], ...
    "Frequency",30e9,"Power",3,"BitRate",10);
gaussianAntenna(gs1Tx1,"DishDiameter",0.6,"ApertureEfficiency",0.5881);
gaussianAntenna(gs1Tx2,"DishDiameter",0.6,"ApertureEfficiency",0.5881);

g1Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);

gs2Rx1 = receiver(g1Gs2,"Name","GS2Rx1","MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",18.5,"RequiredEbNo",0);
gs2Rx2 = receiver(g2Gs2,"Name","GS2Rx2","MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",18.5,"RequiredEbNo",0);
gaussianAntenna(gs2Rx1,"DishDiameter",0.6,"ApertureEfficiency",0.5872);
gaussianAntenna(gs2Rx2,"DishDiameter",0.6,"ApertureEfficiency",0.5872);

% Set the tracking of the targets for Gimbals
% GS1Tx1 points at Sat1Rx
% GS1Tx2 points at Sat2Rx
% Sat1Tx points at GS2Rx1
% Sat2Tx points at GS2Rx2
pointAt(g1Gs1,sat1);
pointAt(g2Gs1,sat2);
pointAt(g2Sat1,gs1);
pointAt(g2Sat2,gs1);
pointAt(g1Gs2,sat1);
pointAt(g2Gs2,sat2);
pointAt(g1Sat1,gs2);
pointAt(g1Sat2,gs2);

% Add Link Analysis
acc1 = access(gs1,sat1);
acc2 = access(gs1,sat2);
acc3 = access(gs2,sat1);
acc4 = access(gs2,sat2);

accint1 = accessIntervals(acc1);
accint2 = accessIntervals(acc2);
accint3 = accessIntervals(acc3);
accint4 = accessIntervals(acc4);

lnk1 = link(gs1Tx1,sat1Rx,sat1Tx,gs2Rx1);
lnk2 = link(gs1Tx2,sat2Rx,sat2Tx,gs2Rx2);
lnk3 = link(gs1Tx1,sat2Rx,sat2Tx,gs2Rx1);
lnk4 = link(gs1Tx2,sat1Rx,sat1Tx,gs2Rx2);

%   Also make distinct links for each pair
link11 = link(gs1Tx1,sat1Rx);
link12 = link(sat1Tx,gs2Rx1);
link21 = link(gs1Tx2,sat2Rx);
link22 = link(sat2Tx,gs2Rx2);
%   For the final SNR, only the last hop is necessary

lint1 = linkIntervals(lnk1);
lint2 = linkIntervals(lnk2);
lint3 = linkIntervals(lnk3);
lint4 = linkIntervals(lnk4);

% [e1,time1] = ebno(lnk1);
% [e2,time2] = ebno(lnk2);

%   Change the links
[e1,time1] = ebno(link12);
[e2,time2] = ebno(link22);
%   Note that time1 is the same as time2

% Plot Eb/No for each link
% t = intersect(time1,time2);
t = time1(e1 ~= -Inf & e2 ~= -Inf);
ebno_branch1 = e1(e1 ~= -Inf & e2 ~= -Inf);
ebno_branch2 = e2(e1 ~= -Inf & e2 ~= -Inf);

% Latency analysis
[delay, time] = latency(sat1,gs1); 
[latency2,time2] = latency(sat1,gs2);

ac = access(sat1,gs2);
intv = accessIntervals(ac);
ac2 = access(sat2,gs2);
intv2 = accessIntervals(ac2);

figure;hold on;
for i = 1:height(intv2)
    startTime = intv2.StartTime(i);
    endTime = intv2.EndTime(i);
    idx = (e1 ~= -Inf) & (e2 ~= -Inf);
    plot(time(idx), e1(idx),...
        'LineWidth', 1.5, 'DisplayName', sprintf('Pass (SC) %d', i));
    plot(time(idx), e2(idx),...
        'LineWidth', 1.5, 'DisplayName', sprintf('Pass (MRC) %d', i));
end
hold off;
axis tight;

figure;
plot(t,ebno_branch1,'LineWidth',1.5); hold on;
plot(t,ebno_branch2,'LineWidth',1.5); hold off;
grid on; axis tight;
ylabel('EbNo (dB)')
xlabel('Simulated time (datetime)');
legend('EbNo for Link1', 'EbNo for Link2');
title('Received Eb/N_o for each link');

% Combine the Eb/N0 from both links
ebno_sc = 10*log10(max(10.^(ebno_branch1./10),10.^(ebno_branch2./10)));
ebno_mrc = 10*log10(10.^(ebno_branch1./10) + 10.^(ebno_branch2./10));

figure;
plot(t,ebno_sc,t,ebno_mrc,'LineWidth',1.5);
grid on;
legend('SC','MRC');
title('Selection Combining VS Maximal Ratio Combining');
ylabel('EbNo (dB)');
xlabel('Simulated time (datetime)');
axis tight;

% Combined BER Calculation
figure;
ber_Theory = berfading(ebno_mrc,"psk",4,1,3);
% ber_qpsk = 0.5*(2*qfunc(sqrt(2*10.^(ebno_sc./10))) - (qfunc(sqrt(2*10.^(ebno_sc./10)))).^2);
ber_qpsk = berfading(ebno_sc,"psk",4,1,3);
ber_qpsk2 = berfading(ebno_mrc,"psk",4,1,3);
% ber_qpsk2 = 0.5*(2*qfunc(sqrt(2*10.^(ebno_mrc./10))) - (qfunc(sqrt(2*10.^(ebno_mrc./10)))).^2);
semilogy(t,ber_qpsk,'-x','LineWidth',1.5); hold on;
semilogy(t,ber_qpsk2,'->','LineWidth',1.5); hold off;
ylabel('Bit Error Rate');
xlabel('Simulation time (datetime)');
legend('SC', 'MRC');
title('Combined SNR Bit-Error Rate (QPSK)');
grid on;

figure; hold on;
avg_pass_duration = 0;
number_of_passes = 0;
avg_latency = 0;

for i = 1:height(intv)
    startTime = intv.StartTime(i);
    endTime = intv.EndTime(i);
    avg_pass_duration = avg_pass_duration + (endTime - startTime);
    number_of_passes = number_of_passes + 1;
    idx = (time >= startTime) & (time <= endTime) & (e1 ~= -Inf) & (e2 ~= -Inf);
    avg_latency = avg_latency + sum(delay(idx));
    plot(time(idx), 10*log10(max(10.^(e1(idx)./10),10.^(e2(idx)./10))),...
        'LineWidth', 1.5, 'DisplayName', sprintf('Pass (SC) %d', i));
    plot(time(idx), 10*log10(10.^(e1(idx)./10)+10.^(e2(idx)./10)),...
        'LineWidth', 1.5, 'DisplayName', sprintf('Pass (MRC) %d', i));
end

for i = 1:height(intv2)
    startTime = intv2.StartTime(i);
    endTime = intv2.EndTime(i);
    idx = (time >= startTime) & (time <= endTime);
    avg_latency = avg_latency + sum(latency2(idx));
end

avpd = avg_pass_duration/number_of_passes;
avl = avg_latency/number_of_passes;

xlabel('Time');
ylabel('Combined Eb/N_0 (dB)');
title('SC and MRC Eb/N_0 per Pass');
grid on; legend('show'); hold off;