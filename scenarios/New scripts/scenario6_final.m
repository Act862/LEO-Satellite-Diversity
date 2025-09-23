clear; clc;
%%  Configurations Section
%   Transceivers specs for the terminals
terminalTxConfig = struct;
terminalTxConfig.Frequency = 3e10; % Hz
terminalTxConfig.TxPower = 3; % dBW
terminalTxConfig.DishDiameter = 0.6; % meters
terminalTxConfig.AntennaGain = 43.2; %  dBi
terminalTxConfig.TransmitSystemLoss = 0; % dB
terminalTxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
    terminalTxConfig.AntennaGain,...
    terminalTxConfig.DishDiameter,...
    physconst('LightSpeed')/terminalTxConfig.Frequency); % ratio
terminalTxConfig.AntennaTemp = 150; % Kelvin
terminalTxConfig.NoiseFigure = 1.2; % dB

%   On the downlink, the receiver works on 20GHz
terminalRxConfig = struct;
terminalRxConfig.Frequency = 2e10;  % Hz
terminalRxConfig.RequiredEbNo = 0;  % dB
terminalRxConfig.MaxGainToNoiseTemperatureRatio = 18.5; % dB/K
terminalRxConfig.ReceiverSystemLoss = 0; % dB
terminalRxConfig.AntennaGain = 39.7; % dBi
terminalRxConfig.DishDiameter = 0.6; % meters
terminalRxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
    terminalRxConfig.AntennaGain,...
    terminalRxConfig.DishDiameter,...
    physconst('LightSpeed')/terminalRxConfig.Frequency); % ratio
terminalRxConfig.AntennaTemp = 150; % kelvin
terminalRxConfig.NoiseFigure = 1.2; % dB
terminalRxConfig.BitRate = 10; % Mbps

%   Transceivers for Satellites
%   The satellites work on the downlink
satTxConfig = struct;
satTxConfig.Frequency = 2e10; % Hz
satTxConfig.DishDiameter = 0.5; % meters
satTxConfig.AntennaGain = 38.5; % dBi
satTxConfig.TxPower = 20; % dB
satTxConfig.EIRPDensity = 4; % dBW/MHz
satTxConfig.ApertureEfficiency = gain2apertureEfficiency(...
    satTxConfig.AntennaGain,...
    satTxConfig.DishDiameter,...
    physconst('LightSpeed')/satTxConfig.Frequency); % ratio
satTxConfig.halfBeamwidth = 1.7647; % degrees
satTxConfig.satelliteBeamDiameter = 20e3; % meters

satRxConfig = struct;
satRxConfig.Frequency = 3e10;
satRxConfig.DishDiameter = 0.33;
satRxConfig.RequiredEbNo = 0;
satRxConfig.MaxGainToNoiseTemperatureRatio = 13;
satRxConfig.AntennaGain = 38.5;
satRxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
    satRxConfig.AntennaGain,...
    satRxConfig.DishDiameter,...
    physconst('LightSpeed')/satRxConfig.Frequency);

%%  Scenario Section
%   Scenario Simulation Time
startTime = datetime(2020,1,11,14,0,0);
stopTime = startTime + days(1);
sampleTime = 10;

%   Scenario Generation
sc = satelliteScenario(startTime,stopTime,sampleTime);

%%  Satellites Section
%   Add satellites to the scenario
sat1 = satellite(sc,earthRadius+600e3,0,60,0,0,0, ...
    "Name","Sat1","OrbitPropagator","two-body-keplerian");
sat2 = satellite(sc,earthRadius+600e3,0,50,0,0,0, ...
    "Name","Sat2","OrbitPropagator","two-body-keplerian");

% Add the gimbals on the satellites
g1Sat1 = gimbal(sat1,"MountingLocation",[0;1;2]);
g2Sat1 = gimbal(sat1,"MountingLocation",[0;-1;2]);
g1Sat2 = gimbal(sat2,"MountingLocation",[0;1;2]);
g2Sat2 = gimbal(sat2,"MountingLocation",[0;-1;2]);

%   Add the transmitters of satellites
sat1Tx = transmitter(g1Sat1,...
    "MountingLocation",[0;0;1], ...
    Frequency=satTxConfig.Frequency,...
    Power=satTxConfig.TxPower);

sat2Tx = transmitter(g1Sat2,...
    "MountingLocation",[0;0;1], ...
    Frequency=satTxConfig.Frequency,...
    Power=satTxConfig.TxPower);
%   Use gaussian Antennas for both transmitters
gaussianAntenna(sat1Tx,"DishDiameter",satTxConfig.DishDiameter,...
    "ApertureEfficiency",satTxConfig.ApertureEfficiency);
gaussianAntenna(sat2Tx,"DishDiameter",satTxConfig.DishDiameter,...
    "ApertureEfficiency",satTxConfig.ApertureEfficiency);

%   Add receivers to the gimbals of the ground station 2
sat1Rx = receiver(g2Sat1,"MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=satRxConfig.MaxGainToNoiseTemperatureRatio,...
    RequiredEbNo=satRxConfig.RequiredEbNo);
sat2Rx = receiver(g2Sat2,"MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=satRxConfig.MaxGainToNoiseTemperatureRatio,...
    RequiredEbNo=satRxConfig.RequiredEbNo);
gaussianAntenna(sat1Rx,"DishDiameter",satRxConfig.DishDiameter,...
    "ApertureEfficiency",satRxConfig.ApertureEfficiency);
gaussianAntenna(sat2Rx,"DishDiameter",satRxConfig.DishDiameter,...
    "ApertureEfficiency",satRxConfig.ApertureEfficiency);

%%  Ground Stations Section
%   Add the ground stations
name = "Lamia Base Station";
lat = 38.875941904578895;
lon = 22.437526584061686;
gs1 = groundStation(sc,Name=name,Latitude=lat,Longitude=lon,...
    MinElevationAngle=5);

latitude = 52.2294963;                                             
longitude = 0.1487094;                                             
gs2 = groundStation(sc,latitude,longitude,Name="Ground Station 2",...
    MinElevationAngle=5);

% Add two gimbals on each ground station
g1Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);
g1Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);

%   Add transmitters on ground station 1 (GS1)
%   Add receivers on ground station 2 (GS2)
gs1Tx1 = transmitter(g1Gs1,...
    "Name","GS1Tx1",...
    "MountingLocation",[0;0;1], ...
    "Frequency",terminalTxConfig.Frequency,...
    "Power",terminalTxConfig.TxPower,...
    "BitRate",10);

gs1Tx2 = transmitter(g2Gs1,"Name","GS1Tx2","MountingLocation",[0;0;1], ...
    "Frequency",30e9,"Power",3,"BitRate",10);
gaussianAntenna(gs1Tx1,"DishDiameter",terminalTxConfig.DishDiameter,...
    "ApertureEfficiency",terminalTxConfig.ApertureEfficiency);
gaussianAntenna(gs1Tx2,"DishDiameter",terminalTxConfig.DishDiameter,...
    "ApertureEfficiency",terminalTxConfig.ApertureEfficiency);

gs2Rx1 = receiver(g1Gs2,"Name","GS2Rx1",...
    "MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",terminalRxConfig.MaxGainToNoiseTemperatureRatio,...
    "RequiredEbNo",terminalRxConfig.RequiredEbNo);

gs2Rx2 = receiver(g2Gs2,"Name","GS2Rx2",...
    "MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",terminalRxConfig.MaxGainToNoiseTemperatureRatio,...
    "RequiredEbNo",terminalRxConfig.RequiredEbNo);
gaussianAntenna(gs2Rx1,"DishDiameter",terminalRxConfig.DishDiameter,...
    "ApertureEfficiency",terminalRxConfig.ApertureEfficiency);
gaussianAntenna(gs2Rx2,"DishDiameter",terminalRxConfig.DishDiameter,...
    "ApertureEfficiency",terminalRxConfig.ApertureEfficiency);

%   Set the tracking of the targets from the gimbals
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

%%  The Access and Link Objects can be created
%   Access objects refer to the visibility of the satellites
access1 = access(gs1,sat1);
access2 = access(gs1,sat2);
access3 = access(gs2,sat1);
access4 = access(gs2,sat2);
accessInterval1 = accessIntervals(access1);
accessInterval2 = accessIntervals(access2);
accessInterval3 = accessIntervals(access3);
accessInterval4 = accessIntervals(access4);

%   Link objects model the link between two transponders
%   Create all the link objects
%   Part1: Links between Ground Station Txs and Satellites (UL)
link1 = link(gs1Tx1,sat1Rx);
link2 = link(gs1Tx1,sat2Rx);
link3 = link(gs1Tx2,sat1Rx);
link4 = link(gs1Tx2,sat2Rx);
%   Part2: Links between Satellites and Ground Station Receivers (DL)
link5 = link(sat1Tx,gs2Rx1);
link6 = link(sat1Tx,gs2Rx2);
link7 = link(sat2Tx,gs2Rx1);
link8 = link(sat2Tx,gs2Rx2);

[ebnr1,time1] = ebno(link5);
[ebnr2,time2] = ebno(link8);

%   Display EbNo from both links
figure;
plot(time1,ebnr1,'LineWidth',1.5);
hold on;
plot(time2,ebnr2,'LineWidth',1.5);
hold off;
grid on;
title("Received E_b/N_o from both satellites");
ylabel('E_b/N_o [dB]');
xlabel('Simulation Time (datetime)');
axis tight;

%%  Get geometrical data for p618 losses
ebnr1_lossy = zeros(1,length(ebnr1));
for n=1:length(time1)
    % Calculate geometry in respect of the ground station
    [~,el,r] = aer(gs2Rx1,sat1Tx,time1(n));
    if el < 5 || el > 175
        ebnr1_lossy(n) = ebnr1(n);
        continue;
    end
    p618cfg = p618Config(AntennaDiameter=satTxConfig.DishDiameter,...
        AntennaEfficiency=satTxConfig.ApertureEfficiency,...
        Frequency=satTxConfig.Frequency,...
        Latitude=gs2.Latitude,...
        Longitude=gs2.Longitude,...
        ElevationAngle=el);
    p618loss = p618PropagationLosses(p618cfg);
    atmoLoss = p618loss.Ac + p618loss.Ag + p618loss.Ar + ...
        p618loss.As + p618loss.At;
    ebnr1_lossy(n) = ebnr1(n) - atmoLoss;
end

ebnr2_lossy = zeros(1,length(ebnr2));
for n=1:length(time2)
    % Calculate geometry in respect of the ground station
    [~,el,r] = aer(gs2Rx2,sat2Tx,time2(n));
    if el < 5 || el > 175
        ebnr2_lossy(n) = ebnr2(n);
        continue;
    end
    p618cfg = p618Config(AntennaDiameter=satTxConfig.DishDiameter,...
        AntennaEfficiency=satTxConfig.ApertureEfficiency,...
        Frequency=satTxConfig.Frequency,...
        Latitude=gs2.Latitude,...
        Longitude=gs2.Longitude,...
        ElevationAngle=el);
    p618loss = p618PropagationLosses(p618cfg);
    atmoLoss = p618loss.Ac + p618loss.Ag + p618loss.Ar + ...
        p618loss.As + p618loss.At;
    ebnr2_lossy(n) = ebnr2(n) - atmoLoss;
end

figure;
plot(time1,ebnr1_lossy,"LineWidth",1.5);
hold on;
plot(time2,ebnr2_lossy,"LineWidth",1.5);
grid on;
axis tight;
title("Received E_b/N_o from both satellites + p618 losses");
ylabel('E_b/N_o [dB]');
xlabel('Simulation Time (datetime)');

%%  Combining at the Receiver
%   Maximal Ratio Combining
%   Selection Combining
ebnr_mrc = 10*log10(10.^(ebnr1_lossy./10) + 10.^(ebnr2_lossy./10));
ebnr_sc = 10*log10(max(10.^(ebnr1_lossy./10),10.^(ebnr2_lossy./10)));

figure;
plot(time1,ebnr_mrc);hold on;
plot(time1,ebnr_sc);hold off;
grid on;
legend('MRC','SC');
title('Selection Combining VS Maximal Ratio Combining');
ylabel('EbNo (dB)');
xlabel('Simulated time (datetime)');
axis tight;

%%  Helper Functions
%   turn antenna gain into the aperture efficiency
function ea = gain2apertureEfficiency(gain,dishDiameter,lambda)
    %   gain = (pi*diameter/lambda)^2*ea
    ea = gain/(pi*dishDiameter/lambda)^2;
end