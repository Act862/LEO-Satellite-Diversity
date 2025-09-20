clean; clc;
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
    3e8/terminalTxConfig.Frequency); % ratio
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
    3e8/terminalRxConfig.Frequency); % ratio
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
satTxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
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
    Power=satTxConfig.txPower);
%   Use gaussian Antennas for both transmitters
gaussianAntenna(sat1Tx,"DishDiameter",satTxConfig.DishDiameter,...
    "ApertureEfficiency",satTxConfig.ApertureEfficiency);
gaussianAntenna(sat2Tx,"DishDiameter",satTxConfig.DishDiameter,...
    "ApertureEfficiency",satTxConfig.ApertureEffiency);

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
gs1 = groundStation(sc,Name=name,Latitude=lat,Longitude=lon);

latitude = 52.2294963;                                             
longitude = 0.1487094;                                             
gs2 = groundStation(sc,latitude,longitude,Name="Ground Station 2");

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
gaussianAntenna(gs1Tx1,"DishDiameter",0.6,"ApertureEfficiency",0.5881);
gaussianAntenna(gs1Tx2,"DishDiameter",0.6,"ApertureEfficiency",0.5881);

%   turn antenna gain into the aperture efficiency
function ea = gain2apertureEfficiency(gain,lambda,dishDiameter)
    %   gain = (pi*diameter/lambda)^2*ea
    ea = gain/(pi*dishDiameter/lambda)^2;
end