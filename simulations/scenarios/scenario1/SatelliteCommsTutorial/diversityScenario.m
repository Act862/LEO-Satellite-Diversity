%% Clean setup
clear; clc;
%% Simulation scenario
startTime = datetime(2020,1,11,14,50,0);
stopTime = startTime + days(3);
sampleTime = 60;
sc = satelliteScenario(startTime, stopTime, sampleTime);
%% Add satellites
earthRadius = 6378100; % Define earthRadius if not already defined in your environment
semiMajorAxis = earthRadius + 600e3;
eccentricity = 0;
RAAN = 0;
AOP = 0;
trueAnomaly = 0;
% Satellite 1
sat1 = satellite(sc, semiMajorAxis, eccentricity, 60, RAAN, AOP, trueAnomaly, ...
    "OrbitPropagator", "sgp4", "Name", "sat1");
% Satellite 2
sat2 = satellite(sc, semiMajorAxis, eccentricity, 10, RAAN, AOP, trueAnomaly, ...
    "OrbitPropagator", "sgp4", "Name", "sat2");
%% Add gimbals and payloads to satellites
% Satellite 1
gimbalRxSat1 = gimbal(sat1, "Name", "Sat1 Rx Gimbal");
gimbalTxSat1 = gimbal(sat1, "Name", "Sat1 Tx Gimbal");
rxSat1 = receiver(gimbalRxSat1, "Name", "Satellite 1 Rx", "GainToNoiseTemperatureRatio", 5, "SystemLoss", 3);
txSat1 = transmitter(gimbalTxSat1, "Name", "Satellite 1 Tx", "Frequency", 30e9, "Power", 20, "BitRate", 20, "SystemLoss", 3);
% Satellite 2
gimbalRxSat2 = gimbal(sat2, "Name", "Sat2 Rx Gimbal");
gimbalTxSat2 = gimbal(sat2, "Name", "Sat2 Tx Gimbal");
rxSat2 = receiver(gimbalRxSat2, "Name", "Satellite 2 Rx", "GainToNoiseTemperatureRatio", 5, "SystemLoss", 3);
txSat2 = transmitter(gimbalTxSat2, "Name", "Satellite 2 Tx", "Frequency", 30e9, "Power", 20, "BitRate", 20, "SystemLoss", 3);
% Antennas
dishDiameter = 0.5;
apertureEfficiency = 0.5;
gaussianAntenna(txSat1, "DishDiameter", dishDiameter, "ApertureEfficiency", apertureEfficiency);
gaussianAntenna(txSat2, "DishDiameter", dishDiameter, "ApertureEfficiency", apertureEfficiency);
gaussianAntenna(rxSat1, "DishDiameter", dishDiameter, "ApertureEfficiency", apertureEfficiency);
gaussianAntenna(rxSat2, "DishDiameter", dishDiameter, "ApertureEfficiency", apertureEfficiency);
%% Ground stations
% Destination GS (receives from both satellites)
gs1 = groundStation(sc, 38.875941904578895, 22.437526584061686, "Name", "Lamia");
gimbalGs1 = gimbal(gs1, "Name", "Lamia GS Gimbal");
rxGs1 = receiver(gimbalGs1, "Name", "Ground Station Receiver", "RequiredEbNo", 14);
gaussianAntenna(rxGs1, "DishDiameter", 5);
% Source GS (sends to both satellites)
gs2 = groundStation(sc, 35.3113, 25.3130, "Name", "Gouves");
gimbalGs2_1 = gimbal(gs2, "Name", "Gimbal 1", "MountingLocation", [-1; 0; 0]);
gimbalGs2_2 = gimbal(gs2, "Name", "Gimbal 2", "MountingLocation", [1; 0; 0]);
txGs2_1 = transmitter(gimbalGs2_1, "Name", "GS2 Tx 1", "Power", 40, "BitRate", 20);
txGs2_2 = transmitter(gimbalGs2_2, "Name", "GS2 Tx 2", "Power", 40, "BitRate", 20);
gaussianAntenna(txGs2_1, "DishDiameter", 5);
gaussianAntenna(txGs2_2, "DishDiameter", 5);
%% Pointing
% GS2 Tx 1 points to Satellite 1
pointAt(gimbalGs2_1, sat1);
% Satellite 1 Rx points to GS2
pointAt(gimbalRxSat1, gs2);

% GS2 Tx 2 points to Satellite 2
pointAt(gimbalGs2_2, sat2);
% Satellite 2 Rx points to GS2
pointAt(gimbalRxSat2, gs2);

% Satellite 1 Tx points to GS1
pointAt(gimbalTxSat1, gs1);
% Satellite 2 Tx points to GS1
pointAt(gimbalTxSat2, gs1);

% GS1 Rx points to Satellite 1 (if receiving from Sat1) OR
% GS1 Rx points to Satellite 2 (if receiving from Sat2)
% Since it's receiving from both, the gimbal will point to whichever is visible
% and has a stronger signal if both are visible, or simply point to one
% if the link is defined as such.
% For simplicity in individual links, we'll assume the GS1 receiver is
% configured to receive from either sat1 or sat2 depending on the specific link.
% In reality, a single gimbal might track a primary target or require multiple receivers/gimbals.
% For the purpose of individual links, `link` function handles the "visibility".
pointAt(gimbalGs1, [sat1 sat2]); % GS1 gimbal points to both satellites to receive

%% Define Links (Single-hop for clarity as per your last instruction)
% Link from Gouves Tx1 to Satellite 1
lnk1_GS2_to_Sat1 = link(txGs2_1, rxSat1);

% Link from Satellite 1 Tx to Lamia GS
lnk1_Sat1_to_GS1 = link(txSat1, rxGs1);

% Link from Gouves Tx2 to Satellite 2
lnk2_GS2_to_Sat2 = link(txGs2_2, rxSat2);

% Link from Satellite 2 Tx to Lamia GS
lnk2_Sat2_to_GS1 = link(txSat2, rxGs1);

%% Access Intervals
disp('Link: Gouves Tx1 to Satellite 1');
disp(linkIntervals(lnk1_GS2_to_Sat1));
disp('Link: Satellite 1 Tx to Lamia GS');
disp(linkIntervals(lnk1_Sat1_to_GS1));
disp('Link: Gouves Tx2 to Satellite 2');
disp(linkIntervals(lnk2_GS2_to_Sat2));
disp('Link: Satellite 2 Tx to Lamia GS');
disp(linkIntervals(lnk2_Sat2_to_GS1));


%% Visualize
v = satelliteScenarioViewer(sc);
% Add access lines for better visualization of individual links
addAccess(lnk1_GS2_to_Sat1);
addAccess(lnk1_Sat1_to_GS1);
addAccess(lnk2_GS2_to_Sat2);
addAccess(lnk2_Sat2_to_GS1);

play(sc);