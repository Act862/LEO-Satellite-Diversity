%   Reinforcement Learning Model for Switching between Diversity
%   Multiple Satellites -- One GS
%   One Antenna on each satellite -- Multiple Antennas on one GS
clear;clc;
%   No simulation like, only sets
%   Each satellite is a struct
sat1 = struct;
sat1.Frequency = 20e9; %    Downlink Frequency of 20 GHz
sat1.transmitPower = 3; % 33dBm based on 38.811, but 3 dBW
sat1.nofAntenna = 1;
sat1.antennaDiameter = 0.6; % 60 cm dish
sat1.antennaTxGain = 43.2; % dBi
sat1.apertureEfficiency = 0.5881;

%   Satellite 2 configuration
sat2 = struct;
sat2.Frequency = 20e9; %    Downlink Frequency of 20 GHz
sat2.transmitPower = 3; % 33dBm based on 38.811, but 3 dBW
sat2.nofAntenna = 1;
sat2.antennaDiameter = 0.6; % 60 cm dish
sat2.antennaTxGain = 43.2; % dBi
sat2.apertureEfficiency = 0.5881;

%   Ground Station Configuration
gs = struct;
gs.GT = 18.5; % Gain-to-Noise-Temperature Ratio
gs.BitRate = 10e6; % bps

%   The link budget will be calculated from the satelliteCNR function
%   for each satellite separately
%   Path 1
path1 = satelliteCNRConfig;
path1.TransmitterPower = sat1.transmitPower;
%   To be calculated
path1.TransmitterSystemLoss = 0;
path1.TransmitterAntennaGain = sat1.antennaTxGain;
path1.Distance = 1200;
path1.Frequency = sat1.Frequency/1e9;
path1.BitRate = gs.BitRate/1e6;
%   To be calculated
path1.MiscellaneousLoss = 0;
path1.ReceiverSystemLoss = 0;
path1.GainToNoiseTemperatureRatio = gs.GT;

%   Path 2
path2 = satelliteCNRConfig;
path2.TransmitterPower = sat2.transmitPower;
%   To be calculated
path2.TransmitterSystemLoss = 0;
path2.TransmitterAntennaGain = sat2.antennaTxGain;
path2.Distance = 1100;
path2.Frequency = sat2.Frequency/1e9;
path2.BitRate = gs.BitRate/1e6;
%   To be calculated
path2.MiscellaneousLoss = 0;
path2.ReceiverSystemLoss = 0;
path2.GainToNoiseTemperatureRatio = gs.GT;

[cn1, info1] = satelliteCNR(path1);
[cn2, info2] = satelliteCNR(path2);