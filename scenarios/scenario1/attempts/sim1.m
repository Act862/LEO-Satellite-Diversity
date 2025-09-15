%   WSN on the ground. Stationary UE and moving LEO satellite.
%   Assuming instance of LEO being present in the elevation angle.
%   Not working for coverage from other consecutively following sats.
%   land mobile-satellite (LMS) channel is considered RuralWooded
clear;clc;
%   Setting up commonParams
commonParams = struct;
commonParams.CarrierFrequency = 20e9;              % In Hz
commonParams.ElevationAngle = 0;                 % In degrees
commonParams.SatelliteAltitude = 600000;          % In m
commonParams.MobileAltitude = 0;                  % In m
commonParams.MobileSpeed = 0;                     % In m/s
commonParams.SampleRate = 7680000;                % In Hz
% Set the random stream and seed, for reproducibility
commonParams.RandomStream = "mt19937ar with seed";
commonParams.Seed = 73;
% Set the number of sinusoids used in generation of Doppler spread
commonParams.NumSinusoids = 48;
% Calculate the Doppler shift due to satellite movement
satelliteDopplerShift = dopplerShiftCircularOrbit(...
    commonParams.ElevationAngle,commonParams.SatelliteAltitude,...
    commonParams.MobileAltitude,commonParams.CarrierFrequency);
% Calculate the maximum Doppler shift due to mobile movement
c = physconst('lightspeed');
mobileMaxDoppler = commonParams.MobileSpeed*commonParams.CarrierFrequency/c;

%   Initialize the NTN flat fading narrowband channel
ntnNarrowbandChan = p681LMSChannel;
ntnNarrowbandChan.SampleRate = commonParams.SampleRate;
ntnNarrowbandChan.CarrierFrequency = commonParams.CarrierFrequency;
ntnNarrowbandChan.ElevationAngle = commonParams.ElevationAngle;
ntnNarrowbandChan.MobileSpeed = commonParams.MobileSpeed;
ntnNarrowbandChan.SatelliteDopplerShift = satelliteDopplerShift;
ntnNarrowbandChan.RandomStream = commonParams.RandomStream;
ntnNarrowbandChan.Seed = commonParams.Seed;
ntnNarrowbandChan.Environment = 'Rural';
ntnNarrowbandChan.AzimuthOrientation = 0;
ntnNarrowbandChan.FadingTechnique = "Sum of sinusoids";
ntnNarrowbandChan.NumSinusoids = commonParams.NumSinusoids;

%   Generate random input
rng(commonParams.Seed);
in = randn(commonParams.SampleRate,1,'like',1i);

[narrowbandOut,narrowbandPathGains,narrowbandSampleTimes]=...
    ntnNarrowbandChan(in);

%   Plot the received spectrum of the faded signal from the NTN flat fading
%   narrowband channel.
ntnNarrowbandAnalyzer = spectrumAnalyzer(...
    SampleRate=ntnNarrowbandChan.SampleRate);
ntnNarrowbandAnalyzer.Title = "Received Signal Spectrum " ...
    + "NTN narrowband with " + string(ntnNarrowbandChan.Environment) + " environment";
ntnNarrowbandAnalyzer.ShowLegend = true;
ntnNarrowbandAnalyzer.ChannelNames = "Rx Antenna 1";
ntnNarrowbandAnalyzer(narrowbandOut)