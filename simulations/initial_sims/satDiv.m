% Satellite Diversity Animation in MATLAB

% Clear environment
clear; clc;

% Create satellite scenario
startTime = datetime(2021,12,10,18,27,57);
stopTime = startTime + hours(3);
sampleTime = 60;

sc = satelliteScenario(startTime, stopTime, sampleTime);
% Define Keplerian elements manually
semiMajorAxis = 6800e3; % meters
eccentricity = 0.001;
inclination = 98; % degrees
raan = 0; % Right Ascension of Ascending Node
argPeriapsis = 0;
trueAnomaly = 0;

sat1 = satellite(sc, semiMajorAxis, eccentricity, inclination, ...
    raan, argPeriapsis, trueAnomaly, "Name", "Sat-A");
sat2 = satellite(sc, semiMajorAxis, eccentricity, inclination + 10, ...
    raan + 60, argPeriapsis, trueAnomaly, "Name", "Sat-B");

% Define a ground station
gs = groundStation(sc, 37.7749, -122.4194); % San Francisco (as an example)

% Access analysis to check which satellite is visible when
access1 = access(sat1, gs);
access2 = access(sat2, gs);

% Setup visualization
v = satelliteScenarioViewer(sc);
v.ShowDetails = true;
v.PlaybackSpeedMultiplier = 10;  % Speed up playback

% Plot access lines dynamically
disp("Playing satellite diversity scenario...");

play(sc); % This will show live animation of both satellites and their coverage to the ground station
