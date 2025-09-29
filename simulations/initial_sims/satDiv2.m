% Complex Satellite Diversity Visualization

clear; clc;
% Create satellite scenario
startTime = datetime(2025,1,1,0,0,0);
stopTime = datetime(2025,1,1,0,30,0);
sampleTime = 10;

% Create the scenario
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Adjust satellite inclination to match the ground station (~38° lat)
a = 6800e3;       % Semi-major axis in meters
e = 0.001;        % Small eccentricity
i = 40;           % Inclination (better for SF latitude)
raan1 = 0; 
raan2 = 5;        % Slight RAAN shift for diversity
argPer = 0;
trueAnom = 0;

% Create two closely spaced LEO satellites
sat1 = satellite(sc, a, e, i, raan1, argPer, trueAnom, "Name", "Sat-1");
sat2 = satellite(sc, a, e, i + 1, raan2, argPer, trueAnom, "Name", "Sat-2");

% Ground station in San Francisco
lat = 37.7749; lon = -122.4194; alt = 0;
gs_blocked = groundStation(sc, lat, lon, "Altitude", alt, ...
    "Name", "GS with Obstacle", "ElevationAngleThreshold", 15);

gs_clear = groundStation(sc, lat, lon, alt, ...
    "Name", "GS clear", "ElevationAngleThreshold", 5);

% Access setup
access1 = access(sat1, gs_blocked);  % Sat-1 affected by elevation mask
access2 = access(sat2, gs_clear);    % Sat-2 more available

% Viewer
v = satelliteScenarioViewer(sc);
v.ShowDetails = true;
v.PlaybackSpeedMultiplier = 10;

disp("Fixed: Satellites adjusted to better match ground station latitude.");
play(sc);
