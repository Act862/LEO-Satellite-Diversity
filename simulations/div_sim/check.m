% Create satellite scenario
startTime = datetime(2025,4,25,0,0,0);
stopTime = startTime + hours(1);
sampleTime = 10;  % seconds

sc = satelliteScenario(startTime, stopTime, sampleTime);

% Add ground station
gs = groundStation(sc, 51.5, -0.12, Name="London");  % London location

% Add two satellites (LEO ~700 km altitude)
earthRadius = 6371e3;
altitude1 = 700e3;
altitude2 = 705e3;

semiMajor1 = earthRadius + altitude1;
semiMajor2 = earthRadius + altitude2;
ecc = 0;   % circular
incl = 98.6;  % degrees

sat1 = satellite(sc, semiMajor1, ecc, incl, 0, 0, 0, ...
    Name="Sat-1", OrbitPropagator="two-body");

sat2 = satellite(sc, semiMajor2, ecc, incl, 180, 0, 0, ...
    Name="Sat-2", OrbitPropagator="two-body");

% Add access analysis (line-of-sight)
access1 = access(sat1, gs);
access2 = access(sat2, gs);

% Visualize
play(sc);  % Launch 3D orbit viewer


% Get access status over time
t = access1.Times;
status1 = access1.Status;
status2 = access2.Status;

dualVisible = status1 & status2;

figure;
plot(t, status1, 'r', t, status2, 'b', t, dualVisible, 'g--', 'LineWidth', 2);
legend('Satellite 1', 'Satellite 2', 'Dual Visibility');
xlabel('Time'); ylabel('Visibility (1=Yes)');
title('Satellite Diversity Visibility Timeline');
grid on;
