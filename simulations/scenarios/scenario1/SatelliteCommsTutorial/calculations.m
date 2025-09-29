% Create a satellite scenario
sc = satelliteScenario;

% Add a satellite (example: circular orbit at 700 km altitude)
sat = satellite(sc, 7000e3, 0, 0);  % 7000 km from Earth's center

% Get the satellite's position over time
positions = states(sat, 'CoordinateFrame', 'geographic'); % [lat; lon; alt]

% Extract the altitude (in meters)
altitude = positions(3, :);  % Altitude is the third row

% Plot altitude over time
time = sc.StartTime + seconds(0:length(altitude)-1);  % Time vector
plot(time, altitude / 1000);  % Altitude in km
xlabel('Time');
ylabel('Altitude (km)');
title('Satellite Altitude vs Time');
grid on;
