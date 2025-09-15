% Filename: animate_polarization.m
% Description: Animates the tip of the electric field vector for different polarizations

clear; clc; close all;

% Parameters
f = 1e9;              % Frequency (1 GHz)
lambda = 3e8 / f;     % Wavelength
omega = 2*pi*f;       % Angular frequency
t = linspace(0, 2*pi, 300);  % Time steps (simulate over one cycle)

% Choose polarization type
% Options: 'linear', 'circular', 'elliptical'
polarization_type = 'circular';  % change this to test different types

% Set amplitudes and phase shift
switch polarization_type
    case 'linear'
        E0x = 1;
        E0y = 1;
        delta = 0;  % in-phase
    case 'circular'
        E0x = 1;
        E0y = 1;
        delta = pi/2;  % 90° out of phase
    case 'elliptical'
        E0x = 1;
        E0y = 0.5;
        delta = pi/3;  % arbitrary phase offset
    otherwise
        error('Unknown polarization type');
end

% Prepare figure
figure('Color', 'w');
axis equal;
axis([-1.2 1.2 -1.2 1.2]);
grid on;
xlabel('E_x'); ylabel('E_y');
title(['Polarization Animation: ', polarization_type], 'FontSize', 14);
hold on;

% Plot zero marker
plot(0, 0, 'ko');

% Create animated vector and trace
vec = plot([0, 0], [0, 0], 'r', 'LineWidth', 2);     % field vector
trace = animatedline('Color', 'b', 'LineWidth', 1.5); % path trace

% Animation loop
for i = 1:length(t)
    Ex = E0x * cos(omega * t(i));
    Ey = E0y * cos(omega * t(i) + delta);
    
    % Update vector
    set(vec, 'XData', [0, Ex], 'YData', [0, Ey]);
    
    % Add point to trace
    addpoints(trace, Ex, Ey);
    
    drawnow;
    pause(0.01);
end
