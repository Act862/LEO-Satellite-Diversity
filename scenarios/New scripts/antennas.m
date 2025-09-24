% Radiation patterns of isotropic, omni-directional, and directional antennas (3D)
clear; clc; close all;

% Define spherical angles
theta = linspace(0, 2*pi, 200);   % azimuth
phi   = linspace(0, pi, 100);     % elevation
[theta, phi] = meshgrid(theta, phi);

%% Radiation patterns
% Isotropic
P_iso = ones(size(theta));

% Omni-directional (dipole-like: depends only on elevation)
P_omni = abs(sin(phi));

% Directional (cos^n(theta) type pattern, beam at theta=0)
n = 6; % larger n = narrower beam
P_dir = cos(phi).^n;
P_dir(P_dir < 0) = 0;

%% Convert spherical to Cartesian
[X_iso, Y_iso, Z_iso] = sph2cart(theta, pi/2 - phi, P_iso);
[X_omni, Y_omni, Z_omni] = sph2cart(theta, pi/2 - phi, P_omni);
[X_dir, Y_dir, Z_dir] = sph2cart(theta, pi/2 - phi, P_dir);

%% Plot
figure;

subplot(1,3,1);
surf(X_iso, Y_iso, Z_iso, P_iso, 'FaceAlpha',0.8); 
title('Isotropic'); axis equal; shading interp; colormap jet; 
xlabel('x'); ylabel('y'); zlabel('z'); grid on;

subplot(1,3,2);
surf(X_omni, Y_omni, Z_omni, P_omni, 'FaceAlpha',0.8);
title('Omni-directional'); axis equal; shading interp; colormap jet;
xlabel('x'); ylabel('y'); zlabel('z'); grid on;

subplot(1,3,3);
surf(X_dir, Y_dir, Z_dir, P_dir, 'FaceAlpha',0.8);
title('Directional'); axis equal; shading interp; colormap jet;
xlabel('x'); ylabel('y'); zlabel('z'); grid on;

sgtitle('3D Radiation Patterns');
