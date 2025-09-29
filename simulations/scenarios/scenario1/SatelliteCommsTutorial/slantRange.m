Re = earthRadius;
el = 0:0.5:90; %    degrees (°)
el_rad = deg2rad(el);
%   LEO satellite elevation
h = 600e3; % 600 km


%   slant range calculation
sR = sqrt(Re^2 + (Re+h)^2 - 2*Re*(Re+h)*sin(el_rad+asin(Re*cos(el_rad)/(Re+h))));

subplot(1,3,1);
plot(el_rad, sR, 'LineWidth',2);
grid on;
axis tight;
title('Slant Range over elevation angle');
xlabel('Elevation Angle (radians)');
ylabel('Slant Range (meters)');
Fc = 30; % GHz

FSPL = 92.45 + 20*log10(sR) + 20*log10(Fc);
subplot(1,3,2);
semilogy(el_rad,FSPL, 'LineWidth',2);
grid on;
axis tight;
title('Pathloss over elevation angle');
xlabel('Elevation Angle (radians)');
ylabel('Pathloss (dBW)');

subplot(1,3,3);
c = 3e8;
tprop = sR/c;
plot(el_rad,tprop, 'LineWidth',2);
grid on;
axis tight;
title('Propagation delay over elevation angle');
xlabel('Elevation Angle (radians)');
ylabel('Propagation Delay (seconds)');
