% The range of Doppler frequencies caused 
% by multipath components, i.e., the 
% difference between the maximum and 
% minimum Doppler shifts of different 
% signal paths received at the same time.

% simulation implicitly assumes that the 
% VSAT (ground terminal) is located directly
% beneath the satellite’s orbital path, 
% i.e., at the nadir point of the satellite's
% pass, also called the sub-satellite point.

%   elevation
theta = 0:0.01:360;  % degrees
fo = 30e9; % 30 GHz (Ka-band)
G = 6.67430e-11; % m^3*kg^-1*s^-2
M_E = 5.972e24; % kg
R_E = 6371e3; % m
h = 600e3; % m
w_sat = sqrt(G*M_E/(R_E+h)^3);
c = 3*10^8; % m

fd = (fo/c)*w_sat*R_E*cosd(theta);
figure(1);
plot(theta,fd,'b','LineWidth',2);
grid on;
hold on;
max_doppler = max(fd); 
theta_max = theta(max_doppler==fd);
plot(theta_max,max_doppler,'b*','LineWidth',2);
min_doppler = min(fd);
theta_min = theta(min_doppler==fd);
plot(theta_min,min_doppler,'c*','LineWidth',2);

fo = 2e9; % 2 GHz (S-band)

fd = (fo/c)*w_sat*R_E*cosd(theta);
plot(theta,fd,'r','LineWidth',2);
max_doppler = max(fd); 
theta_max = theta(max_doppler==fd);
plot(theta_max,max_doppler,'r*','LineWidth',2);
min_doppler = min(fd);
theta_min = theta(min_doppler==fd);
plot(theta_min,min_doppler,'m*','LineWidth',2);
hold off;
axis tight;
legend('Ka-band','Ka-band Max', 'Ka-band Min','S-band','S-band Max','S-band Min','Location','Best');
