clear;clc;
h_sat = 600; % km
v_sat = 7.5; % km/s
fc = 2*10^9; % 1/s
R = 6378; % km
c = 3*10^5; % km/s
v_ter = 1000; % km
h_ter = 0;
t = 0:0.01:40;
g_sat = getGamma(R,h_sat);
g_ter = getGamma(R,h_ter);

angleU_sat = getUang(v_sat, t, R, h_sat);
angleU_ter = getUang(v_ter, t, R, h_ter);

fd_sat = (fc*v_sat/c)*sin(angleU_sat)./sqrt(1+g_sat^2-2*g_sat*cos(angleU_sat));
fd_ter = (fc*v_ter/c)*sin(angleU_ter)./sqrt(1+g_ter^2-2*g_ter*cos(angleU_ter));

plot(t,fd_sat+fd_ter,'LineWidth',2);
grid on;
title('Doppler Shift for LEO Satellite (2GHz)');
xlabel('time (s)');
ylabel('doppler shift (Hz)');


