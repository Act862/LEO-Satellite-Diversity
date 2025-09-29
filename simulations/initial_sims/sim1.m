Fs = 1e6;
%   sampling rate
fd = 50e3;
%   dippler shift
t = (0:1/Fs:1e-3);
%   time vector
x = cos(2*pi*1e5*t);
%   data signal
x_doppler = x.*exp(1j*2*pi*fd*t);
%   application of doppler
%   shift in frequency is multiplication in time
%   by the exponent term
subplot(1,2,1);
plot(t,x);
grid on;
subplot(1,2,2);
plot(t,x_doppler);
grid on;