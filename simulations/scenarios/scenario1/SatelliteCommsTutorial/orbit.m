clear;clc;
Elem = struct;
Re = earthRadius;
h = 10000e3;
Elem.a = Re+h;
Elem.e = 0.5;
Elem.Inc = 30;
Elem.omega = 45;
Elem.RAAN = 15;
Elem.TA = 90;

Param = struct;
Param.mu = 3.9860e14;
Param.Periods = 1;
Param.dt = 10;
Param.T = 2*pi*sqrt(Elem.a^3/Param.mu);
Param.t = 0:Param.dt:Param.T;
OrbitB = GenOrbit(Param,Elem);
