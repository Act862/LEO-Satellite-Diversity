%   Simple channel simulation
%   Satellite is the transmitter and base station is the receiver
clear;clc;

startTime = datetime(2025,8,30,16,2,0);
endTime = startTime + days(1);
sampleTime = 30;                    % seconds

% startTime = datetime(2025,8,30,16,2,0);
% endTime = startTime + minutes(12);
t = startTime:seconds(sampleTime):endTime;

%   Initialize the satellite scenario
sc = satelliteScenario(startTime,endTime,sampleTime);

%   Add satellites
sat = satellite(sc,"leoSatelliteConstellation.tle");
%   add tx on each satellite
satTxs = [];
gTxs = [];
systemLoss = 3;
dishDiameter = 0.5;
apertureEfficiency = 0.5;

for s=sc.Satellites
    gimbalTx = gimbal(s);
    gTxs = [gTxs gimbalTx];
    tx = transmitter(gimbalTx,Name="TxSat"+s.ID,Frequency=30e9,...
        Power=20,BitRate=20,SystemLoss=systemLoss);
    gaussianAntenna(tx,DishDiameter=dishDiameter,...
        ApertureEfficiency=apertureEfficiency);
    satTxs = [satTxs tx];
end

%   Add ground stations
gsLat = 37.9838;     % deg
gsLon = 23.7275;     % deg
gsAlt = 0.1;         % km
gs = groundStation(sc,gsLat,gsLon,"Altitude",gsAlt,"MinElevationAngle",...
    5,"Name","Base Station");
gimbalRx = gimbal(gs);
gainToNoiseTemperatureRatio = 5;
rx = receiver(gimbalRx,Name="Rx_GS", ...
    GainToNoiseTemperatureRatio=gainToNoiseTemperatureRatio,...
    SystemLoss = systemLoss);

%   make all the satellites point at the base station
for s=1:length(sc.Satellites)
    pointAt(gTxs(s), gs);
end

%   generate all the links
links = [];
for s=1:length(sc.Satellites)
    %   for each satellite
    %   make the base station point at it
    pointAt(gimbalRx,sat(s));
    %   get the link data for the whole simulation
    lnk = link(satTxs(s),rx);
    links = [links lnk];
end

%   check during each time slot k, the satellites that are visible.
%   the satellites that are visible are the ones that have closed the link!
%   each link is a channel!
%   The Combining should be done inside this simulation.
average_ebnr = zeros(1,length(t));
average_rp = zeros(1,length(t));

for k=1:length(t)
    ebnr_sum = 0;
    pr_sum = 0;
    count_non_zeros_ebnr = 0;
    count_non_zeros_pr = 0;
    numSat_underService = 0;
    for s=1:length(sc.Satellites)
        ebnr = ebno(links(s),t(k));
        [~,pr] = sigstrength(links(s),t(k));
        if ebnr == -Inf
            ebnr = 0;
        else
            count_non_zeros_ebnr = count_non_zeros_ebnr + 1;
            numSat_underService = numSat_underService + 1;
        end
        
        if pr == -Inf
            pr = 0;
        else
            count_non_zeros_pr = count_non_zeros_pr + 1;
        end
        pr_sum = pr_sum + pr;
        ebnr_sum = ebnr_sum + ebnr;
    end
    disp("The number of satellites servicing the base station at time "+ string(t(k)) +" is: "+numSat_underService);
    average_ebnr(k) = ebnr_sum/count_non_zeros_ebnr;
    average_rp(k) = pr_sum/count_non_zeros_pr;
end

%   display data
%   on linear scale
figure;
subplot(1,2,1);
plot(t,10.^(average_ebnr./10),'LineWidth',2);
title('average Eb/N_0 against time');
ylabel('Eb/N_0 (linear scale) Watts');
xlabel('simulation time (datetime)');
grid on;

subplot(1,2,2);
plot(t,10.^(average_rp./10),'LineWidth',2);
title('average receiver power against time');
ylabel('P (linear scale) Watts');
xlabel('simulation time (datetime)');
grid on;

disp("The number of satellites servicing the base station is: "+numSat_underService/length(t));