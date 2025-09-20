%   This script checks alternative ways to simulate the scenario
clear;
clc;
%%  Scenario Setup
duration = minutes(1);
%   generate the scenario timeline
startTime = datetime(2025,9,11,9,16,1);
endTime = startTime + duration;
sampleTime = 10;

sc = satelliteScenario(startTime,endTime,sampleTime);

%   add the satellites in the scenario
sat = satellite(sc, "largeConstellation.tle");
% sat = satellite(sc,"leoSatelliteConstellation.tle");

%   add two ground stations
%   Minimum elevation between 10° and 20° is specified
%   for LEO applications for maximum visibility.
minElevationAngle = 10;
%   Calculated based on earth horizon inclination
minSatElevation = -15;
%%  Some antenna calculations
%   Calculated based on the VSAT specification of TR38.811
satDishDiameter = 0.5;
gsDishDiameter = 0.6;
%   calculated based on TR38.811
transmitterApertureEfficiency = 0.5881;
receiverApertureEfficiency = 0.85;
%   Antenna Gain is 43.2 dBi for Tx
%   Antenna Gain is 39.7 dB for Rx
rho_a = transmitterApertureEfficiency;
d = satDishDiameter;
c = 3e8;
f = 30e9;
lambda = c/f;
boresightGain = rho_a * ((pi * d / lambda)^2);
beamwidth_3dB = 70 * lambda / d;
boresightGain_dB = 10 * log10(boresightGain);
%%  Ground Station Configurations
name="MDSCC";
lat = 40.43139;
lon = -4.24806;
gsSource = groundStation(sc,...
    Name=name,Latitude=lat,Longitude=lon,...
    MinElevationAngle=minElevationAngle);

%   add another ground station
name2 = "Lamia Base Station";
lat2=38.875941904578895;
lon2=22.437526584061686;
gsTarget = groundStation(sc, Name=name2, Latitude=lat2, Longitude=lon2,...
    MinElevationAngle=minElevationAngle);

%   denote the transmitter configuration as a struct
%   this is an on satellite transmitter
txConfig = struct;
%   Frequencies as specified in TR38.811 Deployment D4
txConfig.ULFrequency = 30e9;   % Hz
txConfig.DLFrequency = 20e9;   % Hz
txConfig.satPower = 20;        % dBW  : based on matlab examples
txConfig.gsPower = 3;           % dBW : 33 dBm on TR38.811
txConfig.BitRate = 10;       % Mbps
txConfig.SystemLoss = 0;    % dB
txConfig.ChannelBandwidth = 16e6;  % Hz
% txConfig.ULWavelength = (3*10^8)/txConfig.ULFrequency;
% txConfig.DLWavelength = (3*10^8)/txConfig.DLFrequency;

%   denote the receiver configuration as a struct
rxConfig = struct;
%   maximum gain-to-noise-temperature in dB/K
%   In TR38.811 the G/T is specified as 18.5 dB/K
%   It will be changed later
rxConfig.MaxGByT = 18.5;
rxConfig.SystemLoss = 0;
rxConfig.PreReceiverLoss = 0;
%   required bit energy to noise power spectral density ratio in dB
rxConfig.RequiredEbNo = 0;

%   add the transmitters and receivers
%   1. configure an isotropic antenna receptor
%   isotropic = arrayConfig(Size=[1 1]);
%   adding gimbals to the base stations
gTxSource = gimbal(gsSource);
gRxTarget = gimbal(gsTarget);

%   Based on TR38.811 we have Ptx = 33 dBm or 3 dBW.
%   The txSource is on the uplink, meaning it uses 30 GHz (Ka-band)
txSource = transmitter(gTxSource,...
    Frequency=txConfig.ULFrequency,...
    Power=txConfig.gsPower,...
    SystemLoss=txConfig.SystemLoss,...
    BitRate=txConfig.BitRate);
gaussianAntenna(txSource,"DishDiameter",gsDishDiameter,"ApertureEfficiency",transmitterApertureEfficiency);

rxTarget = receiver(gRxTarget,...
    SystemLoss=rxConfig.SystemLoss,...
    PreReceiverLoss=rxConfig.PreReceiverLoss,...
    GainToNoiseTemperatureRatio=rxConfig.MaxGByT,...
    RequiredEbNo=rxConfig.RequiredEbNo);
gaussianAntenna(rxTarget,"DishDiameter",gsDishDiameter,"ApertureEfficiency",receiverApertureEfficiency);

%   add transmitters on the satellites
%   add two gimbals on each satellite
gSatTx = gimbal(sat);
gSatRx = gimbal(sat);

%   The initial transmission frequency is chosen as 30 GHz
%   The frequencies are used as:
%   Uplink: 30 GHz   (GroundStation-to-Satellite)
%   Downlink: 20 GHz (Satellite-to-GroundStation)
%   ISL: 30 GHz      (Satellite-to-Sattelite)
tx = transmitter(gSatTx,...
    Frequency=txConfig.ULFrequency,...
    Power=txConfig.satPower,...
    SystemLoss=txConfig.SystemLoss,...
    BitRate=txConfig.BitRate);
gaussianAntenna(tx,"ApertureEfficiency",transmitterApertureEfficiency,"DishDiameter",satDishDiameter);

%   add receivers on all satellites
rx = receiver(gSatRx,...
    SystemLoss=rxConfig.SystemLoss,...
    PreReceiverLoss=rxConfig.PreReceiverLoss,...
    GainToNoiseTemperatureRatio=rxConfig.MaxGByT,...
    RequiredEbNo=rxConfig.RequiredEbNo);
gaussianAntenna(rx,"ApertureEfficiency",receiverApertureEfficiency,"DishDiameter",satDishDiameter);


%% Carrier-to-noise Ratio object creation
cfg = satelliteCNRConfig;
cfg.TransmitterSystemLoss = txConfig.SystemLoss;
% cfg distance is going to be calculated during each step, based on the
% slant range
%   The frequency will be changing based on the current
cfg.Frequency = txConfig.ULFrequency/1e9;
cfg.MiscellaneousLoss = 0;
cfg.TransmitterAntennaGain = 43.2;  % dBi
cfg.GainToNoiseTemperatureRatio = rxConfig.MaxGByT;
cfg.ReceiverSystemLoss = rxConfig.SystemLoss;
cfg.BitRate = txConfig.BitRate;
cfg.Bandwidth = 16e6;

%%  begin the simulation
timeVec = startTime:seconds(sampleTime):endTime;
ebnr = zeros(1,length(timeVec));
%   iterate in all time instances of the simulation
links = [];
infos = [];
for t=1:length(timeVec)
    %   check visibility, we need aer information
    [~,sourceElevation] = aer(gsSource,sat,timeVec(t));
    [~,targetElevation] = aer(gsTarget,sat,timeVec(t));
    %   check visibility
    elSourceVisible = (sourceElevation>=minElevationAngle);
    elTargetVisible = (targetElevation>=minElevationAngle);
    %   for all time instances, determine the best
    %   satellite for initial access to constellation
    trueId = find(elSourceVisible == true);
    if isempty(trueId)
        disp("No route at "+string(timeVec(t)));
        continue;
    end
    %   determine the ranges of these satellites
    [~,~,range] = aer(sat(trueId),gsTarget,timeVec(t));
    %   determine the index of the element in range bearing the minimum
    %   value
    [~,minRangeId] = min(range);
    bestId = trueId(minRangeId);
    %   this is the index of the best satellite for initial access to the
    %   constellation. This will be the 1st hop in the path.
    %   The 'node' variable stores the first two nodes of the routing.
    %   First node: {gsSource, bestSat}
    %   Point the base station gimbal to the best satellite
    pointAt(gTxSource,sat(bestId));
    pointAt(gSatRx(bestId),gsSource);
    %   Update the node list
    nodes = {gsSource sat(bestId)};
    %   Update the link node list
    link_nodes = {txSource rx(bestId)};
    %   Now we can determine the remaining nodes in the path
    pathFound = false;
    while ~pathFound
        %   the index of current node is
        idCurrent = bestId;
        % disp("Satellite"+idCurrent);
        if elTargetVisible(idCurrent)
            %   This is the last satellite before the groundStation, 
            %   The frequency of this satellite must be changed. 
            %   Return all frequencies to 30GHz at the end of the loop
            %   Before the next iteration of the execution
            tx(idCurrent).Frequency = txConfig.DLFrequency;
            %   calculate losses
            cfg.TransmitterPower = txConfig.satPower;
            cfg.TransmitterAntennaGain = 43.2;
            cfg.Distance = range(minRangeId)/1000;
            cfg.Frequency = txConfig.DLFrequency/1e9;
            [cn, info] = satelliteCNR(cfg);
            infos = [infos info];
            nodes = {nodes{:} gsTarget};
            link_nodes = {link_nodes{:} tx(idCurrent) rxTarget};
            pointAt(gSatTx(idCurrent),gsTarget);
            pointAt(gRxTarget,sat(idCurrent));
            pathFound = true;
            lnk = link(link_nodes{:});
            linkStatus(lnk,timeVec(t));
            links = [links lnk];
            ebnr(t) = ebno(lnk,timeVec(t));
            %   change the frequency back before continuing
            tx(idCurrent).Frequency = txConfig.ULFrequency;
            continue;
        end
        [~,elevations] = aer(sat(idCurrent),sat,timeVec(t));
        elevations(idCurrent) = -90;
        %   find the best visible satellite
        s = elevations >= minSatElevation;
        trueId = find(s == true);
        [~,~,range] = aer(sat(trueId),gsTarget,timeVec(t));
        [~,minRangeId] = min(range);
        bestId = trueId(minRangeId);
        pointAt(gSatTx(idCurrent),sat(bestId));
        pointAt(gSatRx(bestId),sat(idCurrent));
        tempLink = link(tx(idCurrent),rx(bestId));
        nodes = {nodes{:} sat(bestId)};
        link_nodes = {link_nodes{:} tx(idCurrent) rx(bestId)};
    end
    disp("Route at "+string(timeVec(t))+": ");
    for i=nodes
        fprintf('%s--', i{:}.Name);
    end
    fprintf("\n");
end
receivedEbNo = zeros(1,length(infos));
for i=1:length(infos)
    receivedEbNo(i) = infos(i).ReceivedEbNo;
end
%   since ebnr is the final ebnr at the receiver, we present the final cn
figure;
plot(timeVec,ebnr,'LineWidth',1.5);
hold on;
plot(timeVec,receivedEbNo,'LineWidth',1.5);
title('E_b/N_0');
xlabel('Simulation time (s)');
ylabel('E_b/N_0 [dB]');
legend('ebno(link)','satelliteCNR(cfg)');
grid on;

%%  BER calculations
%   Using 8-QAM, 16-QAM and QPSK
figure;
%   8-QAM
M = 8;
k = log2(M);
K = 10;
hLOS = exp(1i*2*pi*rand(1,length(ebnr))); % Thank you Emil Bjornson!
hNLOS = sqrt(1/2)*(randn(1, length(ebnr)) + 1i*randn(1, length(ebnr)));
h = sqrt(K/(K + 1))*hLOS + sqrt(1/(K + 1))*hNLOS;
test_ebnr = 10*log10((abs(h).^2).*(10.^((ebnr-50)./10)));
numSymPerFrame = 100;

snrdB = convertSNR(test_ebnr, "ebno", "snr", BitsPerSymbol=k);
berEst = zeros(size(test_ebnr));

for n = 1:length(snrdB)
    numErrs = 0;
    numBits = 0;
    while numErrs < 200 && numBits < 1e6
        dataIn = randi([0 1],numSymPerFrame*k,1);
        dataSym = bit2int(dataIn, k);
        txSig = qammod(dataSym, M);
        rxSig = awgn(txSig,snrdB(n),'measured')./h(n);
        rxSym = qamdemod(rxSig,M);
        dataOut = int2bit(rxSym,k);
        nErrors = biterr(dataIn, dataOut);
        numErrs = numErrs + nErrors;
        numBits = numBits + numSymPerFrame*k;
    end
    berEst(n) = numErrs/numBits;
end
EbNoVec_theory = -10:30;
berTheory = berfading(EbNoVec_theory,"qam",M,1,10);


semilogy(test_ebnr,berEst,'*'); hold on;
semilogy(EbNoVec_theory,berTheory,'LineWidth',1.5);

%   QPSK
M = 4;
k = log2(M);

snrdB = convertSNR(test_ebnr, "ebno", "snr", BitsPerSymbol=k);
berEst = zeros(size(test_ebnr));

for n = 1:length(snrdB)
    numErrs = 0;
    numBits = 0;
    while numErrs < 200 && numBits < 1e6
        dataIn = randi([0 1],numSymPerFrame*k,1);
        dataSym = bit2int(dataIn, k);
        txSig = pskmod(dataSym, M, pi/M);
        rxSig = awgn(txSig,snrdB(n),'measured')./h(n);
        rxSym = pskdemod(rxSig,M,pi/M);
        dataOut = int2bit(rxSym,k);
        nErrors = biterr(dataIn, dataOut);
        numErrs = numErrs + nErrors;
        numBits = numBits + numSymPerFrame*k;
    end
    berEst(n) = numErrs/numBits;
end
EbNoVec_theory = -10:30;
berTheory = berfading(EbNoVec_theory,"psk",M,1,10);


semilogy(test_ebnr,berEst,'*'); hold on;
semilogy(EbNoVec_theory,berTheory,'LineWidth',1.5);

% 16-QAM
M = 16;
k = log2(M);

snrdB = convertSNR(test_ebnr, "ebno", "snr", BitsPerSymbol=k);
berEst = zeros(size(test_ebnr));

for n = 1:length(snrdB)
    numErrs = 0;
    numBits = 0;
    while numErrs < 200 && numBits < 1e6
        dataIn = randi([0 1],numSymPerFrame*k,1);
        dataSym = bit2int(dataIn, k);
        txSig = qammod(dataSym, M);
        rxSig = awgn(txSig,snrdB(n),'measured')./h(n);
        rxSym = qamdemod(rxSig,M);
        dataOut = int2bit(rxSym,k);
        nErrors = biterr(dataIn, dataOut);
        numErrs = numErrs + nErrors;
        numBits = numBits + numSymPerFrame*k;
    end
    berEst(n) = numErrs/numBits;
end
EbNoVec_theory = -10:30;
berTheory = berfading(EbNoVec_theory,"qam",M,1,10);


semilogy(test_ebnr,berEst,'*'); hold on;
semilogy(EbNoVec_theory,berTheory,'LineWidth',1.5); hold off;

grid on;
legend('Estimated BER', 'Theoretical BER');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate');
text = sprintf("%d-QAM BER Analysis",M);
title(text);
legend('8-QAM (Simulation)','8-QAM (Theoretical)',...
    'QPSK (Simulation)', 'QPSK (Theoretical)',...
    '16-QAM (Simulation)','16-QAM (Theoretical)');