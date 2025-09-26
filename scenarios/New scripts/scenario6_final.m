clear; clc;
%%  Configurations Section
%   Transceivers specs for the terminals
terminalTxConfig = struct;
terminalTxConfig.Frequency = 3e10; % Hz
terminalTxConfig.TxPower = 3; % dBW
terminalTxConfig.DishDiameter = 0.6; % meters
terminalTxConfig.AntennaGain = 43.2; %  dBi
terminalTxConfig.SystemLoss = 3; % dB
terminalTxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
    terminalTxConfig.AntennaGain,...
    terminalTxConfig.DishDiameter,...
    physconst('LightSpeed')/terminalTxConfig.Frequency); % ratio
terminalTxConfig.AntennaTemp = 150; % Kelvin
terminalTxConfig.NoiseFigure = 1.2; % dB
terminalTxConfig.BitRate = 10; % Mbps

%   On the downlink, the receiver works on 20GHz
terminalRxConfig = struct;
terminalRxConfig.Frequency = 2e10;  % Hz
terminalRxConfig.RequiredEbNo = 0;  % dB
terminalRxConfig.MaxGainToNoiseTemperatureRatio = 18.5; % dB/K
terminalRxConfig.SystemLoss = 3; % dB
terminalRxConfig.AntennaGain = 39.7; % dBi
terminalRxConfig.DishDiameter = 0.6; % meters
terminalRxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
    terminalRxConfig.AntennaGain,...
    terminalRxConfig.DishDiameter,...
    physconst('LightSpeed')/terminalRxConfig.Frequency); % ratio
terminalRxConfig.AntennaTemp = 150; % kelvin
terminalRxConfig.NoiseFigure = 1.2; % dB
terminalRxConfig.BitRate = 10; % Mbps

%   Transceivers for Satellites
%   The satellites work on the downlink
satTxConfig = struct;
satTxConfig.Frequency = 2e10; % Hz
satTxConfig.DishDiameter = 0.5; % meters
satTxConfig.AntennaGain = 38.5; % dBi
satTxConfig.TxPower = 7; % dB
satTxConfig.SystemLoss = 3; % dB
satTxConfig.EIRPDensity = 4; % dBW/MHz
satTxConfig.ApertureEfficiency = gain2apertureEfficiency(...
    satTxConfig.AntennaGain,...
    satTxConfig.DishDiameter,...
    physconst('LightSpeed')/satTxConfig.Frequency); % ratio
satTxConfig.halfBeamwidth = 1.7647; % degrees
satTxConfig.satelliteBeamDiameter = 20e3; % meters
satTxConfig.BitRate = 10; % Mbps

satRxConfig = struct;
satRxConfig.Frequency = 3e10;
satRxConfig.DishDiameter = 0.33;
satRxConfig.RequiredEbNo = 0;
satRxConfig.SystemLoss = 3; % dB
satRxConfig.MaxGainToNoiseTemperatureRatio = 13;
satRxConfig.AntennaGain = 38.5;
satRxConfig.ApertureEfficiency = gain2apertureEfficiency(...,
    satRxConfig.AntennaGain,...
    satRxConfig.DishDiameter,...
    physconst('LightSpeed')/satRxConfig.Frequency);

%%  Scenario Section
%   Scenario Simulation Time
startTime = datetime(2020,1,11,14,0,0);
stopTime = startTime + days(1);
sampleTime = 10;

%   Scenario Generation
sc = satelliteScenario(startTime,stopTime,sampleTime);

%%  Satellites Section
%   Add satellites to the scenario
sat1 = satellite(sc,earthRadius+600e3,0,60,0,0,0, ...
    "Name","Sat1","OrbitPropagator","two-body-keplerian");
sat2 = satellite(sc,earthRadius+600e3,0,50,0,0,0, ...
    "Name","Sat2","OrbitPropagator","two-body-keplerian");

% Add the gimbals on the satellites
g1Sat1 = gimbal(sat1,"MountingLocation",[0;1;2]);
g2Sat1 = gimbal(sat1,"MountingLocation",[0;-1;2]);
g1Sat2 = gimbal(sat2,"MountingLocation",[0;1;2]);
g2Sat2 = gimbal(sat2,"MountingLocation",[0;-1;2]);

%   Add the transmitters of satellites
sat1Tx = transmitter(g1Sat1,...
    "MountingLocation",[0;0;1], ...
    Frequency=satTxConfig.Frequency,...
    Power=satTxConfig.TxPower,...
    SystemLoss=satTxConfig.SystemLoss);

sat2Tx = transmitter(g1Sat2,...
    "MountingLocation",[0;0;1], ...
    Frequency=satTxConfig.Frequency,...
    Power=satTxConfig.TxPower,...
    SystemLoss=satTxConfig.SystemLoss);

%   Use gaussian Antennas for both transmitters
gaussianAntenna(sat1Tx,...
    "DishDiameter",satTxConfig.DishDiameter,...
    "ApertureEfficiency",satTxConfig.ApertureEfficiency);

gaussianAntenna(sat2Tx,...
    "DishDiameter",satTxConfig.DishDiameter,...
    "ApertureEfficiency",satTxConfig.ApertureEfficiency);

%   Add receivers to the gimbals of the ground station 2
sat1Rx = receiver(g2Sat1,...
    "MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=satRxConfig.MaxGainToNoiseTemperatureRatio,...
    RequiredEbNo=satRxConfig.RequiredEbNo,...
    SystemLoss=satRxConfig.SystemLoss);

sat2Rx = receiver(g2Sat2,...
    "MountingLocation",[0;0;1], ...
    GainToNoiseTemperatureRatio=satRxConfig.MaxGainToNoiseTemperatureRatio,...
    RequiredEbNo=satRxConfig.RequiredEbNo,...
    SystemLoss=satRxConfig.SystemLoss);

gaussianAntenna(sat1Rx,...
    "DishDiameter",satRxConfig.DishDiameter,...
    "ApertureEfficiency",satRxConfig.ApertureEfficiency);

gaussianAntenna(sat2Rx,...
    "DishDiameter",satRxConfig.DishDiameter,...
    "ApertureEfficiency",satRxConfig.ApertureEfficiency);

%%  Ground Stations Section
%   Add the ground stations
name = "Lamia Base Station";
lat = 38.875941904578895;
lon = 22.437526584061686;
%   Both VSATs placed at 10 meters
alt = 0;
minElevationAngle = 0;
gs1 = groundStation(sc,Name=name,Latitude=lat,Longitude=lon,...
    Altitude=alt,...
    MinElevationAngle=minElevationAngle);

latitude = 52.2294963;                                             
longitude = 0.1487094;                                             
gs2 = groundStation(sc,latitude,longitude,Altitude=alt,...
    Name="Cambridge VSAT",...
    MinElevationAngle=minElevationAngle);

% Add two gimbals on each ground station
g1Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs1 = gimbal(gs1,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);
g1Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;1;-5]);
g2Gs2 = gimbal(gs2,"MountingAngles",[0;180;0],"MountingLocation",[0;-1;-5]);

%   Add transmitters on ground station 1 (GS1)
%   Add receivers on ground station 2 (GS2)
gs1Tx1 = transmitter(g1Gs1,...
    "Name","GS1Tx1",...
    "MountingLocation",[0;0;1], ...
    "Frequency",terminalTxConfig.Frequency,...
    "Power",terminalTxConfig.TxPower,...
    "SystemLoss",terminalTxConfig.SystemLoss,...
    "BitRate",terminalTxConfig.BitRate);

gs1Tx2 = transmitter(g2Gs1,"Name","GS1Tx2",...
    "MountingLocation",[0;0;1], ...
    "Frequency",terminalTxConfig.Frequency,...
    "Power",terminalTxConfig.TxPower,...
    "BitRate",terminalTxConfig.BitRate,...
    "SystemLoss",terminalTxConfig.SystemLoss);

gaussianAntenna(gs1Tx1,...
    "DishDiameter",terminalTxConfig.DishDiameter,...
    "ApertureEfficiency",terminalTxConfig.ApertureEfficiency);

gaussianAntenna(gs1Tx2,...
    "DishDiameter",terminalTxConfig.DishDiameter,...
    "ApertureEfficiency",terminalTxConfig.ApertureEfficiency);

%   Ground Station Receivers
gs2Rx1 = receiver(g1Gs2,...
    "Name","GS2Rx1",...
    "MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",terminalRxConfig.MaxGainToNoiseTemperatureRatio,...
    "RequiredEbNo",terminalRxConfig.RequiredEbNo,...
    "SystemLoss",terminalRxConfig.SystemLoss);

gs2Rx2 = receiver(g2Gs2,...
    "Name","GS2Rx2",...
    "MountingLocation",[0;0;1], ...
    "GainToNoiseTemperatureRatio",terminalRxConfig.MaxGainToNoiseTemperatureRatio,...
    "RequiredEbNo",terminalRxConfig.RequiredEbNo,...
    "SystemLoss",terminalRxConfig.SystemLoss);

gaussianAntenna(gs2Rx1,...
    "DishDiameter",terminalRxConfig.DishDiameter,...
    "ApertureEfficiency",terminalRxConfig.ApertureEfficiency);
gaussianAntenna(gs2Rx2,"DishDiameter",terminalRxConfig.DishDiameter,...
    "ApertureEfficiency",terminalRxConfig.ApertureEfficiency);

%   Set the tracking of the targets from the gimbals
% GS1Tx1 points at Sat1Rx
% GS1Tx2 points at Sat2Rx
% Sat1Tx points at GS2Rx1
% Sat2Tx points at GS2Rx2
pointAt(g1Gs1,sat1);
pointAt(g2Gs1,sat2);
pointAt(g2Sat1,gs1);
pointAt(g2Sat2,gs1);
pointAt(g1Gs2,sat1);
pointAt(g2Gs2,sat2);
pointAt(g1Sat1,gs2);
pointAt(g1Sat2,gs2);

%%  The Access and Link Objects can be created
%   Access objects refer to the visibility of the satellites
access1 = access(gs1,sat1);
access2 = access(gs1,sat2);
access3 = access(gs2,sat1);
access4 = access(gs2,sat2);
accessInterval1 = accessIntervals(access1);
accessInterval2 = accessIntervals(access2);
accessInterval3 = accessIntervals(access3);
accessInterval4 = accessIntervals(access4);

%   Link objects model the link between two transponders
%   Create all the link objects
%   Part1: Links between Ground Station Txs and Satellites (UL)
link1 = link(gs1Tx1,sat1Rx);
link2 = link(gs1Tx1,sat2Rx); % no link
link3 = link(gs1Tx2,sat1Rx); % no link
link4 = link(gs1Tx2,sat2Rx);
%   Part2: Links between Satellites and Ground Station Receivers (DL)
link5 = link(sat1Tx,gs2Rx1);
link6 = link(sat1Tx,gs2Rx2); % no link
link7 = link(sat2Tx,gs2Rx1); % no link
link8 = link(sat2Tx,gs2Rx2);

%   overall path
path1 = link(gs1Tx1,sat1Rx,sat1Tx,gs2Rx1);
path2 = link(gs1Tx2,sat2Rx,sat2Tx,gs2Rx2);

%   Downlink Eb/No calculation
[ebnr1,time1] = ebno(link5);
[ebnr2,time2] = ebno(link8);
ebnr1(ebnr1 == -Inf | ebnr2 == -Inf) = -Inf;
ebnr2(ebnr1 == -Inf | ebnr2 == -Inf) = -Inf;
%   Path Eb/No calculation
[path1_ebnr,path_time_1] = ebno(path1);
[path2_ebnr,path_time_2] = ebno(path2);
path1_ebnr(path1_ebnr == -Inf | path2_ebnr == -Inf) = -Inf;
path2_ebnr(path1_ebnr == -Inf | path2_ebnr == -Inf) = -Inf;
%   Display EbNo from both DL links
figure;
plot(time1,ebnr1,'-','LineWidth',1.5);
hold on;
plot(time2,ebnr2,'-','LineWidth',1.5);
hold off;
grid on;
title("Received E_b/N_o from both satellites");
legend('Path 1','Path 2');
ylabel('E_b/N_o [dB]');
xlabel('Simulation Time (datetime)');
axis tight;

%   Display EbNo from both paths
figure;
plot(path_time_1,path1_ebnr,'-','LineWidth',1.5);
hold on;
plot(path_time_2,path2_ebnr,'-','LineWidth',1.5);
hold off;
grid on;
title('E_b/N_o of both paths');
legend('Path 1', 'Path 2');
ylabel('E_b/N_o [dB]');
xlabel('Simulation Time (datetime)');
axis tight;

% [ebnr1_ul,time1_ul] = ebno(link1);
%   Plot multi-hop ebno to dl ebno
figure;
plot(path_time_1,path1_ebnr,'LineWidth',1.5);
hold on;
plot(time1,ebnr1,'-','LineWidth',1.5);
% plot(time1_ul,ebnr1_ul,'-.','LineWidth',1.5);
hold off;
grid on;
axis tight;
legend('Path EbNo','DL EbNo');
title('E_b/N_0 for path 1');
xlabel('Simulation Time');
ylabel('E_b/N_0 [dB]');

%%  Get geometrical data for p618 losses
ebnr1_lossy = zeros(1,length(ebnr1));
atmo_losses1 = zeros(1,length(ebnr1));
atmo_losses2 = zeros(1,length(ebnr2));

for n=1:length(time1)
    % Calculate geometry in respect of the ground station
    [~,el,r] = aer(gs2Rx1,sat1Tx,time1(n));
    if el < 5 || el > 175
        ebnr1_lossy(n) = -Inf;
        atmo_losses1(n) = -Inf;
        continue;
    end
    p618cfg = p618Config(AntennaDiameter=satTxConfig.DishDiameter,...
        AntennaEfficiency=satTxConfig.ApertureEfficiency,...
        Frequency=satTxConfig.Frequency,...
        Latitude=gs2.Latitude,...
        Longitude=gs2.Longitude,...
        ElevationAngle=el);
    p618loss = p618PropagationLosses(p618cfg);
    atmoLoss = p618loss.Ac + p618loss.Ag + p618loss.Ar + ...
        p618loss.As + p618loss.At;
    atmo_losses1(n) = atmoLoss;
    ebnr1_lossy(n) = ebnr1(n) - atmoLoss;
end

ebnr2_lossy = zeros(1,length(ebnr2));
for n=1:length(time2)
    % Calculate geometry in respect of the ground station
    [~,el,r] = aer(gs2Rx2,sat2Tx,time2(n));
    if el < 5 || el > 175
        ebnr2_lossy(n) = -Inf;
        atmo_losses2(n) = -Inf;
        continue;
    end
    p618cfg = p618Config(AntennaDiameter=satTxConfig.DishDiameter,...
        AntennaEfficiency=satTxConfig.ApertureEfficiency,...
        Frequency=satTxConfig.Frequency,...
        Latitude=gs2.Latitude,...
        Longitude=gs2.Longitude,...
        ElevationAngle=el);
    p618loss = p618PropagationLosses(p618cfg);
    atmoLoss = p618loss.Ac + p618loss.Ag + p618loss.Ar + ...
        p618loss.As + p618loss.At;
    atmo_losses2(n) = atmoLoss;
    ebnr2_lossy(n) = ebnr2(n) - atmoLoss;
end

%   Display Received EbNo with atmosphairic losses
figure;
plot(time1,ebnr1_lossy,"LineWidth",1.5);
hold on;
plot(time2,ebnr2_lossy,"LineWidth",1.5);
grid on;
axis tight;
title("Received E_b/N_o from both satellites + p618 losses");
ylabel('E_b/N_o [dB]');
xlabel('Simulation Time (datetime)');

%   Display path EbNo with atmosphairic losses
figure;
plot(path_time_1,path1_ebnr-atmo_losses1,"LineWidth",1.5);
hold on;
plot(path_time_2,path2_ebnr-atmo_losses2,"LineWidth",1.5);
grid on;
axis tight;
title("Received E_b/N_o from both satellites + p618 losses");
ylabel('E_b/N_o [dB]');
xlabel('Simulation Time (datetime)');

%   Display the Atmosphairic Losses
figure;
plot(time1,atmo_losses1,time2,atmo_losses2);
legend('branch 1', 'branch 2');
grid on;
axis tight;
title('Atmosphairic Losses (ITU618)');
xlabel('Simulation time (datetime)');
ylabel('Losses [dB]');

%%  Combining at the Receiver
%   Maximal Ratio Combining
%   Selection Combining
ebnr_mrc = 10*log10(10.^(ebnr1_lossy./10) + 10.^(ebnr2_lossy./10));
ebnr_sc = 10*log10(max(10.^(ebnr1_lossy./10),10.^(ebnr2_lossy./10)));
%   Display the Combined EbNo with both combining techniques
figure;
plot(time1,ebnr_mrc,'-','LineWidth',1.5); hold on;
plot(time1,ebnr_sc,'-','LineWidth',1.5); hold off;
grid on;
legend('MRC','SC');
title('Selection Combining VS Maximal Ratio Combining');
ylabel('EbNo (dB)');
xlabel('Simulated time (datetime)');
axis tight;
%   Display the Margin between the two techniques
figure;
plot(time1,abs(ebnr_mrc-ebnr_sc),'LineWidth',1.5);
grid on;
title('Selection Combining VS Maximal Ratio Combining');
ylabel('EbNo (dB)');
xlabel('Simulated time (datetime)');
axis tight;

%%  Combining at the Receiver for full path
%   Maximal Ratio Combining
%   Selection Combining
path1_ebnr_lossy = path1_ebnr - atmo_losses1;
path2_ebnr_lossy = path2_ebnr - atmo_losses2;
path_ebnr_mrc = 10*log10(10.^(path1_ebnr_lossy./10) + 10.^(path2_ebnr_lossy./10));
path_ebnr_sc = 10*log10(max(10.^(path1_ebnr_lossy./10),10.^(path2_ebnr_lossy./10)));
%   Display the Combined EbNo with both combining techniques
figure;
plot(path_time_1,path_ebnr_mrc,'-','LineWidth',1.5); hold on;
plot(path_time_1,path_ebnr_sc,'-','LineWidth',1.5); hold off;
grid on;
legend('MRC','SC');
title('Selection Combining VS Maximal Ratio Combining');
ylabel('EbNo (dB)');
xlabel('Simulated time (datetime)');
axis tight;

%   Display the Margin between the two techniques
figure;
plot(time1,abs(path_ebnr_mrc-path_ebnr_sc),'LineWidth',1.5);
grid on;
title('Selection Combining VS Maximal Ratio Combining');
ylabel('EbNo (dB)');
xlabel('Simulated time (datetime)');
axis tight;

%%  BER Calculation
%   Using 8-QAM, 16-QAM and QPSK
symbolRate = terminalRxConfig.BitRate/6;
%   Set the sample rate to at least 4-times the symbol rate
sampleRateChannel = max(1e6,4*symbolRate);
figure;
%   8-QAM
M = 8;
k = log2(M);
K = 10;
ebnr_mrc_clean = path_ebnr_mrc(~isnan(path_ebnr_mrc) & path_ebnr_mrc ~= -Inf);
hLOS = exp(1i*2*pi*rand(1,length(ebnr_mrc_clean))); % Thank you Emil Bjornson!
hNLOS = sqrt(1/2)*(randn(1, ...
    length(ebnr_mrc_clean)) + 1i*randn(1, length(ebnr_mrc_clean)));
h = sqrt(K/(K + 1))*hLOS + sqrt(1/(K + 1))*hNLOS;
test_ebnr = 10*log10((abs(h).^2).*(10.^((ebnr_mrc_clean)./10)));
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




%%  Helper Functions
%   turn antenna gain into the aperture efficiency
function ea = gain2apertureEfficiency(gain,dishDiameter,lambda)
    %   gain = (pi*diameter/lambda)^2*ea
    gain_linear = 10^(gain/10);
    ea = gain_linear/(pi*dishDiameter/lambda)^2;
end