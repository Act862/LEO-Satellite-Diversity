rng('default');
%%  Initial Setup
%   Radiofrequency parameters
fc = 38e9;  % carrier frequency (Hz) : 38 GHz
c = physconst('LightSpeed'); % speed of light (m/s)
bw = 100e6; % bandwidth (Hz)
sampleRate = bw; % sample rate (Hz)

%   Tx and Rx parameters
Pt = 0.1; % Peak Power (W)
Gtx = 20; % Antenna Gain dBi
Grx = 20; % Antenna Gain dBi
Nf = 2.9; % Noise figure (dB)

%   Transmitter
antenna = phased.IsotropicAntennaElement('BackBaffled',false);
transmitter = phased.Transmitter("Gain",Gtx, "PeakPower",Pt);
radiator = phased.Radiator('Sensor',antenna, 'OperatingFrequency',fc);
collector = phased.Collector('Sensor',antenna,'OperatingFrequency',fc);

%   Receiver
receiver = phased.Receiver('AddInputNoise',true,'Gain',Grx,...
    'NoiseFigure',Nf, 'SampleRate',sampleRate);

%   Create a one-way free-space propagation channel
channel = phased.FreeSpace('PropagationSpeed',c, 'OperatingFrequency',...
    fc, 'SampleRate',sampleRate,'TwoWayPropagation',false);

%   Device platform
tgtpos = [8;4;0];   % device position (m)
tgtvel = [1; -1; 0.5]; % device velocity (m/s)
tgtplatform = phased.Platform('InitialPosition',tgtpos, 'Velocity',tgtvel);

%   Anchor platform
anchorpos = [0 30 10 20 15;...
    10 -5 -20 15 10;...
    4 -10 10 5 -20];
numAnchor = size(anchorpos,2);
anchorvel = zeros(3, numAnchor);
anchorplatform = phased.Platform('InitialPosition',...
    anchorpos,'Velocity',anchorvel);

%%  Waveform generation
%   Generate the waveform
N = 1024; % number of subcarriers
M = 8; % Number of OFDM symbols (channel samples)
freqSpacing = bw/N; % Frequency spacing (Hz)
tsym = 1/freqSpacing; % Symbol Duration
maxDelay = 200e-9; % maximum delay
rmax = maxDelay*c;
tcp = range2time(rmax);
Ncp = ceil(sampleRate*tcp);
tcp = Ncp/sampleRate;
tWave = tsym + tcp;
Ns = N + Ncp;

%% Channel Estimation

X = cell(1, numAnchor);

for idxAnchor = 1:numAnchor
    bpskSymbol = randi([0,1], [N M])*2-1;
    %   Generate OFDM modulated signal with cyclic prefix
    sigmod = ofdmmod(bpskSymbol, N, Ncp);
    %   Reshape the OFDM modulated signal
    %   Ns: Total waveform samples
    %   M: Number of OFDM symbols per subcarrier
    sig = reshape(sigmod, Ns, M);
    
    %   Power normalization OFDM signal
    sig = sig/max(abs(sig),[],'all');
    %   Initialize estimated channel
    x = complex(zeros(size(sig)));

    %   For each OFDM symbol
    for m=1:M
        [tx_pos, tx_vel] = anchorplatform(tWave);
        [rx_pos, rx_vel] = tgtplatform(tWave);

        %   Calculate the transmit angle
        [~, txang] = rangeangle(rx_pos, tx_pos(:,idxAnchor));

        %   Form Transmitted Signal
        txSig = transmitter(sig);

        %   Radiate the transmitted signal
        radtxsig = radiator(txSig(:,m),txang);

        %   Propagate the signal
        chansig = channel(radtxsig, tx_pos(:,idxAnchor), rx_pos,...
            tx_vel(:,idxAnchor),rx_vel);

        %   Calculate the receive angle
        [~,rxang] = rangeangle(tx_pos(:,idxAnchor),rx_pos);

        %   Collect signal at the receive antenna
        rxSig = collector(chansig, rxang);

        %   Receive signal at the receiver
        x(:,m) = receiver(rxSig);
    end

    xvec = reshape(x, Ns*M, 1);
    xdemod = ofdmdemod(xvec, N, Ncp, Ncp);

    X{idxAnchor} = xdemod./bpskSymbol;

    reset(anchorplatform);
    reset(tgtplatform);
end

%%  TOA

spectrumMethod = 'FFT';

toaEstimator = phased.TOAEstimator('PropagationSpeed',c,...
    'Measurement','TOA','SpectrumMethod',spectrumMethod,...
    'VarianceOutputPort',true,'DelayOffsetInputPort',true);

delayoffset = 0;
[Y, variance] = toaEstimator(X,freqSpacing,delayoffset);

figure
plotTOASpectrum(toaEstimator,freqSpacing,'AnchorIndex',1,'MaxDelay',maxDelay);

% Obtain TOA position estimate
tgtposest = toaposest(Y,variance,anchorpos);

% View TOA position estimate
helperPlotTOADevicePositions(c,Y,tgtposest,anchorpos,tgtpos);

% RMSE of the TOA position estimate
RMSE = rmse(tgtposest,tgtpos);
disp(['RMS Localization error = ', num2str(RMSE), ' meters.'])

function helperPlotTOADevicePositions(c,toaEst,tgtPosEst,anchorPos,tgtPos)
%helperPlotTOADevicePositions Plots 2D positions of all anchors and
%targets, estimated positions of targets, and trilateration circles

% Range estimation
rngEst = toaEst*c;

% Plot positions of all anchors, target, and estimated position of the target
figure
plot(anchorPos(1,:),anchorPos(2,:),'b^','LineWidth',2,'MarkerSize',10),hold on
plot(tgtPos(1,:),tgtPos(2,:),'rx','LineWidth',2,'MarkerSize',10),hold on
plot(tgtPosEst(1),tgtPosEst(2),'go','LineWidth',2,'MarkerSize',12),hold on

% Get distance from anchors to target to keep consistent plot axes
radiustx = sqrt(anchorPos(1,:).^2 + anchorPos(2,:).^2);
range = max(abs(anchorPos+repmat(radiustx,[size(anchorPos,1) 1])),[],'all');
xlim([-range range]),ylim([-range range]);
axis equal
grid on;
xlabel('x-axis (meters)'),ylabel('y-axis (meters)');

% Plot trilateration circles
numAnchor = length(rngEst);
angles = 0:2*pi/720:2*pi;
for anchorIdx = 1:numAnchor
    x = rngEst(anchorIdx) * cos(angles) + anchorPos(1,anchorIdx);
    y = rngEst(anchorIdx) * sin(angles) + anchorPos(2,anchorIdx);
    plot(x,y,'c--','LineWidth',1);
    hold on; grid on
end
legend({'Anchor Positions','Target Position','TOA Position Estimate','Trilateration Circles'},'Location','best','FontSize',10)
end