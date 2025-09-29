channel = nrTDLChannel;
channel.DelayProfile="TDL-A";
channel.DelaySpread = 30e-9;
channel.TransmissionDirection = "Uplink";
channel.MIMOCorrelation = "Low";
channel.Polarization = "Co-Polar";
channel.SampleRate = 7680000;
channel.MaximumDopplerShift = 0;
channel.SatelliteDopplerShift = dopplerShiftCircularOrbit(...
    50,600000,...
    0,2e9); 
channel.RandomStream = "mt19937ar with seed";
channel.Seed = 73;

if ~strcmpi(channel.MIMOCorrelation, 'Custom')
    channel.NumTransmitAntennas = 1;
    channel.NumReceiveAntennas = 2;
else
    if any(strcmpi(channel.Polarization, {'Co-Polar','Cross-Polar'}))
        channel.TransmitCorrelationMatrix = 1;
        channel.ReceiveCorrelationMatrix = [1 0;0 1];
    end
    if strcmpi(channel.Polarization,'Cross-Polar')
        channel.TransmitPolarizationAngles = [45 -45];
        channel.ReceivePolarizationAngles = [90 0];
        channel.XPR = 10;
    end
    if strcmpi(channel.Polarization,'Custom')
        channel.SpatialCorrelationMatrix = [1 0; 0 1];
    end
end

channelInfo = info(channel);
rng(73);
in = randn(channel.SampleRate,channelInfo.NumTransmitAntennas,'like',1i);
[tdlOut,tdlPathGains,tdlSampleTimes] = channel(in);

ntnTDLAnalyzer = spectrumAnalyzer(SampleRate = channel.SampleRate);
ntnTDLAnalyzer.Title = "Received Signal Spectrum " ...
    + "NTN TDL with " + string(channel.DelayProfile) + " delay profile";
ntnTDLAnalyzer.ShowLegend = true;
for nRx = 1:size(tdlOut,2)
    ntnTDLAnalyzer.ChannelNames{nRx} = "Rx Antenna " + nRx;
end
ntnTDLAnalyzer(tdlOut)