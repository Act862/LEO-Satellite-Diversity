%   uplink
cfg = satelliteCNRConfig;
cfg.TransmitterPower = 3;
cfg.TransmitterSystemLoss= 0;
cfg.TransmitterAntennaGain = 43.2;
cfg.Distance = 1230;
cfg.Frequency = 30;
cfg.Bandwidth = 16;

cfg.MiscellaneousLoss = 3;
cfg.GainToNoiseTemperatureRatio = 18.5;
cfg.ReceiverSystemLoss = 0;
cfg.BitRate = 10;
disp(cfg)
[cn, info] = satelliteCNR(cfg)