%% Editable study inputs

config = struct();
config.sBand_Hz = [2e9 4e9];
config.xBand_Hz = [8e9 12e9];
config.detectorWavelengths_m = [0.76e-6 2.3e-6];
config.writeOutputs = true;
results = runSpectrumStudy(config);
