%% Editable study inputs

config = struct();
config.launchLatitude_deg = 34.7;
config.launchLongitude_deg = -120.6;
config.minimumAzimuth_deg = 170;
config.maximumAzimuth_deg = 300;
config.boundaryRadius_deg = 8;
config.arcRadius_deg = 5.5;
config.arcSampleCount = 80;
config.launchAzimuth_deg = 190;
config.distances_km = linspace(0,900,90);
config.outputDirectory = fullfile(fileparts(mfilename('fullpath')),'..','Latex_Code','6.Lanzadores');
config.writeOutputs = true;
results = runAzimuthStudy(config);
