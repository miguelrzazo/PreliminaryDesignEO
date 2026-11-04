%% Editable phasing-manoeuvre inputs
config = struct();
config.earthGravitationalParameter_m3_s2 = 3.986004418e14;
config.earthRadius_m = 6378.137e3;
config.finalAltitude_m = 630e3;
config.phaseApogeeAltitude_m = 830e3;
config.timeStep_s = 120;
config.outputDirectory = fullfile(fileparts(mfilename('fullpath')),'..','Latex_Code','6.Lanzadores');
config.writeOutputs = true;
results = runPhasingStudy(config);
