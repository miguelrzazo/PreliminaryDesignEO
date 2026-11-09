%% Editable study inputs
% Historical rate assumptions: ITU-R P.838-3 at 8.1 GHz, k=0.00425, alpha=1.32,
% with 60-degree elevation and 2 km rain height in the original study notes.
% Seasonal probabilities were attributed to NOAA Climate Normals 1991-2020,
% Fairbanks Intl Airport (USW00026411). Preserve these inputs without recalibration.
% Seasonal rows: spring, summer, autumn, winter; condition columns follow names.
config = struct();
config.dayCount = 365;
config.randomSeed = 42;
config.conditionNames = {'Cond. favorables','Lluvia moderada','Lluvia intensa','Nubes densas'};
config.rateRanges_Mbps = [130 150;120 140;90 120;125 145];
config.seasonProbabilities_percent = [76 12 3 9;52 36 7 5;64 22 5 9;82 1 0 17];
config.outputDirectory = fullfile(fileparts(mfilename('fullpath')),'..','Latex_Code','7.Segmento_Tierra');
config.writeOutputs = true;
results = runWeatherLinkStudy(config);
