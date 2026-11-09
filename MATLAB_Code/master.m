%% Earth observation study — editable inputs
% Edit this file to change the study. runAnalysis coordinates the calculations.
% Wavelength roles remain CO2 band 1, CO2 band 2, and O2 A-band, respectively.
config = struct();

%% Execution and exports
config.writeOutputs = true;
config.outputRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Compiled data','Pipeline');

%% Mission requirements and design axes
config.groundSampleDistance_m = 80;
config.wavelengths_m = [1.61e-6, 2.01e-6, 0.76e-6];
config.altitudes_km = 200:10:1000;
config.pupilDiameters_mm = 20:1:200;
config.swaths_km = 20:10:500;
config.snrRequirement = 400;
config.mtfRequirement = 0.25;
config.revisitRequirement_days = 7;
config.maxDetectorCount = 3;
config.overlapFraction = 0.05;
config.cloudFraction = 1/6;
config.localTimeAscendingNode_h = 6;
config.maintenanceLowerFactor = 0.98;
config.maintenanceUpperFactor = 1.02;

% Each row is [satellites, telescopes per satellite].
config.satelliteConfigurations = [1 1; 2 1; 2 2; 1 2; 3 1];

%% CO2 detectors — shared by the first two wavelength bands
% Property vectors follow the detector names in the same order.
config.co2Detectors.names = {'CAPYORK', 'H2RG', 'SATURN VISIR'};
config.co2Detectors.quantumEfficiency = [0.85, 0.70, 0.60];
config.co2Detectors.detectorMtf = [0.45, 0.50, 0.48];
config.co2Detectors.pixelCount = [1200, 2048, 1000];
config.co2Detectors.pixelSize_m = [15, 18, 30] * 1e-6;

%% O2 detectors — third wavelength band
config.o2Detectors.names = {'CCD98-50', 'H2RG', 'SATURN VISIR'};
config.o2Detectors.quantumEfficiency = [0.65, 0.70, 0.60];
config.o2Detectors.detectorMtf = [0.50, 0.50, 0.45];
config.o2Detectors.pixelCount = [6000, 2048, 1000];
config.o2Detectors.pixelSize_m = [13, 18, 30] * 1e-6;

%% Telescopes
% The existing selection model treats the first entry as refractive (200 mm cap).
config.telescopeNames = {'Refractivo', 'Korsch', 'Cassegrain', 'TMA'};
config.telescopeAlignmentMtf = [0.90, 0.85, 0.80, 0.70];
config.telescopeTransmittance = [0.80, 0.65, 0.70, 0.60];
config.telescopeFov_deg = [10, 3, 3, 8];
config.obstructionRatios = [0, 0.2, 0.2, 0];

%% Optional command-line overrides
% Unset these environment variables when using only the inputs above.
gsdOverride = getenv('TFG_GSD');
if ~isempty(gsdOverride)
    config.groundSampleDistance_m = str2double(gsdOverride);
    fprintf('TFG_GSD override: %g m\n', config.groundSampleDistance_m);
end
outputOverride = getenv('TFG_OUTPUT_ROOT');
if ~isempty(outputOverride), config.outputRoot = outputOverride; end

%% Run the analytical study
% An unsuccessful run must not leave previous results looking current.
results = [];
results = runAnalysis(config);
