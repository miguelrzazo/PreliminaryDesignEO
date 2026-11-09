function results = runAnalysis(config, revisitSolver)
%RUNANALYSIS Execute the analytical study using explicit, validated inputs.
% master.m owns the editable inputs. This function does not read environment
% overrides, clear the caller workspace, or close the caller's figures.
% The optional solver handle supports deterministic coverage integration tests.
analysisDirectory = fileparts(mfilename('fullpath'));
plotDefaultsCleanup = configureThesisPlotDefaults();
previousPath = path;
pathCleanup = onCleanup(@() path(previousPath));
addpath(fullfile(analysisDirectory,'lib'));
validateStudyConfig(config);
if nargin < 2, revisitSolver = @RevisitCalc; end
validateattributes(revisitSolver,{'function_handle'},{'scalar'});

% Relative paths are confined to compatibility exporters. Restore the caller's
% directory on return or error; compute-only runs never change directory.
if config.writeOutputs && ~isempty(config.outputRoot)
    if ~exist(config.outputRoot,'dir'), mkdir(config.outputRoot); end
    previousDirectory = pwd;
    directoryCleanup = onCleanup(@() cd(previousDirectory));
    cd(config.outputRoot);
end
if config.writeOutputs
    outputDirectories = {'MTF','SNR','coverage','HvsDmin','Masa_seca','Masa_total','OptimumConfigs'};
    for directoryIndex = 1:numel(outputDirectories)
        if ~exist(outputDirectories{directoryIndex},'dir'), mkdir(outputDirectories{directoryIndex}); end
    end
end

% Optical cells: {band, detector, telescope}; coverage cells: {configuration,
% detector, telescope}. Each cell holds one grid with explicit physical axes.
results.config = config;
results.executionDirectory = pwd;
results.revisitSolver = func2str(revisitSolver);
results.axes.wavelengths_m = config.wavelengths_m(:).';
results.axes.altitudes_km = config.altitudes_km(:);
results.axes.pupilDiameters_mm = config.pupilDiameters_mm(:).';
results.axes.swaths_km = config.swaths_km(:).';
results.axes.telescopeNames = config.telescopeNames;
results.axes.satelliteConfigurations = config.satelliteConfigurations;
telescopeCount = numel(config.telescopeNames);
configurationCount = size(config.satelliteConfigurations,1);
results.optics.mtf = cell(3,3,telescopeCount);
results.optics.snr = cell(3,3,telescopeCount);
results.coverage.days = cell(configurationCount,3,telescopeCount);
results.coverage.diagnostics = cell(size(results.coverage.days));
results.coverage.stageErrors = cell(size(results.coverage.days));

% Flatten hardware combinations; each iteration owns one altitude/diameter grid.
opticalCases = analysisCases(numel(config.wavelengths_m),3,telescopeCount);
for caseIndex = 1:size(opticalCases,1)
    bandIndex = opticalCases(caseIndex,1);
    detectorIndex = opticalCases(caseIndex,2);
    telescopeIndex = opticalCases(caseIndex,3);
    if bandIndex <= 2
        detectors = config.co2Detectors;
        detectorExportOffset = 0;
    else
        detectors = config.o2Detectors;
        detectorExportOffset = 3;
    end
    wavelength_m = config.wavelengths_m(bandIndex);
    mtfPrefix = sprintf('MTF/MTF_Lambda%d_Detector%d_Telescopio%d', ...
        bandIndex,detectorIndex+detectorExportOffset,telescopeIndex);
    snrPrefix = sprintf('SNR/SNR_Lambda%d_Detector%d_Telescopio%d', ...
        bandIndex,detectorIndex+detectorExportOffset,telescopeIndex);
    results.optics.mtf{bandIndex,detectorIndex,telescopeIndex} = MTFfunction( ...
        wavelength_m,detectors.pixelSize_m(detectorIndex), ...
        detectors.detectorMtf(detectorIndex),config.telescopeAlignmentMtf(telescopeIndex), ...
        config.groundSampleDistance_m,config.obstructionRatios(telescopeIndex), ...
        config.altitudes_km,config.pupilDiameters_mm,mtfPrefix, ...
        config.telescopeNames{telescopeIndex},detectorIndex, ...
        config.mtfRequirement,config.writeOutputs);
    results.optics.snr{bandIndex,detectorIndex,telescopeIndex} = SNRfunction( ...
        wavelength_m,detectors.pixelSize_m(detectorIndex), ...
        detectors.quantumEfficiency(detectorIndex),config.telescopeTransmittance(telescopeIndex), ...
        config.groundSampleDistance_m,config.obstructionRatios(telescopeIndex), ...
        config.altitudes_km,config.pupilDiameters_mm,snrPrefix, ...
        config.telescopeNames{telescopeIndex},detectorIndex, ...
        config.snrRequirement,config.writeOutputs);
end

coverageCache = containers.Map('KeyType','char','ValueType','double');
coverageSolver = @(varargin) cachedRevisit(coverageCache,revisitSolver,varargin{:});
coverageCases = analysisCases(configurationCount,3,telescopeCount);
for caseIndex = 1:size(coverageCases,1)
    configurationIndex = coverageCases(caseIndex,1);
    detectorIndex = coverageCases(caseIndex,2);
    telescopeIndex = coverageCases(caseIndex,3);
    satelliteCount = config.satelliteConfigurations(configurationIndex,1);
    telescopesPerSatellite = config.satelliteConfigurations(configurationIndex,2);
    [coverageGrid_days,diagnostics] = computeCoverageGrid( ...
        config.groundSampleDistance_m,config.altitudes_km,config.swaths_km, ...
        config.revisitRequirement_days,config.co2Detectors.pixelCount(detectorIndex), ...
        config.telescopeFov_deg(telescopeIndex),satelliteCount,telescopesPerSatellite, ...
        config.overlapFraction,config.cloudFraction,config.maxDetectorCount,coverageSolver);
    results.coverage.days{configurationIndex,detectorIndex,telescopeIndex} = coverageGrid_days;
    results.coverage.diagnostics{configurationIndex,detectorIndex,telescopeIndex} = diagnostics;
    if config.writeOutputs
        configurationName = sprintf('%d Satelite(s); %d Telescopio(s): %s; Detector %d', ...
            satelliteCount,telescopesPerSatellite, ...
            strrep(config.telescopeNames{telescopeIndex},' ',''),detectorIndex);
        writeCoverageResult(coverageGrid_days,config.altitudes_km,config.swaths_km, ...
            config.revisitRequirement_days,diagnostics,configurationName);
    end
end
coverageSolverFailed = cellfun(@(entry) ...
    any(~cellfun(@isempty,entry.solverErrors(:))) || any(entry.nonfiniteSolverResult(:)), ...
    results.coverage.diagnostics(:));
results.coverage.executionSucceeded = ~any(coverageSolverFailed);

results.selection = CrossDataFunction(config.altitudes_km,config.swaths_km, ...
    config.pupilDiameters_mm,config.telescopeNames,config.telescopeFov_deg, ...
    config.satelliteConfigurations,config.mtfRequirement,config.snrRequirement, ...
    config.revisitRequirement_days,results,config.writeOutputs);
results.configuration = ConfigAnalysis(legacyAnalysisParameters(config,coverageSolver),results,config.writeOutputs);
results.coverage.uniqueSolverCalls = coverageCache.Count;
if ~results.coverage.executionSucceeded
    warning('eo:IncompleteAnalysis', ...
        'Coverage solver failures occurred. Inspect results.coverage.diagnostics.');
end
end

function params = legacyAnalysisParameters(config,revisitSolver)
% Keep the older mass/selection interface behind one explicit naming adapter.
params.satellite_configs = config.satelliteConfigurations;
params.telescope_names = config.telescopeNames;
params.alturas_orbitales = config.altitudes_km;
params.swaths_km = config.swaths_km;
params.diametros_pupila = config.pupilDiameters_mm;
params.Cov_Requirement = config.revisitRequirement_days;
params.MTF_Requirement = config.mtfRequirement;
params.SNR_Requirement = config.snrRequirement;
params.GSD = config.groundSampleDistance_m;
params.Npix_detectors = config.co2Detectors.pixelCount;
params.fov_limit = config.telescopeFov_deg;
params.LTAN_hour = config.localTimeAscendingNode_h;
params.Pixel_size = config.co2Detectors.pixelSize_m;
params.maintenance_lower_factor = config.maintenanceLowerFactor;
params.maintenance_upper_factor = config.maintenanceUpperFactor;
params.solapamiento = config.overlapFraction;
params.cobertura_nubes = config.cloudFraction;
params.regionalRevisitSolver = revisitSolver;
end

function revisit = cachedRevisit(cache,solver,varargin)
% Revisit depends on the orbit and viewing geometry, not detector identity.
key = sprintf('%.17g_', [varargin{:}]);
if isKey(cache,key)
    revisit = cache(key);
else
    revisit = solver(varargin{:});
    cache(key) = revisit;
end
end
