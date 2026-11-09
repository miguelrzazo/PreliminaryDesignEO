function selection = CrossDataFunction(altitudes_km, swaths_km, pupilDiameters_mm, ...
    telescopeNames, telescopeFov_deg, satelliteConfigurations, mtfRequirement, ...
    snrRequirement, revisitRequirement_days, analysisResults, writeOutputs)
%CROSSDATAFUNCTION Compatibility adapter for selection, import and presentation.
if nargin < 10, analysisResults = []; end
if nargin < 11, writeOutputs = true; end
if isempty(analysisResults)
    analysisResults = importAnalysisResults(altitudes_km,pupilDiameters_mm,swaths_km,telescopeNames,satelliteConfigurations);
else
    validateAnalysisAxes(analysisResults,altitudes_km,pupilDiameters_mm,swaths_km,telescopeNames,satelliteConfigurations);
end
selection = computeSelectionResults(analysisResults,telescopeFov_deg,mtfRequirement,snrRequirement,revisitRequirement_days);
if writeOutputs
    requirements = struct('mtfRequirement',mtfRequirement,'snrRequirement',snrRequirement, ...
        'revisitRequirement_days',revisitRequirement_days);
    writeSelectionResults(selection,altitudes_km,telescopeNames,satelliteConfigurations,requirements);
end
end
