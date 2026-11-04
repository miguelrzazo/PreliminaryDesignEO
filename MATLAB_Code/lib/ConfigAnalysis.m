function summary = ConfigAnalysis(params, analysisResults, writeOutputs)
%CONFIGANALYSIS Compatibility adapter for input, mass calculations and reports.
if nargin < 2, analysisResults = []; end
if nargin < 3, writeOutputs = true; end
if isempty(analysisResults)
    analysisResults = importAnalysisResults(params.alturas_orbitales,params.diametros_pupila, ...
        params.swaths_km,params.telescope_names,params.satellite_configs);
    analysisResults.selection = computeSelectionResults(analysisResults,params.fov_limit, ...
        params.MTF_Requirement,params.SNR_Requirement,params.Cov_Requirement);
end
validateAnalysisAxes(analysisResults,params.alturas_orbitales,params.diametros_pupila, ...
    params.swaths_km,params.telescope_names,params.satellite_configs);
cases = analysisCases(size(params.satellite_configs,1),3,numel(params.telescope_names));
for caseIndex = 1:size(cases,1)
    c = cases(caseIndex,1); d = cases(caseIndex,2); t = cases(caseIndex,3);
    validateAnalysisMatrices(analysisResults.optics.mtf{2,d,t},analysisResults.optics.snr{3,d,t}, ...
        analysisResults.coverage.days{c,d,t},params.alturas_orbitales,params.diametros_pupila,params.swaths_km);
end
summary = computeConfigurationSummary(params,analysisResults);
if writeOutputs, writeConfigurationSummary(summary); end
end
