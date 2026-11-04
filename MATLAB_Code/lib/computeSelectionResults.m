function selection = computeSelectionResults(data, telescopeFov_deg, mtfRequirement, snrRequirement, revisitRequirement_days)
%COMPUTESELECTIONRESULTS Select designs from explicit matrices without file access.
studyAxes = data.axes;
configurationCount = size(studyAxes.satelliteConfigurations,1);
telescopeCount = numel(studyAxes.telescopeNames);
selection = cell(configurationCount,1);
for configurationIndex = 1:configurationCount
    selection{configurationIndex}.viableDiameters_mm = nan(numel(studyAxes.altitudes_km),telescopeCount,3);
    selection{configurationIndex}.theoreticalDiameters_mm = selection{configurationIndex}.viableDiameters_mm;
end
cases = analysisCases(configurationCount,3,telescopeCount);
for caseIndex = 1:size(cases,1)
    configurationIndex = cases(caseIndex,1); detectorIndex = cases(caseIndex,2); telescopeIndex = cases(caseIndex,3);
    [selection{configurationIndex}.viableDiameters_mm(:,telescopeIndex,detectorIndex), ...
        selection{configurationIndex}.theoreticalDiameters_mm(:,telescopeIndex,detectorIndex)] = ...
        selectMinimumDiameters(data.optics.mtf{2,detectorIndex,telescopeIndex}, ...
        data.optics.snr{3,detectorIndex,telescopeIndex},data.coverage.days{configurationIndex,detectorIndex,telescopeIndex}, ...
        studyAxes.altitudes_km,studyAxes.pupilDiameters_mm,studyAxes.swaths_km,telescopeFov_deg(telescopeIndex), ...
        studyAxes.satelliteConfigurations(configurationIndex,2),telescopeIndex,mtfRequirement,snrRequirement,revisitRequirement_days);
end
end
