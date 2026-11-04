function [viableDiameters_mm, theoreticalDiameters_mm] = selectMinimumDiameters( ...
    mtfGrid, snrGrid, coverageGrid_days, altitudes_km, pupilDiameters_mm, ...
    swaths_km, telescopeFov_deg, telescopeCount, telescopeIndex, ...
    mtfRequirement, snrRequirement, revisitRequirement_days)
%SELECTMINIMUMDIAMETERS Preserve the original first-feasible input-order policy.
validateAnalysisMatrices(mtfGrid,snrGrid,coverageGrid_days,altitudes_km,pupilDiameters_mm,swaths_km);
[hasDiameter,diameterIndices] = max(mtfGrid >= mtfRequirement & snrGrid >= snrRequirement,[],2);
pupilDiameters_mm = pupilDiameters_mm(:);
theoreticalDiameters_mm = pupilDiameters_mm(diameterIndices);
theoreticalDiameters_mm(~hasDiameter) = NaN;
[hasSwath,swathIndices] = max(coverageGrid_days <= revisitRequirement_days & ~isnan(coverageGrid_days),[],2);
swaths_km = swaths_km(:);
selectedSwaths_km = swaths_km(swathIndices);
effectiveFov_deg = telescopeFov_deg * (1 + (telescopeCount == 2));
actualFov_deg = 2 * atand(selectedSwaths_km ./ (2 * altitudes_km(:)));
viableDiameters_mm = theoreticalDiameters_mm;
viableDiameters_mm(~hasSwath | actualFov_deg > effectiveFov_deg) = NaN;
if telescopeIndex == 1, viableDiameters_mm(viableDiameters_mm > 200) = NaN; end
end
