function validateAnalysisAxes(analysisResults, altitudes_km, pupilDiameters_mm, ...
    swaths_km, telescopeNames, satelliteConfigurations)
%VALIDATEANALYSISAXES Reject data from a different study or axis ordering.
% Equal dimensions alone cannot detect reordered coordinates or configurations.
expectedAxes = struct('altitudes_km',altitudes_km(:), ...
    'pupilDiameters_mm',pupilDiameters_mm(:).', 'swaths_km',swaths_km(:).', ...
    'telescopeNames',{telescopeNames}, 'satelliteConfigurations',satelliteConfigurations);
axisNames = fieldnames(expectedAxes);
for axisIndex = 1:numel(axisNames)
    axisName = axisNames{axisIndex};
    if ~isfield(analysisResults,'axes') || ~isfield(analysisResults.axes,axisName) || ...
            ~isequal(analysisResults.axes.(axisName),expectedAxes.(axisName))
        error('eo:ResultAxisMismatch', ...
            'Current-run results do not match the requested %s axis.',axisName);
    end
end
end
