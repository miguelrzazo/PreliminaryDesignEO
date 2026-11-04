function validateAnalysisMatrices(mtfGrid, snrGrid, coverageGrid_days, ...
    altitudes_km, pupilDiameters_mm, swaths_km)
%VALIDATEANALYSISMATRICES Enforce the shapes associated with declared axes.
% NaN cells are allowed because they represent rejection or solver failures.
validateattributes(mtfGrid,{'numeric'},{'2d','real'});
validateattributes(snrGrid,{'numeric'},{'2d','real'});
validateattributes(coverageGrid_days,{'numeric'},{'2d','real'});
opticalShape = [numel(altitudes_km),numel(pupilDiameters_mm)];
coverageShape = [numel(altitudes_km),numel(swaths_km)];
assert(isequal(size(mtfGrid),opticalShape) && isequal(size(snrGrid),opticalShape) && ...
    isequal(size(coverageGrid_days),coverageShape), 'eo:ResultDimensions', ...
    'Current-run matrices do not match their declared axes.');
end
