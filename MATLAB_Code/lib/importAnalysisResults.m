function data = importAnalysisResults(altitudes_km, pupilDiameters_mm, swaths_km, telescopeNames, satelliteConfigurations)
%IMPORTANALYSISRESULTS Explicit standalone import boundary for historical CSVs.
% Optical labels are checked against caller coordinates. Historical coverage CSVs
% have no labels, so their row/column axes are supplied explicitly by the caller.
data.axes = struct('altitudes_km',altitudes_km(:),'pupilDiameters_mm',pupilDiameters_mm(:).', ...
    'swaths_km',swaths_km(:).','telescopeNames',{telescopeNames},'satelliteConfigurations',satelliteConfigurations);
data.optics.mtf = cell(3,3,numel(telescopeNames));
data.optics.snr = cell(size(data.optics.mtf));
data.coverage.days = cell(size(satelliteConfigurations,1),3,numel(telescopeNames));
% Optical results are shared by constellations; read each file once.
opticalCases = analysisCases(1,3,numel(telescopeNames));
for caseIndex = 1:size(opticalCases,1)
    detectorIndex = opticalCases(caseIndex,2);
    telescopeIndex = opticalCases(caseIndex,3);
    mtfFilename = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv',detectorIndex,telescopeIndex);
    snrFilename = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv',detectorIndex+3,telescopeIndex);
    data.optics.mtf{2,detectorIndex,telescopeIndex} = readOpticalGrid(mtfFilename,altitudes_km,pupilDiameters_mm);
    data.optics.snr{3,detectorIndex,telescopeIndex} = readOpticalGrid(snrFilename,altitudes_km,pupilDiameters_mm);
end
coverageCases = analysisCases(size(satelliteConfigurations,1),3,numel(telescopeNames));
for caseIndex = 1:size(coverageCases,1)
    configurationIndex = coverageCases(caseIndex,1);
    detectorIndex = coverageCases(caseIndex,2);
    telescopeIndex = coverageCases(caseIndex,3);
    filename = sprintf('coverage/coverage_%d Satelite(s); %d Telescopio(s): %s; Detector %d.csv', ...
        satelliteConfigurations(configurationIndex,1),satelliteConfigurations(configurationIndex,2), ...
        strrep(telescopeNames{telescopeIndex},' ',''),detectorIndex);
    data.coverage.days{configurationIndex,detectorIndex,telescopeIndex} = readmatrix(filename);
    validateAnalysisMatrices(data.optics.mtf{2,detectorIndex,telescopeIndex}, ...
        data.optics.snr{3,detectorIndex,telescopeIndex}, ...
        data.coverage.days{configurationIndex,detectorIndex,telescopeIndex},altitudes_km,pupilDiameters_mm,swaths_km);
end
end

function values = readOpticalGrid(filename,altitudes_km,pupilDiameters_mm)
opticalTable = readtable(filename,'ReadRowNames',true);
expectedRows = cellstr(compose('Alt_%d',altitudes_km(:)));
expectedColumns = cellstr(compose('Diam_%d',pupilDiameters_mm(:)));
assert(isequal(opticalTable.Properties.RowNames(:),expectedRows(:)) && ...
    isequal(opticalTable.Properties.VariableNames(:),expectedColumns(:)), ...
    'eo:ImportedAxisMismatch','Optical CSV labels do not match the requested coordinate order: %s',filename);
values = table2array(opticalTable);
end
