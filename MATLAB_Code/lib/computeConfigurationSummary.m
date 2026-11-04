function summary = computeConfigurationSummary(parameters, analysisResults)
%COMPUTECONFIGURATIONSUMMARY Preserve candidate equations, row order and mass ties.
% Case enumeration removes nested metadata loops; cached integrations are reused.
if ~isfield(parameters,'maintenance_lower_factor'), parameters.maintenance_lower_factor = 0.98; end
if ~isfield(parameters,'maintenance_upper_factor'), parameters.maintenance_upper_factor = 1.02; end
if ~isfield(parameters,'solapamiento'), parameters.solapamiento = 0.05; end
if ~isfield(parameters,'cobertura_nubes'), parameters.cobertura_nubes = 1/6; end
if ~isfield(parameters,'regionalRevisitSolver'), parameters.regionalRevisitSolver = @RevisitCalc; end
regionalCache = containers.Map('KeyType','char','ValueType','double');
cases = analysisCases(size(parameters.satellite_configs,1),3,numel(parameters.telescope_names),numel(parameters.alturas_orbitales));
rows = cell(2*size(cases,1),numel(configurationHeaders()));
rowCount = 0; viableCases = false(size(cases,1),1);
viableKeys = containers.Map('KeyType','char','ValueType','logical');
massCache = containers.Map('KeyType','char','ValueType','any');
for caseIndex = 1:size(cases,1)
    indices = cases(caseIndex,:);
    diameter_mm = analysisResults.selection{indices(1)}.viableDiameters_mm(indices(4),indices(3),indices(2));
    if isnan(diameter_mm), continue; end
    row = viableRow(parameters,analysisResults,indices,diameter_mm,massCache,regionalCache);
    if isempty(row), continue; end
    rowCount = rowCount+1; rows(rowCount,:) = row; viableCases(caseIndex) = true;
    viableKeys(candidateKey(parameters,indices)) = true;
end
for caseIndex = 1:size(cases,1)
    if viableCases(caseIndex), continue; end
    indices = cases(caseIndex,:);
    if isKey(viableKeys,candidateKey(parameters,indices)), continue; end
    diameter_mm = analysisResults.selection{indices(1)}.theoreticalDiameters_mm(indices(4),indices(3),indices(2));
    if isnan(diameter_mm), continue; end
    row = theoreticalRow(parameters,indices,diameter_mm,massCache);
    if isempty(row), continue; end
    rowCount = rowCount+1; rows(rowCount,:) = row;
end
summary.rows = rows(1:rowCount,:);
summary.best = struct('mass',Inf);
if any(viableCases)
    [summary.best.mass,bestIndex] = min(cell2mat(summary.rows(1:sum(viableCases),29)));
    summary.best.data = summary.rows(bestIndex,:);
end
summary.table = cell2table(summary.rows,'VariableNames',configurationHeaders());
summary.maintenance = computeMaintenanceComparison(summary.rows);
end
function row = viableRow(parameters, analysisResults, indices, pupilDiameter_mm, massCache, regionalCache)
configurationIndex = indices(1); detectorIndex = indices(2); telescopeIndex = indices(3); altitudeIndex = indices(4);
satelliteCount = parameters.satellite_configs(configurationIndex,1); telescopeCount = parameters.satellite_configs(configurationIndex,2);
altitude_km = parameters.alturas_orbitales(altitudeIndex);
coverageRow_days = analysisResults.coverage.days{configurationIndex,detectorIndex,telescopeIndex}(altitudeIndex,:);
validSwathIndices = find(coverageRow_days <= parameters.Cov_Requirement & ~isnan(coverageRow_days));
if isempty(validSwathIndices), row = {}; return; end
groundSampleDistance_m = parameters.GSD; pixelSizes_m = parameters.Pixel_size; localTimeAscendingNode_h = parameters.LTAN_hour;
swaths_km = parameters.swaths_km; telescopeNames = parameters.telescope_names;
maintenanceLowerFactor = parameters.maintenance_lower_factor; maintenanceUpperFactor = parameters.maintenance_upper_factor;
earthRadius_m = 6378.137 * 1000; earthGravitationalParameter_m3_s2 = 3.986004418e14; earthJ2 = 1.082629989052e-3;
[selectedSwath_km,revisit_days] = selectRegionalSwath( ...
    altitude_km,satelliteCount,telescopeCount,swaths_km(validSwathIndices), ...
    parameters,regionalCache);
if isnan(selectedSwath_km), row = {}; return; end
detectorPixelCount = parameters.Npix_detectors(detectorIndex);
requiredDetectorCount = max(telescopeCount,ceil(selectedSwath_km / ...
    (groundSampleDistance_m * detectorPixelCount / 1000)));
pixelsPerTelescope = ceil(selectedSwath_km*1000/(groundSampleDistance_m*telescopeCount));
acquiredSwath_km = telescopeCount*pixelsPerTelescope*groundSampleDistance_m/1000;
actualSwath_km = min(requiredDetectorCount*groundSampleDistance_m*detectorPixelCount/1000, ...
    2*altitude_km*tand(telescopeCount*parameters.fov_limit(telescopeIndex)/2));
if acquiredSwath_km > actualSwath_km, row = {}; return; end
excessSwath_km = actualSwath_km-selectedSwath_km;
excessSwath_percent = 100*excessSwath_km/selectedSwath_km;
diameterIndex = find(parameters.diametros_pupila == pupilDiameter_mm,1);
mtfValue = analysisResults.optics.mtf{2,detectorIndex,telescopeIndex}(altitudeIndex,diameterIndex);
snrValue = analysisResults.optics.snr{3,detectorIndex,telescopeIndex}(altitudeIndex,diameterIndex);
altitude_m = altitude_km * 1000;
pixelSize_m = pixelSizes_m(detectorIndex);
pupilDiameter_m = pupilDiameter_mm / 1000;
focalLength_m = (pixelSize_m * altitude_m) / groundSampleDistance_m;
focalRatio = focalLength_m / pupilDiameter_m;
massProperties = computeDryMass(pupilDiameter_mm,telescopeCount);
[totalMass_kg, fuelMass_kg, impulseCount, totalDeltaV_m_s] = cachedTotalMass(massCache, altitude_km, massProperties.dryMass_kg, massProperties.satelliteArea_m2, maintenanceLowerFactor, maintenanceUpperFactor);
if isnan(totalMass_kg), row = {}; return; end
constellationMass_kg = totalMass_kg * satelliteCount;
semimajorAxis_m = (altitude_km * 1000) + earthRadius_m;
sunSynchronousRate_rad_s = 2 * pi / (365.2421897 * 86400);
inclinationCosine = - (2 * sunSynchronousRate_rad_s * semimajorAxis_m^(3.5)) / (3 * earthJ2 * earthRadius_m^2 * sqrt(earthGravitationalParameter_m3_s2));
inclinationCosine = min(max(inclinationCosine, -1), 1);
inclination_deg = acosd(inclinationCosine);
relativeAnomaly_deg = 360 / satelliteCount;
row = {true, ...
    satelliteCount, telescopeCount, requiredDetectorCount, detectorIndex, telescopeNames{telescopeIndex}, altitude_km, pupilDiameter_mm, ...
    selectedSwath_km, revisit_days, ...
    actualSwath_km, excessSwath_km, excessSwath_percent, ...
    inclination_deg, localTimeAscendingNode_h, relativeAnomaly_deg, ...
    mtfValue, snrValue, ...
    focalLength_m, focalRatio, ...
    massProperties.instrumentMass_kg, massProperties.instrumentLength_m, massProperties.instrumentArea_m2, massProperties.instrumentVolume_m3, massProperties.instrumentPower_W, ...
    massProperties.dryMass_kg, fuelMass_kg, totalMass_kg, constellationMass_kg, ...
    impulseCount, totalDeltaV_m_s, ...
    massProperties.thematicMapperDryMass_kg, massProperties.seosatDryMass_kg, ...
    massProperties.satelliteLength_m, massProperties.satelliteArea_m2, massProperties.satelliteVolume_m3, massProperties.cubesatUnits, acquiredSwath_km
};
end
function row = theoreticalRow(parameters, indices, pupilDiameter_mm, massCache)
configurationIndex = indices(1); detectorIndex = indices(2); telescopeIndex = indices(3);
satelliteCount = parameters.satellite_configs(configurationIndex,1); telescopeCount = parameters.satellite_configs(configurationIndex,2);
altitude_km = parameters.alturas_orbitales(indices(4));
groundSampleDistance_m = parameters.GSD; pixelSizes_m = parameters.Pixel_size; localTimeAscendingNode_h = parameters.LTAN_hour;
swaths_km = parameters.swaths_km; telescopeNames = parameters.telescope_names;
maintenanceLowerFactor = parameters.maintenance_lower_factor; maintenanceUpperFactor = parameters.maintenance_upper_factor;
earthRadius_m = 6378.137 * 1000; earthGravitationalParameter_m3_s2 = 3.986004418e14; earthJ2 = 1.082629989052e-3;
massProperties = computeDryMass(pupilDiameter_mm,telescopeCount);
[totalMass_kg, fuelMass_kg, impulseCount, totalDeltaV_m_s] = cachedTotalMass(massCache, altitude_km, massProperties.dryMass_kg, massProperties.satelliteArea_m2, maintenanceLowerFactor, maintenanceUpperFactor);
if isnan(totalMass_kg), row = {}; return; end
relativeAnomaly_deg = 360 / satelliteCount;
row = {false, satelliteCount, telescopeCount, NaN, detectorIndex, telescopeNames{telescopeIndex}, altitude_km, pupilDiameter_mm, ...
    NaN, NaN, NaN, NaN, NaN, ...
    NaN, localTimeAscendingNode_h, relativeAnomaly_deg, ...
    NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, NaN, NaN, ...
    massProperties.dryMass_kg, fuelMass_kg, totalMass_kg, totalMass_kg * satelliteCount, ...
    impulseCount, totalDeltaV_m_s, massProperties.thematicMapperDryMass_kg, massProperties.seosatDryMass_kg, ...
    massProperties.satelliteLength_m, massProperties.satelliteArea_m2, massProperties.satelliteVolume_m3, massProperties.cubesatUnits, NaN};
end
function [totalMass_kg, fuelMass_kg, impulseCount, totalDeltaV_m_s] = cachedTotalMass(cache, altitude_km, dryMass_kg, dragArea_m2, maintenanceLowerFactor, maintenanceUpperFactor)
    key = sprintf('%.6f_%.9f_%.9f_%.6f_%.6f', altitude_km, dryMass_kg, dragArea_m2, maintenanceLowerFactor, maintenanceUpperFactor);
    if isKey(cache, key)
        values = cache(key);
    else
        mass = computeTotalMass(altitude_km,dryMass_kg,dragArea_m2,maintenanceLowerFactor,maintenanceUpperFactor);
        values = [mass.totalMass_kg,mass.fuelMass_kg,mass.impulseCounts,mass.totalDeltaV_m_s];
        cache(key) = values;
    end
    totalMass_kg = values(1);
    fuelMass_kg = values(2);
    impulseCount = values(3);
    totalDeltaV_m_s = values(4);
end

function key = candidateKey(parameters,indices)
% Preserve the historical key and duplicate-configuration handling.
key = sprintf('%d_%s_%d_%d_%d',indices(2),parameters.telescope_names{indices(3)}, ...
    parameters.alturas_orbitales(indices(4)),parameters.satellite_configs(indices(1),1), ...
    parameters.satellite_configs(indices(1),2));
end
