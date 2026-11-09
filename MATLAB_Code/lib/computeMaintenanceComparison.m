function comparison = computeMaintenanceComparison(resultRows)
%COMPUTEMAINTENANCECOMPARISON Preserve the two original maintenance-policy cases.
if isempty(resultRows), comparison = table(); return; end
viableMask = cell2mat(resultRows(:, 1));
candidateRows = resultRows(viableMask, :);
if isempty(candidateRows), comparison = table(); return; end
candidateRows = unique(cell2mat(candidateRows(:, [2, 3, 5, 7, 8, 26, 35])), 'rows');
policies = [0.98, 1.00; 0.98, 1.02];
comparisonRows = zeros(size(candidateRows, 1) * 2, 16);
cases = analysisCases(size(candidateRows,1),size(policies,1),1);
for caseIndex = 1:size(cases,1)
    candidateIndex = cases(caseIndex,1); policyIndex = cases(caseIndex,2);
    satelliteCount = candidateRows(candidateIndex, 1);
    telescopeCount = candidateRows(candidateIndex, 2);
    detectorIndex = candidateRows(candidateIndex, 3);
    altitude_km = candidateRows(candidateIndex, 4);
    pupilDiameter_mm = candidateRows(candidateIndex, 5);
    dryMass_kg = candidateRows(candidateIndex, 6);
    dragArea_m2 = candidateRows(candidateIndex, 7);
    mass = computeTotalMass(altitude_km,dryMass_kg,dragArea_m2,policies(policyIndex,1),policies(policyIndex,2));
    comparisonRows(caseIndex,:) = [satelliteCount,telescopeCount,detectorIndex,altitude_km, ...
        pupilDiameter_mm,policies(policyIndex,:),mass.cycleTimes_s/86400,mass.impulseCounts, ...
        mass.cycleDeltaV_m_s,mass.totalDeltaV_m_s,mass.unmarginedFuelMass_kg,mass.fuelMass_kg, ...
        dryMass_kg,dragArea_m2,mass.totalMass_kg];
end
headers = {'Num_Satelites', 'Num_Telescopios', 'ID_Detector', 'Altura_nominal_km', ...
    'Diametro_Pupila_mm', 'Factor_inferior', 'Factor_reset', 'T_ciclo_dias', ...
    'N_impulsos', 'DeltaV_ciclo_ms', 'DeltaV_total_ms', 'Propelente_sin_margen_kg', ...
    'Propelente_final_kg', 'Masa_seca_kg', 'Area_frontal_m2', 'Masa_total_satelite_kg'};
comparison = array2table(comparisonRows,'VariableNames',headers);
end
