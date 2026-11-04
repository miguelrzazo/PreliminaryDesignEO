function [totalMass_kg,fuelMass_kg,impulseCounts,totalDeltaV_m_s, ...
    cycleTimes_s,cycleDeltaV_m_s,unmarginedFuelMass_kg] = calcularMasaTotal( ...
    altitudes_km,dryMass_kg,dragAreas_m2,maintenanceLowerFactor,maintenanceUpperFactor)
%CALCULARMASATOTAL Preserve the historical positional output interface.
if nargin < 4, maintenanceLowerFactor = 0.98; end
if nargin < 5, maintenanceUpperFactor = 1.02; end
mass = computeTotalMass(altitudes_km,dryMass_kg,dragAreas_m2,maintenanceLowerFactor,maintenanceUpperFactor);
totalMass_kg = mass.totalMass_kg; fuelMass_kg = mass.fuelMass_kg;
impulseCounts = mass.impulseCounts; totalDeltaV_m_s = mass.totalDeltaV_m_s;
cycleTimes_s = mass.cycleTimes_s; cycleDeltaV_m_s = mass.cycleDeltaV_m_s;
unmarginedFuelMass_kg = mass.unmarginedFuelMass_kg;
end
