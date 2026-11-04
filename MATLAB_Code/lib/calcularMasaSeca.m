function [dryMass_kg,thematicMapperDryMass_kg,seosatDryMass_kg, ...
    instrumentLength_m,instrumentArea_m2,instrumentVolume_m3, ...
    instrumentMass_kg,instrumentPower_W,satelliteVolume_m3, ...
    satelliteLength_m,satelliteArea_m2,cubesatUnits] = calcularMasaSeca(pupilDiameters_mm,telescopeCount)
%CALCULARMASASECA Preserve the historical positional output interface.
mass = computeDryMass(pupilDiameters_mm,telescopeCount);
dryMass_kg = mass.dryMass_kg;
thematicMapperDryMass_kg = mass.thematicMapperDryMass_kg;
seosatDryMass_kg = mass.seosatDryMass_kg;
instrumentLength_m = mass.instrumentLength_m;
instrumentArea_m2 = mass.instrumentArea_m2;
instrumentVolume_m3 = mass.instrumentVolume_m3;
instrumentMass_kg = mass.instrumentMass_kg;
instrumentPower_W = mass.instrumentPower_W;
satelliteVolume_m3 = mass.satelliteVolume_m3;
satelliteLength_m = mass.satelliteLength_m;
satelliteArea_m2 = mass.satelliteArea_m2;
cubesatUnits = mass.cubesatUnits;
end
