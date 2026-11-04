function mass = computeDryMass(pupilDiameters_mm,telescopeCount)
%COMPUTEDRYMASS Evaluate the original instrument-scaling and dry-mass laws.
% Output arrays preserve the diameter shape. Two telescopes increase instrument
% mass by 1.5; reference estimates and remaining dimensions retain their formulas.
validateattributes(pupilDiameters_mm,{'numeric'},{'vector','nonempty','real','finite','positive'});
validateattributes(telescopeCount,{'numeric'},{'scalar','integer','>=',1,'<=',2});
thematicMapperAperture_mm = 400;
thematicMapperLength_m = 1;
thematicMapperMass_kg = 240;
thematicMapperPower_W = 280;
seosatAperture_mm = 250;
seosatLength_m = 1;
seosatMass_kg = 100;
seosatPower_W = 100;
[tmLength_m, tmArea_m2, tmVolume_m3, tmMass_kg, tmPower_W] = scaleInstrument(pupilDiameters_mm, thematicMapperAperture_mm, thematicMapperLength_m, thematicMapperMass_kg, thematicMapperPower_W);
[seosatLengthScaled_m, seosatArea_m2, seosatVolume_m3, seosatMassScaled_kg, seosatPowerScaled_W] = scaleInstrument(pupilDiameters_mm, seosatAperture_mm, seosatLength_m, seosatMass_kg, seosatPower_W);
instrumentLength_m = (tmLength_m + seosatLengthScaled_m) / 2;
instrumentArea_m2 = (tmArea_m2 + seosatArea_m2) / 2;
instrumentVolume_m3 = (tmVolume_m3 + seosatVolume_m3) / 2;
instrumentMass_kg = (tmMass_kg + seosatMassScaled_kg) / 2;
instrumentPower_W = (tmPower_W + seosatPowerScaled_W) / 2;
if telescopeCount == 2
    instrumentMass_kg = instrumentMass_kg * 1.5;
end
thematicMapperDryMass_kg = 4 * tmMass_kg;
seosatDryMass_kg = 4 * seosatMassScaled_kg;
dryMass_kg = 4 * instrumentMass_kg;
satelliteDensity_kg_m3 = 79;
satelliteVolume_m3 = dryMass_kg ./ satelliteDensity_kg_m3;
satelliteLength_m = satelliteVolume_m3.^(1/3);
satelliteArea_m2 = satelliteLength_m.^2;
cubesatUnits = satelliteVolume_m3 / (0.1^3);
mass = struct('dryMass_kg',dryMass_kg,'thematicMapperDryMass_kg',thematicMapperDryMass_kg, ...
    'seosatDryMass_kg',seosatDryMass_kg,'instrumentLength_m',instrumentLength_m, ...
    'instrumentArea_m2',instrumentArea_m2,'instrumentVolume_m3',instrumentVolume_m3, ...
    'instrumentMass_kg',instrumentMass_kg,'instrumentPower_W',instrumentPower_W, ...
    'satelliteVolume_m3',satelliteVolume_m3,'satelliteLength_m',satelliteLength_m, ...
    'satelliteArea_m2',satelliteArea_m2,'cubesatUnits',cubesatUnits);
end
function [length_m, area_m2, volume_m3, mass_kg, power_W] = scaleInstrument(aperture_mm, referenceAperture_mm, referenceLength_m, referenceMass_kg, referencePower_W)
    apertureRatio = aperture_mm / referenceAperture_mm;
    smallApertureFactor = ones(size(apertureRatio));
    smallApertureFactor(apertureRatio <= 0.5) = 2;
    length_m = apertureRatio .* referenceLength_m;
    area_m2 = length_m.^2;
    volume_m3 = length_m.^3;
    mass_kg = (apertureRatio.^3) .* referenceMass_kg .* smallApertureFactor;
    power_W = (apertureRatio.^3) .* referencePower_W .* smallApertureFactor;
end
