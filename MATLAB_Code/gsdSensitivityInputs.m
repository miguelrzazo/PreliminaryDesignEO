function config = gsdSensitivityInputs(design)
%GSDSENSITIVITYINPUTS Hold the selected H2RG/refractive architecture fixed.
assert(design.ID_Detector==2 && strcmp(design.Tipo_Telescopio,'Refractivo'), ...
    'eo:GsdArchitecture','This sensitivity study uses the adopted H2RG refractive design.');
config.groundSampleDistances_m = 30:10:120;
config.altitude_km = design.Altura_km;
config.satelliteCount = design.Num_Satelites;
config.telescopesPerSatellite = design.Num_Telescopios;
config.swath_km = design.Swath_Requerido_km;
config.pupilDiameters_mm = 20:200;
config.pixelSize_m = 18e-6;
config.detectorMtf = 0.50;
config.quantumEfficiency = 0.70;
config.detectorPixelCount = 2048;
config.alignmentMtf = 0.90;
config.transmittance = 0.80;
config.telescopeFov_deg = 10;
config.mtfRequirement = 0.25;
config.snrRequirement = 400;
config.maxDetectorCount = 3;
config.writeOutputs = true;
end
