function [snrGrid, model] = computeSnrGrid(wavelength_m, pixelSize_m, quantumEfficiency, ...
    transmittance, groundSampleDistance_m, obstructionRatio, ...
    altitudes_km, pupilDiameters_mm, snrRequirement)
%COMPUTESNRGRID Evaluate the historical SNR model without figures or files.
% Rows correspond to altitudes_km; columns correspond to pupilDiameters_mm.
% Values below snrRequirement become NaN. The threshold itself is accepted.
validateOpticalInputs(wavelength_m, pixelSize_m, groundSampleDistance_m, ...
    obstructionRatio, altitudes_km, pupilDiameters_mm);
validateattributes(quantumEfficiency, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(transmittance, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(snrRequirement, {'numeric'}, {'scalar','real','finite','nonnegative'});

% Keep axis-dependent quantities as a column or row until they are combined.
altitudes_m = altitudes_km(:) * 1e3;
pupilDiameters_m = pupilDiameters_mm(:).' * 1e-3;
focalLengths_m = pixelSize_m * altitudes_m / groundSampleDistance_m;
orbitalSpeeds_m_s = sqrt(3.986e14 ./ (6371e3 + altitudes_m));
integrationTimes_s = groundSampleDistance_m ./ orbitalSpeeds_m_s;

% Preserve the original model constants. Radiance is specified per micrometre,
% so the historical bandwidth value 20e-3 is expressed here in micrometres.
model.spectralBandwidth_um = 20e-3;
model.referenceRadiance_W_m2_sr_um = 100;
model.tdiStages = 1;
model.snrMargin = 0.9;
model.noiseRms_e = struct('dark',50, 'read',100, 'preamplifier',5, ...
    'video',10, 'jitter',5, 'electromagnetic',5, 'quantization',2, 'nonlinearity',2);
planckConstant_J_s = 6.626e-34;
speedOfLight_m_s = 3e8;

% Implicit expansion combines an altitude column with a diameter row.
effectiveDiameters_m = sqrt(pupilDiameters_m.^2 * (1 - obstructionRatio^2));
focalRatios = focalLengths_m ./ effectiveDiameters_m;
irradiance_W_m2 = (pi * transmittance * model.spectralBandwidth_um * ...
    model.referenceRadiance_W_m2_sr_um) ./ (1 + 4 * focalRatios.^2);
pixelArea_m2 = pixelSize_m^2;
photoelectrons = (irradiance_W_m2 * pixelArea_m2 * quantumEfficiency * ...
    wavelength_m * model.tdiStages .* integrationTimes_s) / ...
    (planckConstant_J_s * speedOfLight_m_s);

% Independent noise variances add; photoelectrons supply shot-noise variance.
noise = model.noiseRms_e;
fixedNoiseVariance_e2 = noise.dark^2 + noise.read^2 + noise.preamplifier^2 + ...
    noise.video^2 + noise.jitter^2 + noise.electromagnetic^2 + ...
    noise.quantization^2 + noise.nonlinearity^2;
totalNoiseRms_e = sqrt(fixedNoiseVariance_e2 + photoelectrons);
snrGrid = model.snrMargin * photoelectrons ./ totalNoiseRms_e;
snrGrid(snrGrid < snrRequirement) = NaN;
end
