function mtfGrid = computeMtfGrid(wavelength_m, pixelSize_m, detectorMtf, ...
    alignmentMtf, groundSampleDistance_m, obstructionRatio, ...
    altitudes_km, pupilDiameters_mm, mtfRequirement)
%COMPUTEMTFGRID Evaluate the historical MTF model without figures or files.
% Rows correspond to altitudes_km; columns correspond to pupilDiameters_mm.
% Values below mtfRequirement become NaN. The threshold itself is accepted.
validateOpticalInputs(wavelength_m, pixelSize_m, groundSampleDistance_m, ...
    obstructionRatio, altitudes_km, pupilDiameters_mm);
validateattributes(detectorMtf, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(alignmentMtf, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(mtfRequirement, {'numeric'}, {'scalar','real','finite','nonnegative'});

% Altitude-dependent quantities are columns; aperture-dependent ones are rows.
altitudes_m = altitudes_km(:) * 1e3;
pupilDiameters_m = pupilDiameters_mm(:).' * 1e-3;
focalLengths_m = (altitudes_m * pixelSize_m) / groundSampleDistance_m;
cutoffFrequencies = pupilDiameters_m / wavelength_m;
nyquistFrequency_per_m = 1 / (2 * pixelSize_m);
sampledFrequencies = nyquistFrequency_per_m * focalLengths_m;

% Preserve the original saturation at the diffraction cutoff.
normalizedFrequencies = min(sampledFrequencies, cutoffFrequencies) ./ cutoffFrequencies;
diffractionMtf = computeDiffractionMtf(normalizedFrequencies, obstructionRatio);

% Keep individual budget factors visible and preserve their original order.
aberrationMtf = 0.95;
manufacturingMtf = 0.98;
vibrationMtf = 0.99;
thermoelasticMtf = 0.95;
marginMtf = 0.9;
otherMtf = marginMtf * vibrationMtf * manufacturingMtf * thermoelasticMtf * ...
    aberrationMtf * detectorMtf;
mtfGrid = diffractionMtf * alignmentMtf * otherMtf;
mtfGrid(mtfGrid < mtfRequirement) = NaN;
end

function diffractionMtf = computeDiffractionMtf(normalizedFrequencies, obstructionRatio)
% Terms A, B and C follow the original annular-aperture equations 6-10 to 6-15.
% Evaluate each expression only on its valid domain, avoiding invalid roots.
apertureTermA = zeros(size(normalizedFrequencies));
obstructionTermB = zeros(size(normalizedFrequencies));
overlapTermC = zeros(size(normalizedFrequencies));

apertureMask = normalizedFrequencies >= 0 & normalizedFrequencies <= 1;
frequencies = normalizedFrequencies(apertureMask);
apertureTermA(apertureMask) = (2/pi) * ...
    (acos(frequencies) - frequencies .* sqrt(1 - frequencies.^2));

% An unobscured aperture has B = C = 0. Skip divisions by obstructionRatio.
if obstructionRatio > 0
    scaledFrequencies = normalizedFrequencies / obstructionRatio;
    angleCosines = (1 + obstructionRatio^2 - 4 * normalizedFrequencies.^2) / ...
        (2 * obstructionRatio);
    overlapAngles_rad = zeros(size(normalizedFrequencies));
    validAngleMask = angleCosines >= -1 & angleCosines <= 1;
    overlapAngles_rad(validAngleMask) = acos(angleCosines(validAngleMask));

    obstructionMask = scaledFrequencies >= 0 & scaledFrequencies <= 1;
    frequencies = scaledFrequencies(obstructionMask);
    obstructionTermB(obstructionMask) = (2 * obstructionRatio^2/pi) * ...
        (acos(frequencies) - frequencies .* sqrt(1 - frequencies.^2));

    fullOverlapMask = normalizedFrequencies > 0 & ...
        normalizedFrequencies <= (1 - obstructionRatio)/2;
    overlapTermC(fullOverlapMask) = -2 * obstructionRatio^2;

    partialOverlapMask = normalizedFrequencies > (1 - obstructionRatio)/2 & ...
        normalizedFrequencies < (1 + obstructionRatio)/2;
    angles_rad = overlapAngles_rad(partialOverlapMask);
    overlapTermC(partialOverlapMask) = (2 * obstructionRatio/pi) * sin(angles_rad) + ...
        ((1 + obstructionRatio^2)/pi) * angles_rad - ...
        ((2 * (1 - obstructionRatio^2))/pi) * ...
        atan((1 + obstructionRatio)/(1 - obstructionRatio) * tan(angles_rad/2)) - ...
        2 * obstructionRatio^2;
end

diffractionMtf = (apertureTermA + obstructionTermB + overlapTermC) / ...
    (1 - obstructionRatio^2);
diffractionMtf(isnan(diffractionMtf) | diffractionMtf < 0) = 0;
end
