function validateOpticalInputs(wavelength_m, pixelSize_m, groundSampleDistance_m, ...
    obstructionRatio, altitudes_km, pupilDiameters_mm)
%VALIDATEOPTICALINPUTS Enforce the shared SNR/MTF input contract.
validateattributes(wavelength_m, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(pixelSize_m, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(groundSampleDistance_m, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(obstructionRatio, {'numeric'}, {'scalar','real','finite','>=',0,'<',1});
validateDesignAxes(altitudes_km, pupilDiameters_mm);
end
