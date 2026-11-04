function validateDesignAxes(altitudes_km, pupilDiameters_mm)
%VALIDATEDESIGNAXES Accept row or column axes and reject invalid coordinates.
validateattributes(altitudes_km, {'numeric'}, {'real','finite','vector','nonempty','positive'});
validateattributes(pupilDiameters_mm, {'numeric'}, {'real','finite','vector','nonempty','positive'});
end
