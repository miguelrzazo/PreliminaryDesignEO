function mass = computeTotalMass(altitudes_km, dryMass_kg, dragAreas_m2, ...
    maintenanceLowerFactor, maintenanceUpperFactor)
%COMPUTETOTALMASS Preserve the eight-year drag/maintenance/propellant model.
% Each counted cycle is a complete two-burn Hohmann manoeuvre.
% Arrays follow the altitude input shape. Scalar mass/area inputs are reused.
% Atmospheric constants, integration tolerances, margins and cutoff are unchanged.
    earthGravitationalParameter_m3_s2 = 3.986004418e14;
    earthRadius_m = 6378e3;
    dragCoefficient = 2.5;
    specificImpulse_s = 220;
    standardGravity_m_s2 = 9.80665;
    missionDuration_s = 8 * 365.25 * 24 * 3600;
    if nargin < 4, maintenanceLowerFactor = 0.98; end
    if nargin < 5, maintenanceUpperFactor = 1.02; end
    validateattributes(altitudes_km,{'numeric'},{'vector','nonempty','real','finite','positive'});
    validateattributes(dryMass_kg,{'numeric'},{'vector','real','finite','positive'});
    validateattributes(dragAreas_m2,{'numeric'},{'vector','real','finite','nonnegative'});
    assert((isscalar(dryMass_kg) || numel(dryMass_kg)==numel(altitudes_km)) && ...
        (isscalar(dragAreas_m2) || numel(dragAreas_m2)==numel(altitudes_km)), ...
        'eo:MassInputSize','Mass and drag area must be scalar or match the altitude count.');
    validateattributes(maintenanceLowerFactor,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(maintenanceUpperFactor,{'numeric'},{'scalar','real','finite','positive'});
    impulseCounts = zeros(size(altitudes_km));
    totalDeltaV_m_s = zeros(size(altitudes_km));
    fuelMass_kg = zeros(size(altitudes_km));
    totalMass_kg = zeros(size(altitudes_km));
    cycleTimes_s = zeros(size(altitudes_km));
    cycleDeltaV_m_s = zeros(size(altitudes_km));
    unmarginedFuelMass_kg = zeros(size(altitudes_km));
    % Each case requires separate piecewise numerical quadrature.
    for altitudeIndex = 1:numel(altitudes_km)
        altitude_km = altitudes_km(altitudeIndex);
        nominalAltitude_m = altitude_km * 1e3;
        lowerAltitude_m = maintenanceLowerFactor * nominalAltitude_m;
        upperAltitude_m = maintenanceUpperFactor * nominalAltitude_m;
        caseDryMass_kg = getIndexedValue(dryMass_kg, altitudeIndex);
        dragAreaMargin = 1.2;
        dragArea_m2 = getIndexedValue(dragAreas_m2, altitudeIndex) * dragAreaMargin;
        cycleTime_s = decayCycleTime(lowerAltitude_m, upperAltitude_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, caseDryMass_kg);
        if isfinite(cycleTime_s) && cycleTime_s > 0
            impulseCount = floor((missionDuration_s - eps(missionDuration_s)) / cycleTime_s);
        else
            impulseCount = 0;
        end
        lowerRadius_m = earthRadius_m + lowerAltitude_m;
        upperRadius_m = earthRadius_m + upperAltitude_m;
        deltaV_m_s = hohmannDeltaV(lowerRadius_m, upperRadius_m, earthGravitationalParameter_m3_s2);
        caseTotalDeltaV_m_s = impulseCount * deltaV_m_s;
        impulseCounts(altitudeIndex) = impulseCount;
        totalDeltaV_m_s(altitudeIndex) = caseTotalDeltaV_m_s;
        unmarginedFuel_kg = caseDryMass_kg * (exp(caseTotalDeltaV_m_s / (specificImpulse_s * standardGravity_m_s2)) - 1);
        unmarginedFuelMass_kg(altitudeIndex) = unmarginedFuel_kg;
        fuelMass_kg(altitudeIndex) = unmarginedFuel_kg * 1.1;
        totalMass_kg(altitudeIndex) = fuelMass_kg(altitudeIndex) + caseDryMass_kg;
        cycleTimes_s(altitudeIndex) = cycleTime_s;
        cycleDeltaV_m_s(altitudeIndex) = deltaV_m_s;
        if totalMass_kg(altitudeIndex) > 100000
            totalMass_kg(altitudeIndex) = NaN;
            fuelMass_kg(altitudeIndex) = NaN;
            unmarginedFuelMass_kg(altitudeIndex) = NaN;
        end
    end
    mass = struct('totalMass_kg',totalMass_kg,'fuelMass_kg',fuelMass_kg, ...
        'impulseCounts',impulseCounts,'totalDeltaV_m_s',totalDeltaV_m_s, ...
        'cycleTimes_s',cycleTimes_s,'cycleDeltaV_m_s',cycleDeltaV_m_s, ...
        'unmarginedFuelMass_kg',unmarginedFuelMass_kg);

end
    function value = getIndexedValue(values, valueIndex)
        if isscalar(values)
            value = values;
        else
            value = values(valueIndex);
        end
    end
    function decayTime_s = decayCycleTime(lowerAltitude_m, upperAltitude_m, earthRadius_m, earthGravitationalParameter_m3_s2, dragCoefficient, dragArea_m2, satelliteMass_kg)
        if upperAltitude_m <= lowerAltitude_m || dragCoefficient <= 0 || dragArea_m2 <= 0 || satelliteMass_kg <= 0
            decayTime_s = Inf;
            return;
        end
        [baseAltitudes_km, ~, ~] = atmosphereTable();
        lowerAltitude_km = lowerAltitude_m / 1e3;
        upperAltitude_km = upperAltitude_m / 1e3;
        splitAltitudes_km = baseAltitudes_km(baseAltitudes_km > lowerAltitude_km & baseAltitudes_km < upperAltitude_km);
        intervalBounds_m = [lowerAltitude_m, splitAltitudes_km * 1e3, upperAltitude_m];
        decayTime_s = 0;
        % Preserve separate integrations at atmospheric discontinuities.
        for intervalIndex = 1:(numel(intervalBounds_m) - 1)
            intervalStart_m = intervalBounds_m(intervalIndex);
            intervalEnd_m = intervalBounds_m(intervalIndex + 1);
            if intervalEnd_m <= intervalStart_m
                continue;
            end
            integrand = @(altitude_m) satelliteMass_kg ./ (dragCoefficient .* dragArea_m2 .* getAtmosphericDensity(altitude_m) .* sqrt(earthGravitationalParameter_m3_s2 .* (earthRadius_m + altitude_m)));
            decayTime_s = decayTime_s + integral(integrand, intervalStart_m, intervalEnd_m, ...
                'RelTol', 1e-8, 'AbsTol', 1e-3, 'ArrayValued', true);
        end
    end
    function deltaV_m_s = hohmannDeltaV(lowerRadius_m, upperRadius_m, earthGravitationalParameter_m3_s2)
        transferSemimajorAxis_m = (lowerRadius_m + upperRadius_m) / 2;
        lowerCircularSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 / lowerRadius_m);
        transferPerigeeSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 * (2/lowerRadius_m - 1/transferSemimajorAxis_m));
        firstDeltaV_m_s = transferPerigeeSpeed_m_s - lowerCircularSpeed_m_s;
        upperCircularSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 / upperRadius_m);
        transferApogeeSpeed_m_s = sqrt(earthGravitationalParameter_m3_s2 * (2/upperRadius_m - 1/transferSemimajorAxis_m));
        secondDeltaV_m_s = upperCircularSpeed_m_s - transferApogeeSpeed_m_s;
        deltaV_m_s = firstDeltaV_m_s + secondDeltaV_m_s;
    end
    function density_kg_m3 = getAtmosphericDensity(altitude_m)
        altitude_km = max(altitude_m / 1e3, 0);
        [baseAltitudes_km, baseDensities_kg_m3, scaleHeights_km] = atmosphereTable();
        % Count lower breakpoints to reproduce find(...,1,'last') for every sample.
        % Column indexing avoids accidental row/column expansion of table values.
        originalShape = size(altitude_km);
        samples_km = altitude_km(:);
        layerIndices = max(sum(samples_km >= baseAltitudes_km,2),1);
        baseDensities_kg_m3 = baseDensities_kg_m3(:);
        baseAltitudes_km = baseAltitudes_km(:);
        scaleHeights_km = scaleHeights_km(:);
        density_kg_m3 = baseDensities_kg_m3(layerIndices) .* ...
            exp(-(samples_km-baseAltitudes_km(layerIndices)) ./ scaleHeights_km(layerIndices));
        density_kg_m3 = reshape(density_kg_m3,originalShape);
    end
    function [baseAltitudes_km, baseDensities_kg_m3, scaleHeights_km] = atmosphereTable()
        % Original Vallado/Curtis piecewise exponential table, consistent with USSA76.
        baseAltitudes_km = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, ...
                     130, 140, 150, 180, 200, 250, 300, 350, 400, 450, ...
                     500, 600, 700, 800, 900, 1000];
        baseDensities_kg_m3 = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, ...
                    3.206e-4, 8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7, ...
                    9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9, ...
                    5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, ...
                    9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13, ...
                    1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15];
        scaleHeights_km = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, ...
                5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149, ...
                22.523, 29.740, 37.105, 45.546, 53.628, 53.298, ...
                58.515, 60.828, 63.822, 71.835, 88.667, 124.64, ...
                181.05, 268.00];
    end
