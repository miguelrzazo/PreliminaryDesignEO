function [coverageGrid_days, diagnostics] = computeCoverageGrid( ...
    groundSampleDistance_m, altitudes_km, swaths_km, ...
    revisitRequirement_days, detectorPixelCount, telescopeFov_deg, ...
    satelliteCount, telescopeCount, overlapFraction, ...
    cloudFraction, maxDetectorCount, revisitSolver)
%COMPUTECOVERAGEGRID Evaluate analytical revisit coverage without file exports.
% Rows correspond to altitudes_km; columns correspond to swaths_km.
% Geometry is vectorized. The external revisit solver is called per valid cell.
% Rejected cells and solver errors remain NaN; diagnostics distinguish them.
% An optional solver handle permits tests without running the external toolkit.
if nargin < 12, revisitSolver = @RevisitCalc; end
validateCoverageInputs(groundSampleDistance_m, altitudes_km, swaths_km, ...
    revisitRequirement_days, detectorPixelCount, telescopeFov_deg, ...
    satelliteCount, telescopeCount, overlapFraction, cloudFraction, maxDetectorCount);
validateattributes(revisitSolver, {'function_handle'}, {'scalar'});

% J2 nodal precession equals the apparent solar motion.
earthRadius_m = 6378.137e3;
earthGravitationalParameter_m3_s2 = 3.986004418e14;
earthJ2 = 1.082629989052e-3;
sunSynchronousNodalRate_rad_s = 360/(365.2421897*24*3600) * (pi/180);
targetLatitude_deg = 25;
solverDayLimits_days = [1, revisitRequirement_days];
% Named adapter arguments preserve the historical external-solver setup.
useFovConstraint = 1;
useElevationConstraint = 0;
orbitalPlaneCount = 1;
relativePlaneSpacing = 0;
includeDescendingPasses = 1;
useOneSidedViewing = 0;

% Orbital quantities depend only on altitude, not on the selected swath.
altitudes_km = altitudes_km(:);
swaths_km = swaths_km(:).';
semimajorAxes_m = earthRadius_m + altitudes_km * 1000;
meanMotions_rad_s = sqrt(earthGravitationalParameter_m3_s2 ./ semimajorAxes_m.^3);
inclinationCosines = -sunSynchronousNodalRate_rad_s ./ ...
    (1.5 * meanMotions_rad_s * earthJ2 .* (earthRadius_m ./ semimajorAxes_m).^2);
sunSynchronousMask = abs(inclinationCosines) <= 1;
inclinations_rad = nan(size(altitudes_km));
inclinations_rad(sunSynchronousMask) = acos(inclinationCosines(sunSynchronousMask));

% Keep the historical two-telescope FOV rule and detector-width constraint.
effectiveFov_deg = telescopeFov_deg;
if telescopeCount == 2, effectiveFov_deg = 2 * telescopeFov_deg; end
maxSwathByFov_km = 2 * altitudes_km * tand(effectiveFov_deg/2);
detectorSwathLimits_km = (1:maxDetectorCount) * ...
    groundSampleDistance_m * detectorPixelCount / 1000;
% Preserve the original operation order for the admissibility threshold.
maxSwathByDetector_km = maxDetectorCount * groundSampleDistance_m * detectorPixelCount / 1000;
geometryMask = swaths_km <= maxSwathByFov_km & ...
    swaths_km <= maxSwathByDetector_km & sunSynchronousMask;
effectiveSwaths_km = swaths_km * (1 - overlapFraction);
halfFovAngles_deg = atand(effectiveSwaths_km ./ (2 * altitudes_km));

coverageGrid_days = nan(size(geometryMask));
diagnostics.geometryFeasible = geometryMask;
diagnostics.sunSynchronousFeasible = sunSynchronousMask;
diagnostics.maxSwathByFov_km = maxSwathByFov_km;
diagnostics.detectorSwathLimits_km = detectorSwathLimits_km;
diagnostics.requirementExceeded = false(size(geometryMask));
diagnostics.nonfiniteSolverResult = false(size(geometryMask));
diagnostics.solverErrors = cell(size(geometryMask));
diagnostics.solverErrorIdentifiers = cell(size(geometryMask));

% Transpose before find to preserve the original altitude-then-swath order.
% Each iteration adapts one valid pair to the scalar external solver API.
[swathIndices, altitudeIndices] = find(geometryMask.');
for caseIndex = 1:numel(altitudeIndices)
    altitudeIndex = altitudeIndices(caseIndex);
    swathIndex = swathIndices(caseIndex);
    orbitalElements = [semimajorAxes_m(altitudeIndex), 0, ...
        inclinations_rad(altitudeIndex), 0, 0]; % [a, eccentricity, i, RAAN, AoP]
    try
        revisit_days = revisitSolver(orbitalElements, targetLatitude_deg, ...
            useFovConstraint, useElevationConstraint, satelliteCount, ...
            orbitalPlaneCount, relativePlaneSpacing, includeDescendingPasses, ...
            useOneSidedViewing, solverDayLimits_days, ...
            halfFovAngles_deg(altitudeIndex, swathIndex));
        validateattributes(revisit_days, {'numeric'}, {'real','scalar'});
    catch solverError
        diagnostics.solverErrors{altitudeIndex, swathIndex} = solverError.message;
        diagnostics.solverErrorIdentifiers{altitudeIndex, swathIndex} = solverError.identifier;
        continue;
    end

    % Apply the cloud correction after the solver call.
    correctedRevisit_days = revisit_days / (1 - cloudFraction);
    diagnostics.nonfiniteSolverResult(altitudeIndex, swathIndex) = ~isfinite(correctedRevisit_days);
    if correctedRevisit_days > revisitRequirement_days
        diagnostics.requirementExceeded(altitudeIndex, swathIndex) = true;
    else
        coverageGrid_days(altitudeIndex, swathIndex) = correctedRevisit_days;
    end
end
end

function validateCoverageInputs(groundSampleDistance_m, altitudes_km, swaths_km, ...
    revisitRequirement_days, detectorPixelCount, telescopeFov_deg, ...
    satelliteCount, telescopeCount, overlapFraction, cloudFraction, maxDetectorCount)
% RevisitCalc advances in whole days and compares its horizon by equality.
validateattributes(revisitRequirement_days, {'numeric'}, {'scalar','real','positive','finite'});
if revisitRequirement_days ~= fix(revisitRequirement_days)
    error('eo:WholeDayHorizon', 'The current revisit solver requires a whole-day horizon.');
end
validateattributes(groundSampleDistance_m, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(altitudes_km, {'numeric'}, {'vector','nonempty','real','positive','finite'});
validateattributes(swaths_km, {'numeric'}, {'vector','nonempty','real','positive','finite'});
validateattributes(detectorPixelCount, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(telescopeFov_deg, {'numeric'}, {'scalar','real','positive','finite','<',90});
validateattributes(satelliteCount, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(telescopeCount, {'numeric'}, {'scalar','integer','>=',1,'<=',2});
validateattributes(overlapFraction, {'numeric'}, {'scalar','real','finite','>=',0,'<',1});
validateattributes(cloudFraction, {'numeric'}, {'scalar','real','finite','>=',0,'<',1});
validateattributes(maxDetectorCount, {'numeric'}, {'scalar','integer','positive','finite'});
end
