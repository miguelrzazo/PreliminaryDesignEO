function results = runContactStudy(config)
%RUNCONTACTSTUDY Size raw acquisition, SGP4 contacts and onboard memory.
% The station equipment rate is documented; service availability is assumed.
previousPath=path; pathCleanup=onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
if exist('satelliteScenario','file')==0 && exist('satelliteScenario','class')==0
    error('eo:MissingSatelliteToolbox','This study requires MATLAB satelliteScenario support.');
end
startTime=config.startTime; stopTime=startTime+days(config.duration_days); sampleTime_s=config.sampleTime_s;
groundStationName=config.groundStationName; groundStationLatitude_deg=config.groundStationLatitude_deg; groundStationLongitude_deg=config.groundStationLongitude_deg;
minimumElevation_deg=config.minimumElevation_deg; satelliteCount=config.satelliteCount; altitude_km=config.altitude_km;
earthRadius_km=config.earthRadius_km; semimajorAxis_m=(earthRadius_km+altitude_km)*1000;
eccentricity=config.eccentricity; inclination_deg=config.inclination_deg; argumentOfPerigee_deg=config.argumentOfPerigee_deg;
earthGravitationalParameter_m3_s2=config.earthGravitationalParameter_m3_s2; orbitalPeriod_s=2*pi*sqrt(semimajorAxis_m^3/earthGravitationalParameter_m3_s2);
localTimeAscendingNode_h=config.localTimeAscendingNode_h; groundSampleDistance_m=config.groundSampleDistance_m;
nominalSwath_km=config.nominalSwath_km; bandCount=config.bandCount; bitsPerBand=config.bitsPerBand;
memoryCapacity_GB=config.memoryCapacity_GB; regionBounds=config.regionBounds;
referenceRate_Mbps=config.referenceRate_Mbps;
bitsPerPixel = bandCount * bitsPerBand;
activePixelsPerTelescope = ceil(nominalSwath_km*1000 / ...
    (groundSampleDistance_m*config.telescopesPerSatellite));
pixelsPerLine = config.telescopesPerSatellite*activePixelsPerTelescope;
acquiredSwath_km = pixelsPerLine*groundSampleDistance_m/1000;
scenario = satelliteScenario(startTime, stopTime, sampleTime_s);
groundStationObject = groundStation(scenario, 'Name', groundStationName, 'Latitude', groundStationLatitude_deg, ...
    'Longitude', groundStationLongitude_deg, 'MinElevationAngle', minimumElevation_deg);
fprintf('Calculando RAAN para LTAN = %d:00...\n', localTimeAscendingNode_h);
try
    posSun = planetEphemeris(juliandate(startTime), 'Earth', 'Sun');
    sunRightAscension_deg = atan2d(posSun(2), posSun(1));
catch
    sunRightAscension_deg = approxSunRightAscensionDeg(startTime);
end
ascendingNode_deg = sunRightAscension_deg + (localTimeAscendingNode_h - 12) * 15;
satellites = [];
anomalySeparation_deg = 360 / satelliteCount;
colors=lines(max(2,satelliteCount)); colors(1:2,:)=[[0.8500,0.3250,0.0980];[0.0,0.4470,0.7410]];
for satelliteIndex = 1:satelliteCount
    trueAnomaly_deg = (satelliteIndex - 1) * anomalySeparation_deg;
    satellites = [satellites, satellite(scenario, semimajorAxis_m, eccentricity, inclination_deg, ...
        ascendingNode_deg, argumentOfPerigee_deg, trueAnomaly_deg, "Name", sprintf('Sat %d', satelliteIndex), ...
        "OrbitPropagator", config.propagator)];
end
timeVector = startTime:seconds(sampleTime_s):stopTime;
accessIntervalsAllSats = cell(satelliteCount, 1);
contactTime_s = zeros(satelliteCount, 1);
dataVolume_GB = zeros(satelliteCount, 1);
trackLength_km = zeros(satelliteCount, 1);
segmentDataBits = cell(satelliteCount, 1);
accessMask = false(satelliteCount, numel(timeVector));
trackLatitudes_deg = cell(satelliteCount, 1);
trackLongitudes_deg = cell(satelliteCount, 1);
fprintf('Calculando contactos y trazas dentro de CONUS...\n');
for satelliteIndex = 1:satelliteCount
    accessObject = access(satellites(satelliteIndex), groundStationObject);
    accessIntervalsAllSats{satelliteIndex} = accessIntervals(accessObject);
    if ~isempty(accessIntervalsAllSats{satelliteIndex}) && height(accessIntervalsAllSats{satelliteIndex}) > 0
        durations = accessIntervalsAllSats{satelliteIndex}.EndTime - accessIntervalsAllSats{satelliteIndex}.StartTime;
        contactTime_s(satelliteIndex) = sum(seconds(durations));
        accessMask(satelliteIndex, :) = accessMaskFromIntervals(accessIntervalsAllSats{satelliteIndex}, timeVector);
    end
    [trackLength_km(satelliteIndex), segmentLengths_km, trackLatitudes_deg{satelliteIndex}, trackLongitudes_deg{satelliteIndex}] = ...
        groundTrackLengthInMask(satellites(satelliteIndex), timeVector, regionBounds);
    segmentDataBits{satelliteIndex} = segmentLengths_km * 1000 / groundSampleDistance_m * pixelsPerLine * bitsPerPixel;
    dataVolume_GB(satelliteIndex) = sum(segmentDataBits{satelliteIndex}) / 8e9;
end
totalData_GB = sum(dataVolume_GB);
minimumRate_Mbps = dataVolume_GB * 8e3 ./ contactTime_s;
missionMinimumRate_Mbps = max(minimumRate_Mbps);
fprintf('\n--- Resultados por satélite ---\n');
for satelliteIndex = 1:satelliteCount
    fprintf('RESULT Sat %d: trace_km=%.3f, data_GB=%.6f, contact_s=%.3f, Rmin_Mbps=%.6f\n', ...
        satelliteIndex, trackLength_km(satelliteIndex), dataVolume_GB(satelliteIndex), contactTime_s(satelliteIndex), minimumRate_Mbps(satelliteIndex));
end
fprintf('RESULT Constellation: data_GB=%.6f, contact_s=%.3f\n', ...
    totalData_GB, sum(contactTime_s));
fprintf('RESULT Link: Rmin_Mbps=%.6f, S_reference_Mbps=%.3f, ratio=%.3f\n', ...
    missionMinimumRate_Mbps, referenceRate_Mbps, ...
    referenceRate_Mbps / missionMinimumRate_Mbps);
if referenceRate_Mbps >= missionMinimumRate_Mbps
    fprintf('RESULT Link viability: VIABLE in S-band.\n');
else
    fprintf('RESULT Link viability: NOT VIABLE in S-band reference case.\n');
end
downloadPerSample_GB = referenceRate_Mbps * 1e6 * sampleTime_s / 8e9;
memory=computeOnboardMemory(segmentDataBits,accessMask,memoryCapacity_GB,downloadPerSample_GB);
memoryState=memory.bounded_GB; unboundedMemoryState=memory.required_GB; overflow_GB=memory.overflow_GB;
for satelliteIndex=1:satelliteCount
    fprintf('RESULT Memory Sat %d: peak_required_GB=%.6f\n',satelliteIndex,max(unboundedMemoryState(satelliteIndex,:)));
end
if any(overflow_GB > 0)
    fprintf('RESULT Memory viability: NOT VIABLE.\n');
else
    fprintf('RESULT Memory viability: VIABLE.\n');
end
results.summary = table((1:satelliteCount).',trackLength_km,dataVolume_GB,contactTime_s, ...
    minimumRate_Mbps,max(memory.required_GB,[],2), ...
    'VariableNames',{'Satellite','Trace_km','Data_GB','Contact_s','MinimumRate_Mbps','PeakMemory_GB'});
results.summary.NominalAltitude_km = repmat(altitude_km,satelliteCount,1);
results.summary.AcquiredSwath_km = repmat(acquiredSwath_km,satelliteCount,1);
results.summary.ActivePixelsPerTelescope = repmat(activePixelsPerTelescope,satelliteCount,1);
results.summary.Propagator = repmat(string(config.propagator),satelliteCount,1);
assert(all(isfinite(results.summary{:,2:6}),'all'),'eo:ContactResult','Non-finite contact result.');
results.config=config; results.time=timeVector; results.contactTime_s=contactTime_s;
results.dataVolume_GB=dataVolume_GB; results.minimumRate_Mbps=minimumRate_Mbps;
results.trackLength_km=trackLength_km; results.memory=memory; results.accessMask=accessMask;
results.trackLatitudes_deg=trackLatitudes_deg; results.trackLongitudes_deg=trackLongitudes_deg;
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
writetable(results.summary,fullfile(config.outputDirectory,['contacts_' char(config.propagator) '.csv']));
scriptDirectory = fileparts(mfilename('fullpath'));
outputDirectory = config.outputDirectory;
elapsedDays = days(timeVector - startTime);
figure('Name', 'Estado de memoria de satélites', 'Visible', 'off');
hold on;
yline(memoryCapacity_GB, '--k', 'LineWidth', 1.2, 'DisplayName', 'Límite de memoria');
for satelliteIndex = 1:satelliteCount
    plot(elapsedDays, memoryState(satelliteIndex, :), 'Color', colors(satelliteIndex, :), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Satélite %d', satelliteIndex));
end
hold off; grid on;
xlabel('Tiempo [días]');
ylabel('Memoria ocupada [GB]', 'Interpreter', 'latex');
title('Estado de memoria a bordo durante una semana');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
exportgraphics(gcf, fullfile(outputDirectory, 'memoria.jpg'), 'Resolution', 300);
close(gcf);
if exist('worldmap', 'file') == 2
    try
        load coastlines;
    catch
        coastlat = NaN; coastlon = NaN;
    end
    figureHandle = figure('Name', 'Traza semanal', 'Visible', 'off');
    worldmap('north america');
    geoshow(coastlat, coastlon, 'Color', 'black');
    hold on;
    for satelliteIndex = 1:satelliteCount
        geoshow(trackLatitudes_deg{satelliteIndex}, trackLongitudes_deg{satelliteIndex}, 'DisplayType', 'line', ...
            'Color', colors(satelliteIndex, :), 'LineWidth', 1.2);
    end
    geoshow(groundStationObject.Latitude, groundStationObject.Longitude, 'DisplayType', 'point', 'Marker', 'o', ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    title('Trazas orbitales durante una semana');
    gridm on; mlabel on; plabel on;
    exportgraphics(figureHandle, fullfile(outputDirectory, 'Traza 2D 1 Semana.jpg'), 'Resolution', 300);
    close(figureHandle);
end
end

function [trackLength_km, segmentLengths_km, latitudes_deg, longitudes_deg] = groundTrackLengthInMask(satelliteObject, timeVector, bounds)
    position = states(satelliteObject,'CoordinateFrame','geographic');
    assert(size(position,2)==numel(timeVector),'eo:ContactTimeGrid','Unexpected state time grid.');
    latitudes_deg = position(1, :);
    longitudes_deg = position(2, :);
    segmentLengths_km = zeros(1, numel(timeVector) - 1);
    for sampleIndex = 1:numel(segmentLengths_km)
        segmentLengths_km(sampleIndex) = clippedSegmentLength(latitudes_deg(sampleIndex), longitudes_deg(sampleIndex), ...
            latitudes_deg(sampleIndex + 1), longitudes_deg(sampleIndex + 1), bounds);
    end
    trackLength_km = sum(segmentLengths_km);
end
function trackLength_km = clippedSegmentLength(lat1, lon1, lat2, lon2, bounds)
    if any(~isfinite([lat1, lon1, lat2, lon2])) || abs(lon2-lon1)>180
        trackLength_km = 0;
        return;
    end
    dx = lon2 - lon1;
    dy = lat2 - lat1;
    tEnter = 0;
    tExit = 1;
    for sampleIndex = 1:2
        if sampleIndex == 1
            d = dx; lower = bounds.lonMin; upper = bounds.lonMax; x = lon1;
        else
            d = dy; lower = bounds.latMin; upper = bounds.latMax; x = lat1;
        end
        if d == 0
            if x < lower || x > upper
                trackLength_km = 0;
                return;
            end
        else
            t1 = (lower - x) / d;
            t2 = (upper - x) / d;
            if t1 > t2
                tmp = t1; t1 = t2; t2 = tmp;
            end
            tEnter = max(tEnter, t1);
            tExit = min(tExit, t2);
        end
    end
    if tEnter > tExit
        trackLength_km = 0;
        return;
    end
    latA = lat1 + tEnter * dy; lonA = lon1 + tEnter * dx;
    latB = lat1 + tExit * dy; lonB = lon1 + tExit * dx;
    trackLength_km = haversineKm(latA, lonA, latB, lonB);
end
function distance_km = haversineKm(lat1, lon1, lat2, lon2)
    radius_km = 6371;
    dLat = deg2rad(lat2 - lat1);
    dLon = deg2rad(lon2 - lon1);
    a = sin(dLat / 2).^2 + cosd(lat1) .* cosd(lat2) .* sin(dLon / 2).^2;
    distance_km = 2 * radius_km * atan2(sqrt(a), sqrt(max(0, 1 - a)));
end
function mask = accessMaskFromIntervals(intervals, timeVector)
    mask = false(size(timeVector));
    for sampleIndex = 1:height(intervals)
        mask = mask | (timeVector >= intervals.StartTime(sampleIndex) & ...
            timeVector <= intervals.EndTime(sampleIndex));
    end
end
function sunRightAscension_deg = approxSunRightAscensionDeg(timeIndex)
    jd = juliandate(timeIndex);
    n = jd - 2451545.0;
    L = mod(280.460 + 0.9856474 * n, 360);
    g = mod(357.528 + 0.9856003 * n, 360);
    lambda = L + 1.915 * sind(g) + 0.020 * sind(2 * g);
    epsilon = 23.439 - 0.0000004 * n;
    sunRightAscension_deg = mod(atan2d(cosd(epsilon) * sind(lambda), cosd(lambda)), 360);
end
