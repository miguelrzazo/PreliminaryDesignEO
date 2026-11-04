function result = analisiscontacto(design, outputDir, propagator)
%% Segmento Tierra: volumen de datos, contactos y enlace S
% La generación de datos se calcula por satélite a partir de la longitud de
% traza dentro de CONUS. La longitud se obtiene con segmentos geodésicos
% recortados en la frontera de la región y una muestra temporal de 10 s.

if nargin < 1
    design = struct('Num_Satelites', 2, 'Num_Telescopios', 2, ...
        'Altura_km', 630, 'Swath_Requerido_km', 190, 'GSD', 80);
end
if nargin < 2, outputDir = fullfile(pwd, 'review-results', 'Contacts'); end
if nargin < 3, propagator = 'sgp4'; end
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

startTime = datetime('28-Apr-2024 00:00:00', 'TimeZone', 'UTC');
stopTime = startTime + days(7);
sampleTime = 10; % [s]

gsName = 'Fairbanks';
gsLat = 64.84;
gsLon = -147.712;
minElevation = 10;

% Entradas de misión fijadas en los capítulos anteriores.
N = design.Num_Satelites;
altitude = design.Altura_km; % [km]
earthRadius = 6378.137; % [km]
semiMajorAxis = (earthRadius + altitude) * 1000;
J2 = 1.082629989052e-3; J3 = -2.5327e-6;
omegaSun = 2*pi/(365.2421897*86400);
mu = 3.986004418e14;
inclination = acosd(-2*semiMajorAxis^(7/2)*omegaSun / ...
    (3*J2*(earthRadius*1000)^2*sqrt(mu)));
eccentricity = -J3/(2*J2)*(earthRadius*1000/semiMajorAxis)*sind(inclination);
argOfPeriapsis = 90;
LTAN_target = 6; % [h]

% Datos brutos de un flujo combinado por satélite.
GSD = design.GSD; % [m]
nominalSwath_km = design.Swath_Requerido_km; % swath nominal combinado de los dos telescopios [km].
% La lectura recorta cada telescopio al ancho adquirido; el solape genera datos.
numBands = 3;
bitsPerBand = 12;
memoryPerSatellite_GB = 2;

% Aproximación reproducible de la región continental de EE. UU.
conusBounds = struct('latMin', 25, 'latMax', 49, ...
    'lonMin', -125, 'lonMax', -66);
assert(clippedSegmentLength(30, 179.9, 30, -179.9, conusBounds) == 0, ...
    'A date-line crossing must not create a CONUS acquisition.');

% Tasa adoptada de 10 Mbps. NASA NENUG Rev. 4, tabla 7-2, documenta
% 1 kbps--10 Mbps (FM/PCM) y hasta 20 Mbps (BPSK/OQPSK) para ASF AS1.
% Es una capacidad de equipo; la disponibilidad del servicio se asume.
sBandReferenceCapacity_Mbps = 10;

bitsPerPixel = numBands * bitsPerBand;
activePixelsPerTelescope = ceil(nominalSwath_km * 1000 / ...
    (GSD * design.Num_Telescopios));
pixelsPerLine = design.Num_Telescopios * activePixelsPerTelescope;
acquiredSwath_km = pixelsPerLine * GSD / 1000;

%% Escenario, estación y satélites
sc = satelliteScenario(startTime, stopTime, sampleTime);
gs = groundStation(sc, 'Name', gsName, 'Latitude', gsLat, ...
    'Longitude', gsLon, 'MinElevationAngle', minElevation);

fprintf('Calculando RAAN para LTAN = %d:00...\n', LTAN_target);
try
    posSun = planetEphemeris(juliandate(startTime), 'Earth', 'Sun');
    alphaSun = atan2d(posSun(2), posSun(1));
catch
    alphaSun = approxSunRightAscensionDeg(startTime);
end
raan = alphaSun + (LTAN_target - 12) * 15;

sats = [];
trueAnomalySeparation = 360 / N;
colors = lines(N);
for i = 1:N
    trueAnomaly = (i - 1) * trueAnomalySeparation;
    sats = [sats, satellite(sc, semiMajorAxis, eccentricity, inclination, ...
        raan, argOfPeriapsis, trueAnomaly, "Name", sprintf('Sat %d', i), ...
        "OrbitPropagator", propagator)];
end

%% Contactos y generación de datos por satélite
timeVector = startTime:seconds(sampleTime):stopTime;
accessIntervalsAllSats = cell(N, 1);
contactTime_s = zeros(N, 1);
dataVolume_GB = zeros(N, 1);
trackLength_km = zeros(N, 1);
segmentDataBits = cell(N, 1);
accessMask = false(N, numel(timeVector));
trackLat = cell(N, 1);
trackLon = cell(N, 1);

fprintf('Calculando contactos y trazas dentro de CONUS...\n');
for i = 1:N
    accessObject = access(sats(i), gs);
    accessIntervalsAllSats{i} = accessIntervals(accessObject);
    if ~isempty(accessIntervalsAllSats{i}) && height(accessIntervalsAllSats{i}) > 0
        durations = accessIntervalsAllSats{i}.EndTime - accessIntervalsAllSats{i}.StartTime;
        contactTime_s(i) = sum(seconds(durations));
        accessMask(i, :) = accessMaskFromIntervals(accessIntervalsAllSats{i}, timeVector);
    end

    [trackLength_km(i), segmentLengths_km, trackLat{i}, trackLon{i}] = ...
        groundTrackLengthInMask(sats(i), timeVector, conusBounds);
    segmentDataBits{i} = segmentLengths_km * 1000 / GSD * pixelsPerLine * bitsPerPixel;
    dataVolume_GB(i) = sum(segmentDataBits{i}) / 8e9;
end

totalData_GB = sum(dataVolume_GB);
minimumRate_Mbps = dataVolume_GB * 8e3 ./ contactTime_s;
missionMinimumRate_Mbps = max(minimumRate_Mbps);

fprintf('\n--- Resultados por satélite ---\n');
for i = 1:N
    fprintf('RESULT Sat %d: trace_km=%.3f, data_GB=%.6f, contact_s=%.3f, Rmin_Mbps=%.6f\n', ...
        i, trackLength_km(i), dataVolume_GB(i), contactTime_s(i), minimumRate_Mbps(i));
end
fprintf('RESULT Constellation: data_GB=%.6f, contact_s=%.3f\n', ...
    totalData_GB, sum(contactTime_s));
fprintf('RESULT Link: Rmin_Mbps=%.6f, S_reference_Mbps=%.3f, ratio=%.3f\n', ...
    missionMinimumRate_Mbps, sBandReferenceCapacity_Mbps, ...
    sBandReferenceCapacity_Mbps / missionMinimumRate_Mbps);
if sBandReferenceCapacity_Mbps >= missionMinimumRate_Mbps
    fprintf('RESULT Link viability: VIABLE in S-band.\n');
else
    fprintf('RESULT Link viability: NOT VIABLE in S-band reference case.\n');
end

%% Memoria de a bordo
memoryState = zeros(N, numel(timeVector));
unboundedMemoryState = zeros(N, numel(timeVector));
overflow_GB = zeros(N, 1);
downloadPerSample_GB = sBandReferenceCapacity_Mbps * 1e6 * sampleTime / 8e9;
for t = 2:numel(timeVector)
    for i = 1:N
        generated_GB = segmentDataBits{i}(t - 1) / 8e9;
        downloaded_GB = accessMask(i, t) * downloadPerSample_GB;
        requestedState_GB = max(0, unboundedMemoryState(i, t - 1) + generated_GB - downloaded_GB);
        unboundedMemoryState(i, t) = requestedState_GB;
        overflow_GB(i) = overflow_GB(i) + max(0, requestedState_GB - memoryPerSatellite_GB);
        memoryState(i, t) = min(memoryPerSatellite_GB, requestedState_GB);
    end
end

for i = 1:N
    fprintf('RESULT Memory Sat %d: capacity_GB=%.3f, peak_required_GB=%.6f\n', ...
        i, memoryPerSatellite_GB, max(unboundedMemoryState(i,:)));
end
if any(overflow_GB > 0)
    fprintf('RESULT Memory viability: NOT VIABLE.\n');
else
    fprintf('RESULT Memory viability: VIABLE.\n');
end

%% Figura mínima de memoria y traza semanal
result = table((1:N).', dataVolume_GB, contactTime_s, minimumRate_Mbps, ...
    max(unboundedMemoryState, [], 2), ...
    'VariableNames', {'Satellite','Data_GB','Contact_s','MinimumRate_Mbps','PeakMemory_GB'});
result.AcquiredSwath_km = repmat(acquiredSwath_km, N, 1);
result.Propagator = repmat(string(propagator), N, 1);
assert(all(isfinite(result{:,2:6}), 'all'), 'Non-finite contact result.');
writetable(result, fullfile(outputDir, ['contacts_', char(propagator), '.csv']));
elapsedDays = days(timeVector - startTime);
figure('Name', 'Estado de memoria de satélites', 'Visible', 'off');
hold on;
yline(memoryPerSatellite_GB, '--k', 'LineWidth', 1.2, 'DisplayName', 'Límite de memoria');
for i = 1:N
    plot(elapsedDays, memoryState(i, :), 'Color', colors(i, :), ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Satélite %d', i));
end
hold off; grid on;
xlabel('Tiempo [días]');
ylabel('Memoria ocupada [GB]', 'Interpreter', 'latex');
title('Estado de memoria a bordo durante una semana');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
exportgraphics(gcf, fullfile(outputDir, 'memoria.jpg'), 'Resolution', 300);
close(gcf);

if exist('worldmap', 'file') == 2
    try
        load coastlines;
    catch
        coastlat = NaN; coastlon = NaN;
    end
    fig = figure('Name', 'Traza semanal', 'Visible', 'off');
    worldmap('north america');
    geoshow(coastlat, coastlon, 'Color', 'black');
    hold on;
    for i = 1:N
        geoshow(trackLat{i}, trackLon{i}, 'DisplayType', 'line', ...
            'Color', colors(i, :), 'LineWidth', 1.2);
    end
    geoshow(gs.Latitude, gs.Longitude, 'DisplayType', 'point', 'Marker', 'o', ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    title('Trazas orbitales durante una semana');
    gridm on; mlabel on; plabel on;
    exportgraphics(fig, fullfile(outputDir, 'Traza 2D 1 Semana.jpg'), 'Resolution', 300);
    close(fig);
end

end

%% Funciones auxiliares
function [length_km, segmentLengths_km, lat, lon] = groundTrackLengthInMask(sat, timeVector, bounds)
    position = states(sat, 'CoordinateFrame', 'geographic');
    assert(size(position, 2) == numel(timeVector), 'Unexpected state time grid.');
    lat = position(1, :);
    lon = position(2, :);
    segmentLengths_km = zeros(1, numel(timeVector) - 1);
    for k = 1:numel(segmentLengths_km)
        segmentLengths_km(k) = clippedSegmentLength(lat(k), lon(k), ...
            lat(k + 1), lon(k + 1), bounds);
    end
    length_km = sum(segmentLengths_km);
end

function length_km = clippedSegmentLength(lat1, lon1, lat2, lon2, bounds)
    if any(~isfinite([lat1, lon1, lat2, lon2])) || abs(lon2-lon1) > 180
        length_km = 0;
        return;
    end
    dx = lon2 - lon1;
    dy = lat2 - lat1;
    tEnter = 0;
    tExit = 1;
    for k = 1:2
        if k == 1
            d = dx; lower = bounds.lonMin; upper = bounds.lonMax; x = lon1;
        else
            d = dy; lower = bounds.latMin; upper = bounds.latMax; x = lat1;
        end
        if d == 0
            if x < lower || x > upper
                length_km = 0;
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
        length_km = 0;
        return;
    end
    latA = lat1 + tEnter * dy; lonA = lon1 + tEnter * dx;
    latB = lat1 + tExit * dy; lonB = lon1 + tExit * dx;
    length_km = haversineKm(latA, lonA, latB, lonB);
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
    for k = 1:height(intervals)
        mask = mask | (timeVector >= intervals.StartTime(k) & ...
            timeVector <= intervals.EndTime(k));
    end
end

function alpha_sun_deg = approxSunRightAscensionDeg(t)
    jd = juliandate(t);
    n = jd - 2451545.0;
    L = mod(280.460 + 0.9856474 * n, 360);
    g = mod(357.528 + 0.9856003 * n, 360);
    lambda = L + 1.915 * sind(g) + 0.020 * sind(2 * g);
    epsilon = 23.439 - 0.0000004 * n;
    alpha_sun_deg = mod(atan2d(cosd(epsilon) * sind(lambda), cosd(lambda)), 360);
end
