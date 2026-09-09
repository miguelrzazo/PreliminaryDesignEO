%% 1. DEFINICIÓN DE PARÁMETROS DE SIMULACIÓN Y MISIÓN
startTime = datetime('28-Apr-2024 00:00:00', 'TimeZone', 'UTC');
stopTime = startTime + days(7);
sampleTime = 60; % [s]

gsName = 'Fairbanks'; gsLat = 64.84; gsLon = -147.712; minElevation = 10;

% --- Parámetros de la Constelación y Satélite (SSO) ---
N = 2; altitude = 630; earthRadius = 6371;
semiMajorAxis = (earthRadius + altitude) * 1000;
eccentricity = 0.001; inclination = 97.91;
argOfPeriapsis = 0;
mu = 3.986004418e14;
orbitalPeriod_s = 2 * pi * sqrt(semiMajorAxis^3 / mu);
LTAN_target = 6; % [horas] LTAN de las 6:00 AM

% Parámetros del detector, de la zona de cobertura y de la memoria
GSD = 80; swath_km = 222; numBands = 3; bitsPerBand = 12;
coverageLat_deg = [25, 49];
coverageLon_deg = [-125, -66];
dataSampleTime = 10; % [s], resolución temporal para integrar las trazas
memoryPerSatellite_GB = 2;

%% 2. CREACIÓN DEL ESCENARIO DE SIMULACIÓN
sc = satelliteScenario(startTime, stopTime, sampleTime);

%% 3. DEFINICIÓN DE LA ESTACIÓN DE TIERRA Y SATÉLITES
gs = groundStation(sc, 'Name', gsName, 'Latitude', gsLat, 'Longitude', gsLon, 'MinElevationAngle', minElevation);

fprintf('Calculando RAAN para un LTAN objetivo de %d:00...\n', LTAN_target);
try
    pos_sun_eci = planetEphemeris(juliandate(startTime),'Earth','Sun');
    [alpha_sun_rad, ~, ~] = cart2sph(pos_sun_eci(1), pos_sun_eci(2), pos_sun_eci(3));
    alpha_sun_deg = rad2deg(alpha_sun_rad);
catch
    alpha_sun_deg = approxSunRightAscensionDeg(startTime);
end
raan_offset = (LTAN_target - 12) * 15;
raan_calculated = alpha_sun_deg + raan_offset;
fprintf('RAAN calculado: %.2f grados.\n', raan_calculated);

sats = [];
trueAnomalySeparation = 360 / N;
colors = [[0.8500, 0.3250, 0.0980]; [0.0, 0.4470, 0.7410]; [0.4660, 0.6740, 0.1880]; [0.4940, 0.1840, 0.5560]];
for i = 1:N
    trueAnomaly = (i-1) * trueAnomalySeparation;
    satName = sprintf('Sat %d', i);
    sats = [sats, satellite(sc, semiMajorAxis, eccentricity, inclination, raan_calculated, ...
        argOfPeriapsis, trueAnomaly, ...
        "Name", satName, ...
        "OrbitPropagator", "two-body-keplerian")];
end

%% 4. CÁLCULO DE INTERVALOS DE ACCESO
fprintf('Pre-calculando intervalos de acceso para optimizar la simulación...\n');
accessIntervalsAllSats = cell(N, 1);
for i = 1:N
    accessObject = access(sats(i), gs);
    accessIntervalsAllSats{i} = accessIntervals(accessObject);
end
fprintf('Cálculo de acceso completado.\n');

%% 5. CÁLCULO DEL VOLUMEN DE DATOS A PARTIR DE LAS TRAZAS
fprintf('Calculando el volumen de datos por satélite a partir de las trazas...\n');
[dataBitsBySatellite, traceLengthBySatellite_km] = calculateTraceData( ...
    sats, startTime, stopTime, dataSampleTime, GSD, swath_km, ...
    coverageLat_deg, coverageLon_deg, numBands, bitsPerBand);
totalDataBits = sum(dataBitsBySatellite);
totalDataToMap_GB = totalDataBits / (8 * 1e9);

for i = 1:N
    fprintf('%s: longitud de traza en cobertura = %.1f km, datos = %.3f GB\n', ...
        sats(i).Name, traceLengthBySatellite_km(i), dataBitsBySatellite(i)/(8*1e9));
end
fprintf('Datos totales de la constelación en una semana: %.3f GB\n', totalDataToMap_GB);


%% 6. TASA NECESARIA Y ANÁLISIS DE VIABILIDAD
fprintf('\n--- Análisis de Viabilidad basado en Tiempo de Descarga ---\n');

% --- Cálculo del tiempo de descarga total obtenido por los satélites ---
tiempoDescargaObtenido_s = 0;
contactSecondsBySatellite = zeros(1, N);
for i = 1:N
    if ~isempty(accessIntervalsAllSats{i}) && height(accessIntervalsAllSats{i}) > 0
        % La duración de cada acceso es EndTime - StartTime
        durations = accessIntervalsAllSats{i}.EndTime - accessIntervalsAllSats{i}.StartTime;
        % Se suma la duración total (en segundos) de los accesos para toda la constelación
        contactSecondsBySatellite(i) = sum(seconds(durations));
    end
end
tiempoDescargaObtenido_s = sum(contactSecondsBySatellite);
tiempoDescargaObtenido_h = tiempoDescargaObtenido_s / 3600;
fprintf('Tiempo de descarga OBTENIDO (total constelación): %.2f horas (%.0f segundos).\n', tiempoDescargaObtenido_h, tiempoDescargaObtenido_s);

if any(contactSecondsBySatellite <= 0)
    error('Al menos un satélite no tiene contactos con la estación de Tierra.');
end
if tiempoDescargaObtenido_s <= 0
    error('No se han encontrado contactos con la estación de Tierra.');
end
minimumRateBySatellite_Mbps = dataBitsBySatellite ./ (contactSecondsBySatellite * 1e6);
minimumRate_Mbps = max(minimumRateBySatellite_Mbps);
operationalMargin = 1.25;
designRate_Mbps = minimumRate_Mbps * operationalMargin;
effectiveRate_Mbps = 10; % Tasa efectiva de referencia adoptada para el enlace S [Mbps]
downloadRate_Mbps = effectiveRate_Mbps;
fprintf('Tasa mínima por satélite [Mbps]: %s\n', mat2str(minimumRateBySatellite_Mbps, 4));
fprintf('Tasa mínima de diseño sin margen: %.3f Mbps\n', minimumRate_Mbps);
fprintf('Tasa de diseño con margen operativo del %.0f%%: %.3f Mbps\n', ...
    (operationalMargin - 1) * 100, designRate_Mbps);
fprintf('Tasa efectiva adoptada para el enlace S: %.1f Mbps\n', effectiveRate_Mbps);

totalDataToMap_Mb = totalDataToMap_GB * 8 * 1000;
tiempoDescargaNecesario_s = totalDataToMap_Mb / downloadRate_Mbps;
tiempoDescargaNecesario_h = tiempoDescargaNecesario_s / 3600;
tiempoDescargaPorSat_s = dataBitsBySatellite ./ (downloadRate_Mbps * 1e6);
fprintf('Datos totales a descargar: %.2f GB\n', totalDataToMap_GB);
fprintf('Tiempo de descarga de diseño: %.2f horas (%.0f segundos).\n', ...
    tiempoDescargaNecesario_h, tiempoDescargaNecesario_s);
fprintf('Tiempo de descarga por satélite [h]: %s\n', mat2str(tiempoDescargaPorSat_s / 3600, 4));

% --- Comprobación de viabilidad basada en el tiempo ---
if all(tiempoDescargaPorSat_s <= contactSecondsBySatellite)
    fprintf('RESULTADO (por tiempo): MISIÓN VIABLE. El tiempo de acceso total es suficiente.\n\n');
else
    deficit_s = tiempoDescargaNecesario_s - tiempoDescargaObtenido_s;
    deficit_h = deficit_s / 3600;
    fprintf('RESULTADO (por tiempo): MISIÓN NO VIABLE. Se necesita más tiempo de acceso.\n');
    fprintf('Déficit de tiempo: %.2f horas (%.0f segundos).\n\n', deficit_h, deficit_s);
end



%% 7. SIMULACIÓN DE MEMORIA
fprintf('Iniciando simulación de memoria en GB...\n');
timeVector = startTime:seconds(sampleTime):stopTime;
memoryState = zeros(N, numel(timeVector));
totalSimSeconds = seconds(stopTime - startTime);
dataGeneratedPerSample_GB = (dataBitsBySatellite / (8 * 1e9) / totalSimSeconds) * sampleTime;
downloadPerSample_GB = (downloadRate_Mbps * 1e6 * sampleTime) / (8 * 1e9);

for t = 2:numel(timeVector)
    currentTime = timeVector(t);
    for i = 1:N
        memoryState(i, t) = memoryState(i, t-1);
        memoryState(i, t) = memoryState(i, t) + dataGeneratedPerSample_GB(i);
        isAccess = any(currentTime >= accessIntervalsAllSats{i}.StartTime & currentTime <= accessIntervalsAllSats{i}.EndTime);
        if isAccess
            memoryState(i, t) = max(0, memoryState(i, t) - downloadPerSample_GB);
        end
    end
end

%% 8. GENERACIÓN DE GRÁFICOS Y COMPROBACIÓN FINAL
elapsedDays = days(timeVector - startTime);
figure('Name', 'Estado de Memoria de Satélites');
hold on;
yline(memoryPerSatellite_GB, '--k', 'LineWidth', 1.5, 'DisplayName', 'Limite de memoria');

for i = 1:N
    plot(elapsedDays, memoryState(i,:), 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Satelite %d', i));
end
hold off;
grid on;
ylim([0, max(memoryPerSatellite_GB * 1.08, max(memoryState, [], 'all') * 1.08)]);
xlabel('$t$ [dias]', 'Interpreter', 'latex');
ylabel('Memoria ocupada [GB]', 'Interpreter', 'latex');
title(sprintf('Estado de memoria a bordo ($C_{max}=%.0f$ GB)', memoryPerSatellite_GB), 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
set(gca, 'TickLabelInterpreter', 'latex');
script_dir = fileparts(mfilename('fullpath'));
exportgraphics(gcf, 'Estado_Memoria_Satelites.png', 'Resolution', 300);
exportgraphics(gcf, fullfile(script_dir, '..', 'Latex_Code', '7.Segmento_Tierra', 'memoria.jpg'), 'Resolution', 300);

% Comprobación de viabilidad basada en memoria (del código original)
fprintf('\n--- Análisis de Viabilidad basado en Llenado de Memoria ---\n');
maxMemoryBySatellite_GB = max(memoryState, [], 2);
fprintf('Memoria máxima por satélite [GB]: %s\n', mat2str(maxMemoryBySatellite_GB, 6));
if all(maxMemoryBySatellite_GB <= memoryPerSatellite_GB)
    fprintf('RESULTADO (por memoria): VIABLE. La memoria de los satélites nunca se llena.\n');
else
    fprintf('RESULTADO (por memoria): NO VIABLE. Al menos un satélite ha llenado su memoria.\n');
end


%% 9. GENERACION DE VISUALIZACIONES 3D DEL ESCENARIO
fprintf('\n--- Generando y exportando visualizaciones 3D del escenario ---\n');
generate3DPlot(sc, sats, gs, orbitalPeriod_s, 'Escenario_3D_1_Orbita', colors);
%generate3DPlot(sc, sats, gs, 86400, 'Escenario_3D_1_Dia', colors);
%generate3DPlot(sc, sats, gs, 7*86400, 'Escenario_3D_1_Semana', colors);
%fprintf('Se han generado 3 archivos PNG con las visualizaciones 3D.\n');

%% 10. GENERACION DE TRAZAS 2D CON WORLDMAP
fprintf('\n--- Generando y exportando visualizaciones 2D (worldmap) ---\n');
if exist('worldmap', 'file') ~= 2
    fprintf('Mapping Toolbox no disponible; se omiten las trazas 2D.\n');
else
    try
        load coastlines;
    catch
        coastlat = NaN;
        coastlon = NaN;
    end
    generateWorldmapPlot(sats, gs, startTime, orbitalPeriod_s, 'Traza 2D 1 Orbita', colors, coastlat, coastlon, sampleTime);
    generateWorldmapPlot(sats, gs, startTime, 86400, 'Traza 2D 1 Dia', colors, coastlat, coastlon, sampleTime);
    generateWorldmapPlot(sats, gs, startTime, 7*86400, 'Traza 2D 1 Semana', colors, coastlat, coastlon, sampleTime);
    fprintf('Se han generado 3 archivos PNG con las trazas en 2D.\n');
end

%% --- FUNCIONES AUXILIARES ---
function [dataBitsBySatellite, traceLengthBySatellite_km] = calculateTraceData( ...
    sats, startTime, stopTime, sampleTime, GSD, swath_km, ...
    coverageLat_deg, coverageLon_deg, numBands, bitsPerBand)
    dataBitsBySatellite = zeros(1, numel(sats));
    traceLengthBySatellite_km = zeros(1, numel(sats));
    timeVector = startTime:seconds(sampleTime):stopTime;
    swath_m = swath_km * 1000;
    bitsPerPixel = numBands * bitsPerBand;

    for i = 1:numel(sats)
        lat = zeros(1, numel(timeVector));
        lon = zeros(1, numel(timeVector));
        for k = 1:numel(timeVector)
            position = states(sats(i), timeVector(k), 'CoordinateFrame', 'geographic');
            lat(k) = position(1);
            lon(k) = position(2);
        end

        segmentLength_km = greatCircleSegments(lat(1:end-1), lon(1:end-1), ...
            lat(2:end), lon(2:end));
        midpointLat = 0.5 * (lat(1:end-1) + lat(2:end));
        midpointLon = mod(0.5 * (lon(1:end-1) + lon(2:end)) + 180, 360) - 180;
        insideCoverage = midpointLat >= coverageLat_deg(1) & ...
            midpointLat <= coverageLat_deg(2) & ...
            midpointLon >= coverageLon_deg(1) & midpointLon <= coverageLon_deg(2);

        traceLengthBySatellite_km(i) = sum(segmentLength_km(insideCoverage));
        numberOfPixels = traceLengthBySatellite_km(i) * 1000 * swath_m / GSD^2;
        dataBitsBySatellite(i) = numberOfPixels * bitsPerPixel;
    end
end

function segmentLength_km = greatCircleSegments(lat1_deg, lon1_deg, lat2_deg, lon2_deg)
    earthRadius_km = 6371;
    lat1 = deg2rad(lat1_deg);
    lat2 = deg2rad(lat2_deg);
    deltaLat = lat2 - lat1;
    deltaLon = deg2rad(mod(lon2_deg - lon1_deg + 180, 360) - 180);
    haversine = sin(deltaLat/2).^2 + ...
        cos(lat1) .* cos(lat2) .* sin(deltaLon/2).^2;
    segmentLength_km = 2 * earthRadius_km * asin(sqrt(max(0, haversine)));
end

function generateWorldmapPlot(sats, gs, startTime, duration_s, title_str, colors, coastlat, coastlon, sampleTime)
    fig = figure('Name', title_str, 'NumberTitle', 'off', 'Visible', 'off');
    
    worldmap('north america'); 
    geoshow(coastlat, coastlon, 'Color', 'black');
    hold on;
    
    timeVecPlot = startTime:seconds(sampleTime):(startTime + seconds(duration_s));
    
    plotHandles = gobjects(1, numel(sats));
    legendNames = cell(1, numel(sats));

    for i = 1:numel(sats)
        fprintf('Calculando traza 2D para %s (%s)...\n', sats(i).Name, strrep(title_str, '_', ' '));
        lat = zeros(1, numel(timeVecPlot));
        lon = zeros(1, numel(timeVecPlot));
        for t_idx = 1:numel(timeVecPlot)
            pos_geo_point = states(sats(i), timeVecPlot(t_idx), 'CoordinateFrame', 'geographic');
            lat(t_idx) = pos_geo_point(1);
            lon(t_idx) = pos_geo_point(2);
        end
        plotHandles(i) = geoshow(lat, lon, 'DisplayType', 'line', 'Color', colors(mod(i-1, size(colors,1))+1,:), 'LineWidth', 1.5);
        legendNames{i} = strrep(sats(i).Name, '_', '\_');
    end
    
    geoshow(gs.Latitude, gs.Longitude, 'DisplayType', 'point', 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    
    % 1. Titulo con interprete LaTeX. Se escapan los guiones bajos.
    title(strrep(title_str, '_', '\_'), 'Interpreter', 'latex');
    
    gridm on; mlabel on; plabel on;
    
    % 2. Se obtiene el manejador de los ejes actuales para modificar las etiquetas.
    ax = gca;
    
    % 3. Se establece el interprete LaTeX para las etiquetas de los meridianos y paralelos.
    ax.TickLabelInterpreter = 'latex';
    
    legend(plotHandles, legendNames, 'Interpreter', 'latex', 'Location', 'best');
    hold off;
    
    filename = [title_str, '.png'];
    exportgraphics(fig, filename, 'Resolution', 300);
    close(fig);
end



function generate3DPlot(sc, sats, gs, duration_s, title_str, colors)
    fig = uifigure('Name', strrep(title_str, '_', ' '), 'NumberTitle', 'off', 'Visible', 'off');
    v = satelliteScenarioViewer(sc, 'Parent', fig);
    for i = 1:numel(sats)
        groundTrack(sats(i), 'LeadTime', duration_s, 'TrailTime', 0, ...
            'LeadLineColor', colors(mod(i-1, size(colors,1))+1,:), 'LineWidth', 1.5);
    end
    %v.Globe.CameraTarget = [gs.Latitude, gs.Longitude, 0];
    %v.Globe.CameraPosition = [gs.Latitude, gs.Longitude, 20000e3];
    %v.Globe.CameraUpVector = [0 1 0];
    gs.ShowLabel = true;
    %gs.LabelFontSize = 14;
    %drawnow;
    %filename = [title_str, '.png'];
    %exportgraphics(fig, filename, 'Resolution', 300);
    %close(fig);
end

function alpha_sun_deg = approxSunRightAscensionDeg(t)
    jd = juliandate(t);
    n = jd - 2451545.0;
    L = mod(280.460 + 0.9856474 * n, 360);
    g = mod(357.528 + 0.9856003 * n, 360);
    lambda = L + 1.915 * sind(g) + 0.020 * sind(2*g);
    epsilon = 23.439 - 0.0000004 * n;
    alpha_sun_deg = mod(atan2d(cosd(epsilon) * sind(lambda), cosd(lambda)), 360);
end
