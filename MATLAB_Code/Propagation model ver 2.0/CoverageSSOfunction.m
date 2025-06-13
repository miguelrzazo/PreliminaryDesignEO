function CoverageSSOfunction(GSD, altitudes, swaths, Cov_Requirement, ...
    N_pix, fov_limit, N_sat, N_telescopes, idx_detector, idx_telescopio, telescope_name, ...
    overlapFactor, clearSkyFraction, max_detectors,LTAN_hour)

% Crear array USA
[latGrid, lonGrid, maskUSA] = createUSAGrid(0.1);

% Calcular swath máximo permitido por detectores y telescopio
swath_max_detector = max_detectors * GSD * N_pix / 1000; % km
swath_max_telescope = 2 * fov_limit * (mean(altitudes) * 1000) * pi / 180 / 1000; % km
swath_limit = min(swath_max_detector, swath_max_telescope);

%% Inicializar matriz de resultados
revisitMap = NaN(length(altitudes), length(swaths));

for i = 1:length(altitudes)
    h = altitudes(i);
    for j = 1:length(swaths)
        sw = swaths(j);

        % Verificar si el swath cumple los requisitos
        if sw > swath_limit
            continue;
        end

        fprintf('  Simulando h = %d km, swath = %d km...\n', h, sw);

        % Simulación de cobertura
        daysToCover = simulateCoverage(h, sw, latGrid, lonGrid, maskUSA, ...
                                           overlapFactor, clearSkyFraction, Cov_Requirement,LTAN_hour);

        if daysToCover <= Cov_Requirement
            revisitMap(i,j) = daysToCover;
        else
            revisitMap(i,j) = NaN; % Para que quede blanco en el heatmap
        end
    end
end

%% Guardar CSV
filename_base = sprintf('Detector%d_Telescopio%d_%dsat_%dtel_%dpx_%dkmGSD', ...
    idx_detector, idx_telescopio, N_sat, N_telescopes, N_pix, GSD);

% Crear tabla de cobertura con nombres de filas y columnas
% Crear nombres para filas (alturas orbitales)
row_names = cell(length(altitudes), 1);
for i = 1:length(altitudes)
    row_names{i} = sprintf('Alt_%d', round(altitudes(i)));
end

% Crear nombres para columnas (swaths en km)
col_names = cell(1, length(swaths));
for j = 1:length(swaths)
    col_names{j} = sprintf('Swath_%dkm', round(swaths(j)));
end

% Crear tabla coverage con nombres de filas y columnas
coverage = array2table(revisitMap, 'RowNames', row_names, 'VariableNames', col_names);

% Guardar tabla como CSV con nombres de filas
writetable(coverage, fullfile('Coverage', [filename_base '_resultados.csv']), 'WriteRowNames', true);

fid = fopen(filename_base, 'w');
fprintf(fid, 'Resultados cobertura\n');
fprintf(fid, 'Detector: %d\nTelescopio: %d (%s)\n', idx_detector, idx_telescopio, telescope_name);
fprintf(fid, 'Nº Satélites: %d\nNº Telescopios: %d\n', N_sat, N_telescopes);
fprintf(fid, 'Requisito cobertura: %d días\n', Cov_Requirement);
fprintf(fid, 'Swath máx por detector: %.2f km\n', swath_max_detector);
fprintf(fid, 'Swath máx por telescopio: %.2f km\n', swath_max_telescope);
fclose(fid);

%% Gráficas
figure('Visible', 'off');
imagesc(swaths, altitudes, revisitMap);
set(gca, 'YDir', 'normal');
xlabel('Swath [km]', 'Interpreter', 'latex');
ylabel('Altura orbital [km]', 'Interpreter', 'latex');
title(sprintf('Tiempo de cobertura (Detector %d, Telescopio %d, %d sat, %d tel)', ...
    idx_detector, idx_telescopio, N_sat, N_telescopes), 'Interpreter', 'latex');
colormap(parula);

hold on;

% Límites swath para 2 y 3 detectores (km)
detector_limits = [2, 3] * GSD * N_pix / 1000;

for i = 1:length(detector_limits)
    xline(detector_limits(i), '--w', sprintf('%d detectores', [2,3], (i)), 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'center', 'FontSize', 12);
end

hold off;

h = imagesc(swaths, altitudes, revisitMap);
clim([1 Cov_Requirement]);
colormap(parula);
set(h, 'AlphaData', ~isnan(revisitMap));
colorbar;

colorbar;
clim([1 Cov_Requirement]);
ylabel(colorbar, 'Dias', 'Interpreter', 'latex');

% Colocar en blanco los NaNs (fuera del requerimiento)
ax = gca;
ax.CLim = [0 Cov_Requirement];
set(gca,'Color',[1 1 1]);

saveas(gcf, fullfile('Coverage', [filename_base, '.png']));
close;
end

function [latGrid, lonGrid, maskUSA] = createUSAGrid(resolution_deg)
    lat_range = [24, 50];
    lon_range = [-125, -66];
    [latGrid, lonGrid] = meshgrid(lat_range(1):resolution_deg:lat_range(2), ...
                                  lon_range(1):resolution_deg:lon_range(2));
    maskUSA = true(size(latGrid)); % Máscara uniforme para simplificar
end

function daysToCover = simulateCoverage(h, swath_km, latGrid, lonGrid, mask, overlap, clearFrac, maxDays,LTAN_hour)
    % Inicialización
    earthRadius = 6371; % km
    a = earthRadius + h;
    mu = 398600.4418; % km^3/s^2
    nOrbitsPerDay = floor(86400 / (2*pi*sqrt(a^3 / mu))); 

    % Paso temporal (una órbita completa)
    dt = 86400 / nOrbitsPerDay;
    timeVec = 0:dt:maxDays*86400;

    % Crear malla de puntos a cubrir
    covered = zeros(size(mask));
    totalPoints = sum(mask(:));
    
    % Simulación
    for t = timeVec
        % Calcular posición subsatélite (simplificada como ground track SSO)
        lon0 = mod(-0.9856 * (t/3600 - LTAN_hour*15), 360); % ° longitud subsatélite
        latTrack = linspace(-latGrid(end), latGrid(end), 500);
        lonTrack = mod(lon0 + zeros(size(latTrack)), 360);
        
        % Verificar cobertura de cada punto por swath
        for k = 1:numel(latGrid)
            if ~mask(k), continue; end
            dist = distance(latGrid(k), lonGrid(k), latTrack, lonTrack);
            if any(dist <= swath_km/2 * (1 - overlap))
                if rand() < clearFrac
                    covered(k) = 1;
                end
            end
        end

        % ¿Ya se cubrió todo?
        if sum(covered(:)) == totalPoints
            daysToCover = t / 86400;
            return;
        end
    end

    % Si no se cubrió en maxDays
    daysToCover = NaN;
end