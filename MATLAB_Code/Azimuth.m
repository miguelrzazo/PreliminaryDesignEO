% Parametros de lanzamiento
lat_vafb = 34.7;   % Latitud de Vandenberg
lon_vafb = -120.6; % Longitud de Vandenberg

% Limites de azimut WTR (170 a 300 grados)
azimut_min = 170;
azimut_max = 300;
radio_linea_deg = 8;

az_min_rad = deg2rad(azimut_min);
lat_limite_min = lat_vafb + radio_linea_deg * cos(az_min_rad);
lon_limite_min = lon_vafb + radio_linea_deg * sin(az_min_rad) / cosd(lat_vafb);

az_max_rad = deg2rad(azimut_max);
lat_limite_max = lat_vafb + radio_linea_deg * cos(az_max_rad);
lon_limite_max = lon_vafb + radio_linea_deg * sin(az_max_rad) / cosd(lat_vafb);

azimuts_arco = linspace(azimut_min, azimut_max, 80);
radio_arco_deg = 5.5;
lats_arco = lat_vafb + radio_arco_deg * cosd(azimuts_arco);
lons_arco = lon_vafb + radio_arco_deg * sind(azimuts_arco) ./ cosd(lat_vafb);

% Trayectoria de lanzamiento para la inclinacion de diseno.
azimut = 190;
distancia_km = linspace(0, 900, 90);
latitudes = lat_vafb + (distancia_km .* cosd(azimut)) / 111;
longitudes = lon_vafb + (distancia_km .* sind(azimut)) ./ (111 * cosd(lat_vafb));

fig = figure('Color', 'white', 'Position', [100 100 900 650]);
gx = geoaxes(fig);
hold(gx, 'on');
geolimits(gx, [27 42], [-128 -115]);

try
    geobasemap(gx, 'grayland');
catch
    geobasemap(gx, 'streets-light');
end

gx.Grid = 'on';
gx.FontSize = 11;

geoplot(gx, [lat_vafb lat_limite_min], [lon_vafb lon_limite_min], ...
    'k--', 'LineWidth', 1.8, 'DisplayName', 'Limites WTR');
geoplot(gx, [lat_vafb lat_limite_max], [lon_vafb lon_limite_max], ...
    'k--', 'LineWidth', 1.8, 'HandleVisibility', 'off');
geoplot(gx, lats_arco, lons_arco, 'k-', 'LineWidth', 1.6, ...
    'DisplayName', 'Rango permitido');
geoplot(gx, latitudes, longitudes, 'r-', 'LineWidth', 3.0, ...
    'DisplayName', 'Trayectoria 190$^\circ$');
geoscatter(gx, lat_vafb, lon_vafb, 120, [0 0.65 0.15], 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 1.2, 'DisplayName', 'Vandenberg SFB');

text(gx, lat_vafb + 0.35, lon_vafb + 0.2, 'Vandenberg SFB', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black', 'Interpreter', 'latex');
text(gx, lat_limite_min - 0.15, lon_limite_min + 0.15, '$170^\circ$', ...
    'FontSize', 11, 'Color', 'black', 'Interpreter', 'latex');
text(gx, lat_limite_max + 0.2, lon_limite_max - 0.85, '$300^\circ$', ...
    'FontSize', 11, 'Color', 'black', 'Interpreter', 'latex');
text(gx, latitudes(42), longitudes(42) + 0.15, '$190^\circ$', ...
    'FontSize', 12, 'Color', 'red', 'Interpreter', 'latex', 'FontWeight', 'bold');

title(gx, 'Azimut de lanzamiento desde Vandenberg SFB', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Color', 'black');
lgd = legend(gx, 'Location', 'northeast', 'Interpreter', 'latex');
lgd.Color = 'white';
lgd.TextColor = 'black';

script_dir = fileparts(mfilename('fullpath'));
exportgraphics(fig, fullfile(script_dir, '..', 'Latex_Code', '6.Lanzadores', 'azimuth.jpg'), ...
    'Resolution', 300);
close(fig);
