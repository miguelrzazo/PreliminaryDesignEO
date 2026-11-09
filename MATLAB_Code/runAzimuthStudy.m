function results = runAzimuthStudy(config)
% Explicit study inputs; keep the caller workspace and MATLAB path intact.
plotDefaultsCleanup = configureThesisPlotDefaults();
previousPath = path;
pathCleanup = onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
results.config = config;
% Preserve the original local latitude/longitude launch-corridor approximation.
launchLatitude_deg=config.launchLatitude_deg; launchLongitude_deg=config.launchLongitude_deg;
minimumAzimuth_deg=config.minimumAzimuth_deg; maximumAzimuth_deg=config.maximumAzimuth_deg; boundaryRadius_deg=config.boundaryRadius_deg;
minimumAzimuth_rad = deg2rad(minimumAzimuth_deg);
minimumBoundaryLatitude_deg = launchLatitude_deg + boundaryRadius_deg * cos(minimumAzimuth_rad);
minimumBoundaryLongitude_deg = launchLongitude_deg + boundaryRadius_deg * sin(minimumAzimuth_rad) / cosd(launchLatitude_deg);
maximumAzimuth_rad = deg2rad(maximumAzimuth_deg);
maximumBoundaryLatitude_deg = launchLatitude_deg + boundaryRadius_deg * cos(maximumAzimuth_rad);
maximumBoundaryLongitude_deg = launchLongitude_deg + boundaryRadius_deg * sin(maximumAzimuth_rad) / cosd(launchLatitude_deg);
arcAzimuths_deg = linspace(minimumAzimuth_deg, maximumAzimuth_deg, config.arcSampleCount);
arcRadius_deg = config.arcRadius_deg;
arcLatitudes_deg = launchLatitude_deg + arcRadius_deg * cosd(arcAzimuths_deg);
arcLongitudes_deg = launchLongitude_deg + arcRadius_deg * sind(arcAzimuths_deg) ./ cosd(launchLatitude_deg);
launchAzimuth_deg = config.launchAzimuth_deg;
distances_km = config.distances_km;
trackLatitudes_deg = launchLatitude_deg + (distances_km .* cosd(launchAzimuth_deg)) / 111;
trackLongitudes_deg = launchLongitude_deg + (distances_km .* sind(launchAzimuth_deg)) ./ (111 * cosd(launchLatitude_deg));
results.trackLatitudes_deg = trackLatitudes_deg;
results.trackLongitudes_deg = trackLongitudes_deg;
results.arcLatitudes_deg = arcLatitudes_deg;
results.arcLongitudes_deg = arcLongitudes_deg;
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
figureHandle = figure('Color', 'white', 'Position', [100 100 900 650]);
mapAxes = geoaxes(figureHandle);
hold(mapAxes, 'on');
geolimits(mapAxes, [27 42], [-128 -115]);
try
    geobasemap(mapAxes, 'grayland');
catch
    geobasemap(mapAxes, 'streets-light');
end
mapAxes.Grid = 'on';
mapAxes.FontSize = 11;
geoplot(mapAxes, [launchLatitude_deg minimumBoundaryLatitude_deg], [launchLongitude_deg minimumBoundaryLongitude_deg], ...
    'k--', 'LineWidth', 1.8, 'DisplayName', 'Limites WTR');
geoplot(mapAxes, [launchLatitude_deg maximumBoundaryLatitude_deg], [launchLongitude_deg maximumBoundaryLongitude_deg], ...
    'k--', 'LineWidth', 1.8, 'HandleVisibility', 'off');
geoplot(mapAxes, arcLatitudes_deg, arcLongitudes_deg, 'k-', 'LineWidth', 1.6, ...
    'DisplayName', 'Rango permitido');
geoplot(mapAxes, trackLatitudes_deg, trackLongitudes_deg, 'r-', 'LineWidth', 3.0, ...
    'DisplayName',sprintf('Trayectoria %g deg',config.launchAzimuth_deg));
geoscatter(mapAxes, launchLatitude_deg, launchLongitude_deg, 120, [0 0.65 0.15], 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 1.2, 'DisplayName', 'Vandenberg SFB');
text(mapAxes, launchLatitude_deg + 0.35, launchLongitude_deg + 0.2, 'Vandenberg SFB', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black', 'Interpreter', 'latex');
text(mapAxes, minimumBoundaryLatitude_deg - 0.15, minimumBoundaryLongitude_deg + 0.15, sprintf('$%g^\\circ$',config.minimumAzimuth_deg), ...
    'FontSize', 11, 'Color', 'black', 'Interpreter', 'latex');
text(mapAxes, maximumBoundaryLatitude_deg + 0.2, maximumBoundaryLongitude_deg - 0.85, sprintf('$%g^\\circ$',config.maximumAzimuth_deg), ...
    'FontSize', 11, 'Color', 'black', 'Interpreter', 'latex');
text(mapAxes, trackLatitudes_deg(min(42,numel(trackLatitudes_deg))), trackLongitudes_deg(min(42,numel(trackLongitudes_deg))) + 0.15, sprintf('$%g^\\circ$',config.launchAzimuth_deg), ...
    'FontSize', 12, 'Color', 'red', 'Interpreter', 'latex', 'FontWeight', 'bold');
title(mapAxes, 'Azimut de lanzamiento desde Vandenberg SFB', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Color', 'black');
legendHandle = legend(mapAxes, 'Location', 'northeast', 'Interpreter', 'latex');
legendHandle.Color = 'white';
legendHandle.TextColor = 'black';
exportgraphics(figureHandle,fullfile(config.outputDirectory,'azimuth.jpg'),'Resolution',300);
close(figureHandle);
end
