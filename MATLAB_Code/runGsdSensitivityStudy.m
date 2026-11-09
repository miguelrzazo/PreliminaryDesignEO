function results = runGsdSensitivityStudy(config)
%RUNGSDSENSITIVITYSTUDY Resize the adopted optical architecture at fixed altitude.
previousPath=path; pathCleanup=onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
outputDirectory=config.outputDirectory;
gsd = sort(config.groundSampleDistances_m(:));
validateattributes(gsd,{'numeric'},{'vector','finite','positive','nonempty'});
rows = cell(numel(gsd),15);
for sampleIndex = 1:numel(gsd)
    rows(sampleIndex,:) = sensitivitySample(config,gsd(sampleIndex));
end
sensitivityTable = cell2table(rows,'VariableNames',{'GSD_m','Num_Satelites', ...
    'Num_Telescopios','Detectores_Necesarios','Altura_km','Diametro_Pupila_mm', ...
    'Swath_Requerido_km','Masa_Seca_Satelite_kg','Masa_Combustible_kg', ...
    'Masa_Total_Constelacion_kg','Tipo_Telescopio','MTF','SNR', ...
    'Swath_Disponible_km','Cumple_Limite_Detectores_y_FoV'});
results.config=config; results.summary=sensitivityTable;
if ~config.writeOutputs, return; end
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
writetable(sensitivityTable,fullfile(outputDirectory,'GSD_sensitivity_summary.csv'));
groundSampleDistances_m = sensitivityTable.GSD_m;
constellationDryMass_kg=config.satelliteCount*sensitivityTable.Masa_Seca_Satelite_kg;
constellationFuelMass_kg=config.satelliteCount*sensitivityTable.Masa_Combustible_kg;
pupilDiameters_mm = sensitivityTable.Diametro_Pupila_mm;
figureHandle = figure('Visible', 'off', 'Position', [100, 100, 900, 520]);
yyaxis left;
massBars = bar(groundSampleDistances_m, [constellationDryMass_kg, constellationFuelMass_kg], 0.72, 'stacked', 'LineWidth', 0.4);
massBars(1).FaceColor = [0.25, 0.45, 0.70];
massBars(2).FaceColor = [0.90, 0.55, 0.20];
massBars(1).EdgeColor = [0.15, 0.25, 0.40];
massBars(2).EdgeColor = [0.55, 0.30, 0.05];
hold on;
nominalLine = xline(80, ':', 'GSD nominal', 'LabelVerticalAlignment', 'bottom', ...
    'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);
ylabel('Masa de seleccion de la constelacion [kg]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');
yyaxis right;
diameterLine = plot(groundSampleDistances_m, pupilDiameters_mm, 'k-o', 'LineWidth', 2, 'MarkerSize', 7, ...
    'MarkerFaceColor', 'k');
ylabel('Diametro pupila [mm]', 'Interpreter', 'latex');
set(gca, 'YColor', 'k');
xlabel('GSD [m]', 'Interpreter', 'latex');
title(sprintf('Sensibilidad al GSD a altura fija de %g km', sensitivityTable.Altura_km(1)), ...
    'Interpreter', 'latex', 'Color', 'k');
legend([massBars(1), massBars(2), nominalLine, diameterLine], ...
    {'Masa seca', 'Combustible de mantenimiento', 'GSD nominal', 'Diametro pupila'}, ...
    'Location', 'northeast', 'Interpreter', 'latex', 'TextColor', 'black', 'Color', 'white');
grid on;
xticks(groundSampleDistances_m);
xlim([min(groundSampleDistances_m)-5, max(groundSampleDistances_m)+5]);
figureHandle.Color = 'white';
set(findall(figureHandle, 'Type', 'Axes'), 'Color', 'white', 'XColor', 'black');
figureDirectory=config.figureDirectory;
if ~exist(figureDirectory,'dir'), mkdir(figureDirectory); end
exportgraphics(figureHandle, fullfile(figureDirectory, 'GSD_sensitivity.pdf'), 'ContentType', 'vector');
exportgraphics(figureHandle, fullfile(outputDirectory, 'GSD_sensitivity.pdf'), 'ContentType', 'vector');
exportgraphics(figureHandle, fullfile(outputDirectory, 'GSD_sensitivity.jpg'), 'Resolution', 300);
close(figureHandle);
fprintf('GSD sensitivity data and plot saved at fixed altitude %.0f km.\n', sensitivityTable.Altura_km(1));
end

function row = sensitivitySample(config,gsd)
% The limiting bands are CO2 at 2.01 um (MTF) and O2 at 0.76 um (SNR).
diameters = config.pupilDiameters_mm;
mtf = computeMtfGrid(2.01e-6,config.pixelSize_m,config.detectorMtf, ...
    config.alignmentMtf,gsd,0,config.altitude_km,diameters,config.mtfRequirement);
snr = computeSnrGrid(0.76e-6,config.pixelSize_m,config.quantumEfficiency, ...
    config.transmittance,gsd,0,config.altitude_km,diameters,config.snrRequirement);
index = find(isfinite(mtf) & isfinite(snr),1);
assert(~isempty(index),'eo:GsdOptics','No aperture satisfies MTF and SNR for GSD %g m.',gsd);
diameter = diameters(index);
dry = computeDryMass(diameter,config.telescopesPerSatellite);
mass = computeTotalMass(config.altitude_km,dry.dryMass_kg,dry.satelliteArea_m2);
detectors = max(config.telescopesPerSatellite, ...
    ceil(config.swath_km*1000/(gsd*config.detectorPixelCount)));
availableSwath = min(detectors*gsd*config.detectorPixelCount/1000, ...
    2*config.altitude_km*tand(config.telescopesPerSatellite*config.telescopeFov_deg/2));
pixels = ceil(config.swath_km*1000/(gsd*config.telescopesPerSatellite));
acquiredSwath = pixels*gsd*config.telescopesPerSatellite/1000;
geometryFits = detectors<=config.maxDetectorCount && acquiredSwath<=availableSwath;
row = {gsd,config.satelliteCount,config.telescopesPerSatellite,detectors, ...
    config.altitude_km,diameter,config.swath_km,dry.dryMass_kg,mass.fuelMass_kg, ...
    config.satelliteCount*mass.totalMass_kg,'Refractivo',mtf(index),snr(index), ...
    availableSwath,geometryFits};
end
