function results = runGsdSensitivityStudy(config)
%RUNGSDSENSITIVITYSTUDY Recompute the historical fixed-altitude samples.
% Preserve sample-specific constellation masses and the explicitly configured
% two-satellite display assumption. The input table is not overwritten by default.
previousPath=path; pathCleanup=onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
outputDirectory=config.outputDirectory;
sensitivityTable = readtable(config.inputFile);
sensitivityTable = sortrows(sensitivityTable, 'GSD_m');
for sampleIndex = 1:height(sensitivityTable)
    massProperties = computeDryMass(sensitivityTable.Diametro_Pupila_mm(sampleIndex), ...
        sensitivityTable.Num_Telescopios(sampleIndex));
    mass = computeTotalMass(sensitivityTable.Altura_km(sampleIndex),massProperties.dryMass_kg,massProperties.satelliteArea_m2);
    sensitivityTable.Masa_Seca_Satelite_kg(sampleIndex) = massProperties.dryMass_kg;
    sensitivityTable.Masa_Combustible_kg(sampleIndex) = mass.fuelMass_kg;
    sensitivityTable.Masa_Total_Constelacion_kg(sampleIndex) = sensitivityTable.Num_Satelites(sampleIndex) * mass.totalMass_kg;
end
results.config=config; results.summary=sensitivityTable;
if ~config.writeOutputs, return; end
if ~exist(config.outputDirectory,'dir'), mkdir(config.outputDirectory); end
writetable(sensitivityTable,fullfile(config.outputDirectory,'GSD_sensitivity_summary_recomputed.csv'));
if config.updateInputSummary, writetable(sensitivityTable,config.inputFile); end
groundSampleDistances_m = sensitivityTable.GSD_m;
constellationDryMass_kg=config.displaySatelliteCount*sensitivityTable.Masa_Seca_Satelite_kg;
constellationFuelMass_kg=config.displaySatelliteCount*sensitivityTable.Masa_Combustible_kg;
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
    'Location', 'northwest', 'Interpreter', 'latex', 'TextColor', 'black', 'Color', 'white');
grid on;
xlim([25, 125]);
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
