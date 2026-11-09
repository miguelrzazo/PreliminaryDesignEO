function writeMtfResult(mtfGrid, wavelength_m, pixelSize_m, ...
    detectorMtf, alignmentMtf, groundSampleDistance_m, ...
    obstructionRatio, altitudes_km, pupilDiameters_mm, ...
    outputPrefix, telescopeName, detectorIndex, ...
    mtfRequirement)
assert(isequal(size(mtfGrid), [numel(altitudes_km), numel(pupilDiameters_mm)]), ...
    'MTF matrix size does not match altitude/diameter axes.');
% Keep plotting and file operations separate from numerical evaluation.
figureHandle = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
figureHandleCleanup = onCleanup(@() closeStudyFigure(figureHandle));
imageHandle = imagesc(pupilDiameters_mm, altitudes_km, mtfGrid);

set(imageHandle, 'AlphaData', ~isnan(mtfGrid));
set(gca, 'Color', [1 1 1]);

colormap(parula);
clim([mtfRequirement max(0.35,mtfRequirement+0.1)]);

axis xy;

ylabel('Altura Orbital (km)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('Di\''ametro de Pupila (mm)', 'FontSize', 12, 'Interpreter', 'latex');
title(sprintf('MTF - %s - Detector %d\n($\\lambda$ = %.2f $\\mu$m, Blanco: MTF $< %g$)', ...
    telescopeName, detectorIndex, wavelength_m*1e6, mtfRequirement), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

colorbarHandle = colorbar;
colorbarHandle.Label.String = 'MTF';
colorbarHandle.Label.FontSize = 12;
colorbarHandle.Label.Interpreter = 'latex';

diameterTicks_mm = pupilDiameters_mm(1:50:end);
altitudeTicks_km = altitudes_km(1:5:end);
set(gca, 'XTick', diameterTicks_mm, 'YTick', altitudeTicks_km);

saveas(figureHandle, [outputPrefix '_heatmap.jpg']);
exportgraphics(figureHandle, [outputPrefix '_heatmap.pdf'], 'ContentType', 'image', 'Resolution', 300);
close(figureHandle);

% Preserve historical CSV axis labels and coordinate ordering.
rowNames = cellstr(compose('Alt_%d', altitudes_km(:)));
columnNames = cellstr(compose('Diam_%d', pupilDiameters_mm(:).'));

outputTable = array2table(mtfGrid, 'RowNames', rowNames, 'VariableNames', columnNames);
writetable(outputTable, [outputPrefix '_resultados.csv'], 'WriteRowNames', true);

fileId = fopen([outputPrefix '_estadisticas.txt'], 'w');
assert(fileId>=0,'eo:OutputFile','Cannot open the requested text output.');
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '=== ANÁLISIS MTF ===\n\n');
fprintf(fileId, 'Telescopio: %s\n', telescopeName);
fprintf(fileId, 'Detector: %d\n', detectorIndex);
fprintf(fileId, 'Longitud de onda: %.2f μm\n', wavelength_m*1e6);
fprintf(fileId, 'GSD: %d m\n', groundSampleDistance_m);
fprintf(fileId, 'Requirement MTF: ≥ %g\n\n', mtfRequirement);

validMtfValues = mtfGrid(~isnan(mtfGrid));
if ~isempty(validMtfValues)
    fprintf(fileId, 'Estadísticas de valores válidos (MTF ≥ %g):\n', mtfRequirement);
    fprintf(fileId, 'Valores válidos: %d de %d (%.1f%%)\n', length(validMtfValues), numel(mtfGrid), 100*length(validMtfValues)/numel(mtfGrid));
    fprintf(fileId, 'MTF mínimo: %.4f\n', min(validMtfValues));
    fprintf(fileId, 'MTF máximo: %.4f\n', max(validMtfValues));
    fprintf(fileId, 'MTF promedio: %.4f\n', mean(validMtfValues));
    fprintf(fileId, 'MTF mediana: %.4f\n', median(validMtfValues));
else
    fprintf(fileId, 'No hay valores válidos que cumplan el requirement MTF ≥ %g\n', mtfRequirement);
end

fprintf(fileId, '\nParámetros del análisis:\n');
fprintf(fileId, 'Rango alturas: %d - %d km\n', min(altitudes_km), max(altitudes_km));
fprintf(fileId, 'Rango diámetros: %d - %d mm\n', min(pupilDiameters_mm), max(pupilDiameters_mm));
fprintf(fileId, 'Obscuración central: %.1f\n', obstructionRatio);
fprintf(fileId, 'MTF alineamiento: %.2f\n', alignmentMtf);
fprintf(fileId, 'MTF detector: %.2f\n', detectorMtf);
fprintf(fileId, 'Tamaño pixel: %.1f μm\n', pixelSize_m*1e6);

clear fileCleanup;

fprintf('MTF calculado para %s, Detector %d\n', telescopeName, detectorIndex);

end
