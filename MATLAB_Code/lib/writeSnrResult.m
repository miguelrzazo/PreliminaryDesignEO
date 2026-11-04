function writeSnrResult(snrGrid, wavelength_m, pixelSize_m, ...
    quantumEfficiency, transmittance, groundSampleDistance_m, ...
    obstructionRatio, altitudes_km, pupilDiameters_mm, ...
    outputPrefix, telescopeName, detectorIndex, ...
    snrRequirement, model)
assert(isequal(size(snrGrid), [numel(altitudes_km), numel(pupilDiameters_mm)]), ...
    'SNR matrix size does not match altitude/diameter axes.');
% Keep plotting and file operations separate from numerical evaluation.
figureHandle = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
figureHandleCleanup = onCleanup(@() closeStudyFigure(figureHandle));
imageHandle = imagesc(pupilDiameters_mm, altitudes_km, snrGrid);

set(imageHandle, 'AlphaData', ~isnan(snrGrid));
set(gca, 'Color', [1 1 1]); 

colormap(parula);
clim([snrRequirement max(2000,snrRequirement+1)]);

axis xy;
xlabel('Di\''ametro de Pupila (mm)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Altura Orbital (km)', 'FontSize', 12, 'Interpreter', 'latex');
title(sprintf('SNR - %s - Detector %d\n($\\lambda$ = %.2f $\\mu$m, Blanco: SNR $< %g$)', ...
    telescopeName, detectorIndex, wavelength_m*1e6, snrRequirement), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

colorbarHandle = colorbar;
colorbarHandle.Label.String = 'SNR';
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

outputTable = array2table(snrGrid, 'RowNames', rowNames, 'VariableNames', columnNames);
writetable(outputTable, [outputPrefix '_resultados.csv'], 'WriteRowNames', true);

fileId = fopen([outputPrefix '_estadisticas.txt'], 'w');
assert(fileId>=0,'eo:OutputFile','Cannot open the requested text output.');
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '=== ANÁLISIS SNR ===\n\n');
fprintf(fileId, 'Telescopio: %s\n', telescopeName);
fprintf(fileId, 'Detector: %d\n', detectorIndex);
fprintf(fileId, 'Longitud de onda: %.2f μm\n', wavelength_m*1e6);
fprintf(fileId, 'GSD: %d m\n', groundSampleDistance_m);
fprintf(fileId, 'Requirement SNR: ≥ %g\n\n', snrRequirement);

validSnrValues = snrGrid(~isnan(snrGrid));
if ~isempty(validSnrValues)
    fprintf(fileId, 'Estadísticas de valores válidos (SNR ≥ %g):\n', snrRequirement);
    fprintf(fileId, 'Valores válidos: %d de %d (%.1f%%)\n', length(validSnrValues), numel(snrGrid), 100*length(validSnrValues)/numel(snrGrid));
    fprintf(fileId, 'SNR mínimo: %.2f\n', min(validSnrValues));
    fprintf(fileId, 'SNR máximo: %.2f\n', max(validSnrValues));
    fprintf(fileId, 'SNR promedio: %.2f\n', mean(validSnrValues));
    fprintf(fileId, 'SNR mediana: %.2f\n', median(validSnrValues));
else
    fprintf(fileId, 'No hay valores válidos que cumplan el requirement SNR ≥ %g\n', snrRequirement);
end

fprintf(fileId, '\nParámetros del análisis:\n');
fprintf(fileId, 'Rango alturas: %d - %d km\n', min(altitudes_km), max(altitudes_km));
fprintf(fileId, 'Rango diámetros: %d - %d mm\n', min(pupilDiameters_mm), max(pupilDiameters_mm));
fprintf(fileId, 'Obscuración central: %.1f\n', obstructionRatio);
fprintf(fileId, 'Transmitancia: %.2f\n', transmittance);
fprintf(fileId, 'Eficiencia cuántica: %.2f\n', quantumEfficiency);
fprintf(fileId, 'Tamaño pixel: %.1f μm\n', pixelSize_m*1e6);
fprintf(fileId, 'Radiancia referencia: %d W/(m²·sr·μm)\n', model.referenceRadiance_W_m2_sr_um);
fprintf(fileId, 'Ancho de banda: %.1f nm\n', model.spectralBandwidth_um*1e3);

fprintf(fileId, '\nRuidos del sistema (e⁻ RMS):\n');
fprintf(fileId, 'Ruido térmico: %d\n', model.noiseRms_e.dark);
fprintf(fileId, 'Ruido de lectura: %d\n', model.noiseRms_e.read);
fprintf(fileId, 'Ruido preamplificador: %d\n', model.noiseRms_e.preamplifier);
fprintf(fileId, 'Ruido línea video: %d\n', model.noiseRms_e.video);
fprintf(fileId, 'Ruido jitter: %d\n', model.noiseRms_e.jitter);
fprintf(fileId, 'Ruido EMC: %d\n', model.noiseRms_e.electromagnetic);
fprintf(fileId, 'Ruido cuantización: %d\n', model.noiseRms_e.quantization);
fprintf(fileId, 'Ruido no linealidad: %d\n', model.noiseRms_e.nonlinearity);

clear fileCleanup;

fprintf('SNR calculado para %s, Detector %d\n', telescopeName, detectorIndex);

end
