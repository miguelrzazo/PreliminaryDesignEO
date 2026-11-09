function results = runSpectrumStudy(config)
% Explicit study inputs; keep the caller workspace and MATLAB path intact.
previousPath = path;
pathCleanup = onCleanup(@() path(previousPath));
addpath(fullfile(fileparts(mfilename('fullpath')),'lib'));
results.config = config;
% Frequencies retain the original speed-of-light conversion and band ranges.
sBandLower_Hz = config.sBand_Hz(1);
sBandUpper_Hz = config.sBand_Hz(2);
xBandLower_Hz = config.xBand_Hz(1);
xBandUpper_Hz = config.xBand_Hz(2);
speedOfLight_m_s = 3e8;
detectorLower_Hz = speedOfLight_m_s / config.detectorWavelengths_m(2);
detectorUpper_Hz = speedOfLight_m_s / config.detectorWavelengths_m(1);
detectorPosition = 1;
xBandPosition = 2;
sBandPosition = 3;
bandHeight = 0.8;
results.sBand_Hz = [sBandLower_Hz sBandUpper_Hz]; results.xBand_Hz = [xBandLower_Hz xBandUpper_Hz];
results.detectorBand_Hz = [detectorLower_Hz detectorUpper_Hz];
if ~config.writeOutputs, return; end
graphicsCleanup = studyGraphicsDefaults();
figure('Position',[100 100 1200 400]);
semilogx([1e9, 1e15], [0, 0], 'HandleVisibility', 'off');
hold on;
rectangle('Position', [sBandLower_Hz, sBandPosition, sBandUpper_Hz-sBandLower_Hz, bandHeight], ...
          'FaceColor', [0.12, 0.47, 0.71], 'EdgeColor', 'none');
rectangle('Position', [xBandLower_Hz, xBandPosition, xBandUpper_Hz-xBandLower_Hz, bandHeight], ...
          'FaceColor', [1.0, 0.5, 0.05], 'EdgeColor', 'none');
rectangle('Position', [detectorLower_Hz, detectorPosition, detectorUpper_Hz-detectorLower_Hz, bandHeight], ...
          'FaceColor', [0.17, 0.63, 0.17], 'EdgeColor', 'none');
sBandHandle = plot(NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', [0.12, 0.47, 0.71], ...
          'MarkerEdgeColor', 'none');
xBandHandle = plot(NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', [1.0, 0.5, 0.05], ...
          'MarkerEdgeColor', 'none');
detectorHandle = plot(NaN, NaN, 's', 'MarkerSize', 15, 'MarkerFaceColor', [0.17, 0.63, 0.17], ...
          'MarkerEdgeColor', 'none');
xlabel('Frecuencia (Hz, escala logaritmica)');
ylabel('Componentes');
yticks([detectorPosition + bandHeight/2, xBandPosition + bandHeight/2, sBandPosition + bandHeight/2]);
yticklabels({'Detector', 'Banda X', 'Banda S'});
title('Ubicacion de las Bandas S, X y Detector en el Espectro de Frecuencias');
xlim([1e9, 1e15]);
ylim([0.5, 4]);
grid on;
legend([sBandHandle, xBandHandle, detectorHandle], {sprintf('Banda S (%g-%g GHz)',config.sBand_Hz/1e9),sprintf('Banda X (%g-%g GHz)',config.xBand_Hz/1e9),sprintf('Detector (%g-%g um)',config.detectorWavelengths_m/1e-6)}, ...
       'Location', 'best');
hold off;
end
