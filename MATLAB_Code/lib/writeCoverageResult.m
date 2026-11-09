function writeCoverageResult(coverageGrid_days, altitudes_km, swaths_km, ...
    revisitRequirement_days, diagnostics, configurationName)
%WRITECOVERAGERESULT Plot and export an already computed coverage matrix.
% Retain the historical coverage directory and filenames for downstream stages.
assert(isequal(size(coverageGrid_days), [numel(altitudes_km), numel(swaths_km)]), ...
    'eo:CoverageAxes', 'Coverage matrix dimensions do not match its axes.');
outputDirectory = 'coverage';
if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
figureHandle = figure('Visible', 'off');
figureCleanup = onCleanup(@() closeIfOpen(figureHandle));
imageHandle = imagesc(swaths_km, altitudes_km, coverageGrid_days);
set(imageHandle, 'AlphaData', ~isnan(coverageGrid_days));
set(gca, 'Color', [1 1 1], 'YDir', 'normal');
colormap([1 1 1; parula(256)]);
hold on;
plot(diagnostics.maxSwathByFov_km, altitudes_km(:), 'r--', ...
    'LineWidth', 1.5, 'DisplayName', 'Límite FOV');

% These few plotting calls correspond to detector annotations, not grid cells.
for detectorCount = 1:numel(diagnostics.detectorSwathLimits_km)
    swathLimit_km = diagnostics.detectorSwathLimits_km(detectorCount);
    if swathLimit_km >= min(swaths_km) && swathLimit_km <= max(swaths_km)
        xline(swathLimit_km, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        text(swathLimit_km + 2, max(altitudes_km) - 15, sprintf('%d Det', detectorCount), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'Color', 'black', 'FontWeight', 'bold');
    end
end
legend('show', 'Location', 'southeast');
hold off;
title(sprintf('Tiempo Revisita - %s\n(Blanco: sin resultado admisible)', ...
    configurationName), 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Swath (km)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 12);
colorbarHandle = colorbar;
ylabel(colorbarHandle, 'Tiempo revisita (dias)', 'Interpreter', 'latex', 'FontSize', 11);
% A one-day horizon still requires distinct lower and upper color limits.
clim([1 max(revisitRequirement_days, 1 + eps(1))]);
axis tight;
imageFilename = fullfile(outputDirectory, sprintf('heatmap_%s.jpg', configurationName));
pdfFilename = fullfile(outputDirectory, sprintf('heatmap_%s.pdf', configurationName));
csvFilename = fullfile(outputDirectory, sprintf('coverage_%s.csv', configurationName));
print(figureHandle, imageFilename, '-djpeg', '-r300');
exportgraphics(figureHandle, pdfFilename, 'ContentType', 'image', 'Resolution', 300);
writematrix(coverageGrid_days, csvFilename);
end

function closeIfOpen(figureHandle)
% Close only this presenter's figure, including when an export fails.
if isgraphics(figureHandle), close(figureHandle); end
end
