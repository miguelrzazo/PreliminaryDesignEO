function writeSelectionResults(selection, altitudes_km, telescopeNames, satelliteConfigurations, requirements)
%WRITESELECTIONRESULTS Export supplied selections without recomputing models.
outputDirectory = 'HvsDmin';
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
for configurationIndex = 1:size(satelliteConfigurations,1)
    writeSelectionCase(selection,altitudes_km,telescopeNames,satelliteConfigurations, ...
        requirements,configurationIndex,outputDirectory);
end
end

function writeSelectionCase(selection,altitudes_km,telescopeNames,satelliteConfigurations, ...
    requirements,configurationIndex,outputDirectory)
    satelliteCount = satelliteConfigurations(configurationIndex,1); telescopeCount = satelliteConfigurations(configurationIndex,2);
    configurationName = sprintf('%dsat_%dtel',satelliteCount,telescopeCount);
    viableDiameters_mm = selection{configurationIndex}.viableDiameters_mm;
    theoreticalDiameters_mm = selection{configurationIndex}.theoreticalDiameters_mm;
    theoreticalFigure = figure('Position', [100, 100, 900, 700], 'Visible', 'off');
    theoreticalFigureCleanup = onCleanup(@() closeStudyFigure(theoreticalFigure));
    configurationTitle = sprintf('Configuracion: %d satelite(s), %d telescopio(s) por satelite', satelliteCount, telescopeCount);
    for detectorIndex = 1:3
        subplot(3, 1, detectorIndex);
        hold on;
        legendEntries = plotSelectionCurves(altitudes_km,theoreticalDiameters_mm(:,:,detectorIndex),telescopeNames);
        hold off;
        title(sprintf('Detector %d: minimo teorico por MTF/SNR - %s', detectorIndex, configurationTitle), ...
              'Interpreter', 'latex', 'FontSize', 10);
        xlabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 9);
        ylabel('$D_{min}$ teorico (mm)', 'Interpreter', 'latex', 'FontSize', 9);
        if ~isempty(legendEntries)
            legend(legendEntries, 'Location', 'northwest', 'Interpreter', 'latex');
        end
        grid on;
        axis tight;
    end
    theoreticalPlotPrefix = sprintf('%s/HvsDmin_Teorico_%s', outputDirectory, configurationName);
    theoreticalFigure.Color = 'white';
    set(findall(theoreticalFigure, 'Type', 'Axes'), 'Color', 'white');
    saveas(theoreticalFigure, [theoreticalPlotPrefix '.jpg']);
    exportgraphics(theoreticalFigure, [theoreticalPlotPrefix '.pdf'], 'ContentType', 'vector');
    close(theoreticalFigure);
    mainFigure = figure('Position', [100, 100, 900, 700], 'Visible', 'off');
    mainFigureCleanup = onCleanup(@() closeStudyFigure(mainFigure));
    for detectorIndex = 1:3
        subplot(3, 1, detectorIndex);
        hold on;
        legendEntries = plotSelectionCurves(altitudes_km,viableDiameters_mm(:,:,detectorIndex),telescopeNames);
        hold off;
        titleText = sprintf('$h$ vs $D_{min}$ viable - Detector %d - %s', detectorIndex, configurationTitle);
        title(titleText, 'Interpreter', 'latex', 'FontSize', 10);
        xlabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 9);
        ylabel('$D_{min}$ (mm)', 'Interpreter', 'latex', 'FontSize', 9);
        if ~isempty(legendEntries)
            legend(legendEntries, 'Location', 'northwest', 'Interpreter', 'latex');
        end
        grid on;
        axis tight;
    end
    mainPlotPrefix = sprintf('%s/HvsDmin_Viable_%s', outputDirectory, configurationName);
    mainPlotFilename = [mainPlotPrefix '.jpg'];
    saveas(mainFigure, mainPlotFilename);
    saveas(mainFigure, sprintf('%s/HvsDmin_%s.jpg', outputDirectory, configurationName));
    mainFigure.Color = 'white';
    set(findall(mainFigure, 'Type', 'Axes'), 'Color', 'white');
    exportgraphics(mainFigure, [mainPlotPrefix '.pdf'], 'ContentType', 'vector');
    close(mainFigure);
    if length(telescopeNames) >= 1
        refractiveFigure = figure('Position', [100, 100, 700, 500], 'Visible', 'off');
        refractiveFigureCleanup = onCleanup(@() closeStudyFigure(refractiveFigure));
        hold on;
        colors = {'r', 'g', 'b'};
        lineStyles = {'--', ':', '-.'};
        refractiveLegend = {};
        for detectorIndex = 1:3
            refractiveDiameters_mm = viableDiameters_mm(:, 1, detectorIndex);
            validRefractiveMask = ~isnan(refractiveDiameters_mm);
            if any(validRefractiveMask)
                plot(altitudes_km(validRefractiveMask), refractiveDiameters_mm(validRefractiveMask), ...
                     [colors{detectorIndex}, lineStyles{detectorIndex}], 'LineWidth', 2);
                refractiveLegend{end+1} = sprintf('Detector %d', detectorIndex);
            end
        end
        hold off;
        if ~isempty(refractiveLegend)
            refractiveTitle = sprintf('$h$ vs $D_{min}$ - %s - %s', telescopeNames{1}, configurationTitle);
            title(refractiveTitle, 'Interpreter', 'latex', 'FontSize', 10);
            xlabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 9);
            ylabel(sprintf('$D_{min}$ %s (mm)', telescopeNames{1}), 'Interpreter', 'latex', 'FontSize', 9);
            legend(refractiveLegend, 'Location', 'best', 'Interpreter', 'latex');
            grid on;
            axis tight;
            refractivePlotPrefix = sprintf('%s/HvsDmin_Viable_%s_Config_%s', outputDirectory, strrep(telescopeNames{1},' ','_'), configurationName);
            refractivePlotFilename = [refractivePlotPrefix '.jpg'];
            saveas(refractiveFigure, refractivePlotFilename);
            saveas(refractiveFigure, sprintf('%s/HvsDmin_%s_Config_%s.jpg', outputDirectory, strrep(telescopeNames{1},' ','_'), configurationName));
            refractiveFigure.Color = 'white';
            set(findall(refractiveFigure, 'Type', 'Axes'), 'Color', 'white');
            exportgraphics(refractiveFigure, [refractivePlotPrefix '.pdf'], 'ContentType', 'vector');
        end
        close(refractiveFigure);
    end
    for detectorIndex = 1:3
        columnNames = [{'Altura_km'}, cellfun(@(x) sprintf('Dmin_%s_mm', strrep(x,' ','_')), telescopeNames, 'UniformOutput', false)];
        outputMatrix = [altitudes_km(:), viableDiameters_mm(:, :, detectorIndex)];
        outputTable = array2table(outputMatrix, 'VariableNames', columnNames);
        csvFilename = sprintf('HvsDmin/HvsDmin_%s_Detector%d.csv', configurationName, detectorIndex);
        writetable(outputTable, csvFilename);
        theoreticalMatrix = [altitudes_km(:), theoreticalDiameters_mm(:, :, detectorIndex)];
        theoreticalTable = array2table(theoreticalMatrix, 'VariableNames', columnNames);
        theoreticalCsvFilename = sprintf('HvsDmin/HvsDmin_Teorico_%s_Detector%d.csv', configurationName, detectorIndex);
        writetable(theoreticalTable, theoreticalCsvFilename);
        textFilename = sprintf('%s/HvsDmin_%s_Detector%d.txt', outputDirectory, configurationName, detectorIndex);
        fileId = fopen(textFilename, 'w');
        assert(fileId>=0,'eo:OutputFile','Cannot open the requested text output.');
        fileCleanup = onCleanup(@() fclose(fileId));
        fprintf(fileId, 'Altura_km');
        fprintf(fileId,'\t%s',columnNames{2:end});
        fprintf(fileId, '\n');
        rowFormat = ['%.1f',repmat('\t%.1f',1,size(outputMatrix,2)-1),'\n'];
        fprintf(fileId,rowFormat,outputMatrix.');
        clear fileCleanup;
    end
    summaryFilename = sprintf('%s/RESUMEN_Config_%s.txt', outputDirectory, configurationName);
    fileId = fopen(summaryFilename, 'w');
    assert(fileId>=0,'eo:OutputFile','Cannot open the requested text output.');
    fileCleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '=== RESUMEN DE ANALISIS CRUZADO ===\n\n');
    fprintf(fileId, 'Configuracion: %s (%d satelite(s), %d telescopio(s) por satelite)\n\n', ...
            configurationName, satelliteCount, telescopeCount);
    fprintf(fileId, 'Criterios aplicados para Dmin:\n');
    fprintf(fileId, '- MTF (Lambda2) >= %g\n',requirements.mtfRequirement);
    fprintf(fileId, '- SNR (Lambda3) >= %g\n',requirements.snrRequirement);
    fprintf(fileId, '- Dias de cobertura <= %g\n',requirements.revisitRequirement_days);
    fprintf(fileId, '- FOV actual <= Limite FOV del telescopio\n');
    fprintf(fileId, '- Para telescopio Refractivo (si idx=1): Diametro <= 200mm\n\n');
    fprintf(fileId, 'Resultados (Dmin en mm) por detector y telescopio:\n\n');
    for caseIndex = 1:3*numel(telescopeNames)
        detectorIndex = ceil(caseIndex/numel(telescopeNames));
        telescopeIndex = mod(caseIndex-1,numel(telescopeNames))+1;
        if telescopeIndex==1, fprintf(fileId,'Detector %d:\n',detectorIndex); end
            summaryDiameters_mm = viableDiameters_mm(:, telescopeIndex, detectorIndex);
            validSummaryMask = ~isnan(summaryDiameters_mm);
            if any(validSummaryMask)
                minimumAltitude_km = min(altitudes_km(validSummaryMask));
                maximumAltitude_km = max(altitudes_km(validSummaryMask));
                minimumDiameter_mm = min(summaryDiameters_mm(validSummaryMask));
                maximumDiameter_mm = max(summaryDiameters_mm(validSummaryMask));
                fprintf(fileId, '  %s: Alturas validas [%.0f-%.0f km], Diametros resultantes [%.1f-%.1f mm]\n', ...
                        telescopeNames{telescopeIndex}, minimumAltitude_km, maximumAltitude_km, minimumDiameter_mm, maximumDiameter_mm);
            else
                fprintf(fileId, '  %s: No se encontraron combinaciones validas que cumplan todos los criterios.\n', telescopeNames{telescopeIndex});
            end
        if telescopeIndex==numel(telescopeNames), fprintf(fileId,'\n'); end
    end
    clear fileCleanup;
end

function labels = plotSelectionCurves(altitudes_km,diameters_mm,telescopeNames)
labels = {};
for telescopeIndex = 1:numel(telescopeNames)
    valid = ~isnan(diameters_mm(:,telescopeIndex));
    if any(valid)
        plot(altitudes_km(valid),diameters_mm(valid,telescopeIndex),'LineWidth',2);
        labels{end+1} = telescopeNames{telescopeIndex};
    end
end
end
