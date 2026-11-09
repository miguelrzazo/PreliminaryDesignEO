function writeMassPlots(resultsTable)
%WRITEMASSPLOTS Plot supplied mass rows; no physical-model evaluations.
if isempty(resultsTable), return; end
dryMassDirectory = 'Masa_seca';
totalMassDirectory = 'Masa_total';
if ~exist(dryMassDirectory, 'dir'), mkdir(dryMassDirectory); end
if ~exist(totalMassDirectory, 'dir'), mkdir(totalMassDirectory); end
fprintf('Generando gráficas de masa para las soluciones encontradas...\n');
configurations = unique(resultsTable(:, {'Num_Satelites', 'Num_Telescopios', 'ID_Detector', 'Tipo_Telescopio'}), 'rows');
for configurationIndex = 1:height(configurations)
    configuration = configurations(configurationIndex, :);
    satelliteCount = configuration.Num_Satelites;
    telescopeCount = configuration.Num_Telescopios;
    detectorIndex = configuration.ID_Detector;
    telescopeNameCell = configuration.Tipo_Telescopio;
    telescopeName = telescopeNameCell{1};
    configurationMask = resultsTable.Num_Satelites == satelliteCount & ...
                  resultsTable.Num_Telescopios == telescopeCount & ...
                  resultsTable.ID_Detector == detectorIndex & ...
                  strcmp(resultsTable.Tipo_Telescopio, telescopeName);
    plotData = resultsTable(configurationMask, :);
    plotData = sortrows(plotData, 'Altura_km');
    if isempty(plotData), continue; end
    dryMassFigure = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
    dryMassFigureCleanup = onCleanup(@() closeStudyFigure(dryMassFigure));
    hold on;
    plot(plotData.Altura_km, plotData.Masa_Seca_TM_kg, '--o', 'LineWidth', 1.5, 'DisplayName', 'Estimacion (Thematic Mapper)');
    plot(plotData.Altura_km, plotData.Masa_Seca_SEOSAT_kg, '--s', 'LineWidth', 1.5, 'DisplayName', 'Estimacion (SEOSAT)');
    plot(plotData.Altura_km, plotData.Masa_Seca_Satelite_kg, '-d', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', sprintf('Propuesta (%d Tel)', telescopeCount));
    if telescopeCount >= 2
        singleTelescopeDryMass_kg = plotData.Masa_Seca_Satelite_kg / 1.5;
        plot(plotData.Altura_km, singleTelescopeDryMass_kg, '--p', 'LineWidth', 1.5, 'DisplayName', 'Promedio (1 Telescopio)');
    end
    hold off;
    titleText = sprintf('Comparativa de Masa Seca vs. Altura\nConfig: %d Sat, %d Tel, Det %d, %s', satelliteCount, telescopeCount, detectorIndex, telescopeName);
    title(titleText, 'Interpreter', 'latex');
    xlabel('Altura Orbital [km]', 'Interpreter', 'latex');
    ylabel('Masa Seca [kg]', 'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex');
    grid on;
    axis tight;
    dryMassFigure.Color = 'white';
    set(findall(dryMassFigure, 'Type', 'Axes'), 'Color', 'white','XColor','k','YColor','k');
    set(findall(dryMassFigure, 'Type', 'Text'), 'Color', 'k');
    set(findall(dryMassFigure, 'Type', 'Legend'), 'Color', 'white','TextColor','k','EdgeColor','k');
    dryMassFilename = sprintf('MasaSeca_Comparativa_%dsat_%dtel_Det%d_%s.jpg', satelliteCount, telescopeCount, detectorIndex, strrep(telescopeName, ' ', '_'));
    saveas(dryMassFigure, fullfile(dryMassDirectory, dryMassFilename));
    exportgraphics(dryMassFigure, fullfile(dryMassDirectory, strrep(dryMassFilename, '.jpg', '.pdf')), 'ContentType', 'vector');
    close(dryMassFigure);
    totalMassFigure = figure('Visible', 'off', 'Position', [100, 100, 900, 600]);
    totalMassFigureCleanup = onCleanup(@() closeStudyFigure(totalMassFigure));
    hold on;
    if ismember('Es_Viable', plotData.Properties.VariableNames)
        viableMask = plotData.Es_Viable == true;
    else
        viableMask = true(size(plotData,1), 1);
    end
    theoreticalMask = ~viableMask;
    if any(theoreticalMask)
        % Smooth display copies only. Selection, exported data and viable markers
        % retain the original discrete model values.
        displayDryMass = exp(smoothdata(log(plotData.Masa_Seca_Satelite_kg),'lowess',9));
        displayFuelMass = exp(smoothdata(log(max(plotData.Masa_Combustible_kg,1e-8)),'lowess',9));
        displayFuelMass(plotData.Masa_Combustible_kg==0) = 0;
        displayTotalMass = displayDryMass + displayFuelMass;
        plot(plotData.Altura_km(theoreticalMask), displayTotalMass(theoreticalMask), ...
            'Color', [0.55 0.55 0.55], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Masa total (teorica)');
        plot(plotData.Altura_km(theoreticalMask), displayDryMass(theoreticalMask), ...
            'Color', [0.65 0.65 0.80], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Masa seca (teorica)');
        plot(plotData.Altura_km(theoreticalMask), displayFuelMass(theoreticalMask), ...
            'Color', [0.60 0.75 0.60], 'LineStyle', ':', 'LineWidth', 1.2, 'DisplayName', 'Masa combustible (teorica)');
    end
    if any(viableMask)
        plot(plotData.Altura_km(viableMask), plotData.Masa_Total_Satelite_kg(viableMask), ...
            'r-d', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'DisplayName', 'Masa total (puntos viables)');
        plot(plotData.Altura_km(viableMask), plotData.Masa_Seca_Satelite_kg(viableMask), ...
            'b--s', 'LineWidth', 1.5, 'DisplayName', 'Masa seca (viable)');
        plot(plotData.Altura_km(viableMask), plotData.Masa_Combustible_kg(viableMask), ...
            'g:^', 'LineWidth', 1.5, 'DisplayName', 'Masa combustible (viable)');
        [minimumMass_kg, minimumIndex] = min(plotData.Masa_Total_Satelite_kg(viableMask));
        viableAltitudes_km = plotData.Altura_km(viableMask);
        minimumAltitude_km = viableAltitudes_km(minimumIndex);
        plot(minimumAltitude_km, minimumMass_kg, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
        annotationText = sprintf('Minimo viable: %.3f kg @ %d km', minimumMass_kg, minimumAltitude_km);
        text(minimumAltitude_km + 25, minimumMass_kg, annotationText, 'VerticalAlignment', 'middle', ...
             'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', 'k','Interpreter','latex');
    end
    hold off;
    titleText = sprintf('Desglose de Masas vs. Altura\nConfig: %d Sat, %d Tel, Det %d, %s', satelliteCount, telescopeCount, detectorIndex, telescopeName);
    title(titleText, 'Interpreter', 'latex');
    xlabel('Altura Orbital [km]', 'Interpreter', 'latex');
    ylabel('Masa [kg]', 'Interpreter', 'latex');
    legend('Location', 'northwest', 'Interpreter', 'latex');
    grid on;
    axis tight;
    if any(viableMask)
        viableAltitudes_km = plotData.Altura_km(viableMask);
        xlim([max(400, min(plotData.Altura_km)), max(viableAltitudes_km) + 10]);
        visibleRows = plotData.Altura_km >= max(400,min(plotData.Altura_km));
        ylim([0,1.1*max(plotData.Masa_Total_Satelite_kg(visibleRows))]);
    end
    totalMassPrefix = sprintf('MasaDesglose_%dsat_%dtel_Det%d_%s', satelliteCount, telescopeCount, detectorIndex, strrep(telescopeName, ' ', '_'));
    totalMassFigure.Color = 'white';
    set(findall(totalMassFigure, 'Type', 'Axes'), 'Color', 'white','XColor','k','YColor','k');
    set(findall(totalMassFigure, 'Type', 'Text'), 'Color', 'k');
    set(findall(totalMassFigure, 'Type', 'Legend'), 'Color', 'white','TextColor','k','EdgeColor','k');
    exportgraphics(totalMassFigure, fullfile(totalMassDirectory, [totalMassPrefix '.pdf']), 'ContentType', 'vector');
    close(totalMassFigure);
end
fprintf('Gráficas de masa guardadas en las carpetas ''%s'' y ''%s''.\n', dryMassDirectory, totalMassDirectory);
end
