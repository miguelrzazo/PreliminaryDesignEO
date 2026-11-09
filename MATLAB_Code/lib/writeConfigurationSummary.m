function writeConfigurationSummary(summary)
if isempty(summary.rows), return; end
%WRITECONFIGURATIONSUMMARY Export precomputed rows and policy results.
outputDirectory = 'OptimumConfigs';
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
resultRows = summary.rows; bestDesign = summary.best; headers = configurationHeaders();
if ~isempty(summary.maintenance)
    writetable(summary.maintenance,fullfile(outputDirectory,'maintenance_policy_comparison.csv'));
end
resultsTable = cell2table(resultRows, 'VariableNames', headers);
resultsTable = sortrows(resultsTable, 'Masa_Total_Constelacion_kg', 'ascend');
csvFilename = fullfile(outputDirectory, 'reporte_completo_configuraciones.csv');
writetable(resultsTable, csvFilename);
fprintf('Reporte CSV completo guardado en: %s\n', csvFilename);
textFilename = fullfile(outputDirectory, 'reporte_resumen_configuraciones.txt');
fileId = fopen(textFilename, 'w');
assert(fileId>=0,'eo:OutputFile','Cannot open the requested text output.');
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '==========================================================\n');
fprintf(fileId, ' INFORME DE ANÁLISIS DE CONFIGURACIONES\n');
fprintf(fileId, '==========================================================\n\n');
if ~isinf(bestDesign.mass)
    fprintf(fileId, '--- ÓPTIMO GLOBAL (Masa Mínima de Constelación) ---\n');
    data = bestDesign.data;
    fprintf(fileId, 'Masa Total de la Constelación: %.2f kg\n\n', bestDesign.mass);
    fprintf(fileId, ' > Configuración: %d Satélite(s), %d Telescopio(s)/sat, %d Detector(es) Necesario(s), ID Detector %d con %s\n', data{2}, data{3}, data{4}, data{5}, data{6});
    fprintf(fileId, ' > Órbita: %d km, %.2f° inc, Desfase Anom. %.1f°\n', data{7}, data{14}, data{16});
    fprintf(fileId, ' > Rendimiento: Swath Requerido=%.1f km (Real=%.1f km, Excedente=%.1f%%) | Cobertura=%.1f días\n', data{9}, data{11}, data{13}, data{10});
    fprintf(fileId, ' > Parámetros Ópticos: D=%.1f mm, F=%.2f m (f/%.1f) | MTF@Nyquist=%.3f, SNR=%.1f\n', data{8}, data{19}, data{20}, data{17}, data{18});
    fprintf(fileId, ' > Prop. Físicas (Instrumento): Masa=%.2f kg, Long=%.2f m, Área=%.2f m^2, Vol=%.6f m^3, Potencia=%.1f W\n', data{21}, data{22}, data{23}, data{24}, data{25});
    fprintf(fileId, ' > Prop. Físicas (Satélite): Long=%.2f m, Área=%.2f m^2, Vol=%.6f m^3, Estándar U=%.2f U\n', data{34}, data{35}, data{36}, data{37});
    fprintf(fileId, ' > Masa Satélite: Total=%.2f kg (Seca=%.2f kg, Combustible=%.2f kg)\n', data{28}, data{26}, data{27});
    fprintf(fileId, ' > Mantenimiento: %d maniobras, Delta-V Total %.1f m/s\n', data{30}, data{31});
end
fprintf(fileId, '\n==========================================================\n');
fprintf(fileId, ' LISTADO DE TODAS LAS CONFIGURACIONES VÁLIDAS\n');
fprintf(fileId, ' (Ordenadas por masa total de la constelación)\n');
fprintf(fileId, '==========================================================\n\n');
viableCount = 0;
for rowIndex = 1:height(resultsTable)
    row = table2cell(resultsTable(rowIndex,:));
    if ~row{1}, continue; end
    viableCount = viableCount + 1;
    fprintf(fileId, '\n--- Solución #%d ---\n', viableCount);
    fprintf(fileId, 'Masa Constelación: %.2f kg | Masa Satélite Total: %.2f kg (Seca: %.2f kg, Comb: %.2f kg)\n', row{29}, row{28}, row{26}, row{27});
    fprintf(fileId, ' > Configuración: %d Sat, %d Tel, %d Dets | ID Detector: %d | Telescopio: %s\n', row{2}, row{3}, row{4}, row{5}, row{6});
    fprintf(fileId, ' > Órbita: %.0f km alt, %.2f° inc (SSO), LTAN %d h, Desfase Anom. %.1f°\n', row{7}, row{14}, row{15}, row{16});
    fprintf(fileId, ' > Rendimiento: Swath Requerido=%.1f km (Real=%.1f km, Excedente=%.1f%%) | Cobertura=%.1f días\n', row{9}, row{11}, row{13}, row{10});
    fprintf(fileId, ' > Parámetros Ópticos: D=%.1f mm, F=%.2f m (f/%.1f) | MTF@Nyquist=%.3f, SNR=%.1f\n', row{8}, row{19}, row{20}, row{17}, row{18});
    fprintf(fileId, ' > Prop. Físicas (Satélite): Long=%.2f m, Área=%.2f m^2, Vol=%.6f m^3, Estándar U=%.2f U\n', row{34}, row{35}, row{36}, row{37});
    fprintf(fileId, '\n');
end
clear fileCleanup;
fprintf('Reporte TXT detallado guardado en: %s\n', textFilename);
fprintf('\nIniciando la generación de gráficas de desglose de masa...\n');
GenerateMassPlots(resultsTable);
fprintf('\nAnálisis finalizado.\n');
end
