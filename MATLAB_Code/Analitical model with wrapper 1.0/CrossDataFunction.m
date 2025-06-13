function CrossDataFunction(alturas_orbitales, swaths_km, diametros_pupila, telescope_names, fov_limit, configs)
% CrossDataFunction: Analiza datos de MTF, SNR y cobertura para encontrar el diámetro
% mínimo de pupila (Dmin) para diversas configuraciones de satélites,
% telescopios, detectores y alturas orbitales.
%
% Inputs:
%   alturas_orbitales - Vector con las alturas orbitales a analizar (km).
%   swaths_km         - Vector con los anchos de barrido (swaths) a considerar (km).
%   diametros_pupila  - Vector con los diámetros de pupila a evaluar (mm).
%   telescope_names   - Cell array con los nombres de los telescopios.
%   fov_limit         - Vector con los límites de campo de visión (FOV) para cada telescopio (grados).
%   configs           - Matriz donde cada fila define una configuración [N_sat, N_telescopes].

% Crear directorio para resultados si no existe
output_dir = 'HvsDmin';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Generar nombres para configuraciones dinámicamente para los nombres de archivo
config_names = cell(1, size(configs, 1));
for idx = 1:size(configs, 1)
    config_names{idx} = sprintf('%dsat_%dtel', configs(idx,1), configs(idx,2));
end

% Procesar cada configuración (combinación de número de satélites y telescopios por satélite)
for config_idx = 1:size(configs, 1)
    N_sat = configs(config_idx, 1);        % Número de satélites en la constelación
    N_telescopes = configs(config_idx, 2); % Número de telescopios por satélite

    % Inicializar matriz para almacenar los resultados de Dmin
    % Dimensiones: (alturas, telescopios, detectores)
    Dmin_results = nan(length(alturas_orbitales), length(telescope_names), 3); % Asumiendo 3 detectores y hasta 4 telescopios (o length(telescope_names))

    % Iterar sobre cada tipo de detector (asumiendo 3 detectores)
    for detector_idx = 1:3
        % Iterar sobre cada tipo de telescopio
        for telescopio_idx = 1:length(telescope_names)
            % Construir nombres de archivo para MTF y SNR
            % Nota: detector_idx+3 para SNR_Lambda3 es específico de la estructura de nombres de archivo del usuario
            mtf_file = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx);
            snr_file = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx + 3, telescopio_idx);

            % Limpiar nombre del telescopio para usar en nombre de archivo de cobertura
            telescope_name_clean = strrep(telescope_names{telescopio_idx}, ' ', '');
            config_name_coverage = sprintf('%dSat_%dTel_%s_Det%d', N_sat, N_telescopes, telescope_name_clean, detector_idx);
            coverage_file = sprintf('coverage/coverage_%s.csv', config_name_coverage);

            % Verificar existencia de los archivos de datos
            if ~exist(mtf_file, 'file') || ~exist(snr_file, 'file') || ~exist(coverage_file, 'file')
                fprintf('Advertencia: Faltan archivos para Detector %d, Telescopio %d (%s), Configuracion %s.\n', ...
                        detector_idx, telescopio_idx, telescope_names{telescopio_idx}, config_names{config_idx});
                fprintf('  MTF file: %s (Existe: %d)\n', mtf_file, exist(mtf_file, 'file'));
                fprintf('  SNR file: %s (Existe: %d)\n', snr_file, exist(snr_file, 'file'));
                fprintf('  Coverage file: %s (Existe: %d)\n', coverage_file, exist(coverage_file, 'file'));
                continue; % Saltar a la siguiente iteración de telescopio
            end

            % Leer archivos CSV
            try
                mtf_data = readtable(mtf_file, 'ReadRowNames', true);
                snr_data = readtable(snr_file, 'ReadRowNames', true);
                coverage_data = readmatrix(coverage_file); % Asume que coverage_data no tiene encabezados de fila/columna o se manejan adecuadamente

                % Convertir tablas a matrices para facilitar el procesamiento
                mtf_matrix = table2array(mtf_data);
                snr_matrix = table2array(snr_data);

                % Iterar sobre cada altura orbital
                for alt_idx = 1:length(alturas_orbitales)
                    altura = alturas_orbitales(alt_idx);

                    % Validar que alt_idx no exceda las dimensiones de las matrices leídas
                    if alt_idx > size(mtf_matrix, 1) || alt_idx > size(snr_matrix, 1) || alt_idx > size(coverage_data, 1)
                        fprintf('Advertencia: alt_idx %d excede las dimensiones de los datos para H=%.1f km, Det %d, Tel %d.\n', ...
                                alt_idx, altura, detector_idx, telescopio_idx);
                        continue;
                    end

                    % Obtener filas de MTF, SNR y cobertura para la altura actual
                    mtf_row = mtf_matrix(alt_idx, :);
                    snr_row = snr_matrix(alt_idx, :);
                    coverage_row = coverage_data(alt_idx, :); % Asume que las columnas corresponden a swaths_km

                    % 1. Encontrar swath óptimo: el mayor swath que cumpla cobertura <= 7 días
                    valid_swath_indices = find(coverage_row <= 7 & ~isnan(coverage_row));
                    if isempty(valid_swath_indices)
                        continue; % No hay swaths válidos para esta altura y configuración de cobertura
                    end
                    % Seleccionar el índice del swath más grande (asumiendo que swaths_km y las columnas de coverage_data están ordenados ascendentemente)
                    selected_swath_col_idx = valid_swath_indices(end); 
                    selected_swath_km = swaths_km(selected_swath_col_idx);

                    % 2. Verificar si el FOV del swath seleccionado excede el límite del telescopio
                    % FOV (full angle) = 2 * atan( (Swath/2) / Altura )
                    fov_actual_deg = 2 * atand(selected_swath_km / (2 * altura));
                    if fov_actual_deg > fov_limit(telescopio_idx)
                        continue; % FOV excede el límite para este telescopio y swath
                    end

                    % 3. Encontrar el diámetro mínimo que cumple MTF >= 0.25 y SNR >= 400 para el swath seleccionado
                    % Se asume que las columnas de mtf_row y snr_row corresponden a diametros_pupila
                    valid_diameter_indices = find(mtf_row >= 0.25 & snr_row >= 400);
                    if isempty(valid_diameter_indices)
                        continue; % No hay diámetros válidos que cumplan MTF y SNR para esta altura
                    end
                    % Seleccionar el índice del diámetro mínimo (asumiendo diametros_pupila y las columnas de mtf/snr están ordenados ascendentemente)
                    min_diam_col_idx = valid_diameter_indices(1);
                    
                    % Asegurar que el índice del diámetro no exceda las dimensiones de diametros_pupila
                    if min_diam_col_idx > length(diametros_pupila)
                        fprintf('Advertencia: Índice de diámetro %d excede la longitud de diametros_pupila (%d).\n', min_diam_col_idx, length(diametros_pupila));
                        continue;
                    end
                    
                    selected_dmin_mm = diametros_pupila(min_diam_col_idx);

                    % 4. Condición específica para telescopio Refractivo (asumiendo telescopio_idx == 1)
                    % Si es Refractivo y el diámetro es > 90mm, no es una solución válida.
                    if telescopio_idx == 1 && selected_dmin_mm > 90
                        continue; 
                    end

                    % Almacenar el Dmin encontrado
                    Dmin_results(alt_idx, telescopio_idx, detector_idx) = selected_dmin_mm;
                end
            catch ME
                fprintf('Error procesando archivos para Detector %d, Telescopio %d (%s), Config %s: %s\n', ...
                        detector_idx, telescopio_idx, telescope_names{telescopio_idx}, config_names{config_idx}, ME.message);
                continue; % Saltar al siguiente telescopio en caso de error de lectura/procesamiento
            end
        end % Fin loop telescopios
    end % Fin loop detectores

    % --- Generación de Gráficos para la configuración actual ---
    main_fig = figure('Position', [100, 100, 900, 700], 'Visible', 'off'); % Crear invisible para no interrumpir
    config_title_tex = sprintf('Configuracion: %d satelite(s), %d telescopio(s) por satelite', N_sat, N_telescopes);

    for detector_idx = 1:3
        subplot(3, 1, detector_idx);
        hold on;
        legend_entries = {};
        for telescopio_idx = 1:length(telescope_names)
            dmin_values = Dmin_results(:, telescopio_idx, detector_idx);
            valid_indices = ~isnan(dmin_values);
            if any(valid_indices) % Graficar solo si hay datos válidos
                plot(alturas_orbitales(valid_indices), dmin_values(valid_indices), 'LineWidth', 2);
                legend_entries{end+1} = telescope_names{telescopio_idx};
            end
        end
        hold off;

        title_str = sprintf('Detector %d: Altura orbital vs $D_{min}$ - %s', detector_idx, config_title_tex);
        title(title_str, 'Interpreter', 'latex', 'FontSize', 10);
        xlabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 9);
        ylabel('$D_{min}$ (mm)', 'Interpreter', 'latex', 'FontSize', 9);
        if ~isempty(legend_entries)
            legend(legend_entries, 'Location', 'northwest', 'Interpreter', 'none'); % Usar 'none' si los nombres tienen caracteres LaTeX especiales
        end
        grid on;
        axis tight; % Ajustar ejes a los datos
    end
    
    % Guardar la figura principal (todos los detectores, todos los telescopios)
    main_plot_filename = sprintf('%s/HvsDmin_%s.png', output_dir, config_names{config_idx});
    saveas(main_fig, main_plot_filename);
    close(main_fig);

    % --- Gráfica Adicional: H vs Dmin para Telescopio Refractivo (Telescopio 1) ---
    % Asumiendo que el primer telescopio (idx=1) es el 'Refractivo'
    if length(telescope_names) >= 1 % Solo si hay al menos un telescopio
        refractivo_fig = figure('Position', [100, 100, 700, 500], 'Visible', 'off');
        hold on;
        colors = {'r', 'g', 'b'};
        line_styles = {'--', ':', '-.'};
        legend_labels_refractivo = {};

        for detector_idx = 1:3
            dmin_refractivo = Dmin_results(:, 1, detector_idx); % Columna 1 para el telescopio refractivo
            valid_indices_refractivo = ~isnan(dmin_refractivo);
            if any(valid_indices_refractivo)
                plot(alturas_orbitales(valid_indices_refractivo), dmin_refractivo(valid_indices_refractivo), ...
                     [colors{detector_idx}, line_styles{detector_idx}], 'LineWidth', 2);
                legend_labels_refractivo{end+1} = sprintf('Detector %d', detector_idx);
            end
        end
        hold off;

        if ~isempty(legend_labels_refractivo) % Solo configurar y guardar si hay datos
            title_str_refractivo = sprintf('$h$ vs $D_{min}$ - %s - Config: %s', telescope_names{1}, config_title_tex);
            title(title_str_refractivo, 'Interpreter', 'latex', 'FontSize', 10);
            xlabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 9);
            ylabel(sprintf('$D_{min}$ %s (mm)', telescope_names{1}), 'Interpreter', 'latex', 'FontSize', 9);
            legend(legend_labels_refractivo, 'Location', 'best', 'Interpreter', 'latex');
            grid on;
            axis tight;
            
            refractivo_plot_filename = sprintf('%s/HvsDmin_%s_Config_%s.png', output_dir, strrep(telescope_names{1},' ','_'), config_names{config_idx});
            saveas(refractivo_fig, refractivo_plot_filename);
        end
        close(refractivo_fig);
    end

    % --- Exportar resultados a CSV y TXT para la configuración actual ---
    for detector_idx = 1:3
        % Crear nombres de variable para la tabla basados en telescope_names
        variable_names_table = [{'Altura_km'}, cellfun(@(x) sprintf('Dmin_%s_mm', strrep(x,' ','_')), telescope_names, 'UniformOutput', false)];
        
        output_data_matrix = [alturas_orbitales', Dmin_results(:, :, detector_idx)];
        output_table = array2table(output_data_matrix, 'VariableNames', variable_names_table);

        % Exportar a CSV
        csv_filename = sprintf('HvsDmin/HvsDmin_%s_Detector%d.csv', config_names{config_idx}, detector_idx);
        writetable(output_table, csv_filename);

        % Exportar a TXT
        txt_filename = sprintf('%s/HvsDmin_%s_Detector%d.txt', output_dir, config_names{config_idx}, detector_idx);
        fileID = fopen(txt_filename, 'w');
        % Escribir cabecera del TXT
        fprintf(fileID, 'Altura_km');
        for tel_name_idx = 1:length(telescope_names)
            fprintf(fileID, '\tDmin_%s_mm', strrep(telescope_names{tel_name_idx},' ','_'));
        end
        fprintf(fileID, '\n');
        
        % Escribir datos al TXT
        for row_idx = 1:size(output_data_matrix, 1)
            fprintf(fileID, '%.1f', output_data_matrix(row_idx, 1)); % Altura
            for col_idx = 2:size(output_data_matrix, 2) % Dmin para cada telescopio
                if isnan(output_data_matrix(row_idx, col_idx))
                    fprintf(fileID, '\tNaN');
                else
                    fprintf(fileID, '\t%.1f', output_data_matrix(row_idx, col_idx));
                end
            end
            fprintf(fileID, '\n');
        end
        fclose(fileID);
    end

    % --- Generar archivo de resumen para la configuración actual ---
    summary_filename = sprintf('%s/RESUMEN_Config_%s.txt', output_dir, config_names{config_idx});
    fileID = fopen(summary_filename, 'w');
    fprintf(fileID, '=== RESUMEN DE ANALISIS CRUZADO ===\n\n');
    fprintf(fileID, 'Configuracion: %s (%d satelite(s), %d telescopio(s) por satelite)\n\n', ...
            config_names{config_idx}, N_sat, N_telescopes);

    fprintf(fileID, 'Criterios aplicados para Dmin:\n');
    fprintf(fileID, '- MTF (Lambda2) >= 0.25\n');
    fprintf(fileID, '- SNR (Lambda3) >= 400\n'); % Asumiendo que SNR > 400 es el criterio
    fprintf(fileID, '- Dias de cobertura <= 7\n');
    fprintf(fileID, '- FOV actual <= Limite FOV del telescopio\n');
    fprintf(fileID, '- Para telescopio Refractivo (si idx=1): Diametro <= 90mm\n\n');

    fprintf(fileID, 'Resultados (Dmin en mm) por detector y telescopio:\n\n');
    for detector_idx = 1:3
        fprintf(fileID, 'Detector %d:\n', detector_idx);
        for telescopio_idx = 1:length(telescope_names)
            dmin_values_for_summary = Dmin_results(:, telescopio_idx, detector_idx);
            valid_summary_indices = ~isnan(dmin_values_for_summary);
            if any(valid_summary_indices)
                min_alt_valid = min(alturas_orbitales(valid_summary_indices));
                max_alt_valid = max(alturas_orbitales(valid_summary_indices));
                min_diam_valid = min(dmin_values_for_summary(valid_summary_indices));
                max_diam_valid = max(dmin_values_for_summary(valid_summary_indices));
                fprintf(fileID, '  %s: Alturas validas [%.0f-%.0f km], Diametros resultantes [%.1f-%.1f mm]\n', ...
                        telescope_names{telescopio_idx}, min_alt_valid, max_alt_valid, min_diam_valid, max_diam_valid);
            else
                fprintf(fileID, '  %s: No se encontraron combinaciones validas que cumplan todos los criterios.\n', telescope_names{telescopio_idx});
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);

end % Fin loop configuraciones

disp('Analisis cruzado completado. Resultados guardados en la carpeta HvsDmin/');

end % Fin de la función
