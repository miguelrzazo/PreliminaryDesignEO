function CrossDataFunction(GSD, alturas_orbitales, swaths_km, N_pix, diametros_pupila, telescope_names, fov_limit, configs)

% Crear directorio para resultados
output_dir = 'HvsDmin';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Generar nombres para configuraciones dinÃ¡micamente
config_names = cell(1, size(configs, 1));
for idx = 1:size(configs, 1)
    config_names{idx} = sprintf('%dsat_ph%d_%dtel', configs(idx,1), configs(idx,2), configs(idx,3));
end

% Para cada configuraciÃ³n
for config_idx = 1:size(configs, 1)
    % Extraer parÃ¡metros de la configuraciÃ³n
    N_sat = configs(config_idx, 1);
    phasing_type = configs(config_idx, 2);
    N_telescopes = configs(config_idx, 3);

    % Inicializar matriz para almacenar resultados
    % Filas: alturas, Columnas: telescopios, PÃ¡ginas: detectores
    Dmin_results = nan(length(alturas_orbitales), 4, 3);

    % Para cada detector
    for detector_idx = 1:3
        % Para cada telescopio
        for telescopio_idx = 1:4
            % Obtener paths de archivos
            mtf_file = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx);
            snr_file = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx+3, telescopio_idx);
            coverage_file = sprintf('Coverage/Sat%d_Detector%d_Telescopio%d_%d_resultados.csv', N_sat, detector_idx, telescopio_idx, N_telescopes);

            % Verificar si existen los archivos
            if ~exist(mtf_file, 'file') || ~exist(snr_file, 'file') || ~exist(coverage_file, 'file')
                warning('Faltan archivos para Detector %d, Telescopio %d, ConfiguraciÃ³n %s', detector_idx, telescopio_idx, config_names{config_idx});
                continue;
            end

            % Leer archivos CSV
            mtf_data = readtable(mtf_file, 'ReadRowNames', true);
            snr_data = readtable(snr_file, 'ReadRowNames', true);
            coverage_data = readtable(coverage_file, 'ReadRowNames', true);

            % Convertir a matrices para facilitar el procesamiento
            mtf_matrix = table2array(mtf_data);
            snr_matrix = table2array(snr_data);
            coverage_matrix = table2array(coverage_data);

            % Para cada altura, encontrar el diÃ¡metro mÃ­nimo que cumple las condiciones
            for alt_idx = 1:length(alturas_orbitales)
                altura = alturas_orbitales(alt_idx);
                row_name = sprintf('Alt_%d', altura);

                % Obtener Ã­ndice de fila para esta altura
                alt_row_idx = find(strcmp(mtf_data.Properties.RowNames, row_name));
                if isempty(alt_row_idx)
                    continue;
                end

                % Obtener valores de MTF, SNR y cobertura para esta altura
                mtf_row = mtf_matrix(alt_row_idx, :);
                snr_row = snr_matrix(alt_row_idx, :);
                coverage_row = coverage_matrix(alt_row_idx, :);

                % Encontrar swath Ã³ptimo para esta altura (el mayor que cumpla cobertura <= 7 dÃ­as)
                valid_swaths = find(coverage_row <= 7 & ~isnan(coverage_row));
                if isempty(valid_swaths)
                    continue; % No hay swaths vÃ¡lidos para esta altura
                end

                % Seleccionar el swath mÃ¡s pequeÃ±o que cumple con la cobertura
                selected_swath_idx = valid_swaths(1);
                selected_swath = swaths_km(selected_swath_idx);

               

                % Verificar si el FOV excede el lÃ­mite del telescopio
                fov_actual_deg = 2 * atand(selected_swath / (2 * altura));
                if fov_actual_deg > fov_limit(telescopio_idx) * N_telescopes
                    continue; % FOV excede el lÃ­mite
                end

                % Encontrar el diÃ¡metro mÃ­nimo que cumple MTF > 0.25 y SNR > 400
                valid_diameters = find(mtf_row >= 0.25 & snr_row >= 400);
                if isempty(valid_diameters)
                    continue; % No hay diÃ¡metros vÃ¡lidos
                end

                % Seleccionar el diÃ¡metro mÃ­nimo
                min_diam_idx = valid_diameters(1);
                
                % Verificar si es telescopio Refractivo y el diÃ¡metro es mayor a 100mm
                if telescopio_idx == 1 && diametros_pupila(min_diam_idx) > 90
                    continue; % DiÃ¡metro demasiado grande para telescopio Refractivo
                end
                
                Dmin_results(alt_idx, telescopio_idx, detector_idx) = diametros_pupila(min_diam_idx);
            end
        end
    end

    % Generar grÃ¡ficos para esta configuraciÃ³n
    figure('Position', [100, 100, 900, 700]);
    
    % Crear tÃ­tulo descriptivo segÃºn la configuraciÃ³n
    N_sat = configs(config_idx, 1);
    phasing_type = configs(config_idx, 2);
    N_telescopes = configs(config_idx, 3);

    % Determinar descripciÃ³n del tipo de desfase
    if phasing_type == 1
        phasing_desc = 'desfase en anomalÃ­a';
    elseif phasing_type == 2
        phasing_desc = 'desfase en RAAN';
    elseif phasing_type == 3
        phasing_desc = 'desfase en anomalÃ­a y RAAN';
    else
        phasing_desc = 'desfase desconocido';
    end

    % Crear tÃ­tulo base para la configuraciÃ³n
    config_title = sprintf('ConfiguraciÃ³n: %d satÃ©lite(s), %s, %d telescopio(s)', N_sat, phasing_desc, N_telescopes);
    
    for detector_idx = 1:3
        subplot(3, 1, detector_idx);
        hold on;
        for telescopio_idx = 1:4
            % Extraer datos para este detector y telescopio
            dmin_values = Dmin_results(:, telescopio_idx, detector_idx);
            % Filtrar valores NaN
            valid_indices = ~isnan(dmin_values);
            if sum(valid_indices) > 0
                % Graficar relaciÃ³n H vs Dmin
                plot(alturas_orbitales(valid_indices), dmin_values(valid_indices), 'LineWidth', 2);
            end
        end
        
        % Configurar grÃ¡fico con tÃ­tulo personalizado
        title(sprintf('Detector %d: Altura vs DiÃ¡metro MÃ­nimo - %s', detector_idx, config_title));
        xlabel('Altura (km)');
        ylabel('DiÃ¡metro pupila mÃ­nimo (mm)');
        legend(telescope_names, 'Location', 'northwest');
        grid on;
        hold off;
    end

    % Guardar figura
    saveas(gcf, sprintf('%s/HvsDmin_%s.fig', output_dir, config_names{config_idx}));
    saveas(gcf, sprintf('%s/HvsDmin_%s.png', output_dir, config_names{config_idx}));

    % Generar archivos CSV y TXT para esta configuraciÃ³n
    for detector_idx = 1:3
        % Crear tabla para exportar
        output_table = array2table(zeros(length(alturas_orbitales), 5), 'VariableNames', ...
            {'Altura_km', 'Dmin_Refractivo', 'Dmin_Korsch', 'Dmin_Cassegrain', 'Dmin_TMA'});
        output_table.Altura_km = alturas_orbitales';
        for telescopio_idx = 1:4
            output_table{:, telescopio_idx + 1} = Dmin_results(:, telescopio_idx, detector_idx);
        end
        % Exportar a CSV
        csv_filename = sprintf('%s/HvsDmin_%s_Detector%d.csv', output_dir, config_names{config_idx}, detector_idx);
        writetable(output_table, csv_filename);
        % Exportar a TXT
        txt_filename = sprintf('%s/HvsDmin_%s_Detector%d.txt', output_dir, config_names{config_idx}, detector_idx);
        fileID = fopen(txt_filename, 'w');
        fprintf(fileID, 'Altura (km)\tDmin Refractivo (mm)\tDmin Korsch (mm)\tDmin Cassegrain (mm)\tDmin TMA (mm)\n');
        for i = 1:length(alturas_orbitales)
            fprintf(fileID, '%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n', ...
                alturas_orbitales(i), ...
                Dmin_results(i, 1, detector_idx), ...
                Dmin_results(i, 2, detector_idx), ...
                Dmin_results(i, 3, detector_idx), ...
                Dmin_results(i, 4, detector_idx));
        end
        fclose(fileID);
    end

    % Generar un archivo de resumen para esta configuraciÃ³n
    summary_filename = sprintf('%s/HvsDmin_%s_resumen.txt', output_dir, config_names{config_idx});
    fileID = fopen(summary_filename, 'w');
    fprintf(fileID, '=== RESUMEN DE ANÃLISIS CRUZADO ===\n\n');
    fprintf(fileID, 'ConfiguraciÃ³n: %s\n', config_names{config_idx});
    fprintf(fileID, 'NÃºmero de satÃ©lites: %d\n', N_sat);
    fprintf(fileID, 'NÃºmero de telescopios por satÃ©lite: %d\n\n', N_telescopes);
    fprintf(fileID, 'Criterios aplicados:\n');
    fprintf(fileID, '- MTF (Lambda2) > 0.25\n');
    fprintf(fileID, '- SNR (Lambda3) > 400\n');
    fprintf(fileID, '- DÃ­as de cobertura <= 7\n');
    fprintf(fileID, '- Para telescopio Refractivo: DiÃ¡metro <= 100mm\n\n');
    fprintf(fileID, 'Resultados por detector y telescopio:\n\n');
    for detector_idx = 1:3
        fprintf(fileID, 'Detector %d:\n', detector_idx);
        for telescopio_idx = 1:4
            dmin_values = Dmin_results(:, telescopio_idx, detector_idx);
            valid_indices = ~isnan(dmin_values);
            if sum(valid_indices) > 0
                min_alt = min(alturas_orbitales(valid_indices));
                max_alt = max(alturas_orbitales(valid_indices));
                min_diam = min(dmin_values(valid_indices));
                max_diam = max(dmin_values(valid_indices));
                fprintf(fileID, ' %s: Alturas vÃ¡lidas %.0f-%.0f km, DiÃ¡metros %.0f-%.0f mm\n', ...
                    telescope_names{telescopio_idx}, min_alt, max_alt, min_diam, max_diam);
            else
                fprintf(fileID, ' %s: No hay combinaciones vÃ¡lidas\n', telescope_names{telescopio_idx});
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end

disp('AnÃ¡lisis cruzado completado. Resultados guardados en la carpeta HvsDmin/');

end