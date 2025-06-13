function CrossDataFunction(GSD, alturas_orbitales, swaths_km, N_pix, diametros_pupila, telescope_names, fov_limit, configs)

% Crear directorio para resultados
output_dir = 'HvsDmin';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Generar nombres para configuraciones dinámicamente
config_names = cell(1, size(configs, 1));
for idx = 1:size(configs, 1)
    config_names{idx} = sprintf('%dsat_%dtel', configs(idx,1), configs(idx,2));
end

% Para cada configuración
for config_idx = 1:size(configs, 1)
    % Extraer parámetros de la configuración
    N_sat = configs(config_idx, 1);
    N_telescopes = configs(config_idx, 2);
    
    % Inicializar matriz para almacenar resultados
    % Filas: alturas, Columnas: telescopios, Páginas: detectores
    Dmin_results = nan(length(alturas_orbitales), 4, 3);
    
    % Para cada detector (1-3)
    for detector_idx = 1:3
        % Para cada telescopio (1-4)
        for telescopio_idx = 1:4
            % Obtener paths de archivos corregidos
            mtf_file = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx);
            snr_file = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx+3, telescopio_idx);
            
            % Generar nombre de archivo de coverage correcto
            telescope_name_clean = strrep(telescope_names{telescopio_idx}, ' ', '');
            config_name_coverage = sprintf('%dSat_%dTel_%s_Det%d', N_sat, N_telescopes, telescope_name_clean, detector_idx);
            coverage_file = sprintf('coverage/coverage_%s.csv', config_name_coverage);
            
            % Verificar si existen los archivos
            if ~exist(mtf_file, 'file') || ~exist(snr_file, 'file') || ~exist(coverage_file, 'file')
                fprintf('Faltan archivos para Detector %d, Telescopio %d, Configuracion %s\n', detector_idx, telescopio_idx, config_names{config_idx});
                fprintf(' MTF file: %s (exists: %d)\n', mtf_file, exist(mtf_file, 'file'));
                fprintf(' SNR file: %s (exists: %d)\n', snr_file, exist(snr_file, 'file'));
                fprintf(' Coverage file: %s (exists: %d)\n', coverage_file, exist(coverage_file, 'file'));
                continue;
            end
            
            % Leer archivos CSV
            try
                mtf_data = readtable(mtf_file, 'ReadRowNames', true);
                snr_data = readtable(snr_file, 'ReadRowNames', true);
                coverage_data = readmatrix(coverage_file);
                
                % Convertir a matrices para facilitar el procesamiento
                mtf_matrix = table2array(mtf_data);
                snr_matrix = table2array(snr_data);
                
                % Para cada altura, encontrar el diámetro mínimo que cumple las condiciones
                for alt_idx = 1:length(alturas_orbitales)
                    altura = alturas_orbitales(alt_idx);
                    
                    % Verificar que el índice no exceda las dimensiones de las matrices
                    if alt_idx > size(mtf_matrix, 1) || alt_idx > size(snr_matrix, 1) || alt_idx > size(coverage_data, 1)
                        continue;
                    end
                    
                    % Obtener valores de MTF, SNR y cobertura para esta altura
                    mtf_row = mtf_matrix(alt_idx, :);
                    snr_row = snr_matrix(alt_idx, :);
                    coverage_row = coverage_data(alt_idx, :);
                    
                    % Encontrar swath óptimo para esta altura (el mayor que cumpla cobertura <= 7 días)
                    valid_swaths = find(coverage_row <= 7 & ~isnan(coverage_row));
                    if isempty(valid_swaths)
                        continue; % No hay swaths válidos para esta altura
                    end
                    
                    % Seleccionar el swath más grande que cumple con la cobertura
                    selected_swath_idx = valid_swaths(end);
                    selected_swath = swaths_km(selected_swath_idx);
                    
                    % Verificar si el FOV excede el límite del telescopio
                    fov_actual_deg = 2 * atand(selected_swath / (2 * altura));
                    if fov_actual_deg > fov_limit(telescopio_idx)
                        continue; % FOV excede el límite
                    end
                    
                    % Encontrar el diámetro mínimo que cumple MTF > 0.25 y SNR > 400
                    valid_diameters = find(mtf_row >= 0.25 & snr_row >= 400);
                    if isempty(valid_diameters)
                        continue; % No hay diámetros válidos
                    end
                    
                    % Seleccionar el diámetro mínimo
                    min_diam_idx = valid_diameters(1);
                    
                    % Verificar que el índice no exceda las dimensiones
                    if min_diam_idx > length(diametros_pupila)
                        continue;
                    end
                    
                    % Verificar si es telescopio Refractivo y el diámetro es mayor a 90mm
                    if telescopio_idx == 1 && diametros_pupila(min_diam_idx) > 90
                        continue; % Diámetro demasiado grande para telescopio Refractivo
                    end
                    
                    Dmin_results(alt_idx, telescopio_idx, detector_idx) = diametros_pupila(min_diam_idx);
                end
                
            catch ME
                fprintf('Error leyendo archivos para Detector %d, Telescopio %d: %s\n', detector_idx, telescopio_idx, ME.message);
                continue;
            end
        end
    end
    
    % Generar gráficos para esta configuración
    figure('Position', [100, 100, 900, 700]);
    
    % Crear título descriptivo según la configuración (sin tildes)
    config_title = sprintf('Configuracion: %d satelite(s), %d telescopio(s)', N_sat, N_telescopes);
    
    for detector_idx = 1:3
        subplot(3, 1, detector_idx);
        hold on;
        
        for telescopio_idx = 1:4
            % Extraer datos para este detector y telescopio
            dmin_values = Dmin_results(:, telescopio_idx, detector_idx);
            
            % Filtrar valores NaN
            valid_indices = ~isnan(dmin_values);
            if sum(valid_indices) > 0
                % Graficar relación H vs Dmin
                plot(alturas_orbitales(valid_indices), dmin_values(valid_indices), 'LineWidth', 2);
            end
        end
        
        % Configurar gráfico con título personalizado usando LaTeX (sin tildes)
        title(sprintf('Detector %d: Altura vs Diametro Minimo - %s', detector_idx, config_title), ...
              'Interpreter', 'latex', 'FontSize', 12);
        xlabel('Altura (km)', 'Interpreter', 'latex', 'FontSize', 11);
        ylabel('Diametro pupila minimo (mm)', 'Interpreter', 'latex', 'FontSize', 11);
        legend(telescope_names, 'Location', 'northwest', 'Interpreter', 'latex');
        grid on;
        hold off;
    end
    
    % Guardar solo PNG (sin .fig)
    saveas(gcf, sprintf('%s/HvsDmin_%s.png', output_dir, config_names{config_idx}));
    
% === GRÁFICAS ADICIONALES: H vs Dmin para Telescopio 1 (Refractivo) ===

fig = figure('Visible', false);

% Definir colores y estilos de línea para cada detector
colors = {'r', 'g', 'b'}; % Rojo, Verde, Azul
line_styles = {'--', ':', '-.'};  % Diferentes tipos de raya

hold on;
legend_labels = {};

for detector_idx = 1:3
    dmin_refractivo = Dmin_results(:, 1, detector_idx); % Telescopio 1 = columna 1
    valid_indices = ~isnan(dmin_refractivo);
    
    if sum(valid_indices) > 0
        plot(alturas_orbitales(valid_indices), dmin_refractivo(valid_indices), ...
             [colors{detector_idx} line_styles{detector_idx}], ...
             'LineWidth', 2);
        legend_labels{end+1} = sprintf('Detector %d', detector_idx);
    end
end

% Configurar la gráfica
grid on;
xlabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Dmin Refractivo (mm)', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('H vs Dmin - Telescopio Refractivo - Config: %d sat, %d tel', ...
    N_sat, N_telescopes), 'Interpreter', 'latex', 'FontSize', 12);

% Añadir leyenda 
if ~isempty(legend_labels)
    legend(legend_labels, 'Location', 'best', 'Interpreter', 'latex');
end

hold off;

% Guardar la figura 
saveas(fig, sprintf('%s/HvsDmin_Refractivo_%s.png', output_dir, config_names{config_idx}));

close(fig);

    
    % Generar archivos CSV y TXT para esta configuración
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
    
    % Generar un archivo de resumen para esta configuración
    summary_filename = sprintf('%s/HvsDmin_%s_resumen.txt', output_dir, config_names{config_idx});
    fileID = fopen(summary_filename, 'w');
    fprintf(fileID, '=== RESUMEN DE ANALISIS CRUZADO ===\n\n');
    fprintf(fileID, 'Configuracion: %s\n', config_names{config_idx});
    fprintf(fileID, 'Numero de satelites: %d\n', N_sat);
    fprintf(fileID, 'Numero de telescopios por satelite: %d\n\n', N_telescopes);
    fprintf(fileID, 'Criterios aplicados:\n');
    fprintf(fileID, '- MTF (Lambda2) > 0.25\n');
    fprintf(fileID, '- SNR (Lambda3) > 400\n');
    fprintf(fileID, '- Dias de cobertura <= 7\n');
    fprintf(fileID, '- Para telescopio Refractivo: Diametro <= 90mm\n\n');
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
                fprintf(fileID, ' %s: Alturas validas %.0f-%.0f km, Diametros %.0f-%.0f mm\n', ...
                    telescope_names{telescopio_idx}, min_alt, max_alt, min_diam, max_diam);
            else
                fprintf(fileID, ' %s: No hay combinaciones validas\n', telescope_names{telescopio_idx});
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end

disp('Analisis cruzado completado. Resultados guardados en la carpeta HvsDmin/');
end

