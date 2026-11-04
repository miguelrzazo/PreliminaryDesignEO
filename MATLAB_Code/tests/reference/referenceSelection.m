function selection = referenceSelection(alturas_orbitales, swaths_km, diametros_pupila, telescope_names, fov_limit, configs,MTF_Requirement,SNR_Requirement,Cov_Requirement, analysisResults, writeOutputs)
%CROSSDATAFUNCTION Select minimum diameters from current-run matrices.
% Optional analysisResults uses the shared axes/optics/coverage contract.
% Calls without it retain the historical CSV-import behavior.
% Returned cells contain altitude-by-telescope-by-detector diameter arrays.
if nargin < 10, analysisResults = []; end
if nargin < 11, writeOutputs = true; end
useStoredFiles = isempty(analysisResults);
if ~useStoredFiles
    validateAnalysisAxes(analysisResults, alturas_orbitales, diametros_pupila, ...
        swaths_km, telescope_names, configs);
end
selection = cell(size(configs,1),1);

% Crear directorio para resultados si no existe
output_dir = 'HvsDmin';
if writeOutputs && ~exist(output_dir, 'dir')
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
    Dmin_results = nan(length(alturas_orbitales), length(telescope_names), 3); % soluciones con cobertura y FOV
    Dmin_theoretical = nan(length(alturas_orbitales), length(telescope_names), 3); % solo MTF y SNR

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
            config_name_coverage = sprintf('%d Satelite(s); %d Telescopio(s): %s; Detector %d', N_sat, N_telescopes, telescope_name_clean, detector_idx);
            coverage_file = sprintf('coverage/coverage_%s.csv', config_name_coverage);

            % Verificar existencia de los archivos de datos
            if useStoredFiles && (~exist(mtf_file, 'file') || ~exist(snr_file, 'file') || ~exist(coverage_file, 'file'))
                fprintf('Advertencia: Faltan archivos para Detector %d, Telescopio %d (%s), Configuracion %s.\n', ...
                        detector_idx, telescopio_idx, telescope_names{telescopio_idx}, config_names{config_idx});
                fprintf('  MTF file: %s (Existe: %d)\n', mtf_file, exist(mtf_file, 'file'));
                fprintf('  SNR file: %s (Existe: %d)\n', snr_file, exist(snr_file, 'file'));
                fprintf('  Coverage file: %s (Existe: %d)\n', coverage_file, exist(coverage_file, 'file'));
                continue; % Saltar a la siguiente iteración de telescopio
            end

            % Leer archivos CSV
            try
                if useStoredFiles
                    mtf_matrix = table2array(readtable(mtf_file, 'ReadRowNames', true));
                    snr_matrix = table2array(readtable(snr_file, 'ReadRowNames', true));
                    coverage_data = readmatrix(coverage_file);
                else
                    mtf_matrix = analysisResults.optics.mtf{2,detector_idx,telescopio_idx};
                    snr_matrix = analysisResults.optics.snr{3,detector_idx,telescopio_idx};
                    coverage_data = analysisResults.coverage.days{config_idx,detector_idx,telescopio_idx};
                    validateAnalysisMatrices(mtf_matrix,snr_matrix,coverage_data, ...
                        alturas_orbitales,diametros_pupila,swaths_km);
                end

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

                    theoretical_diameter_indices = find(mtf_row >= MTF_Requirement & snr_row >= SNR_Requirement);
                    if ~isempty(theoretical_diameter_indices)
                        theoretical_diam_col_idx = theoretical_diameter_indices(1);
                        if theoretical_diam_col_idx <= length(diametros_pupila)
                            Dmin_theoretical(alt_idx, telescopio_idx, detector_idx) = diametros_pupila(theoretical_diam_col_idx);
                        end
                    end

                    % Preserve the existing first-feasible swath in input-axis order.
                    valid_swath_indices = find(coverage_row <= Cov_Requirement & ~isnan(coverage_row));
                    if isempty(valid_swath_indices)
                        continue; % No hay swaths válidos para esta altura y configuración de cobertura
                    end
                    % On ascending axes, the first feasible entry is the smallest swath.
                    selected_swath_col_idx = valid_swath_indices(1);
                    selected_swath_km = swaths_km(selected_swath_col_idx);

                    % 2. Verificar si el FOV del swath seleccionado excede el límite del telescopio
                    % FOV (full angle) = 2 * atan( (Swath/2) / Altura )
                    current_fov_limit = fov_limit(telescopio_idx);
                    if N_telescopes == 2
                        effective_fov_limit = current_fov_limit * 2;
                    else
                        effective_fov_limit = current_fov_limit;
                    end
                    fov_actual_deg = 2 * atand(selected_swath_km / (2 * altura));
                    if fov_actual_deg > effective_fov_limit
                        continue; % FOV excede el límite para este telescopio y swath
                    end

                    % 3. Encontrar el diámetro mínimo que cumple MTF >= 0.25 y SNR >= 400 para el swath seleccionado
                    % Se asume que las columnas de mtf_row y snr_row corresponden a diametros_pupila
                    valid_diameter_indices = find(mtf_row >= MTF_Requirement & snr_row >= SNR_Requirement);
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
                    % Preserve the existing 200 mm limit for the first telescope.
                    if telescopio_idx == 1 && selected_dmin_mm > 200
                        continue; 
                    end

                    % Almacenar el Dmin encontrado
                    Dmin_results(alt_idx, telescopio_idx, detector_idx) = selected_dmin_mm;
                end
            catch ME
                if ~useStoredFiles, rethrow(ME); end
                fprintf('Error procesando archivos para Detector %d, Telescopio %d (%s), Config %s: %s\n', ...
                        detector_idx, telescopio_idx, telescope_names{telescopio_idx}, config_names{config_idx}, ME.message);
                continue; % Saltar al siguiente telescopio en caso de error de lectura/procesamiento
            end
        end % Fin loop telescopios
    end % Fin loop detectores

    selection{config_idx}.viableDiameters_mm = Dmin_results;
    selection{config_idx}.theoreticalDiameters_mm = Dmin_theoretical;
    if ~writeOutputs, continue; end

end % Fin loop configuraciones

disp('Cross-data analysis complete.');

end % Fin de la función
