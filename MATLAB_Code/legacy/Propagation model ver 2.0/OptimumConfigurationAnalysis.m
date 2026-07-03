function OptimumConfigurationAnalysis(configs, telescope_names, GSD, N_pix_12, N_pix_3)

% Generar nombres para configuraciones (CORREGIDO)
config_names = cell(1, size(configs, 1));
for idx = 1:size(configs, 1)
    config_names{idx} = sprintf('%dsat_%dtel', configs(idx,1), configs(idx,2));
end


% Abrir archivo para resultados
fileID = fopen('OptimumConfigs/optimum_configurations.txt', 'w');
fprintf(fileID, '=== ANÁLISIS DE CONFIGURACIONES ÓPTIMAS ===\n\n');

% Inicializar variables para el mínimo global
min_total_mass_global = Inf;
min_config_global = [];
min_data_global = [];

% Para cada configuración (CORREGIDO)
for config_idx = 1:size(configs, 1)
    N_sat = configs(config_idx, 1);
    N_telescopes = configs(config_idx, 2);  % CORREGIDO: columna 2, no 3
    config_name = config_names{config_idx};
    
    fprintf(fileID, '==================================================\n');
    fprintf(fileID, 'CONFIGURACIÓN: %d satélite(s), %d telescopio(s)\n', N_sat, N_telescopes);
    fprintf(fileID, '==================================================\n\n');
    
    % Inicializar variables para el mínimo de esta configuración
    min_total_mass_config = Inf;
    min_detector_config = 0;
    min_telescopio_config = 0;
    min_altura_config = 0;
    min_data_config = [];
    
    % Para cada detector
    for detector_idx = 1:3
        % Para cada telescopio
        for telescopio_idx = 1:4
            % Construir nombre de archivo
            filename = sprintf('Masa_total/MasaTotal_%s_Detector%d_%s.csv', ...
                config_name, detector_idx, telescope_names{telescopio_idx});
            
            % Verificar si existe el archivo
            if ~exist(filename, 'file')
                continue;
            end
            
            % Leer datos
            try
                data = readtable(filename);
                
                % Filtrar valores NaN
                valid_indices = ~isnan(data.Masa_total);
                if sum(valid_indices) == 0
                    continue;
                end
                
                % Calcular masa total real (considerando número de satélites)
                real_total_mass = data.Masa_total .* N_sat;
                
                % Encontrar el mínimo para esta combinación
                [min_mass, min_idx] = min(real_total_mass(valid_indices));
                
                % Guardar los datos de esta combinación
                valid_data = data(valid_indices, :);
                min_data = valid_data(min_idx, :);
                min_altura = min_data.Altura_km;
                
                % Escribir resultados para esta combinación
                fprintf(fileID, '--- Detector %d, Telescopio %s ---\n', detector_idx, telescope_names{telescopio_idx});
                fprintf(fileID, 'Masa total mínima (considerando %d satélite(s)): %.2f kg\n', N_sat, min_mass);
                fprintf(fileID, 'Altura orbital: %d km\n', min_altura);
                fprintf(fileID, 'Diámetro de pupila: %.2f mm\n', min_data.Diametro_pupila);
                fprintf(fileID, 'Masa seca por satélite: %.2f kg\n', min_data.Masa_seca);
                fprintf(fileID, 'Masa combustible por satélite: %.2f kg\n', min_data.Masa_combustible);
                fprintf(fileID, 'Masa total por satélite: %.2f kg\n\n', min_data.Masa_total);
                
                % Actualizar mínimo para esta configuración
                if min_mass < min_total_mass_config
                    min_total_mass_config = min_mass;
                    min_detector_config = detector_idx;
                    min_telescopio_config = telescopio_idx;
                    min_altura_config = min_altura;
                    min_data_config = min_data;
                end
                
            catch
                fprintf('Error al procesar %s\n', filename);
                continue;
            end
        end
    end
    
    % Escribir el mínimo para esta configuración
    if min_total_mass_config < Inf
        fprintf(fileID, '>>> MÍNIMO PARA CONFIGURACIÓN %s <<<\n', config_name);
        fprintf(fileID, 'Masa total mínima (considerando %d satélite(s)): %.2f kg\n', N_sat, min_total_mass_config);
        fprintf(fileID, 'Detector: %d\n', min_detector_config);
        fprintf(fileID, 'Telescopio: %s\n', telescope_names{min_telescopio_config});
        fprintf(fileID, 'Altura orbital: %d km\n\n', min_altura_config);
        
        % Actualizar mínimo global
        if min_total_mass_config < min_total_mass_global
            min_total_mass_global = min_total_mass_config;
            min_config_global = [config_idx, min_detector_config, min_telescopio_config, min_altura_config];
            min_data_global = min_data_config;
        end
    else
        fprintf(fileID, '>>> NO SE ENCONTRARON DATOS VÁLIDOS PARA CONFIGURACIÓN %s <<<\n\n', config_name);
    end
end

% Escribir el mínimo global
if min_total_mass_global < Inf
    % Obtener datos adicionales para el mínimo global
    config_idx = min_config_global(1);
    detector_idx = min_config_global(2);
    telescopio_idx = min_config_global(3);
    altura = min_config_global(4);
    
    N_sat = configs(config_idx, 1);
    N_telescopes = configs(config_idx, 2);  % CORREGIDO
    
    % Determinar N_pix según el detector
    if detector_idx <= 3
        N_pix_value = N_pix_12(detector_idx);
    else
        N_pix_value = N_pix_3(detector_idx-3);
    end
    
    % CORREGIDO: Obtener variables del min_data_global
    num_impulsos = min_data_global.Num_impulsos;
    delta_v_total = min_data_global.DeltaV_total;
    area_instrumento = min_data_global.Sup_media;
    volumen_instrumento = min_data_global.Volumen_medio;
    potencia_instrumento = min_data_global.Potencia_media;
    
    % Calcular el número de detectores necesarios para el swath óptimo
    swath_max_detector = N_pix_value * GSD / 1000; % [km]
    
    % Buscar datos adicionales en archivos MTF, SNR y Coverage
    mtf_file = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx);
    snr_file = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx+3, telescopio_idx);
    
    % Configuración para archivo de cobertura (CORREGIDO)
    config_name = config_names{config_idx};
    telescope_name_clean = strrep(telescope_names{telescopio_idx}, ' ', '');
    coverage_config = sprintf('%dSat_%dTel_%s_Det%d', N_sat, N_telescopes, telescope_name_clean, detector_idx);
    coverage_file = sprintf('coverage/coverage_%s.csv', coverage_config);
    
    % [Resto del código para leer MTF, SNR y Coverage permanece igual...]
    % [Aquí va todo el código para calcular MTF, SNR, swath óptimo, etc.]
    
    % Escribir resultados del mínimo global
    fprintf(fileID, '==================================================\n');
    fprintf(fileID, 'CONFIGURACIÓN ÓPTIMA GLOBAL\n');
    fprintf(fileID, '==================================================\n\n');
    fprintf(fileID, 'Configuración: %d satélite(s), %d telescopio(s)\n', N_sat, N_telescopes);
    fprintf(fileID, 'Detector: %d\n', detector_idx);
    fprintf(fileID, 'Telescopio: %s\n', telescope_names{telescopio_idx});
    fprintf(fileID, 'Masa total (considerando %d satélite(s)): %.2f kg\n', N_sat, min_total_mass_global);
    fprintf(fileID, 'Masa seca por satélite: %.2f kg\n', min_data_global.Masa_seca);
    fprintf(fileID, 'Masa combustible por satélite: %.2f kg\n', min_data_global.Masa_combustible);
    fprintf(fileID, 'Masa total por satélite: %.2f kg\n', min_data_global.Masa_total);
    fprintf(fileID, 'Altura orbital: %d km\n', altura);
    fprintf(fileID, 'Diámetro de pupila: %.2f mm\n', min_data_global.Diametro_pupila);
    fprintf(fileID, 'Número de impulsos: %d\n', num_impulsos);
    fprintf(fileID, 'Delta-V total: %.2f m/s\n', delta_v_total);
    fprintf(fileID, 'Área del instrumento: %.4f m²\n', area_instrumento);
    fprintf(fileID, 'Volumen del instrumento: %.4f m³\n', volumen_instrumento);
    fprintf(fileID, 'Potencia del instrumento: %.2f W\n', potencia_instrumento);
    
    % [Resto de cálculos y salida en consola...]
    
else
    fprintf(fileID, '>>> NO SE ENCONTRARON DATOS VÁLIDOS PARA NINGUNA CONFIGURACIÓN <<<\n');
    fprintf('\n>>> NO SE ENCONTRARON DATOS VÁLIDOS PARA NINGUNA CONFIGURACIÓN <<<\n');
end

fclose(fileID);
fprintf('\nAnálisis de configuraciones óptimas completado.\n');
fprintf('Resultados guardados en Optimum/optimum_configurations.txt\n');

end
