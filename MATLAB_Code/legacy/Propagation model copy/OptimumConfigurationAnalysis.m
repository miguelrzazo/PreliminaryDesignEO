function OptimumConfigurationAnalysis(configs, telescope_names, GSD, lambda2, lambda3, LTAN_hour, N_pix_12, N_pix_3)
    %% AnÃ¡lisis de configuraciones Ã³ptimas
    % Esta funciÃ³n analiza todas las configuraciones de satÃ©lites y encuentra
    % los valores Ã³ptimos (mÃ­nimos) de masa total para cada configuraciÃ³n y
    % el mÃ­nimo absoluto entre todas las configuraciones.
    
    % Generar nombres para configuraciones
    config_names = cell(1, size(configs, 1));
    for idx = 1:size(configs, 1)
        config_names{idx} = sprintf('%dsat_ph%d_%dtel', configs(idx,1), configs(idx,2), configs(idx,3));
    end
    
    % Crear directorio para resultados
    if ~exist('Optimum', 'dir')
        mkdir('Optimum');
    end
    
    % Abrir archivo para resultados
    fileID = fopen('Optimum/optimum_configurations.txt', 'w');
    fprintf(fileID, '=== ANÃLISIS DE CONFIGURACIONES Ã“PTIMAS ===\n\n');
    %fprintf(fileID, 'Fecha de anÃ¡lisis: %s\n\n', datetime(now, 'dd/mm/yyyy HH:MM:SS'));
    
    % Inicializar variables para el mÃ­nimo global
    min_total_mass_global = Inf;
    min_config_global = [];
    
    % Para cada configuraciÃ³n
    for config_idx = 1:size(configs, 1)
        N_sat = configs(config_idx, 1);
        phasing_type = configs(config_idx, 2);
        N_telescopes = configs(config_idx, 3);
        
        config_name = config_names{config_idx};
        
        fprintf(fileID, '==================================================\n');
        fprintf(fileID, 'CONFIGURACIÃ“N: %d satÃ©lite(s), %d telescopio(s)\n', N_sat, N_telescopes);
        fprintf(fileID, '==================================================\n\n');
        
        % Inicializar variables para el mÃ­nimo de esta configuraciÃ³n
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
                filename = sprintf('Masa/MasaTotal_%s_Detector%d_%s.csv', ...
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
                    
                    % Calcular masa total real (considerando nÃºmero de satÃ©lites)
                    real_total_mass = data.Masa_total .* N_sat;
                    
                    % Encontrar el mÃ­nimo para esta combinaciÃ³n
                    [min_mass, min_idx] = min(real_total_mass(valid_indices));
                    
                    % Guardar los datos de esta combinaciÃ³n
                    valid_data = data(valid_indices, :);
                    min_data = valid_data(min_idx, :);
                    min_altura = min_data.Altura_km;
                    
                    % Escribir resultados para esta combinaciÃ³n
                    fprintf(fileID, '--- Detector %d, Telescopio %s ---\n', detector_idx, telescope_names{telescopio_idx});
                    fprintf(fileID, 'Masa total mÃ­nima (considerando %d satÃ©lite(s)): %.2f kg\n', N_sat, min_mass);
                    fprintf(fileID, 'Altura orbital: %d km\n', min_altura);
                    fprintf(fileID, 'DiÃ¡metro de pupila: %.2f mm\n', min_data.Diametro_pupila);
                    fprintf(fileID, 'Masa seca por satÃ©lite: %.2f kg\n', min_data.Masa_seca);
                    fprintf(fileID, 'Masa combustible por satÃ©lite: %.2f kg\n', min_data.Masa_combustible);
                    fprintf(fileID, 'Masa total por satÃ©lite: %.2f kg\n\n', min_data.Masa_total);
                    
                    % Actualizar mÃ­nimo para esta configuraciÃ³n
                    if min_mass < min_total_mass_config
                        min_total_mass_config = min_mass;
                        min_detector_config = detector_idx;
                        min_telescopio_config = telescopio_idx;
                        min_altura_config = min_altura;
                        min_data_config = min_data;
                        num_impulsos = min_data_global.Num_impulsos; % NÃºmero de impulsos
                        delta_v_total = min_data_global.DeltaV_total; % Delta-V total
                        area_instrumento = min_data_global.Sup_media; % Ãrea del instrumento
                        volumen_instrumento = min_data_global.Volumen_medio; % Volumen del instrumento
                        potencia_instrumento = min_data_global.Potencia_media; % Potencia del instrumento
                    end
                    
                catch
                    fprintf('Error al procesar %s\n', filename);
                    continue;
                end
            end
        end
        
        % Escribir el mÃ­nimo para esta configuraciÃ³n
        if min_total_mass_config < Inf
            fprintf(fileID, '>>> MÃNIMO PARA CONFIGURACIÃ“N %s <<<\n', config_name);
            fprintf(fileID, 'Masa total mÃ­nima (considerando %d satÃ©lite(s)): %.2f kg\n', N_sat, min_total_mass_config);
            fprintf(fileID, 'Detector: %d\n', min_detector_config);
            fprintf(fileID, 'Telescopio: %s\n', telescope_names{min_telescopio_config});
            fprintf(fileID, 'Altura orbital: %d km\n\n', min_altura_config);
            
            % Actualizar mÃ­nimo global
            if min_total_mass_config < min_total_mass_global
                min_total_mass_global = min_total_mass_config;
                min_config_global = [config_idx, min_detector_config, min_telescopio_config, min_altura_config];
                min_data_global = min_data_config;
            end
        else
            fprintf(fileID, '>>> NO SE ENCONTRARON DATOS VÃLIDOS PARA CONFIGURACIÃ“N %s <<<\n\n', config_name);
        end
    end
    
    % Escribir el mÃ­nimo global
    if min_total_mass_global < Inf
        % Obtener datos adicionales para el mÃ­nimo global
        config_idx = min_config_global(1);
        detector_idx = min_config_global(2);
        telescopio_idx = min_config_global(3);
        altura = min_config_global(4);
        
        N_sat = configs(config_idx, 1);
        N_telescopes = configs(config_idx, 3);
        
        % Determinar N_pix segÃºn el detector
            if detector_idx <= 3
                N_pix_value = N_pix_12(detector_idx);
            else
                N_pix_value = N_pix_3(detector_idx-3);
            end
        % Calcular el nÃºmero de detectores necesarios para el swath Ã³ptimo
        swath_max_detector = N_pix_value * GSD / 1000; % [km]
        


        % Buscar datos adicionales en archivos MTF, SNR y Coverage
        mtf_file = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx);
        snr_file = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx+3, telescopio_idx);
        
        % ConfiguraciÃ³n para archivo de cobertura
        config_name = config_names{config_idx};
        coverage_file = sprintf('Coverage/Sat%d_Detector%d_Telescopio%d_%d_resultados.csv', ...
            N_sat, detector_idx, telescopio_idx, N_telescopes);
        
        % Leer archivos adicionales
        mtf_value = NaN;
        snr_value = NaN;
        swath = NaN;
        revisit_time = NaN;
        
        % Leer MTF
        if exist(mtf_file, 'file')
            try
                mtf_data = readtable(mtf_file, 'ReadRowNames', true);
                row_name = sprintf('Alt_%d', altura);
                if ismember(row_name, mtf_data.Properties.RowNames)
                    mtf_row = mtf_data(row_name, :);
                    % Buscar el diÃ¡metro mÃ¡s cercano
                    diam_cols = mtf_data.Properties.VariableNames;
                    target_diam = min_data_global.Diametro_pupila;
                    
                    closest_diam = NaN;
                    min_diff = Inf;
                    for i = 1:length(diam_cols)
                        diam_str = diam_cols{i};
                        diam_val = str2double(regexprep(diam_str, 'Diam_', ''));
                        diff = abs(diam_val - target_diam);
                        if diff < min_diff
                            min_diff = diff;
                            closest_diam = diam_str;
                        end
                    end
                    
                    if ~isnan(closest_diam)
                        mtf_value = mtf_row.(closest_diam);
                    end
                end
            catch
                fprintf('Error al leer MTF para el mÃ­nimo global\n');
            end
        end
        
        % Leer SNR
        if exist(snr_file, 'file')
            try
                snr_data = readtable(snr_file, 'ReadRowNames', true);
                row_name = sprintf('Alt_%d', altura);
                if ismember(row_name, snr_data.Properties.RowNames)
                    snr_row = snr_data(row_name, :);
                    % Buscar el diÃ¡metro mÃ¡s cercano
                    diam_cols = snr_data.Properties.VariableNames;
                    target_diam = min_data_global.Diametro_pupila;
                    
                    closest_diam = NaN;
                    min_diff = Inf;
                    for i = 1:length(diam_cols)
                        diam_str = diam_cols{i};
                        diam_val = str2double(regexprep(diam_str, 'Diam_', ''));
                        diff = abs(diam_val - target_diam);
                        if diff < min_diff
                            min_diff = diff;
                            closest_diam = diam_str;
                        end
                    end
                    
                    if ~isnan(closest_diam)
                        snr_value = snr_row.(closest_diam);
                    end
                end
            catch
                fprintf('Error al leer SNR para el mÃ­nimo global\n');
            end
        end
        
        % Leer Coverage
        if exist(coverage_file, 'file')
            try
                coverage_data = readtable(coverage_file, 'ReadRowNames', true);
                row_name = sprintf('Alt_%d', altura);
                if ismember(row_name, coverage_data.Properties.RowNames)
                    coverage_row = coverage_data(row_name, :);
                    
                    % Encontrar el swath Ã³ptimo (menor valor de dÃ­as de revisita)
                    swath_cols = coverage_data.Properties.VariableNames;
                    min_revisit = Inf;
                    best_swath = NaN;
                    
                    for i = 1:length(swath_cols)
                        swath_str = swath_cols{i};
                        revisit_val = coverage_row.(swath_str);
                        
                        if ~isnan(revisit_val) && revisit_val < min_revisit && revisit_val <= 7
                            min_revisit = revisit_val;
                            best_swath = str2double(regexprep(swath_str, 'Swath_', ''));
                        end
                    end
                    
                    swath = best_swath;
                    revisit_time = min_revisit;
                end
            catch
                fprintf('Error al leer Coverage para el mÃ­nimo global\n');
            end
        end
        num_detectores = ceil(swath / swath_max_detector);
        % Calcular parÃ¡metros adicionales
        % InclinaciÃ³n para Ã³rbita heliosÃ­ncrona
        R_tierra = 6371; % Radio medio de la Tierra [km]
        J2 = 1.08263e-3; % Coeficiente J2 de achatamiento terrestre
        arg = -2/3 * (R_tierra/(R_tierra + altura))^(7/2) * 86400/(2*pi) * 2*pi/(365.2422*86400) / J2;
        arg = max(min(arg,1),-1); % Acotar entre -1 y 1
        inc_rad = acos(arg);
        inc_deg = inc_rad * 180/pi;

        
        % Calcular periodo orbital
        mu = 3.986e14; % Constante gravitacional de la Tierra [mÂ³/sÂ²]
        R_altura = (R_tierra + altura) * 1000; % Radio orbital en metros
        periodo_orbital_s = 2 * pi * sqrt(R_altura^3 / mu); % Periodo orbital en segundos
        periodo_orbital_min = periodo_orbital_s / 60; % Periodo orbital en minutos

        % Longitud focal
        pixel_size = 15e-6; % m (valor tÃ­pico para detector 1-3)
        if detector_idx > 3
            pixel_size = 10e-6; % m (valor para detectores 4-6)
        end
        focal_length = (pixel_size * altura * 1000) / GSD; % m
        
        % Escribir resultados del mÃ­nimo global
        fprintf(fileID, '==================================================\n');
        fprintf(fileID, 'CONFIGURACIÃ“N Ã“PTIMA GLOBAL\n');
        fprintf(fileID, '==================================================\n\n');
        
        fprintf(fileID, 'ConfiguraciÃ³n: %d satÃ©lite(s), %d telescopio(s)\n', N_sat, N_telescopes);
        fprintf(fileID, 'Detector: %d\n', detector_idx);
        fprintf(fileID, 'Telescopio: %s\n', telescope_names{telescopio_idx});
        fprintf(fileID, 'Masa total (considerando %d satÃ©lite(s)): %.2f kg\n', N_sat, min_total_mass_global);
        fprintf(fileID, 'Masa seca por satÃ©lite: %.2f kg\n', min_data_global.Masa_seca);
        fprintf(fileID, 'Masa combustible por satÃ©lite: %.2f kg\n', min_data_global.Masa_combustible);
        fprintf(fileID, 'Masa total por satÃ©lite: %.2f kg\n', min_data_global.Masa_total);
        fprintf(fileID, 'Altura orbital: %d km\n', altura);
        fprintf(fileID, 'DiÃ¡metro de pupila: %.2f mm\n', min_data_global.Diametro_pupila);
        fprintf(fileID, 'GSD: %d m\n', GSD);
        fprintf(fileID, 'MTF en Lambda 2 (%.2f Âµm): %.4f\n', lambda2*1e6, mtf_value);
        fprintf(fileID, 'SNR en Lambda 3 (%.2f Âµm): %.2f\n', lambda3*1e6, snr_value);
        fprintf(fileID, 'Swath Ã³ptimo: %.2f km\n', swath);
        fprintf(fileID, 'Tiempo de revisita: %.2f dÃ­as\n', revisit_time);
        fprintf(fileID, 'InclinaciÃ³n orbital: %.2f grados\n', inc_deg);
        fprintf(fileID, 'Longitud focal: %.4f m\n', focal_length);
        fprintf(fileID, 'Longitud media instrumento: %.4f\n', min_data_global.Longitud_media);
        fprintf(fileID, 'Superficie media instrumento: %.4f\n', min_data_global.Sup_media);
        fprintf(fileID, 'Volumen medio instrumento: %.4f\n', min_data_global.Volumen_medio);
        fprintf(fileID, 'Potencia media instrumento: %.4f W\n', min_data_global.Potencia_media);
        fprintf(fileID, 'Periodo orbital: %.2f minutos\n', periodo_orbital_min);
        fprintf(fileID, 'Hora local del nodo ascendente (LTAN): %.2f horas\n', LTAN_hour);
        fprintf(fileID, 'NÃºmero de impulsos: %d\n', num_impulsos);
        fprintf(fileID, 'Delta-V total: %.2f m/s\n', delta_v_total);
        fprintf(fileID, 'Ãrea del instrumento: %.4f mÂ²\n', area_instrumento);
        fprintf(fileID, 'NÃºmero de detectores necesarios: %d\n', num_detectores);


        % Imprimir en consola tambiÃ©n
        fprintf('\n==================================================\n');
        fprintf('CONFIGURACIÃ“N Ã“PTIMA GLOBAL\n');
        fprintf('==================================================\n\n');
        
        fprintf('ConfiguraciÃ³n: %d satÃ©lite(s), %d telescopio(s)\n', N_sat, N_telescopes);
        fprintf('Detector: %d\n', detector_idx);
        fprintf('Telescopio: %s\n', telescope_names{telescopio_idx});
        fprintf('Masa total (considerando %d satÃ©lite(s)): %.2f kg\n', N_sat, min_total_mass_global);
        fprintf('Masa seca por satÃ©lite: %.2f kg\n', min_data_global.Masa_seca);
        fprintf('Masa combustible por satÃ©lite: %.2f kg\n', min_data_global.Masa_combustible);
        fprintf('Masa total por satÃ©lite: %.2f kg\n', min_data_global.Masa_total);
        fprintf('Altura orbital: %d km\n', altura);
        fprintf('DiÃ¡metro de pupila: %.2f mm\n', min_data_global.Diametro_pupila);
        fprintf('GSD: %d m\n', GSD);
        fprintf('MTF en Lambda 2 (%.2f Âµm): %.4f\n', lambda2*1e6, mtf_value);
        fprintf('SNR en Lambda 3 (%.2f Âµm): %.2f\n', lambda3*1e6, snr_value);
        fprintf('Swath Ã³ptimo: %.2f km\n', swath);
        fprintf('Tiempo de revisita: %.2f dÃ­as\n', revisit_time);
        fprintf('InclinaciÃ³n orbital: %.2f grados\n', inc_deg);
        fprintf('Longitud focal: %.4f m\n', focal_length);
        fprintf('Periodo orbital: %.2f minutos\n', periodo_orbital_min);
        fprintf('Hora local del nodo ascendente (LTAN): %.2f horas\n', LTAN_hour);
        fprintf('NÃºmero de impulsos: %d\n', num_impulsos);
        fprintf('Delta-V total: %.2f m/s\n', delta_v_total);
        fprintf('Ãrea del instrumento: %.4f mÂ²\n', area_instrumento);
        fprintf('Volumen del instrumento: %.4f mÂ³\n', volumen_instrumento);
        fprintf('Consumo de potencia del instrumento: %.2f W\n', potencia_instrumento);
        fprintf('NÃºmero de detectores necesarios: %d\n', num_detectores);
    else
        
        fprintf(fileID, '>>> NO SE ENCONTRARON DATOS VÃLIDOS PARA NINGUNA CONFIGURACIÃ“N <<<\n');
        fprintf('\n>>> NO SE ENCONTRARON DATOS VÃLIDOS PARA NINGUNA CONFIGURACIÃ“N <<<\n');
    end
    
    fclose(fileID);
    fprintf('\nAnÃ¡lisis de configuraciones Ã³ptimas completado.\n');
    fprintf('Resultados guardados en Optimum/optimum_configurations.txt\n');
end