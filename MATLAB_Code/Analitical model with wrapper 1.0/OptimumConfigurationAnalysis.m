function OptimumConfigurationAnalysis(configs, telescope_names)
% OptimumConfigurationAnalysis: Analiza los resultados de masa para encontrar la 
% configuración óptima global.
%
% VERSIÓN REFACTORIZADA:
% - Adaptada para funcionar con los archivos CSV generados por los scripts optimizados.
% - Modularizada en funciones locales para mayor claridad y mantenimiento.
% - Eliminada la relectura innecesaria de archivos MTF/SNR/Coverage.
% - Mejor manejo de errores y estructura de datos para los resultados.
% - Se han eliminado los argumentos de entrada no utilizados (GSD, N_pix).

%% --- Configuración Inicial ---
OUTPUT_DIR = 'OptimumConfigs';
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

output_filename = fullfile(OUTPUT_DIR, 'optimum_configurations.txt');
fileID = fopen(output_filename, 'w');
% Asegura que el archivo se cierre al finalizar la función, incluso si hay un error
cleanupObj = onCleanup(@() fclose(fileID));

fprintf(fileID, '=== ANÁLISIS DE CONFIGURACIONES ÓPTIMAS ===\n\n');

%% --- Búsqueda del Óptimo Global ---
global_best.mass = Inf;
global_best.data = [];
global_best.config_idx = 0;
global_best.detector_idx = 0;
global_best.telescopio_idx = 0;

for config_idx = 1:size(configs, 1)
    N_sat = configs(config_idx, 1);
    N_telescopes = configs(config_idx, 2);
    config_name = sprintf('%dsat_%dtel', N_sat, N_telescopes);

    fprintf(fileID, '==================================================\n');
    fprintf(fileID, 'CONFIGURACIÓN: %d satélite(s), %d telescopio(s)\n', N_sat, N_telescopes);
    fprintf(fileID, '==================================================\n\n');

    config_best.mass = Inf;

    for detector_idx = 1:3
        for telescopio_idx = 1:length(telescope_names)
            
            % 1. Encontrar el mínimo en el archivo correspondiente
            filename = sprintf('Masa_total/MasaTotal_%s_Detector%d_%s.csv', ...
                               config_name, detector_idx, telescope_names{telescopio_idx});

            [min_data, min_mass_constellation] = find_minimum_in_file(filename, N_sat);

            if ~isinf(min_mass_constellation)
                % 2. Escribir el resumen de esta combinación
                write_combination_summary(fileID, detector_idx, telescope_names{telescopio_idx}, ...
                                          N_sat, min_mass_constellation, min_data);

                % 3. Actualizar el mejor de la configuración actual
                if min_mass_constellation < config_best.mass
                    config_best.mass = min_mass_constellation;
                    config_best.data = min_data;
                    config_best.detector_idx = detector_idx;
                    config_best.telescopio_idx = telescopio_idx;
                end
            end
        end
    end

    % 4. Escribir el resumen de la configuración
    if ~isinf(config_best.mass)
        write_config_summary(fileID, config_name, N_sat, config_best, telescope_names);
        
        % 5. Actualizar el óptimo global
        if config_best.mass < global_best.mass
            global_best = config_best;
            global_best.config_idx = config_idx;
        end
    else
        fprintf(fileID, '>>> NO SE ENCONTRARON DATOS VÁLIDOS PARA LA CONFIGURACIÓN %s <<<\n\n', config_name);
    end
end

% 6. Escribir el resumen del óptimo global
if ~isinf(global_best.mass)
    write_global_optimum(fileID, global_best, configs, telescope_names);
else
    fprintf(fileID, '\n>>> NO SE ENCONTRARON DATOS VÁLIDOS PARA NINGUNA CONFIGURACIÓN <<<\n');
    fprintf('\n>>> ADVERTENCIA: No se encontraron datos válidos. El análisis no pudo completarse. <<<\n');
end

fprintf('\nAnálisis de configuraciones óptimas completado.\n');
fprintf('Resultados guardados en %s\n', output_filename);
end


%% =========== FUNCIONES AUXILIARES (LOCAL FUNCTIONS) ===========

function [min_row, min_mass_constellation] = find_minimum_in_file(filename, N_sat)
    % Lee un archivo CSV, encuentra el mínimo y devuelve la fila y la masa.
    min_row = [];
    min_mass_constellation = Inf;

    if ~exist(filename, 'file'), return; end

    try
        data = readtable(filename);
        if isempty(data) || ~ismember('Masa_total_kg', data.Properties.VariableNames), return; end

        valid_indices = ~isnan(data.Masa_total_kg);
        if ~any(valid_indices), return; end
        
        valid_data = data(valid_indices, :);
        constellation_mass = valid_data.Masa_total_kg * N_sat;

        [min_mass_constellation, min_idx] = min(constellation_mass);
        min_row = valid_data(min_idx, :);

    catch ME
        fprintf('Error al procesar el archivo %s: %s\n', filename, ME.message);
    end
end

function write_combination_summary(fileID, det_idx, tel_name, N_sat, min_mass, min_data)
    % Escribe en el archivo de texto el resumen para una combinación específica.
    fprintf(fileID, '--- Detector %d, Telescopio %s ---\n', det_idx, tel_name);
    fprintf(fileID, 'Masa total mínima (considerando %d satélite(s)): %.2f kg\n', N_sat, min_mass);
    fprintf(fileID, 'Altura orbital: %d km\n', min_data.Altura_km);
    fprintf(fileID, 'Diámetro de pupila: %.2f mm\n', min_data.Diametro_pupila_mm);
    fprintf(fileID, 'Masa seca por satélite: %.2f kg\n', min_data.Masa_seca_kg);
    fprintf(fileID, 'Masa combustible por satélite: %.2f kg\n', min_data.Masa_combustible_kg);
    fprintf(fileID, 'Masa total por satélite: %.2f kg\n\n', min_data.Masa_total_kg);
end

function write_config_summary(fileID, config_name, N_sat, config_best, telescope_names)
    % Escribe en el archivo de texto el resumen para la mejor opción de una configuración.
    fprintf(fileID, '>>> MÍNIMO PARA CONFIGURACIÓN %s <<<\n', config_name);
    fprintf(fileID, 'Masa total mínima (considerando %d satélite(s)): %.2f kg\n', N_sat, config_best.mass);
    fprintf(fileID, 'Detector: %d\n', config_best.detector_idx);
    fprintf(fileID, 'Telescopio: %s\n', telescope_names{config_best.telescopio_idx});
    fprintf(fileID, 'Altura orbital: %d km\n\n', config_best.data.Altura_km);
end

function write_global_optimum(fileID, global_best, configs, telescope_names)
    % Escribe el informe final detallado de la configuración óptima global.
    % Extrae todos los datos de la estructura 'global_best'.
    
    config_idx = global_best.config_idx;
    N_sat = configs(config_idx, 1);
    N_telescopes = configs(config_idx, 2);
    
    fprintf(fileID, '==================================================\n');
    fprintf(fileID, '           CONFIGURACIÓN ÓPTIMA GLOBAL\n');
    fprintf(fileID, '==================================================\n\n');
    
    fprintf(fileID, 'Configuración: %d satélite(s), %d telescopio(s)\n', N_sat, N_telescopes);
    fprintf(fileID, 'Detector: %d\n', global_best.detector_idx);
    fprintf(fileID, 'Telescopio: %s\n', telescope_names{global_best.telescopio_idx});
    fprintf(fileID, 'Masa total de la constelación: %.2f kg\n\n', global_best.mass);

    fprintf(fileID, '--- Detalles del Satélite Óptimo ---\n');
    fprintf(fileID, 'Masa total por satélite: %.2f kg\n', global_best.data.Masa_total_kg);
    fprintf(fileID, '  - Masa seca por satélite: %.2f kg\n', global_best.data.Masa_seca_kg);
    fprintf(fileID, '  - Masa combustible por satélite: %.2f kg\n', global_best.data.Masa_combustible_kg);
    fprintf(fileID, 'Altura orbital óptima: %d km\n', global_best.data.Altura_km);
    fprintf(fileID, 'Diámetro de pupila requerido: %.2f mm\n\n', global_best.data.Diametro_pupila_mm);

    fprintf(fileID, '--- Parámetros del Instrumento y Misión ---\n');
    fprintf(fileID, 'Número de impulsos de mantenimiento: %d\n', global_best.data.Num_impulsos);
    fprintf(fileID, 'Delta-V total: %.2f m/s\n', global_best.data.DeltaV_total_ms);
    fprintf(fileID, 'Área media del instrumento: %.4f m²\n', global_best.data.Sup_media);
    fprintf(fileID, 'Volumen medio del instrumento: %.4f m³\n', global_best.data.Volumen_medio);
    fprintf(fileID, 'Potencia media del instrumento: %.2f W\n', global_best.data.Potencia_media_W);
end
