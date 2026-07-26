function ConfigAnalysis(params)

% ConfigAnalysis: Realiza un análisis exhaustivo de todas las
% configuraciones viables, integrando rendimiento, coste y parámetros orbitales.

%% ========================================================================
% 1. INICIALIZACIÓN Y DESEMPAQUETADO DE PARÁMETROS
% ========================================================================
fprintf('\n==========================================================\n');
fprintf('INICIANDO ANÁLISIS EXHAUSTIVO DE CONFIGURACIONES...\n');
fprintf('==========================================================\n\n');

% --- Desempaquetar el struct de parámetros
satellite_configs = params.satellite_configs;
telescope_names = params.telescope_names;
alturas_orbitales = params.alturas_orbitales;
swaths_km = params.swaths_km;
Cov_Requirement = params.Cov_Requirement;
GSD = params.GSD;
Pixel_size = params.Pixel_size;
LTAN_hour = params.LTAN_hour;

% --- Constantes físicas y orbitales
R_E_m = 6378.137 * 1000;
mu_E = 3.986004418e14;
J2 = 1.08262668e-3;

% --- Configuración de salida
output_dir = 'OptimumConfigs';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ========================================================================
% 2. PROCESAMIENTO DE SOLUCIONES Y CÁLCULO DE PARÁMETROS
% ========================================================================
results_data = {};
global_best.mass = Inf;

for config_idx = 1:size(satellite_configs, 1)
    N_sat = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    config_name = sprintf('%dsat_%dtel', N_sat, N_telescopes);

    for detector_idx = 1:3
        h_vs_dmin_filename = sprintf('HvsDmin/HvsDmin_%s_Detector%d.csv', config_name, detector_idx);
        if ~exist(h_vs_dmin_filename, 'file'), continue; end
        
        datos_validos_h_dmin = readtable(h_vs_dmin_filename);
        
        for telescopio_idx = 1:length(telescope_names)
            diam_col_name = sprintf('Dmin_%s_mm', strrep(telescope_names{telescopio_idx},' ','_'));
            if ~ismember(diam_col_name, datos_validos_h_dmin.Properties.VariableNames), continue; end
            
            try
                mtf_data_table = readtable(sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx), 'ReadRowNames', true);
                snr_data_table = readtable(sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx + 3, telescopio_idx), 'ReadRowNames', true);
            catch
                fprintf('Advertencia: No se pudieron cargar archivos MTF/SNR para Det %d, Tel %d. Saltando...\n', detector_idx, telescopio_idx);
                continue;
            end

            for sol_idx = 1:height(datos_validos_h_dmin)
                h_km = datos_validos_h_dmin.Altura_km(sol_idx);
                D_mm = datos_validos_h_dmin.(diam_col_name)(sol_idx);

                if isnan(D_mm), continue; end

                coverage_matrix = readmatrix(sprintf('coverage/coverage_%d Satelite(s); %d Telescopio(s): %s; Detector %d.csv', N_sat, N_telescopes, strrep(telescope_names{telescopio_idx}, ' ', ''), detector_idx));
                alt_idx_coverage = find(alturas_orbitales == h_km, 1);
                if isempty(alt_idx_coverage), continue; end
                
                coverage_row = coverage_matrix(alt_idx_coverage, :);
                valid_swath_indices = find(coverage_row <= Cov_Requirement & ~isnan(coverage_row));
                if isempty(valid_swath_indices), continue; end
                
                selected_swath_km = swaths_km(valid_swath_indices(1));
                dias_cobertura = coverage_row(valid_swath_indices(1));
                
                current_Npix = params.Npix_detectors(detector_idx);
                num_detectores_necesarios = ceil(selected_swath_km / (GSD * current_Npix / 1000));

                swath_real_sistema_km = num_detectores_necesarios * GSD * current_Npix / 1000;
                excedente_swath_km = swath_real_sistema_km - selected_swath_km;
                excedente_swath_porc = (excedente_swath_km / selected_swath_km) * 100;
                
                row_name_perf = sprintf('Alt_%d', h_km);
                col_name_perf = sprintf('Diam_%d', round(D_mm));
                
                try
                    mtf_val = mtf_data_table{row_name_perf, col_name_perf};
                    snr_val = snr_data_table{row_name_perf, col_name_perf};
                catch
                    mtf_val = NaN; snr_val = NaN;
                end
                
                h_m = h_km * 1000;
                pixel_size_m = Pixel_size(detector_idx);
                D_m = D_mm / 1000;
                focal_length_m = (pixel_size_m * h_m) / GSD;
                focal_ratio = focal_length_m / D_m;

                [masa_seca, masa_seca_TM, masa_seca_SEOSAT, L_instrumento, S_instrumento, ...
                V_instrumento, W_instrumento, P_instrumento, V_sat, L_sat, S_sat, U_sat] = calcularMasaSeca(D_mm, N_telescopes);

                [masa_total_satelite, masa_comb, num_impulsos, delta_v] = calcularMasaTotal(h_km, masa_seca, S_sat);

                if isnan(masa_total_satelite), continue; end

                masa_total_constelacion = masa_total_satelite * N_sat;
                a_m = (h_km * 1000) + R_E_m;
                SSO_nodal_rate_rad_s = 2 * pi / (365.2422 * 86400);
                cos_i = - (2 * SSO_nodal_rate_rad_s * a_m^(3.5)) / (3 * J2 * R_E_m^2 * sqrt(mu_E));
                cos_i = min(max(cos_i, -1), 1);
                inclinacion_deg = acosd(cos_i);
                desfase_anomalia = 360 / N_sat;
                
                solution_data = {
                    N_sat, N_telescopes, num_detectores_necesarios, detector_idx, telescope_names{telescopio_idx}, h_km, D_mm, ...
                    selected_swath_km, dias_cobertura, ...
                    swath_real_sistema_km, excedente_swath_km, excedente_swath_porc, ...
                    inclinacion_deg, LTAN_hour, desfase_anomalia, ...
                    mtf_val, snr_val, ...
                    focal_length_m, focal_ratio, ...
                    W_instrumento, L_instrumento, S_instrumento, V_instrumento, P_instrumento, ...
                    masa_seca, masa_comb, masa_total_satelite, masa_total_constelacion, ...
                    num_impulsos, delta_v, ...
                    masa_seca_TM, masa_seca_SEOSAT, ...
                    L_sat, S_sat, V_sat, U_sat
                };
                
                results_data(end+1, :) = solution_data;

                if masa_total_constelacion < global_best.mass
                    global_best.mass = masa_total_constelacion;
                    global_best.data = solution_data;
                end
            end
        end
    end
end

fprintf('Análisis detallado completado. Se encontraron %d soluciones viables.\n', size(results_data,1));

%% ========================================================================
% 3. EXPORTACIÓN DE RESULTADOS (CSV Y TXT)
% ========================================================================
if isempty(results_data)
    fprintf('ADVERTENCIA: No se encontraron soluciones válidas.\n');
    return;
end

headers = {
    'Num_Satelites', 'Num_Telescopios', 'Detectores_Necesarios', 'ID_Detector', 'Tipo_Telescopio', 'Altura_km', 'Diametro_Pupila_mm', ...
    'Swath_Requerido_km', 'Dias_Cobertura', ...
    'Swath_Real_Sistema_km', 'Excedente_Swath_km', 'Excedente_Swath_porc', ...
    'Inclinacion_SSO_deg', 'LTAN_h', 'Desfase_Anomalia_Verdadera_deg', ...
    'MTF_a_Nyquist_lambda2', 'SNR_lambda3', ...
    'Longitud_Focal_m', 'Relacion_Focal_f_num', ...
    'Masa_Instrumento_kg', 'Longitud_Instrumento_m', 'Area_Instrumento_m2', 'Volumen_Instrumento_m3', 'Potencia_Electrica_Instrumento_W', ...
    'Masa_Seca_Satelite_kg', 'Masa_Combustible_kg', 'Masa_Total_Satelite_kg', 'Masa_Total_Constelacion_kg', ...
    'Num_Impulsos_Mantenimiento', 'DeltaV_Total_ms', ...
    'Masa_Seca_TM_kg', 'Masa_Seca_SEOSAT_kg', ...
    'Longitud_Satelite_m', 'Area_Satelite_m2', 'Volumen_Satelite_m3', 'Estandar_U'
};

results_table = cell2table(results_data, 'VariableNames', headers);
results_table = sortrows(results_table, 'Masa_Total_Constelacion_kg', 'ascend');

csv_filename = fullfile(output_dir, 'reporte_completo_configuraciones.csv');
writetable(results_table, csv_filename);
fprintf('Reporte CSV completo guardado en: %s\n', csv_filename);

txt_filename = fullfile(output_dir, 'reporte_resumen_configuraciones.txt');
fileID = fopen(txt_filename, 'w');
fprintf(fileID, '==========================================================\n');
fprintf(fileID, ' INFORME DE ANÁLISIS DE CONFIGURACIONES\n');
fprintf(fileID, '==========================================================\n\n');

if ~isinf(global_best.mass)
    fprintf(fileID, '--- ÓPTIMO GLOBAL (Masa Mínima de Constelación) ---\n');
    data = global_best.data;
    fprintf(fileID, 'Masa Total de la Constelación: %.2f kg\n\n', global_best.mass);
    fprintf(fileID, ' > Configuración: %d Satélite(s), %d Telescopio(s)/sat, %d Detector(es) Necesario(s), ID Detector %d con %s\n', data{1}, data{2}, data{3}, data{4}, data{5});
    fprintf(fileID, ' > Órbita: %d km, %.2f° inc, Desfase Anom. %.1f°\n', data{6}, data{13}, data{15});
    fprintf(fileID, ' > Rendimiento: Swath Requerido=%.1f km (Real=%.1f km, Excedente=%.1f%%) | Cobertura=%.1f días\n', data{8}, data{10}, data{12}, data{9});
    fprintf(fileID, ' > Parámetros Ópticos: D=%.1f mm, F=%.2f m (f/%.1f) | MTF@Nyquist=%.3f, SNR=%.1f\n', data{7}, data{18}, data{19}, data{16}, data{17});
    fprintf(fileID, ' > Prop. Físicas (Instrumento): Masa=%.2f kg, Long=%.2f m, Área=%.2f m^2, Vol=%.6f m^3, Potencia=%.1f W\n', data{20}, data{21}, data{22}, data{23}, data{24});
    fprintf(fileID, ' > Prop. Físicas (Satélite): Long=%.2f m, Área=%.2f m^2, Vol=%.6f m^3, Estándar U=%.2f U\n', data{33}, data{34}, data{35}, data{36});
    fprintf(fileID, ' > Masa Satélite: Total=%.2f kg (Seca=%.2f kg, Combustible=%.2f kg)\n', data{27}, data{25}, data{26});
    fprintf(fileID, ' > Mantenimiento: %d impulsos, Delta-V Total %.1f m/s\n\n', data{29}, data{30});
end

fprintf(fileID, '\n==========================================================\n');
fprintf(fileID, ' LISTADO DE TODAS LAS CONFIGURACIONES VÁLIDAS\n');
fprintf(fileID, ' (Ordenadas por masa total de la constelación)\n');
fprintf(fileID, '==========================================================\n\n');

for i = 1:height(results_table)
    row = table2cell(results_table(i,:));
    fprintf(fileID, '--- Solución #%d ---\n', i);
    fprintf(fileID, 'Masa Constelación: %.2f kg | Masa Satélite Total: %.2f kg (Seca: %.2f kg, Comb: %.2f kg)\n', row{28}, row{27}, row{25}, row{26});
    fprintf(fileID, ' > Configuración: %d Sat, %d Tel, %d Dets | ID Detector: %d | Telescopio: %s\n', row{1}, row{2}, row{3}, row{4}, row{5});
    fprintf(fileID, ' > Órbita: %.0f km alt, %.2f° inc (SSO), LTAN %d h, Desfase Anom. %.1f°\n', row{6}, row{13}, row{14}, row{15});
    fprintf(fileID, ' > Rendimiento: Swath Requerido=%.1f km (Real=%.1f km, Excedente=%.1f%%) | Cobertura=%.1f días\n', row{8}, row{10}, row{12}, row{9});
    fprintf(fileID, ' > Parámetros Ópticos: D=%.1f mm, F=%.2f m (f/%.1f) | MTF@Nyquist=%.3f, SNR=%.1f\n', row{7}, row{18}, row{19}, row{16}, row{17});
    fprintf(fileID, ' > Prop. Físicas (Satélite): Long=%.2f m, Área=%.2f m^2, Vol=%.6f m^3, Estándar U=%.2f U\n', row{33}, row{34}, row{35}, row{36});
    fprintf(fileID, '\n');
end

fclose(fileID);
fprintf('Reporte TXT detallado guardado en: %s\n', txt_filename);

%% ========================================================================
% 4. GENERACIÓN DE GRÁFICAS DE MASA
% ========================================================================
fprintf('\nIniciando la generación de gráficas de desglose de masa...\n');
GenerateMassPlots(results_table);
fprintf('\nAnálisis finalizado.\n');
end
