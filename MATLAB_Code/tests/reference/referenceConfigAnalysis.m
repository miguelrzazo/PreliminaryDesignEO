function summary = referenceConfigAnalysis(params, analysisResults, writeOutputs)
%CONFIGANALYSIS Consume current-run matrices and optionally export reports.
% Omitting analysisResults retains the historical CSV-import entry point.
if nargin < 2, analysisResults = []; end
if nargin < 3, writeOutputs = true; end
useStoredFiles = isempty(analysisResults);
if ~useStoredFiles
    validateAnalysisAxes(analysisResults, params.alturas_orbitales, ...
        params.diametros_pupila, params.swaths_km, ...
        params.telescope_names, params.satellite_configs);
end

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
MTF_Requirement = params.MTF_Requirement;
SNR_Requirement = params.SNR_Requirement;
diametros_pupila = params.diametros_pupila;
GSD = params.GSD;
Pixel_size = params.Pixel_size;
LTAN_hour = params.LTAN_hour;
if isfield(params, 'maintenance_lower_factor')
    maintenance_lower_factor = params.maintenance_lower_factor;
else
    maintenance_lower_factor = 0.98;
end
if isfield(params, 'maintenance_upper_factor')
    maintenance_upper_factor = params.maintenance_upper_factor;
else
    maintenance_upper_factor = 1.02;
end

% --- Constantes físicas y orbitales
R_E_m = 6378.137 * 1000;
mu_E = 3.986004418e14;
J2 = 1.08262668e-3;

% --- Configuración de salida
output_dir = 'OptimumConfigs';
if writeOutputs && ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ========================================================================
% 2. PROCESAMIENTO DE SOLUCIONES Y CÁLCULO DE PARÁMETROS
% ========================================================================
results_data = {};
global_best.mass = Inf;
mass_cache = containers.Map('KeyType','char','ValueType','any');

for config_idx = 1:size(satellite_configs, 1)
    N_sat = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    config_name = sprintf('%dsat_%dtel', N_sat, N_telescopes);

    for detector_idx = 1:3
        h_vs_dmin_filename = sprintf('HvsDmin/HvsDmin_%s_Detector%d.csv', config_name, detector_idx);
        if useStoredFiles
            if ~exist(h_vs_dmin_filename, 'file'), continue; end
            datos_validos_h_dmin = readtable(h_vs_dmin_filename);
        else
            % Preserve the report schema without serializing intermediate data.
            columnNames = [{'Altura_km'}, cellfun(@(name) ...
                sprintf('Dmin_%s_mm',strrep(name,' ','_')), ...
                telescope_names, 'UniformOutput',false)];
            selectedDiameters_mm = analysisResults.selection{config_idx}.viableDiameters_mm(:,:,detector_idx);
            datos_validos_h_dmin = array2table([alturas_orbitales(:), selectedDiameters_mm], ...
                'VariableNames',columnNames);
        end

        for telescopio_idx = 1:length(telescope_names)
            diam_col_name = sprintf('Dmin_%s_mm', strrep(telescope_names{telescopio_idx},' ','_'));
            if ~ismember(diam_col_name, datos_validos_h_dmin.Properties.VariableNames), continue; end
            
            try
                if useStoredFiles
                    mtf_data_table = readtable(sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx), 'ReadRowNames', true);
                    snr_data_table = readtable(sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx + 3, telescopio_idx), 'ReadRowNames', true);
                else
                    mtf_matrix = analysisResults.optics.mtf{2,detector_idx,telescopio_idx};
                    snr_matrix = analysisResults.optics.snr{3,detector_idx,telescopio_idx};
                    coverage_matrix = analysisResults.coverage.days{config_idx,detector_idx,telescopio_idx};
                    validateAnalysisMatrices(mtf_matrix,snr_matrix,coverage_matrix, ...
                        alturas_orbitales,diametros_pupila,swaths_km);
                end
            catch dataError
                if ~useStoredFiles, rethrow(dataError); end
                fprintf('Advertencia: No se pudieron cargar archivos MTF/SNR para Det %d, Tel %d. Saltando...\n', detector_idx, telescopio_idx);
                continue;
            end

            for sol_idx = 1:height(datos_validos_h_dmin)
                h_km = datos_validos_h_dmin.Altura_km(sol_idx);
                D_mm = datos_validos_h_dmin.(diam_col_name)(sol_idx);

                if isnan(D_mm), continue; end

                if useStoredFiles
                    coverage_matrix = readmatrix(sprintf('coverage/coverage_%d Satelite(s); %d Telescopio(s): %s; Detector %d.csv', N_sat, N_telescopes, strrep(telescope_names{telescopio_idx}, ' ', ''), detector_idx));
                end
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
                    if useStoredFiles
                        mtf_val = mtf_data_table{row_name_perf, col_name_perf};
                        snr_val = snr_data_table{row_name_perf, col_name_perf};
                    else
                        diameterIndex = find(diametros_pupila == D_mm, 1);
                        assert(~isempty(diameterIndex), 'eo:SelectedDiameter', ...
                            'Selected diameter is absent from the declared axis.');
                        mtf_val = mtf_matrix(alt_idx_coverage,diameterIndex);
                        snr_val = snr_matrix(alt_idx_coverage,diameterIndex);
                    end
                catch dataError
                    if ~useStoredFiles, rethrow(dataError); end
                    mtf_val = NaN; snr_val = NaN;
                end
                
                h_m = h_km * 1000;
                pixel_size_m = Pixel_size(detector_idx);
                D_m = D_mm / 1000;
                focal_length_m = (pixel_size_m * h_m) / GSD;
                focal_ratio = focal_length_m / D_m;

                [masa_seca, masa_seca_TM, masa_seca_SEOSAT, L_instrumento, S_instrumento, ...
                V_instrumento, W_instrumento, P_instrumento, V_sat, L_sat, S_sat, U_sat] = referenceDryMass(D_mm, N_telescopes);

                [masa_total_satelite, masa_comb, num_impulsos, delta_v] = cachedMasaTotal(mass_cache, h_km, masa_seca, S_sat, maintenance_lower_factor, maintenance_upper_factor);

                if isnan(masa_total_satelite), continue; end

                masa_total_constelacion = masa_total_satelite * N_sat;
                a_m = (h_km * 1000) + R_E_m;
                SSO_nodal_rate_rad_s = 2 * pi / (365.2422 * 86400);
                cos_i = - (2 * SSO_nodal_rate_rad_s * a_m^(3.5)) / (3 * J2 * R_E_m^2 * sqrt(mu_E));
                cos_i = min(max(cos_i, -1), 1);
                inclinacion_deg = acosd(cos_i);
                desfase_anomalia = 360 / N_sat;
                
                solution_data = {true, ...
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

% --- Computar masa teorica para todas las alturas (solo MTF+SNR, sin cobertura) ---
viable_keys = containers.Map('KeyType','char','ValueType','logical');
for r = 1:size(results_data,1)
    rd = results_data(r,:);
    viable_keys(sprintf('%d_%s_%d_%d_%d', rd{5}, rd{6}, rd{7}, rd{2}, rd{3})) = true;
end

for config_idx = 1:size(satellite_configs, 1)
    N_sat = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    for detector_idx = 1:3
        for telescopio_idx = 1:length(telescope_names)
            try
                if useStoredFiles
                    mtf_t = readtable(sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d_resultados.csv', detector_idx, telescopio_idx), 'ReadRowNames', true);
                    snr_t = readtable(sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d_resultados.csv', detector_idx + 3, telescopio_idx), 'ReadRowNames', true);
                    mtf_mat = table2array(mtf_t);
                    snr_mat = table2array(snr_t);
                else
                    mtf_mat = analysisResults.optics.mtf{2,detector_idx,telescopio_idx};
                    snr_mat = analysisResults.optics.snr{3,detector_idx,telescopio_idx};
                    coverage_matrix = analysisResults.coverage.days{config_idx,detector_idx,telescopio_idx};
                    validateAnalysisMatrices(mtf_mat,snr_mat,coverage_matrix, ...
                        alturas_orbitales,diametros_pupila,swaths_km);
                end
            catch dataError
                if ~useStoredFiles, rethrow(dataError); end
                continue;
            end
            max_alt_idx = min([length(alturas_orbitales), size(mtf_mat,1), size(snr_mat,1)]);
            max_diam_idx = min([length(diametros_pupila), size(mtf_mat,2), size(snr_mat,2)]);
            for alt_idx = 1:max_alt_idx
                h_km = alturas_orbitales(alt_idx);
                D_mm = NaN;
                for d = 1:max_diam_idx
                    if mtf_mat(alt_idx, d) >= MTF_Requirement && snr_mat(alt_idx, d) >= SNR_Requirement
                        D_mm = diametros_pupila(d); break;
                    end
                end
                if isnan(D_mm), continue; end
                key = sprintf('%d_%s_%d_%d_%d', detector_idx, telescope_names{telescopio_idx}, h_km, N_sat, N_telescopes);
                if isKey(viable_keys, key), continue; end
                [masa_seca, masa_seca_TM, masa_seca_SEOSAT, ~, ~, ~, ~, ~, ...
                 V_sat, L_sat, S_sat, U_sat] = referenceDryMass(D_mm, N_telescopes);
                [masa_total_satelite, masa_comb, num_impulsos, delta_v] = cachedMasaTotal(mass_cache, h_km, masa_seca, S_sat, maintenance_lower_factor, maintenance_upper_factor);
                if isnan(masa_total_satelite), continue; end
                desfase_anomalia = 360 / N_sat;
                solution_data = {false, N_sat, N_telescopes, NaN, detector_idx, telescope_names{telescopio_idx}, h_km, D_mm, ...
                    NaN, NaN, NaN, NaN, NaN, ...
                    NaN, LTAN_hour, desfase_anomalia, ...
                    NaN, NaN, NaN, NaN, ...
                    NaN, NaN, NaN, NaN, NaN, ...
                    masa_seca, masa_comb, masa_total_satelite, masa_total_satelite * N_sat, ...
                    num_impulsos, delta_v, masa_seca_TM, masa_seca_SEOSAT, ...
                    L_sat, S_sat, V_sat, U_sat};
                results_data(end+1, :) = solution_data;
            end
        end
    end
end

summary.rows = results_data;
summary.best = global_best;
if isempty(results_data)
    num_viables = 0;
else
    num_viables = sum(cell2mat(results_data(:, 1)));
end
fprintf('Análisis detallado completado. Se encontraron %d soluciones viables y %d soluciones teóricas.\n', ...
    num_viables, size(results_data,1) - num_viables);

end

function [total_mass, fuel_mass, num_impulses, total_delta_v] = cachedMasaTotal(cache, h_km, masa_seca, area_sat, lower_factor, upper_factor)
    key = sprintf('%.6f_%.9f_%.9f_%.6f_%.6f', h_km, masa_seca, area_sat, lower_factor, upper_factor);
    if isKey(cache, key)
        values = cache(key);
    else
        [total_mass, fuel_mass, num_impulses, total_delta_v] = referenceTotalMass(h_km, masa_seca, area_sat, lower_factor, upper_factor);
        values = [total_mass, fuel_mass, num_impulses, total_delta_v];
        cache(key) = values;
    end
    total_mass = values(1);
    fuel_mass = values(2);
    num_impulses = values(3);
    total_delta_v = values(4);
end

