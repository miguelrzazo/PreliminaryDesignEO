%%
clear; clc; close all;
                
lambda1 = 1.61e-06;                 % Banda 1 Co2
lambda2 = 2.01e-6;                  % Banda 2 CO2
lambda3 = 0.76e-6;                  % Banda Visible O2 A-band
GSD = 80; % GSD

%% Definir rangos
alturas_orbitales = 400:10:1400;    % Alturas orbitales en km
diametros_pupila = 20:1:600;         % Diámetros pupilares en mm
SNR_Requirement = 400;              % SNR minimo
MTF_Requirement = 0.25;             % MTF minimo
Cov_Requirement = 7;                % días
max_detectores = 3;                 % Límite máximo de detectores
swaths_km = 20:1:200;              % Array de swaths a calcular (km)                   
solapamiento = 0.05;                % Factor de solapamiento del swath
cobertura_nubes= 1/6;               % Pasadas no efectivas
LTAN_hour = 10.5;                   % Hora de paso local (10:30)

% Configuraciones a simular
satellite_configs = [
    1, 1; % 1 satélite, 1 telescopio
    2, 1; % 2 satélites, 1 telescopio
    1, 2; % 1 satélites, 2 telescopios
];

%% DETECTOR
% Detector 1 CAPYORK
eta1 = 0.6;
MTF_detector1 = 0.45;
N_pix1 = 1024;
pixel_size1 = 15e-6;

% Detector 2 CHROMA D
eta2 = 0.8;
MTF_detector2 = 0.4;
N_pix2 = 3000;
pixel_size2 = 18e-6;

% Detector 3 SATURN VISIR
eta3 = 0.6;
MTF_detector3 = 0.45;
N_pix3 = 1000;
pixel_size3 = 30e-6;

% Vector detectores 1-3 (para lambda1 y lambda2)
eta_12 = [eta1, eta2, eta3];
MTF_detector_12 = [MTF_detector1, MTF_detector2, MTF_detector3];
N_pix_12 = [N_pix1, N_pix2, N_pix3];
pixel_size_12 = [pixel_size1, pixel_size2, pixel_size3];

%% O2 Band
% Detector 4 SATURN VISIR
eta4 = 0.6;
MTF_detector4 = 0.5;
N_pix4 = 1000;
pixel_size4 = 30e-6;

% Detector 5 Chroma
eta5 = 0.8;
MTF_detector5 = 0.4;
N_pix5 = 3000;
pixel_size5 = 18e-6;

% Detector 6 CMOS
eta6 = 0.35;
MTF_detector6 = 0.36;
N_pix6 = 512;
pixel_size6 = 25e-6;

% Vector detectores 4-6 (para lambda3)
eta_3 = [eta4, eta5, eta6];
MTF_detector_3 = [MTF_detector4, MTF_detector5, MTF_detector6];
N_pix_3 = [N_pix4, N_pix5, N_pix6];
pixel_size_3 = [pixel_size4, pixel_size5, pixel_size6];

%% TELESCOPIOS
% Refractivo 1
MTF_alineamiento1 = 0.9;
tau1 = 0.8;
fov_limit_deg1 = 10;

% Korsch 2
MTF_alineamiento2 = 0.85;
tau2 = 0.65;
fov_limit_deg2 = 3;

% Cassegrein 3
MTF_alineamiento3 = 0.8;
tau3 = 0.7;
fov_limit_deg3 = 3;

% TMA 4
MTF_alineamiento4 = 0.7;
tau4 = 0.6;
fov_limit_deg4 = 8;

% Matriz telescopio
MTF_telescope = [MTF_alineamiento1, MTF_alineamiento2, MTF_alineamiento3, MTF_alineamiento4];
tau_telescope = [tau1, tau2, tau3, tau4];
fov_limit = [fov_limit_deg1, fov_limit_deg2, fov_limit_deg3, fov_limit_deg4];
telescope_names = {'Refractivo', 'Korsch', 'Cassegrain', 'TMA'};
R = [0, 0.2, 0.2, 0]; % Radio 0 para Refractivo y TMA que no son obscurados

%% Crear directorios
dirs = {'MTF', 'SNR', 'Coverage', 'HvsDmin', 'Masa_seca', 'Masa_total', 'OptimumConfigs'};
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i});
    end
end

%% MTF para lambda1 (detectores 1-3)
for i = 1:3 % Bucle detectores 1-3
    for j = 1:4 % Bucle telescopios
        % Configurar parámetros para esta combinación
        current_pixel_size = pixel_size_12(i);
        current_MTF_detector = MTF_detector_12(i);
        current_MTF_alineamiento = MTF_telescope(j);
        R_obs = R(j);
        
        % Nombre de archivo para esta combinación
        filename_prefix = sprintf('MTF/MTF_Lambda1_Detector%d_Telescopio%d', i, j);
        
        % Llamar a la función MTF con los parámetros actuales
        MTFfunction(lambda1, current_pixel_size, current_MTF_detector, ...
            current_MTF_alineamiento, GSD, R_obs, alturas_orbitales, ...
            diametros_pupila, filename_prefix, telescope_names{j}, i);
    end
end

%% MTF para lambda2 (detectores 1-3)
for i = 1:3 % Bucle detectores 1-3
    for j = 1:4 % Bucle telescopios
        % Configurar parámetros para esta combinación
        current_pixel_size = pixel_size_12(i);
        current_MTF_detector = MTF_detector_12(i);
        current_MTF_alineamiento = MTF_telescope(j);
        R_obs = R(j);
        
        % Nombre de archivo para esta combinación
        filename_prefix = sprintf('MTF/MTF_Lambda2_Detector%d_Telescopio%d', i, j);
        
        % Llamar a la función MTF con los parámetros actuales
        MTFfunction(lambda2, current_pixel_size, current_MTF_detector, ...
            current_MTF_alineamiento, GSD, R_obs, alturas_orbitales, ...
            diametros_pupila, filename_prefix, telescope_names{j}, i);
    end
end

%% MTF para lambda3 (detectores 4-6)
for i = 1:3 % Bucle detectores 4-6
    for j = 1:4 % Bucle telescopios
        % Configurar parámetros para esta combinación
        current_pixel_size = pixel_size_3(i);
        current_MTF_detector = MTF_detector_3(i);
        current_MTF_alineamiento = MTF_telescope(j);
        R_obs = R(j);
        
        % Nombre de archivo para esta combinación
        filename_prefix = sprintf('MTF/MTF_Lambda3_Detector%d_Telescopio%d', i+3, j);
        
        % Llamar a la función MTF con los parámetros actuales
        MTFfunction(lambda3, current_pixel_size, current_MTF_detector, ...
            current_MTF_alineamiento, GSD, R_obs, alturas_orbitales, ...
            diametros_pupila, filename_prefix, telescope_names{j}, i);
    end
end

%% SNR para lambda1 (detectores 1-3)
for i = 1:3 % Bucle detectores 1-3
    for j = 1:4 % Bucle telescopios
        % Configurar parámetros para esta combinación
        current_pixel_size = pixel_size_12(i);
        current_eta = eta_12(i);
        current_tau = tau_telescope(j);
        R_obs = R(j);
        
        % Nombre de archivo para esta combinación
        filename_prefix = sprintf('SNR/SNR_Lambda1_Detector%d_Telescopio%d', i, j);
        
        % Llamar a la función SNR con los parámetros actuales
        SNRfunction(lambda1, current_pixel_size, current_eta, ...
            current_tau, GSD, R_obs, alturas_orbitales, ...
            diametros_pupila, filename_prefix, telescope_names{j}, i);
    end
end

%% SNR para lambda2 (detectores 1-3)
for i = 1:3 % Bucle detectores 1-3
    for j = 1:4 % Bucle telescopios
        % Configurar parámetros para esta combinación
        current_pixel_size = pixel_size_12(i);
        current_eta = eta_12(i);
        current_tau = tau_telescope(j);
        R_obs = R(j);
        
        % Nombre de archivo para esta combinación
        filename_prefix = sprintf('SNR/SNR_Lambda2_Detector%d_Telescopio%d', i, j);
        
        % Llamar a la función SNR con los parámetros actuales
        SNRfunction(lambda2, current_pixel_size, current_eta, ...
            current_tau, GSD, R_obs, alturas_orbitales, ...
            diametros_pupila, filename_prefix, telescope_names{j}, i);
    end
end

%% SNR para lambda3 (detectores 4-6)
for i = 1:3 % Bucle detectores 4-6
    for j = 1:4 % Bucle telescopios
        % Configurar parámetros para esta combinación
        current_pixel_size = pixel_size_3(i);
        current_eta = eta_3(i);
        current_tau = tau_telescope(j);
        R_obs = R(j);
        
        % Nombre de archivo para esta combinación
        filename_prefix = sprintf('SNR/SNR_Lambda3_Detector%d_Telescopio%d', i+3, j);
        
        % Llamar a la función SNR con los parámetros actuales
        SNRfunction(lambda3, current_pixel_size, current_eta, ...
            current_tau, GSD, R_obs, alturas_orbitales, ...
            diametros_pupila, filename_prefix, telescope_names{j}, i);
    end
end

%% COVERAGE 
% Bucle para cada configuración
parfor config_idx = 1:size(satellite_configs, 1)
    % Extraer parámetros de la configuración actual
    N_sat = satellite_configs(config_idx, 1); % Número de satélites
    N_telescopes = satellite_configs(config_idx, 2); % Número de telescopios
    
    % Descripción de la configuración para mostrar en consola
    config_desc = sprintf('Configuración %d: %d satélite(s), %d telescopio(s)', ...
        N_sat, N_telescopes);
    disp(['Procesando ' config_desc]);
    
    % Bucle para cada detector
    for i = 1:3
        % Bucle para cada telescopio
        for j = 1:4
            % Configurar parámetros para esta combinación
            current_Npix = N_pix_12(i);
            current_fov_limit = fov_limit(j);
            
            % Llamar a la función de cobertura con los parámetros actuales
            fprintf(' Calculando cobertura para Detector %d, Telescopio %d...\n', i, j);
            
            % LLAMADA
            CoverageSSOAnaliticalfunction(GSD, alturas_orbitales, swaths_km, Cov_Requirement, ...
               current_Npix, current_fov_limit, N_sat, ...
               N_telescopes, i, j, telescope_names{j}, solapamiento, cobertura_nubes, max_detectores);

        end
    end
    disp(['Completada ' config_desc]);
end

disp('Análisis de cobertura completado para todas las configuraciones.');


%% CRUCE DE DATOS
CrossDataFunction(GSD, alturas_orbitales, swaths_km, N_pix_12, diametros_pupila, telescope_names, fov_limit,satellite_configs)

%% CÁLCULO DE MASA SECA Y TOTAL

disp('Calculando masa seca y total para cada configuración...');

for config_idx = 1:size(satellite_configs, 1)
    N_sat = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    config_name = sprintf('%dsat_%dtel', N_sat, N_telescopes);

    for detector_idx = 1:3
        for telescopio_idx = 1:4
            filename = sprintf('HvsDmin/HvsDmin_%s_Detector%d.csv', config_name, detector_idx);
            if ~exist(filename, 'file')
                warning('No se encontró el archivo %s. Saltando...', filename);
                continue;
            end

            datos = readtable(filename);
            alturas_km = datos.Altura_km;
            diam_col_name = sprintf('Dmin_%s', telescope_names{telescopio_idx});
            if ~ismember(diam_col_name, datos.Properties.VariableNames)
                continue;
            end
            diametros_mm = datos.(diam_col_name);

            valid_indices = ~isnan(diametros_mm);
            if sum(valid_indices) == 0
                continue;
            end

            alturas_valid = alturas_km(valid_indices);
            diametros_valid = diametros_mm(valid_indices);

            [masa_seca, tabla_final] = calcularMasaSeca(alturas_valid, diametros_valid, N_telescopes);

            masa_seca_filename = sprintf('%s/MasaSeca_%s_Detector%d_%s.csv', ...
                'Masa_seca', config_name, detector_idx, telescope_names{telescopio_idx});
            writetable(tabla_final, masa_seca_filename);

            [total_mass, fuel_mass, num_impulses, total_delta_v] = calcularMasaTotal(alturas_valid, masa_seca, tabla_final.Sup_media);

            tabla_final.Masa_combustible = fuel_mass;
            tabla_final.Masa_total = total_mass;
            tabla_final.Num_impulsos = num_impulses;
            tabla_final.DeltaV_total = total_delta_v;

            masa_total_filename = sprintf('%s/MasaTotal_%s_Detector%d_%s.csv', ...
                'Masa_total', config_name, detector_idx, telescope_names{telescopio_idx});
            writetable(tabla_final, masa_total_filename);

            txt_filename = sprintf('%s/Resumen_%s_Detector%d_%s.txt', ...
                'Masa_total', config_name, detector_idx, telescope_names{telescopio_idx});
            fid = fopen(txt_filename, 'w');
            fprintf(fid, '=== ANÁLISIS DE MASA ===\n\n');
            fprintf(fid, 'Configuración: %d satélite(s), %d telescopio(s)\n', N_sat, N_telescopes);
            fprintf(fid, 'Detector: %d\n', detector_idx);
            fprintf(fid, 'Telescopio: %s\n\n', telescope_names{telescopio_idx});
            valid_mass_indices = ~isnan(total_mass);
            if sum(valid_mass_indices) > 0
                valid_total_mass = total_mass(valid_mass_indices);
                valid_alturas_mass = alturas_valid(valid_mass_indices);
                [min_mass, min_idx] = min(valid_total_mass);
                min_altura = valid_alturas_mass(min_idx);
                fprintf(fid, 'Estadísticas de masa total:\n');
                fprintf(fid, 'Masa mínima: %.2f kg (Altura: %d km)\n', min_mass, min_altura);
                fprintf(fid, 'Masa máxima: %.2f kg\n', max(valid_total_mass));
                fprintf(fid, 'Masa promedio: %.2f kg\n', mean(valid_total_mass));
                fprintf(fid, 'Rango de alturas válidas: %d - %d km\n', min(valid_alturas_mass), max(valid_alturas_mass));
            else
                fprintf(fid, 'No hay valores válidos de masa total\n');
            end
            fclose(fid);

            fig1 = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
            subplot(3,1,1);
            plot(alturas_valid, num_impulses, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            title(sprintf('N\\''umero de Impulsos - %d sat - %d tel - Detector %d - %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{telescopio_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Altura Orbital (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('N\\''umero de Impulsos', 'Interpreter', 'latex', 'FontSize', 11);
            grid on;
            subplot(3,1,2);
            plot(alturas_valid, total_delta_v, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            title(sprintf('$\\Delta V$ Total - %d sat - %d tel - Detector %d - %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{telescopio_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Altura Orbital (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('$\\Delta V$ Total (m/s)', 'Interpreter', 'latex', 'FontSize', 11);
            grid on;
            subplot(3,1,3);
            plot(alturas_valid, fuel_mass, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            title(sprintf('Masa de Combustible - %d sat - %d tel - Detector %d - %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{telescopio_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Altura Orbital (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('Masa de Combustible (kg)', 'Interpreter', 'latex', 'FontSize', 11);
            grid on;
            saveas(fig1, sprintf('%s/Parametros_Orbitales_%s_Detector%d_%s.png', ...
                masa_total_dir, config_name, detector_idx, telescope_names{telescopio_idx}));
            close(fig1);

            fig2 = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
            hold on;
            plot(alturas_valid, fuel_mass, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Masa Combustible');
            plot(alturas_valid, masa_seca, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Masa Seca');
            plot(alturas_valid, total_mass, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Masa Total');
            valid_mass_indices = ~isnan(total_mass);
            if sum(valid_mass_indices) > 0
                valid_total_mass = total_mass(valid_mass_indices);
                valid_alturas_mass = alturas_valid(valid_mass_indices);
                [min_mass, min_idx] = min(valid_total_mass);
                min_altura = valid_alturas_mass(min_idx);
                plot(min_altura, min_mass, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
                plot([min_altura min_altura], [0 min_mass], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                text(min_altura, min_mass*1.1, sprintf('M\\''in: %.1f kg\\newlineAltura: %d km', min_mass, min_altura), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 10, ...
                    'BackgroundColor', 'white', 'EdgeColor', 'black');
            end
            title(sprintf('Desglose de Masas - %d sat - %d tel - Detector %d - %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{telescopio_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Altura Orbital (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('Masa (kg)', 'Interpreter', 'latex', 'FontSize', 11);
            legend('Location', 'best', 'Interpreter', 'latex');
            grid on;
            hold off;
            saveas(fig2, sprintf('%s/DesgloseMasas_%s_Detector%d_%s.png', ...
                masa_total_dir, config_name, detector_idx, telescope_names{telescopio_idx}));
            close(fig2);
        end
    end
end

disp('Análisis de masa completado para todas las configuraciones.');


%% === ANÁLISIS DE CONFIGURACIONES ÓPTIMAS ===
OptimumConfigurationAnalysis(satellite_configs, telescope_names, GSD, N_pix_12, N_pix_3);
