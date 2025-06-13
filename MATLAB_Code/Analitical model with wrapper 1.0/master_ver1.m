    %%
clear; clc; close all;
                
lambda1 = 1.61e-06;                 % Banda 1 Co2
lambda2 = 2.01e-6;                  % Banda 2 CO2
lambda3 = 0.76e-6;                  % Banda Visible O2 A-band
GSD = 80; % GSD

%% Definir rangos
alturas_orbitales = 200:10:1000;    % Alturas orbitales en km
diametros_pupila = 20:1:600;        % Diámetros pupila en mm
SNR_Requirement = 400;              % SNR minimo
MTF_Requirement = 0.25;             % MTF minimo
Cov_Requirement = 7;                % días revisita
max_detectores = 3;                 % Límite máximo de detectores
swaths_km = 20:10:200;               % Array de swaths a calcular (km)                   
solapamiento = 0.05;                % Factor de solapamiento del swath
cobertura_nubes= 1/6;               % Pasadas no efectivas
LTAN_hour = 6;                   % Hora de paso local Dawn Dusk

% Configuraciones a simular
satellite_configs = [
    %1, 1;
    2, 1;
    1, 2;
    3, 1; % 1 satélite,  1 telescopio
    4, 1;
];

%% DETECTOR
% Detector 1 CAPYORK
eta1 = 0.85;
MTF_detector1 = 0.45;
N_pix1 = 1200;
pixel_size1 = 15e-6;

% Detector 2 H2RG
eta2 = 0.7;
MTF_detector2 = 0.45;
N_pix2 = 2048;
pixel_size2 = 18e-6;

% Detector 3 SATURN VISIR
eta3 = 0.6;
MTF_detector3 = 0.5;
N_pix3 = 1000;
pixel_size3 = 30e-6;

% Vector detectores 1-3 (para lambda1 y lambda2)
eta_12 = [eta1, eta2, eta3];
MTF_detector_12 = [MTF_detector1, MTF_detector2, MTF_detector3];
N_pix_12 = [N_pix1, N_pix2, N_pix3];
pixel_size_12 = [pixel_size1, pixel_size2, pixel_size3];

%% O2 Band
% Detector 4 CMOS
eta4 = 0.35;
MTF_detector4 = 0.36;
N_pix4 = 512;
pixel_size4 = 25e-6;


% Detector 5 H2RG
eta5 = 0.7;
MTF_detector5 = 0.5;
N_pix5 = 2048;
pixel_size5 = 18e-6;

% Detector 6 SATURN VISIR
eta6 = 0.6;
MTF_detector6 = 0.45;
N_pix6 = 1000;
pixel_size6 = 30e-6;


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
            diametros_pupila, filename_prefix, telescope_names{j}, i,MTF_Requirement);
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
            diametros_pupila, filename_prefix, telescope_names{j}, i,MTF_Requirement);
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
            diametros_pupila, filename_prefix, telescope_names{j}, i, MTF_Requirement);
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
            diametros_pupila, filename_prefix, telescope_names{j}, i,SNR_Requirement);
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
            diametros_pupila, filename_prefix, telescope_names{j}, i,SNR_Requirement);
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
            diametros_pupila, filename_prefix, telescope_names{j}, i,SNR_Requirement);
    end
end

%% COVERAGE (USANDO REVISITCALC)
disp('====================================================================');
disp('Iniciando análisis de cobertura con RevisitCalc');
disp('====================================================================');

% Bucle para cada configuración de satélite/telescopio
for config_idx = 1:size(satellite_configs, 1)
    N_sat = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);

    config_desc = sprintf('Config: %d satélite(s), %d telescopio(s)', N_sat, N_telescopes);
    disp(['Procesando ' config_desc]);

    % Bucle para cada tipo de detector (IDs 1-3)
    for detector_id = 1:3
        % Bucle para cada tipo de telescopio (IDs 1-4)
        for telescope_id = 1:4
            
            % Obtener parámetros para la combinación actual
            current_Npix = N_pix_12(detector_id);
            current_fov_limit = fov_limit(telescope_id);
            current_telescope_name = telescope_names{telescope_id};

            fprintf(' -> Calculando para Detector %d, Telescopio %s...\n', detector_id, current_telescope_name);

            % Llamada a la función de cálculo de cobertura
            try
                CoverageRevisitCalc(GSD, alturas_orbitales, swaths_km, ...
                                    Cov_Requirement, current_Npix, ...
                                    current_fov_limit, N_sat, N_telescopes, ...
                                    detector_id, current_telescope_name, ...
                                    solapamiento, cobertura_nubes, max_detectores);
            catch ME
                fprintf('    -> X Error al procesar: %s\n', ME.message);
                continue;
            end
        end
    end
    
    disp(['✓ ' config_desc ' completada.']);
    fprintf('\n');

end

disp('====================================================================');
disp('Análisis de cobertura finalizado para todas las configuraciones.');
disp('====================================================================');


%% CRUCE DE DATOS
CrossDataFunction(alturas_orbitales, swaths_km, diametros_pupila, telescope_names, fov_limit, satellite_configs)

%% ========================================================================
%  ANÁLISIS EXHAUSTIVO FINAL Y GENERACIÓN DE REPORTES
%  (Esta sección reemplaza los cálculos de masa y el análisis de óptimo anteriores)
%  ========================================================================

% --- Empaquetar todos los parámetros de entrada en un único struct
%     para pasarlos a la función de análisis.
params.satellite_configs = satellite_configs;
params.telescope_names   = telescope_names;
params.alturas_orbitales = alturas_orbitales;
params.swaths_km         = swaths_km;
params.diametros_pupila  = diametros_pupila;
params.Cov_Requirement   = Cov_Requirement;
params.GSD               = GSD;
params.N_pix_12          = N_pix_12; % Necesario para calcular N_detectores_swath
params.fov_limit         = fov_limit;
params.LTAN_hour         = LTAN_hour;

% --- Llamar a la función principal de análisis
%     Esta función se encargará de todo el procesamiento final y la
%     generación de los reportes en CSV y TXT.
ComprehensiveAnalysisFunction(params);

disp('====================================================================');
disp('SCRIPT PRINCIPAL FINALIZADO.');
disp('Todos los resultados han sido generados.');
disp('====================================================================');