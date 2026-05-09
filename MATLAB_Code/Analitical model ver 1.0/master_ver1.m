%%
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

%% Wavelengths (m)
lambda = [1.61e-6, 2.01e-6, 0.76e-6];   % [CO2 band 1, CO2 band 2, O2 A-band]

GSD = 50;                                % Ground sampling distance (m)

%% Orbital and mission ranges
alturas_orbitales = 400:10:1400;         % Orbital altitudes (km)
diametros_pupila  = 20:1:600;            % Pupil diameters (mm)
SNR_Requirement   = 400;
MTF_Requirement   = 0.25;
Cov_Requirement   = 7;                   % days
max_detectores    = 3;
swaths_km         = 20:1:200;            % Swath widths to evaluate (km)
solapamiento      = 0.05;                % Swath overlap fraction
cobertura_nubes   = 1/6;                 % Fraction of passes blocked by cloud
LTAN_hour         = 10.5;               % Local Time of Ascending Node (10:30)

satellite_configs = [
    1, 1;   % 1 satellite,  1 telescope
    2, 1;   % 2 satellites, 1 telescope
    1, 2;   % 1 satellite,  2 telescopes
];

%% Detectors — CO2 bands (lambda 1 & 2): detectors 1-3
%                    [CAPYORK, CHROMA D, SATURN VISIR]
eta_12        = [0.60,  0.80,  0.60];
MTF_det_12    = [0.45,  0.40,  0.45];
N_pix_12      = [1024,  3000,  1000];
pixel_size_12 = [15,    18,    30  ] * 1e-6;

%% Detectors — O2 A-band (lambda 3): detectors 4-6
%                    [SATURN VISIR, CHROMA D, CMOS]
eta_3         = [0.60,  0.80,  0.35];
MTF_det_3     = [0.50,  0.40,  0.36];
N_pix_3       = [1000,  3000,  512 ];
pixel_size_3  = [30,    18,    25  ] * 1e-6;

%% Telescopes
%                       [Refractivo, Korsch, Cassegrain, TMA]
MTF_telescope   = [0.90, 0.85, 0.80, 0.70];
tau_telescope   = [0.80, 0.65, 0.70, 0.60];
fov_limit       = [10,   3,    3,    8   ];
telescope_names = {'Refractivo', 'Korsch', 'Cassegrain', 'TMA'};
R               = [0, 0.2, 0.2, 0];     % Secondary mirror obstruction ratio

%% Create output directories
dirs = {'MTF', 'SNR', 'Coverage', 'HvsDmin', 'Masa_seca', 'Masa_total', 'OptimumConfigs'};
for k = 1:numel(dirs)
    if ~exist(dirs{k}, 'dir'), mkdir(dirs{k}); end
end

%% MTF analysis — all bands
% {lambda, pixel_size_vec, MTF_det_vec, detector_id_offset, band_tag}
mtf_bands = {
    lambda(1), pixel_size_12, MTF_det_12, 0, 'Lambda1';
    lambda(2), pixel_size_12, MTF_det_12, 0, 'Lambda2';
    lambda(3), pixel_size_3,  MTF_det_3,  3, 'Lambda3';
};
for b = 1:size(mtf_bands, 1)
    lam    = mtf_bands{b,1};
    psz    = mtf_bands{b,2};
    mtfd   = mtf_bands{b,3};
    offset = mtf_bands{b,4};
    tag    = mtf_bands{b,5};
    for i = 1:3
        for j = 1:4
            fname = sprintf('MTF/MTF_%s_Detector%d_Telescopio%d', tag, i+offset, j);
            MTFfunction(lam, psz(i), mtfd(i), MTF_telescope(j), GSD, R(j), ...
                alturas_orbitales, diametros_pupila, fname, telescope_names{j}, i, MTF_Requirement);
        end
    end
end

%% SNR analysis — all bands
% {lambda, pixel_size_vec, eta_vec, detector_id_offset, band_tag}
snr_bands = {
    lambda(1), pixel_size_12, eta_12, 0, 'Lambda1';
    lambda(2), pixel_size_12, eta_12, 0, 'Lambda2';
    lambda(3), pixel_size_3,  eta_3,  3, 'Lambda3';
};
for b = 1:size(snr_bands, 1)
    lam    = snr_bands{b,1};
    psz    = snr_bands{b,2};
    eta    = snr_bands{b,3};
    offset = snr_bands{b,4};
    tag    = snr_bands{b,5};
    for i = 1:3
        for j = 1:4
            fname = sprintf('SNR/SNR_%s_Detector%d_Telescopio%d', tag, i+offset, j);
            SNRfunction(lam, psz(i), eta(i), tau_telescope(j), GSD, R(j), ...
                alturas_orbitales, diametros_pupila, fname, telescope_names{j}, i, SNR_Requirement);
        end
    end
end

%% Coverage analysis
parfor config_idx = 1:size(satellite_configs, 1)
    N_sat        = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    fprintf('Processing config: %d sat, %d tel\n', N_sat, N_telescopes);
    for i = 1:3
        for j = 1:4
            fprintf('  Coverage: Detector %d, Telescope %d...\n', i, j);
            CoverageSSOAnaliticalfunction(GSD, alturas_orbitales, swaths_km, Cov_Requirement, ...
                N_pix_12(i), fov_limit(j), N_sat, N_telescopes, i, j, telescope_names{j}, ...
                solapamiento, cobertura_nubes, max_detectores);
        end
    end
    fprintf('Done: %d sat, %d tel\n', N_sat, N_telescopes);
end
disp('Coverage analysis complete.');

%% Cross-data analysis
CrossDataFunction(GSD, alturas_orbitales, swaths_km, N_pix_12, diametros_pupila, ...
    telescope_names, fov_limit, satellite_configs);

%% Mass analysis
disp('Calculating dry and total mass for each configuration...');

for config_idx = 1:size(satellite_configs, 1)
    N_sat        = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    config_name  = sprintf('%dsat_%dtel', N_sat, N_telescopes);

    for detector_idx = 1:3
        for tel_idx = 1:4
            filename = sprintf('HvsDmin/HvsDmin_%s_Detector%d.csv', config_name, detector_idx);
            if ~exist(filename, 'file')
                warning('File not found: %s. Skipping.', filename);
                continue;
            end

            datos     = readtable(filename);
            diam_col  = sprintf('Dmin_%s', telescope_names{tel_idx});
            if ~ismember(diam_col, datos.Properties.VariableNames), continue; end

            alturas_km   = datos.Altura_km;
            diametros_mm = datos.(diam_col);
            valid        = ~isnan(diametros_mm);
            if ~any(valid), continue; end

            h_valid = alturas_km(valid);
            d_valid = diametros_mm(valid);

            [masa_seca, tabla_final] = calcularMasaSeca(h_valid, d_valid, N_telescopes);

            writetable(tabla_final, sprintf('Masa_seca/MasaSeca_%s_Detector%d_%s.csv', ...
                config_name, detector_idx, telescope_names{tel_idx}));

            [total_mass, fuel_mass, num_impulses, total_delta_v] = ...
                calcularMasaTotal(h_valid, masa_seca, tabla_final.Sup_media);

            tabla_final.Masa_combustible = fuel_mass;
            tabla_final.Masa_total       = total_mass;
            tabla_final.Num_impulsos     = num_impulses;
            tabla_final.DeltaV_total     = total_delta_v;

            writetable(tabla_final, sprintf('Masa_total/MasaTotal_%s_Detector%d_%s.csv', ...
                config_name, detector_idx, telescope_names{tel_idx}));

            % Summary text report
            fid = fopen(sprintf('Masa_total/Resumen_%s_Detector%d_%s.txt', ...
                config_name, detector_idx, telescope_names{tel_idx}), 'w');
            fprintf(fid, '=== MASS ANALYSIS ===\n\n');
            fprintf(fid, 'Config: %d satellite(s), %d telescope(s)\n', N_sat, N_telescopes);
            fprintf(fid, 'Detector: %d\n', detector_idx);
            fprintf(fid, 'Telescope: %s\n\n', telescope_names{tel_idx});
            valid_mass = ~isnan(total_mass);
            if any(valid_mass)
                vm = total_mass(valid_mass);
                vh = h_valid(valid_mass);
                [min_m, min_i] = min(vm);
                fprintf(fid, 'Total mass statistics:\n');
                fprintf(fid, '  Min: %.2f kg at %d km\n', min_m, vh(min_i));
                fprintf(fid, '  Max: %.2f kg\n', max(vm));
                fprintf(fid, '  Mean: %.2f kg\n', mean(vm));
                fprintf(fid, '  Altitude range: %d - %d km\n', min(vh), max(vh));
            else
                fprintf(fid, 'No valid total mass values.\n');
            end
            fclose(fid);

            % Orbital parameters plot
            fig1 = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
            subplot(3,1,1);
            plot(h_valid, num_impulses, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            title(sprintf('Number of Impulses — %d sat, %d tel, Det %d, %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{tel_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Orbital Altitude (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('Number of Impulses', 'Interpreter', 'latex', 'FontSize', 11); grid on;

            subplot(3,1,2);
            plot(h_valid, total_delta_v, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            title(sprintf('$\\Delta V$ Total — %d sat, %d tel, Det %d, %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{tel_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Orbital Altitude (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('$\Delta V$ Total (m/s)', 'Interpreter', 'latex', 'FontSize', 11); grid on;

            subplot(3,1,3);
            plot(h_valid, fuel_mass, 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            title(sprintf('Fuel Mass — %d sat, %d tel, Det %d, %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{tel_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Orbital Altitude (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('Fuel Mass (kg)', 'Interpreter', 'latex', 'FontSize', 11); grid on;

            saveas(fig1, sprintf('Masa_total/OrbitalParams_%s_Detector%d_%s.png', ...
                config_name, detector_idx, telescope_names{tel_idx}));
            close(fig1);

            % Mass breakdown plot
            fig2 = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
            hold on;
            plot(h_valid, fuel_mass,  'b-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Fuel Mass');
            plot(h_valid, masa_seca,  'g-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Dry Mass');
            plot(h_valid, total_mass, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Total Mass');
            if any(valid_mass)
                vm = total_mass(valid_mass);
                vh = h_valid(valid_mass);
                [min_m, min_i] = min(vm);
                plot(vh(min_i), min_m, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
                plot([vh(min_i) vh(min_i)], [0 min_m], 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                text(vh(min_i), min_m*1.1, sprintf('Min: %.1f kg\nAlt: %d km', min_m, vh(min_i)), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'FontWeight', 'bold', 'Interpreter', 'latex', 'FontSize', 10, ...
                    'BackgroundColor', 'white', 'EdgeColor', 'black');
            end
            title(sprintf('Mass Breakdown — %d sat, %d tel, Det %d, %s', ...
                N_sat, N_telescopes, detector_idx, telescope_names{tel_idx}), ...
                'Interpreter', 'latex', 'FontSize', 12);
            xlabel('Orbital Altitude (km)', 'Interpreter', 'latex', 'FontSize', 11);
            ylabel('Mass (kg)', 'Interpreter', 'latex', 'FontSize', 11);
            legend('Location', 'best', 'Interpreter', 'latex'); grid on; hold off;

            saveas(fig2, sprintf('Masa_total/MassBreakdown_%s_Detector%d_%s.png', ...
                config_name, detector_idx, telescope_names{tel_idx}));
            close(fig2);
        end
    end
end
disp('Mass analysis complete.');

%% Optimum configuration analysis
OptimumConfigurationAnalysis(satellite_configs, telescope_names, GSD, N_pix_12, N_pix_3);
