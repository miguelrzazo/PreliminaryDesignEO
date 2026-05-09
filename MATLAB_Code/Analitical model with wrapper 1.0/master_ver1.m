%%
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'lib'));

%% Wavelengths (m)
lambda = [1.61e-6, 2.01e-6, 0.76e-6];   % [CO2 band 1, CO2 band 2, O2 A-band]

GSD = 80;                                % Ground sampling distance (m)

%% Orbital and mission ranges
alturas_orbitales = 200:10:1000;         % Orbital altitudes (km)
diametros_pupila  = 20:1:600;            % Pupil diameters (mm)
SNR_Requirement   = 400;
MTF_Requirement   = 0.25;
Cov_Requirement   = 7;                   % days
max_detectores    = 3;
swaths_km         = 20:10:200;           % Swath widths to evaluate (km)
solapamiento      = 0.05;                % Swath overlap fraction
cobertura_nubes   = 1/6;                 % Fraction of passes blocked by cloud
LTAN_hour         = 6;                   % Local Time of Ascending Node (Dawn-Dusk)

satellite_configs = [
    %1, 1;
    2, 1;
    1, 2;
    3, 1;
    4, 1;
];

%% Detectors — CO2 bands (lambda 1 & 2): detectors 1-3
%                    [CAPYORK, H2RG, SATURN VISIR]
eta_12        = [0.85,  0.70,  0.60];
MTF_det_12    = [0.45,  0.45,  0.50];
N_pix_12      = [1200,  2048,  1000];
pixel_size_12 = [15,    18,    30  ] * 1e-6;

%% Detectors — O2 A-band (lambda 3): detectors 4-6
%                    [CMOS, H2RG, SATURN VISIR]
eta_3         = [0.35,  0.70,  0.60];
MTF_det_3     = [0.36,  0.50,  0.45];
N_pix_3       = [512,   2048,  1000];
pixel_size_3  = [25,    18,    30  ] * 1e-6;

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
disp('================================================================');
disp('Starting coverage analysis with RevisitCalc');
disp('================================================================');

for config_idx = 1:size(satellite_configs, 1)
    N_sat        = satellite_configs(config_idx, 1);
    N_telescopes = satellite_configs(config_idx, 2);
    fprintf('Processing config: %d sat, %d tel\n', N_sat, N_telescopes);

    for detector_id = 1:3
        for telescope_id = 1:4
            fprintf('  -> Detector %d, Telescope %s...\n', detector_id, telescope_names{telescope_id});
            try
                CoverageRevisitCalc(GSD, alturas_orbitales, swaths_km, Cov_Requirement, ...
                    N_pix_12(detector_id), fov_limit(telescope_id), N_sat, N_telescopes, ...
                    detector_id, telescope_names{telescope_id}, solapamiento, cobertura_nubes, max_detectores);
            catch ME
                fprintf('    -> Error: %s\n', ME.message);
                continue;
            end
        end
    end
    fprintf('Done: %d sat, %d tel\n\n', N_sat, N_telescopes);
end

disp('================================================================');
disp('Coverage analysis complete.');
disp('================================================================');

%% Cross-data analysis
CrossDataFunction(alturas_orbitales, swaths_km, diametros_pupila, telescope_names, fov_limit, satellite_configs);

%% Comprehensive final analysis
params.satellite_configs = satellite_configs;
params.telescope_names   = telescope_names;
params.alturas_orbitales = alturas_orbitales;
params.swaths_km         = swaths_km;
params.diametros_pupila  = diametros_pupila;
params.Cov_Requirement   = Cov_Requirement;
params.GSD               = GSD;
params.N_pix_12          = N_pix_12;
params.fov_limit         = fov_limit;
params.LTAN_hour         = LTAN_hour;

ComprehensiveAnalysisFunction(params);

disp('================================================================');
disp('Script complete. All results generated.');
disp('================================================================');
