SNRCal(2.01e-6, 15e-6,0.6,0.8,30,0.3,500,120)

function SNR_table = SNRCal(lambda_c, pixel_size, eta, tau, GSD, r_obs, altitudes, diameters)

% === Constantes fisicas ===
h = 6.626e-34;
c = 3e8;
bandwidth_m = 20e-3; % 20 nm
TDI = 1;
rad_ref = 100;

% === Ruidos del sistema (en e-) ===
Ndark = 50;
Nread = 100;
Npreamp = 5;
Nvideo = 10;
Njitter = 5;
Nemc = 5;
Nquant = 2;
Nnonlin = 2;

%% === Inicialización ===
SNR_table = zeros(length(altitudes), length(diameters));

for i = 1:length(altitudes)
    h_orb = altitudes(i) * 1e3;
    
    % Focal length basado en GSD, altura y pixel size
    focal_length = (pixel_size * h_orb) / GSD;
    
    % Velocidad sobre el suelo
    v_orb = sqrt(3.986e14 / (6371e3 + h_orb));
    
    % Tiempo de integración
    integration_time = GSD / v_orb;
    
    for j = 1:length(diameters)
        D = diameters(j) / 1000; % Diámetro de la pupila en metros
        
        % Irradiancia en el detector (W/m²)
        irradiance = (pi * tau * bandwidth_m * rad_ref) / (1 + 4 * (focal_length /sqrt(D^2*(1-r_obs^2)))^2)
        
        % Área del píxel
        pixel_area = pixel_size^2;
        
        % Cálculo de Ne (electrones generados)
        Ne = (irradiance * pixel_area * eta * lambda_c * TDI * integration_time) / (h * c)
        
        % Ruido total
        N_total = sqrt(Ndark^2 + Nread^2 + Npreamp^2 + Nvideo^2 + ...
            Njitter^2 + Nemc^2 + Nquant^2 + Nnonlin^2 + Ne)
        
        % SNR
        SNR = Ne / N_total;
        SNR_table(i,j) = SNR
    end
end

%% Crear matriz para heatmap (NaN para valores que no cumplen requisito)
SNR_heatmap = SNR_table;
SNR_heatmap(SNR_table < SNR_Requirement) = NaN;

%% Generar heatmap
fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600]);
h = heatmap(SNR_heatmap);

% Configurar etiquetas
h.XDisplayLabels = repmat({''}, 1, length(diameters));
h.YDisplayLabels = repmat({''}, 1, length(altitudes));

for i = 1:50:length(diameters)
    h.XDisplayLabels{i} = sprintf('%d', diameters(i));
end

for i = 1:100:length(altitudes)
    h.YDisplayLabels{i} = sprintf('%d', altitudes(i));
end

h.Title = sprintf('SNR - Det%d, %s, GSD=%dm, R=%.1f, %s', ...
    detector_num, telescope_name, GSD, r_obs, lambda_name);
h.XLabel = 'Diámetro Pupila (mm)';
h.YLabel = 'Altura Orbital (km)';

h.Colormap = parula(256);
h.ColorLimits = [0, 2000];
h.MissingDataColor = [1, 1, 1]; % Blanco para NaN
h.GridVisible = 'off';

%% Guardar archivos
filename_base = sprintf('SNR_Det%d_%s_GSD%d_R%.1f_%s', ...
    detector_num, strrep(telescope_name, ' ', ''), GSD, r_obs*10, lambda_name);

% Guardar heatmap
saveas(fig, fullfile('SNR', [filename_base '.png']));
close(fig);

% Guardar datos en archivo MAT
SNR_data.SNR_values = SNR_table;
SNR_data.altitudes = altitudes;
SNR_data.diameters = diameters;
SNR_data.detector_num = detector_num;
SNR_data.telescope_name = telescope_name;
SNR_data.lambda_name = lambda_name;
SNR_data.GSD = GSD;
SNR_data.r_obs = r_obs;
SNR_data.SNR_Requirement = SNR_Requirement;
SNR_data.lambda_c = lambda_c;
SNR_data.pixel_size = pixel_size;
SNR_data.eta = eta;
SNR_data.tau = tau;

save(fullfile('SNR', [filename_base '.mat']), 'SNR_data');

end