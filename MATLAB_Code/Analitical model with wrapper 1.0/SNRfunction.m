function SNRfunction(lambda_c, pixel_size, eta, tau, GSD, r_obs, altitudes, diameters, filename_prefix, telescope_name, detector_idx,SNR_req)

%% Constantes físicas
h = 6.626e-34; % Constante de Planck (J·s)
c = 3e8; % Velocidad de la luz (m/s)
bandwidth_m = 20e-3; % Ancho de banda espectral Δλ (20 nm)
TDI = 1; % Número de etapas TDI
rad_ref = 100; % Radiancia de referencia [W/(m²·sr·μm)]

%% Ruidos del sistema (en e⁻ RMS)
Ndark = 50; % Ruido térmico alto por mal enfriamiento
Nread = 100; % Sistema de lectura básico o rápido
Npreamp = 5; % Electrónica con ganancia moderada
Nvideo = 10; % Línea de vídeo sin filtrado
Njitter = 5; % Plataforma con jitter notable
Nemc = 5; % Alta interferencia electromagnética
Nquant = 2; % ADC de baja resolución
Nnonlin = 2; % Mala calibración o no linealidad no compensada

%% Inicialización
SNR_table = zeros(length(altitudes), length(diameters));

%% Cálculo SNR para cada combinación altura-diámetro
for i = 1:length(altitudes)
    h_orb = altitudes(i) * 1e3; % convertir a metros
    
    % Focal length basado en GSD, altura y pixel size
    focal_length = (pixel_size * h_orb) / GSD;
    
    % Velocidad sobre el suelo
    v_orb = sqrt(3.986e14 / (6371e3 + h_orb));
    
    % Tiempo de integración
    integration_time = GSD / v_orb;
    
    for j = 1:length(diameters)
        D = diameters(j) / 1000; % Diámetro de la pupila en metros
        
        % Irradiancia en el detector (W/m²)
        irradiance = (pi * tau * bandwidth_m * rad_ref) / (1 + 4 * (focal_length / sqrt(D^2*(1-r_obs^2)))^2);
        
        % Área del píxel
        pixel_area = pixel_size^2;
        
        % Cálculo de Ne (electrones generados)
        Ne = (irradiance * pixel_area * eta * lambda_c * TDI * integration_time) / (h * c);
        
        % Ruido total
        N_total = sqrt(Ndark^2 + Nread^2 + Npreamp^2 + Nvideo^2 + ...
            Njitter^2 + Nemc^2 + Nquant^2 + Nnonlin^2 + Ne);
        
        % SNR
        SNR_value = Ne / N_total;
        
        % Aplicar requirement: valores < 400 se convierten a NaN
        if SNR_value < SNR_req
            SNR_table(i, j) = NaN;
        else
            SNR_table(i, j) = SNR_value;
        end
    end
end

%% Generar heatmap con paleta parula y límites optimizados
fig = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
h = imagesc(diameters, altitudes, SNR_table);

% Aplicar técnica AlphaData para mostrar NaN como blanco
set(h, 'AlphaData', ~isnan(SNR_table));
set(gca, 'Color', [1 1 1]); % Fondo blanco

% Configurar paleta y límites
colormap(parula);
caxis([400 2000]); % Límites optimizados para SNR

% Configuración de ejes
axis xy;
xlabel('Di\''ametro de Pupila (mm)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('Altura Orbital (km)', 'FontSize', 12, 'Interpreter', 'latex');
title(sprintf('SNR - %s - Detector %d\n($\\lambda$ = %.2f $\\mu$m, Blanco: SNR $< 400$)', ...
    telescope_name, detector_idx, lambda_c*1e6), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Colorbar
cb = colorbar;
cb.Label.String = 'SNR';
cb.Label.FontSize = 12;
cb.Label.Interpreter = 'latex';

% Configurar ticks
xticks_vals = diameters(1:50:end);
yticks_vals = altitudes(1:5:end);
set(gca, 'XTick', xticks_vals, 'YTick', yticks_vals);

% Guardar heatmap
saveas(fig, [filename_prefix '_heatmap.png']);
close(fig);

%% Guardar datos en CSV
% Crear tabla con nombres de filas y columnas
row_names = cell(length(altitudes), 1);
for i = 1:length(altitudes)
    row_names{i} = sprintf('Alt_%d', altitudes(i));
end

col_names = cell(1, length(diameters));
for j = 1:length(diameters)
    col_names{j} = sprintf('Diam_%d', diameters(j));
end

SNR_table_export = array2table(SNR_table, 'RowNames', row_names, 'VariableNames', col_names);
writetable(SNR_table_export, [filename_prefix '_resultados.csv'], 'WriteRowNames', true);

%% Generar archivo de estadísticas TXT
fid = fopen([filename_prefix '_estadisticas.txt'], 'w');
fprintf(fid, '=== ANÁLISIS SNR ===\n\n');
fprintf(fid, 'Telescopio: %s\n', telescope_name);
fprintf(fid, 'Detector: %d\n', detector_idx);
fprintf(fid, 'Longitud de onda: %.2f μm\n', lambda_c*1e6);
fprintf(fid, 'GSD: %d m\n', GSD);
fprintf(fid, 'Requirement SNR: ≥ 400\n\n');

% Estadísticas de valores válidos (no NaN)
valid_snr = SNR_table(~isnan(SNR_table));
if ~isempty(valid_snr)
    fprintf(fid, 'Estadísticas de valores válidos (SNR ≥ 400):\n');
    fprintf(fid, 'Valores válidos: %d de %d (%.1f%%)\n', length(valid_snr), numel(SNR_table), 100*length(valid_snr)/numel(SNR_table));
    fprintf(fid, 'SNR mínimo: %.2f\n', min(valid_snr));
    fprintf(fid, 'SNR máximo: %.2f\n', max(valid_snr));
    fprintf(fid, 'SNR promedio: %.2f\n', mean(valid_snr));
    fprintf(fid, 'SNR mediana: %.2f\n', median(valid_snr));
else
    fprintf(fid, 'No hay valores válidos que cumplan el requirement SNR ≥ 400\n');
end

fprintf(fid, '\nParámetros del análisis:\n');
fprintf(fid, 'Rango alturas: %d - %d km\n', min(altitudes), max(altitudes));
fprintf(fid, 'Rango diámetros: %d - %d mm\n', min(diameters), max(diameters));
fprintf(fid, 'Obscuración central: %.1f\n', r_obs);
fprintf(fid, 'Transmitancia: %.2f\n', tau);
fprintf(fid, 'Eficiencia cuántica: %.2f\n', eta);
fprintf(fid, 'Tamaño pixel: %.1f μm\n', pixel_size*1e6);
fprintf(fid, 'Radiancia referencia: %d W/(m²·sr·μm)\n', rad_ref);
fprintf(fid, 'Ancho de banda: %.1f nm\n', bandwidth_m*1e9);

fprintf(fid, '\nRuidos del sistema (e⁻ RMS):\n');
fprintf(fid, 'Ruido térmico: %d\n', Ndark);
fprintf(fid, 'Ruido de lectura: %d\n', Nread);
fprintf(fid, 'Ruido preamplificador: %d\n', Npreamp);
fprintf(fid, 'Ruido línea video: %d\n', Nvideo);
fprintf(fid, 'Ruido jitter: %d\n', Njitter);
fprintf(fid, 'Ruido EMC: %d\n', Nemc);
fprintf(fid, 'Ruido cuantización: %d\n', Nquant);
fprintf(fid, 'Ruido no linealidad: %d\n', Nnonlin);

fclose(fid);

fprintf('SNR calculado para %s, Detector %d\n', telescope_name, detector_idx);

end