function MTFfunction(lambda, pixel_size, MTF_Detector, MTF_alineamiento, GSD, R, alturas_orbitales, diametros_pupila, filename_prefix, telescope_name, detector_idx,MTF_req)

%% Parámetros MTF adicionales
MTF_aberraciones = 0.95;
MTF_fabricacion = 0.98;
MTF_vibraciones = 0.99;
MTF_Termoelastico = 0.95;
MTF_Margen = 0.9;
MTF_resto = MTF_Margen * MTF_vibraciones * MTF_fabricacion * MTF_Termoelastico * MTF_aberraciones * MTF_Detector;

%% Inicializar matriz MTF
MTF_total = zeros(length(alturas_orbitales), length(diametros_pupila));

%% Cálculo MTF para cada combinación altura-diámetro
for i = 1:length(alturas_orbitales)
    for j = 1:length(diametros_pupila)
        altura_m = alturas_orbitales(i) * 1e3;
        diametro_m = diametros_pupila(j) * 1e-3;
        
        % Distancia focal
        distancia_focal = (altura_m * pixel_size) / GSD;
        
        % Frecuencias de corte
        f_co = diametro_m / lambda;
        f_Nyquist = 1 / (2 * pixel_size);
        f_x = f_Nyquist * distancia_focal;
        
        if f_x > f_co
            f_x = f_co;
        end
        
        % Cálculo para sistema con obscuración central (ecuaciones 6-10 a 6-15)
        X = f_x / f_co;
        Y = X / R;
        
        % Cálculo de alpha según la ecuación 6-12
        if (1 + R^2 - 4*X^2) / (2*R) >= -1 && (1 + R^2 - 4*X^2) / (2*R) <= 1
            alpha = acos((1 + R^2 - 4*X^2) / (2*R));
        else
            alpha = 0;
        end
        
        % Cálculo de A (ecuación 6-13)
        if 0 <= X && X <= 1
            A = (2/pi) * (acos(X) - X*sqrt(1 - X^2));
        else
            A = 0;
        end
        
        % Cálculo de B (ecuación 6-14)
        if 0 <= Y && Y <= 1
            B = (2*R^2/pi) * (acos(Y) - Y*sqrt(1 - Y^2));
        else
            B = 0;
        end
        
        % Cálculo de C (ecuaciones 6-15)
        if 0 < X && X <= (1-R)/2
            C = -2*R^2;
        elseif (1-R)/2 < X && X < (1+R)/2
            C = (2*R/pi)*sin(alpha) + ((1+R^2)/pi)*alpha - ((2*(1-R^2))/pi)*atan((1+R)/(1-R)*tan(alpha/2)) - 2*R^2;
        elseif X >= (1+R)/2
            C = 0;
        else
            C = 0;
        end
        
        % Cálculo del OTF difracción según ecuación 6-10
        MTF_difraccion = (A + B + C) / (1 - R^2);
        
        % Si el valor es negativo o NaN, ajustarlo a 0
        if isnan(MTF_difraccion) || MTF_difraccion < 0
            MTF_difraccion = 0;
        end
        
        % MTF total
        MTF_value = MTF_difraccion * MTF_alineamiento * MTF_resto;
        
        % Aplicar requirement: valores < 0.25 se convierten a NaN
        if MTF_value < MTF_req
            MTF_total(i, j) = NaN;
        else
            MTF_total(i, j) = MTF_value;
        end
    end
end

%% Generar heatmap con paleta parula y límites optimizados
fig = figure('Visible', 'off', 'Position', [100, 100, 900, 700]);
h = imagesc(diametros_pupila, alturas_orbitales, MTF_total);

% Aplicar técnica AlphaData para mostrar NaN como blanco
set(h, 'AlphaData', ~isnan(MTF_total));
set(gca, 'Color', [1 1 1]); % Fondo blanco

% Configurar paleta y límites
colormap(parula);
clim([0.25 0.35]); % Límites optimizados para MTF

% Configuración de ejes
axis xy;

ylabel('Altura Orbital (km)', 'FontSize', 12, 'Interpreter', 'latex');
xlabel('Di\''ametro de Pupila (mm)', 'FontSize', 12, 'Interpreter', 'latex');
title(sprintf('MTF - %s - Detector %d\n($\\lambda$ = %.2f $\\mu$m, Blanco: MTF $< 0.25$)', ...
    telescope_name, detector_idx, lambda*1e6), 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Colorbar
cb = colorbar;
cb.Label.String = 'MTF';
cb.Label.FontSize = 12;
cb.Label.Interpreter = 'latex';


% Configurar ticks
xticks_vals = diametros_pupila(1:50:end);
yticks_vals = alturas_orbitales(1:5:end);
set(gca, 'XTick', xticks_vals, 'YTick', yticks_vals);

% Guardar heatmap
saveas(fig, [filename_prefix '_heatmap.jpg']);
close(fig);

%% Guardar datos en CSV
% Crear tabla con nombres de filas y columnas
row_names = cell(length(alturas_orbitales), 1);
for i = 1:length(alturas_orbitales)
    row_names{i} = sprintf('Alt_%d', alturas_orbitales(i));
end

col_names = cell(1, length(diametros_pupila));
for j = 1:length(diametros_pupila)
    col_names{j} = sprintf('Diam_%d', diametros_pupila(j));
end

MTF_table = array2table(MTF_total, 'RowNames', row_names, 'VariableNames', col_names);
writetable(MTF_table, [filename_prefix '_resultados.csv'], 'WriteRowNames', true);

%% Generar archivo de estadísticas TXT
fid = fopen([filename_prefix '_estadisticas.txt'], 'w');
fprintf(fid, '=== ANÁLISIS MTF ===\n\n');
fprintf(fid, 'Telescopio: %s\n', telescope_name);
fprintf(fid, 'Detector: %d\n', detector_idx);
fprintf(fid, 'Longitud de onda: %.2f μm\n', lambda*1e6);
fprintf(fid, 'GSD: %d m\n', GSD);
fprintf(fid, 'Requirement MTF: ≥ 0.25\n\n');

% Estadísticas de valores válidos (no NaN)
valid_mtf = MTF_total(~isnan(MTF_total));
if ~isempty(valid_mtf)
    fprintf(fid, 'Estadísticas de valores válidos (MTF ≥ 0.25):\n');
    fprintf(fid, 'Valores válidos: %d de %d (%.1f%%)\n', length(valid_mtf), numel(MTF_total), 100*length(valid_mtf)/numel(MTF_total));
    fprintf(fid, 'MTF mínimo: %.4f\n', min(valid_mtf));
    fprintf(fid, 'MTF máximo: %.4f\n', max(valid_mtf));
    fprintf(fid, 'MTF promedio: %.4f\n', mean(valid_mtf));
    fprintf(fid, 'MTF mediana: %.4f\n', median(valid_mtf));
else
    fprintf(fid, 'No hay valores válidos que cumplan el requirement MTF ≥ 0.25\n');
end

fprintf(fid, '\nParámetros del análisis:\n');
fprintf(fid, 'Rango alturas: %d - %d km\n', min(alturas_orbitales), max(alturas_orbitales));
fprintf(fid, 'Rango diámetros: %d - %d mm\n', min(diametros_pupila), max(diametros_pupila));
fprintf(fid, 'Obscuración central: %.1f\n', R);
fprintf(fid, 'MTF alineamiento: %.2f\n', MTF_alineamiento);
fprintf(fid, 'MTF detector: %.2f\n', MTF_Detector);
fprintf(fid, 'Tamaño pixel: %.1f μm\n', pixel_size*1e6);

fclose(fid);

fprintf('MTF calculado para %s, Detector %d\n', telescope_name, detector_idx);

end