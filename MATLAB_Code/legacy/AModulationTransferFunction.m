%% Parámetros de entrada
lambda = 2.01e-6; % Longitud de onda en metros
pixel_size = 18e-6; % Tamaño del píxel en metros
GSD = 50; % GSD en metros
% Parámetro de obscuración central
R = 0; % Radio de obscuración (d_obs/D_o)

MTF_alineamiento = 0.90; % MTF de alineamiento del telescopio (Refractivo)
CTF = 0.6;
MTF_Detector = 0.6;%0.636 * CTF;
MTF_aberraciones = 0.95;
MTF_fabricacion = 0.98;
MTF_vibraciones = 0.99;
MTF_Termoelastico = 0.95;
MTF_Margen = 0.9;

MTF_resto = MTF_Margen*MTF_vibraciones*MTF_fabricacion*MTF_Termoelastico*MTF_aberraciones*MTF_Detector;

%% Definir rangos

alturas_orbitales = 410; % Alturas orbitales en km
diametros_pupila = 25; % Diámetros pupilares en mm

%% Inicializar matriz MTF

MTF_total = zeros(length(alturas_orbitales), length(diametros_pupila));

% Cálculo MTF

for i = 1:length(alturas_orbitales)
    for j = 1:length(diametros_pupila)
        altura_m = alturas_orbitales(i) * 1e3;
        diametro_m = diametros_pupila(j) * 1e-3;
        distancia_focal = (altura_m * pixel_size) / GSD;

        f_co = diametro_m / lambda;
        f_Nyquist = 1 / (2 * pixel_size);
        f_x = f_Nyquist * distancia_focal;

        if f_x > f_co
            f_x = f_co;
        end

        % Cálculo para sistema con obscuración central
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
        
        MTF_total(i, j) = MTF_difraccion * MTF_alineamiento * MTF_resto;
    end
end

%% Preparar heatmap

MTF_visual = MTF_total;

MTF_visual(MTF_visual < 0.25) = NaN; % Valores <0.25 a NaN

% Generar heatmap

figure('Position', [100, 100, 900, 700]);

imagesc(diametros_pupila, alturas_orbitales, MTF_visual);

colormap('parula');

set(gca, 'YDir', 'normal');

set(gca, 'Color', 'white'); % Fondo blanco

% Eliminar valores NaN del colormap

cmap = colormap;

cmap(1,:) = [1 1 1]; % Primer color (blanco) para NaN

colormap(cmap);

% Añadir los valores como texto en blanco, sólo para MTF > 0.25

[X, Y] = meshgrid(diametros_pupila, alturas_orbitales);

for i = 1:length(alturas_orbitales)
    for j = 1:length(diametros_pupila)
        % Solo mostrar texto para valores mayores o iguales a 0.25
        if MTF_total(i, j) >= 0.25
            %text(X(i,j), Y(i,j), sprintf('%.2f', MTF_total(i,j)), ...
            %'HorizontalAlignment', 'center', ...
            %'VerticalAlignment', 'middle', ...
            %'Color', 'white', ...
            %'FontWeight', 'bold');
        end
    end
end

% Etiquetas

title(sprintf('Heatmap de MTF Total - \\lambda = %.2f \\mum - Pixel = %.1f \\mum - GSD = %d m - R = %.2f', lambda*1e6, pixel_size*1e6, GSD, R));

xlabel('Diámetro de pupila [mm]');

ylabel('Altitud orbital [km]');

colorbar;

%% Generar archivo de texto con todos los valores

fileID = fopen('resumen_MTF.txt','w');

fprintf(fileID, '=== RESUMEN DE SIMULACIÓN DE MTF ===\n\n');

fprintf(fileID, '>> Parámetros iniciales:\n');

fprintf(fileID, ' - Tamaño de píxel : %.2f µm\n', pixel_size * 1e6);

fprintf(fileID, ' - GSD : %.1f m\n', GSD);

fprintf(fileID, ' - Longitud de onda : %.3f µm\n', lambda * 1e6);
fprintf(fileID, ' - Radio de obscuración (R) : %.2f\n', R);
fprintf(fileID, ' - MTF alineamiento (Refractivo) : %.2f\n', MTF_alineamiento);

fprintf(fileID, ' - MTF detector : %.4f\n', MTF_Detector);

fprintf(fileID, '>> Tabla de MTF (Altitud [km] en columnas, Diámetro [mm] en filas):\n');

fprintf(fileID, ' D \\ h');

for i = 1:length(alturas_orbitales)
    fprintf(fileID, '%7d', alturas_orbitales(i));
end

fprintf(fileID, '\n');

for j = 1:length(diametros_pupila)
    fprintf(fileID, '%7.1f', diametros_pupila(j));
    for i = 1:length(alturas_orbitales)
        fprintf(fileID, '%7.2f', MTF_total(i, j));
    end
    fprintf(fileID, '\n');
end

fclose(fileID);

disp('Archivo resumen_MTF.txt generado con éxito.');

%% tabla

% Verificar dimensiones

[rows, cols] = size(MTF_total);

fprintf('Matriz MTF: %d filas x %d columnas\n', rows, cols);

fprintf('Alturas: %d elementos\n', length(alturas_orbitales));

fprintf('Diámetros: %d elementos\n', length(diametros_pupila));

% Crear nombres de variables válidos para la tabla
alt_names = cell(1, length(alturas_orbitales));

for i = 1:length(alturas_orbitales)
    alt_names{i} = sprintf('Alt_%d', alturas_orbitales(i));
end

diam_names = cell(1, length(diametros_pupila));

for i = 1:length(diametros_pupila)
    diam_names{i} = sprintf('Diam_%d', diametros_pupila(i));
end

% Crear tabla con los resultados (asegurando que las dimensiones coincidan)

if rows == length(alturas_orbitales) && cols == length(diametros_pupila)
    % La matriz está orientada como [alturas, diámetros]
    MTF_table = array2table(MTF_total, 'VariableNames', diam_names);
    MTF_table.Properties.RowNames = alt_names;
else
    % La matriz podría estar transpuesta, verificar y ajustar
    MTF_table = array2table(MTF_total', 'VariableNames', alt_names);
    MTF_table.Properties.RowNames = diam_names;
end

% Mostrar la tabla en una figura separada

fig_table = figure('Name', 'Tabla de MTF Total', 'Position', [150, 150, 800, 600]);

uit = uitable(fig_table);

uit.Data = MTF_table{:,:};
uit.ColumnName = MTF_table.Properties.VariableNames;

uit.RowName = MTF_table.Properties.RowNames;

uit.Position = [20 20 760 560];

uit.ColumnWidth = {50};

% Exportar tabla a formato Excel

writetable(MTF_table, 'MTF_resultados.xlsx', 'WriteRowNames', true);
writetable(MTF_table, 'MTF_resultados.csv', 'WriteRowNames', true);

disp('Tabla de resultados exportada a MTF_resultados.xlsx');
