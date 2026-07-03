%% === Parámetros configurables ===
pixel_size = 15e-6;         % Tamaño de píxel en metros (d_x = d_y)
GSD = 30;                   % Ground Sampling Distance en metros
r_obs = 0.3;
% === Constantes físicas ===
h = 6.626e-34;              % Constante de Planck (J·s)
c = 3e8;                    % Velocidad de la luz (m/s)
lambda_c = 0.76e-6;         % Longitud de onda central (m) (más restrictiva O2 A-band)
bandwidth_m = 20e-3;        % Ancho de banda espectral Δλ (20 nm)
eta = 0.8;                  % Eficiencia cuántica
tau = 0.8;                  % Transmisión óptica del sistema
TDI = 1;                    % Número de etapas TDI
rad_ref = 100;              % Radiancia de referencia [W/(m²·sr·µm)]

% === Ruidos del sistema   (en e⁻ RMS) ===
Ndark   = 50;   % Ruido térmico alto por mal enfriamiento
Nread   = 100;  % Sistema de lectura básico o rápido
Npreamp = 5;    % Electrónica con ganancia moderada
Nvideo  = 10;   % Línea de vídeo sin filtrado
Njitter = 5;    % Plataforma con jitter notable
Nemc    = 5;    % Alta interferencia electromagnética
Nquant  = 2;    % ADC de baja resolución
Nnonlin = 2;    % Mala calibración o no linealidad no compensada

%% === Rango de variables ===
altitudes = 500;     % en km
diameters = 120;     % en mm


%% === Inicialización ===
SNR_table = zeros(length(altitudes), length(diameters));

for i = 1:length(altitudes)
    h_orb = altitudes(i) * 1e3;  % convertir a metros
    
    % Focal length basado en GSD, altura y pixel size
    focal_length = (pixel_size * h_orb) / GSD;

    % Velocidad sobre el suelo
    v_orb = sqrt(3.986e14 / (6371e3 + h_orb));
    
    % Tiempo de integración
    integration_time = GSD / v_orb;

    for j = 1:length(diameters)
        D = diameters(j) / 1000;     % Diámetro de la pupila en metros

        % Irradiancia en el detector (W/m²)
        irradiance = (pi * tau * bandwidth_m * rad_ref) / (1 + 4 * (focal_length /sqrt(D^2*(1-r_obs^2)))^2);

        % Área del píxel
        pixel_area = pixel_size^2;

        % Cálculo de Ne (electrones generados)
        Ne = (irradiance * pixel_area * eta * lambda_c * TDI * integration_time) / (h * c);

        % Ruido total
        N_total = sqrt(Ndark^2 + Nread^2 + Npreamp^2 + Nvideo^2 + ...
                       Njitter^2 + Nemc^2 + Nquant^2 + Nnonlin^2 + Ne);

        % SNR
        SNR = Ne / N_total;
        SNR_table(i,j) = SNR;
    end
end

%% === Crear copia para visualización ===
SNR_visual = SNR_table;
% Poner NaN en valores menores a 400 para que aparezcan en blanco
SNR_visual(SNR_visual < 400) = NaN;



%% === Crear el heatmap ===
figure('Position', [100, 100, 900, 700]);
h = imagesc(diameters, altitudes, SNR_visual);
colorbar;
title(sprintf('Heatmap de SNR - \\lambda = %.2f \\mum - Pixel = %.1f \\mum - GSD = %d m - Umbral SNR = 400- Topt= %.2f', lambda_c*1e6, pixel_size*1e6, GSD, tau));
xlabel('Diámetro de pupila [mm]');
ylabel('Altura orbital [km]');
colormap('parula');
set(gca, 'YDir', 'normal');
set(gca, 'Color', 'white'); % Fondo blanco para valores NaN

% Eliminar valores NaN del colormap
cmap = colormap;
cmap(1,:) = [1 1 1]; % Primer color (blanco) para NaN
colormap(cmap);

% Añadir los valores como texto en blanco, sólo para SNR ≥ 400
[X, Y] = meshgrid(diameters, altitudes);
for i = 1:length(altitudes)
    for j = 1:length(diameters)
        % Solo mostrar texto para valores mayores o iguales a 400
        if SNR_table(i, j) >= 400
            %text(X(i,j), Y(i,j), sprintf('%.0f', SNR_table(i,j)), ...
             %   'HorizontalAlignment', 'center', ...
              %  'VerticalAlignment', 'middle', ...
               % 'Color', 'white', ...
                %'FontWeight', 'bold');
        end
    end
end

%% === Crear archivo de resumen txt ===
fileID = fopen('resumen_SNR.txt','w');

% === Escribir encabezado con parámetros iniciales ===
fprintf(fileID, '=== RESUMEN DE SIMULACIÓN DE SNR ===\n\n');

fprintf(fileID, '>> Parámetros iniciales:\n');
fprintf(fileID, '  - Tamaño de píxel        : %.2f µm\n', pixel_size * 1e6);
fprintf(fileID, '  - GSD                    : %.1f m\n', GSD);
fprintf(fileID, '  - Longitud de onda       : %.3f µm\n', lambda_c * 1e6);
fprintf(fileID, '  - Ancho de banda         : %.3f µm\n', bandwidth_m);
fprintf(fileID, '  - Eficiencia cuántica    : %.2f\n', eta);
fprintf(fileID, '  - Tipo de telescopio     : Refractivo\n');
fprintf(fileID, '  - Transmisión óptica     : %.2f\n', tau);
fprintf(fileID, '  - TDI                    : %d\n', TDI);
fprintf(fileID, '  - Radiancia de referencia: %.1f W/m²·sr·µm\n', rad_ref);



fprintf(fileID, '>> Ruidos (caso peor):\n');
fprintf(fileID, '  - Ndark   = %d e⁻ RMS\n', Ndark);
fprintf(fileID, '  - Nread   = %d e⁻ RMS\n', Nread);
fprintf(fileID, '  - Npreamp = %d e⁻ RMS\n', Npreamp);
fprintf(fileID, '  - Nvideo  = %d e⁻ RMS\n', Nvideo);
fprintf(fileID, '  - Njitter = %d e⁻ RMS\n', Njitter);
fprintf(fileID, '  - Nemc    = %d e⁻ RMS\n', Nemc);
fprintf(fileID, '  - Nquant  = %d e⁻ RMS\n', Nquant);
fprintf(fileID, '  - Nnonlin = %.1f e⁻ RMS\n\n', Nnonlin);

% === Escribir tabla de SNR ===
fprintf(fileID, '>> Tabla de SNR (Altitud [km] en columnas, Diámetro [mm] en filas):\n');

% Imprimir encabezado con las altitudes como columnas
fprintf(fileID, '  D \\ h');
for i = 1:length(altitudes)
    fprintf(fileID, '%7d', altitudes(i));
end
fprintf(fileID, '\n');

% Imprimir los valores de SNR para cada combinación de diámetro y altitud
for j = 1:length(diameters)
    fprintf(fileID, '%7.1f', diameters(j));
    for i = 1:length(altitudes)
        fprintf(fileID, '%7.0f', SNR_table(i, j));
    end
    fprintf(fileID, '\n');
end

fclose(fileID);

disp('Archivo resumen_SNR.txt generado con éxito.');


%% TABLA

% Verificar dimensiones
[rows, cols] = size(SNR_table);
fprintf('Matriz SNR: %d filas x %d columnas\n', rows, cols);
fprintf('Alturas: %d elementos\n', length(altitudes));
fprintf('Diámetros: %d elementos\n', length(diameters));

% Crear nombres de variables válidos para la tabla
alt_names = cell(1, length(altitudes));
for i = 1:length(altitudes)
    alt_names{i} = sprintf('Alt_%d', altitudes(i));
end

diam_names = cell(1, length(diameters));
for i = 1:length(diameters)
    diam_names{i} = sprintf('Diam_%d', diameters(i));
end

% Crear tabla con los resultados (asegurando que las dimensiones coincidan)
if rows == length(altitudes) && cols == length(diameters)
    % La matriz está orientada como [alturas, diámetros]
    SNR_table_data = array2table(SNR_table, 'VariableNames', diam_names);
    SNR_table_data.Properties.RowNames = alt_names;
else
    % La matriz podría estar transpuesta, verificar y ajustar
    SNR_table_data = array2table(SNR_table', 'VariableNames', alt_names);
    SNR_table_data.Properties.RowNames = diam_names;
end

% Mostrar la tabla en una figura separada
fig_table = figure('Name', 'Tabla de SNR', 'Position', [150, 150, 800, 600]);
uit = uitable(fig_table);
uit.Data = SNR_table_data{:,:};
uit.ColumnName = SNR_table_data.Properties.VariableNames;
uit.RowName = SNR_table_data.Properties.RowNames;
uit.Position = [20 20 760 560];
uit.ColumnWidth = {50};

% Exportar tabla a formato Excel
writetable(SNR_table_data, 'SNR_resultados.xlsx', 'WriteRowNames', true);
disp('Tabla de resultados exportada a SNR_resultados.xlsx');