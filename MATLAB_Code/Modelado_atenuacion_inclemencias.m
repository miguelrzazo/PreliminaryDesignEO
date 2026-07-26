%% Modelado de Atenuación por Inclemencias del Tiempo en Banda X
% Este script simula la velocidad de descarga en banda X bajo diferentes
% condiciones meteorológicas durante un año completo

%% Configuración inicial
clear; clc; close all;

% Establecer semilla para reproducibilidad
rng(42);

% Configurar intérprete LaTeX por defecto
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Parámetros de simulación
dias_simulacion = 365; % Un año completo

% Parámetros de velocidad (Mbps) - [min, max]
% Rangos basados en ITU-R P.838-3 a 8.1 GHz: gamma_R = k*R^alpha, k=0.00425, alpha=1.32
% Atenuación en trayecto oblicuo (elev. 60°, altura lluvia 2 km): L_eff = 2/sin(60°) = 2.31 km
% Lluvia moderada (5 mm/h): ~0.08 dB; Lluvia intensa (30 mm/h): ~0.74 dB
% Margen de enlace disponible: 22.8 dB >> atenuaciones anteriores
VELOCIDADES = containers.Map(...
    {'Cond. favorables', 'Lluvia moderada', 'Lluvia intensa', 'Nubes densas'}, ...
    {[130, 150], [120, 140], [90, 120], [125, 145]});

% Probabilidades de condiciones meteorológicas por estación (%)
% Fuente: NOAA Climate Normals 1991-2020, Fairbanks Intl Airport (USW00026411)
% Días de lluvia real: invierno ~0, primavera ~4, verano ~37, otoño ~9 (total año ~103 precip.)
% Nieve y nubes sin lluvia = "Cond. favorables" a efectos de enlace en banda X
PROBABILIDADES_ESTACION = containers.Map(...
    {'Primavera', 'Verano', 'Otono', 'Invierno'}, ...
    {[76, 12, 3, 9], [52, 36, 7, 5], [64, 22, 5, 9], [82, 1, 0, 17]});

condiciones_nombres = {'Cond. favorables', 'Lluvia moderada', 'Lluvia intensa', 'Nubes densas'};

%% Función para asignar estación
function estacion = asignar_estacion(dia)
    if dia >= 80 && dia <= 172
        estacion = 'Primavera';
    elseif dia >= 173 && dia <= 265
        estacion = 'Verano';
    elseif dia >= 266 && dia <= 355
        estacion = 'Otono';
    else
        estacion = 'Invierno';
    end
end

%% Generar datos meteorológicos diarios
condiciones_diarias = cell(dias_simulacion, 1);
velocidades_diarias = zeros(dias_simulacion, 1);

for dia = 1:dias_simulacion
    estacion = asignar_estacion(dia);
    probabilidades = PROBABILIDADES_ESTACION(estacion);
    
    % Generar condición meteorológica según probabilidades
    rand_val = rand() * 100;
    if rand_val <= probabilidades(1)
        condicion = condiciones_nombres{1};
    elseif rand_val <= sum(probabilidades(1:2))
        condicion = condiciones_nombres{2};
    elseif rand_val <= sum(probabilidades(1:3))
        condicion = condiciones_nombres{3};
    else
        condicion = condiciones_nombres{4};
    end
    
    % Generar velocidad según la condición
    rango_velocidad = VELOCIDADES(condicion);
    velocidad = rango_velocidad(1) + rand() * (rango_velocidad(2) - rango_velocidad(1));
    
    condiciones_diarias{dia} = condicion;
    velocidades_diarias(dia) = velocidad;
end

%% Cálculo de velocidad promedio semanal
num_semanas = floor(dias_simulacion / 7);
semanas = 1:num_semanas;
velocidades_semanales = zeros(1, num_semanas);

for i = 1:num_semanas
    inicio_semana = (i-1) * 7 + 1;
    fin_semana = min(i * 7, dias_simulacion);
    velocidades_semanales(i) = mean(velocidades_diarias(inicio_semana:fin_semana));
end
media_anual = mean(velocidades_semanales);

%% Gráfica de velocidad promedio semanal
figure('Position', [100, 100, 1000, 600]);
hold on;

% Usar colores más visibles y contrastantes
color_principal = [0, 0.4470, 0.7410]; % Azul más visible
color_media = [0.8500, 0.3250, 0.0980]; % Naranja para contraste

plot(semanas, velocidades_semanales, 'o-', 'Color', color_principal, ...
    'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', color_principal);

% Línea de media anual con color contrastante
yline(media_anual, '--', 'Color', color_media, 'LineWidth', 3, ...
    'DisplayName', sprintf('Media anual (%.1f Mbps)', media_anual));

xlim([1 52]);

xlabel('Semana del a\~no', 'FontSize', 12);
ylabel('Velocidad promedio (Mbps)', 'FontSize', 12);
title('Velocidad promedio semanal en banda X durante un a\~no', 'FontSize', 14);
grid on; grid minor;

% Corregir la leyenda - agregar ambas series
legend({'Velocidad promedio semanal', sprintf('Media anual (%.1f Mbps)', media_anual)}, ...
    'Location', 'best', 'FontSize', 11);
hold off;

out_dir = '/Users/miguelrosa/Desktop/UPM/TFG/Preliminary DesignEO github/PreliminaryDesignEO/Latex_Code/7.Segmento_Tierra';
exportgraphics(gcf, fullfile(out_dir, 'velmediasemanal.png'), 'Resolution', 150);


%% Gráfica de velocidad diaria
figure('Position', [150, 150, 1000, 600]);
hold on;

% Definir colores más visibles y diferenciables
colores_visibles = [
    0, 0.4470, 0.7410;      % Azul para cond. favorables
    0.8500, 0.3250, 0.0980; % Naranja para lluvia moderada
    0.6350, 0.0780, 0.1840; % Rojo oscuro para lluvia intensa
    0.4660, 0.6740, 0.1880  % Verde para nubes densas
];

colores_map = containers.Map(condiciones_nombres, ...
    {colores_visibles(1,:), colores_visibles(2,:), ...
     colores_visibles(3,:), colores_visibles(4,:)});

% Graficar puntos por condición con mayor tamaño
for i = 1:length(condiciones_nombres)
    condicion = condiciones_nombres{i};
    indices = strcmp(condiciones_diarias, condicion);
    dias_condicion = find(indices);
    velocidades_condicion = velocidades_diarias(indices);
    
    scatter(dias_condicion, velocidades_condicion, 50, ...
        'MarkerFaceColor', colores_map(condicion), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8, ...
        'DisplayName', condicion);
end

% Línea de media anual en negro para mejor contraste
media_anual_diaria = mean(velocidades_diarias);
yline(media_anual_diaria, '--', 'Color', 'k', 'LineWidth', 3, ...
    'DisplayName', sprintf('Media anual (%.1f Mbps)', media_anual_diaria));
xlim([1 365]);

xlabel('D\''ia del a\~no', 'FontSize', 12);
ylabel('Velocidad de descarga (Mbps)', 'FontSize', 12);
title('Velocidad diaria en banda X durante un a\~no', 'FontSize', 14);
grid on; grid minor;
legend('Location', 'best', 'FontSize', 11);
hold off;

exportgraphics(gcf, fullfile(out_dir, 'veldiariabandax.png'), 'Resolution', 150);

%% Resumen estadístico
fprintf('\nResumen Estadístico:\n');
fprintf('%-20s %-10s %-10s %-12s\n', 'Condición', 'Media', 'Días', 'Porcentaje');
fprintf('%-20s %-10s %-10s %-12s\n', repmat('-', 1, 20), repmat('-', 1, 10), ...
        repmat('-', 1, 10), repmat('-', 1, 12));

for i = 1:length(condiciones_nombres)
    condicion = condiciones_nombres{i};
    indices = strcmp(condiciones_diarias, condicion);
    velocidades_condicion = velocidades_diarias(indices);
    
    media_condicion = mean(velocidades_condicion);
    num_dias = sum(indices);
    porcentaje = (num_dias / dias_simulacion) * 100;
    
    fprintf('%-20s %-10.1f %-10d %-12.1f\n', condicion, media_condicion, ...
            num_dias, porcentaje);
end

fprintf('\nTotal de días simulados: %d\n', dias_simulacion);
fprintf('Velocidad promedio general: %.1f Mbps\n', mean(velocidades_diarias));

