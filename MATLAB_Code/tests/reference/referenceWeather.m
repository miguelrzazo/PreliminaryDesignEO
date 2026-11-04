function results = referenceWeather()
oldRng=rng; cleanup=onCleanup(@() rng(oldRng)); rng(42);
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


results.dailyRates_Mbps=velocidades_diarias;
results.weeklyRates_Mbps=velocidades_semanales;
results.dailyConditions=condiciones_diarias;
end
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
