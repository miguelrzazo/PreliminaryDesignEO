% Script para simular la generación de potencia de paneles solares en un satélite
% durante 8 años, considerando eclipses y dimensionamiento de batería
% Válido para órbitas con inclinaciones entre 50 y 140 grados

%% Parámetros configurables
% Órbita
altitud = 700;            % km, altura orbital
inclinacion = 98;         % grados, inclinación de la órbita (entre 50-140)

% Paneles solares
superficiePaneles = 3.06;   % m^2, superficie total de los paneles solares
eficienciaPaneles = 0.30; % eficiencia de conversión (30%)
tasaDegradacionPaneles = 0.0045; % 0.45% de degradación anual

% Consumo de energía
consumoPromedio = 500;    % W, consumo promedio del satélite

% Batería
capacidadBateria = 592.1;  % Wh, capacidad a evaluar
profundidadDescargaMaxima = 0.2; % 20% de capacidad como mínimo permitido
densidadEnergeticaBateria = 200; % Wh/kg, densidad energética típica de baterías espaciales

% Peso de los paneles solares
densidadPaneles = 2;     % kg/m^2, densidad típica de paneles solares espaciales

%% Constantes físicas
radioTierra = 6371;       % km, radio de la Tierra
mu = 398600;              % km^3/s^2, constante gravitacional de la Tierra
irradiancia = 1366;       % W/m^2, constante solar

%% Cálculos orbitales
radioOrbita = altitud + radioTierra;
periodoOrbital = 2*pi*sqrt(radioOrbita^3/mu); % segundos
periodoDias = periodoOrbital / 86400;  % días
orbitasPorDia = 24*3600 / periodoOrbital;

% Cálculo del tiempo en eclipse
anguloEclipse = 2 * asin(radioTierra / radioOrbita) * 180/pi; % grados
fraccionEclipse = anguloEclipse / 360; % fracción de la órbita en eclipse
tiempoEclipsePorOrbita = periodoOrbital * fraccionEclipse; % segundos

%% Parámetros de simulación
anosSimulacion = 8;       % años a simular
diasSimulacion = anosSimulacion * 365; % días a simular
puntosOrbitasPorDia = 100; % puntos de cálculo por día
puntosTotales = diasSimulacion * puntosOrbitasPorDia;
tiempoSimulacion = linspace(0, diasSimulacion, puntosTotales); % días

%% Potencia máxima teórica
potenciaMaxima = superficiePaneles * eficienciaPaneles * irradiancia; % W

%% Inicialización de vectores
potenciaGenerada = zeros(puntosTotales, 1);
energiaBateria = zeros(puntosTotales, 1);
estadoEnEclipse = zeros(puntosTotales, 1);

% Propiedades de la batería
eficienciaCarga = 0.9;     % eficiencia de carga
eficienciaDescarga = 0.95; % eficiencia de descarga
energiaBateria(1) = capacidadBateria * 0.8; % 80% de carga inicial

% Degradación de la batería
tasaDegradacionAnual = 0.02; % 2% de pérdida de capacidad por año
capacidadBateriaTiempo = @(dias) capacidadBateria * (1 - tasaDegradacionAnual * dias / 365);

% Degradación de los paneles solares
eficienciaPanelesTiempo = @(dias) eficienciaPaneles * (1 - tasaDegradacionPaneles * dias / 365);

%% Simulación de generación y consumo
for i = 2:puntosTotales
    % Posición orbital
    tiempoDias = tiempoSimulacion(i);
    orbitasCompletadas = tiempoDias / periodoDias;
    fraccionOrbita = orbitasCompletadas - floor(orbitasCompletadas);
    
    % Eficiencia actual de los paneles (considerando degradación)
    eficienciaActual = eficienciaPanelesTiempo(tiempoDias);
    
    % Potencia máxima actual (considerando degradación)
    potenciaMaximaActual = superficiePaneles * eficienciaActual * irradiancia;

    % Variación anual de la irradiancia solar
    diaDelAno = mod(tiempoDias, 365);
    variacionAnual = 1 + 0.034 * cos(2*pi*(diaDelAno-4)/365);

    % Determinar si está en eclipse
    anguloOrbita = fraccionOrbita * 360;
    limiteInferiorEclipse = 180 - anguloEclipse/2;
    limiteSuperiorEclipse = 180 + anguloEclipse/2;

    if anguloOrbita >= limiteInferiorEclipse && anguloOrbita <= limiteSuperiorEclipse
        estadoEnEclipse(i) = 1; % En eclipse
        potenciaGenerada(i) = 0; % No hay generación en eclipse
    else
        estadoEnEclipse(i) = 0; % No en eclipse
        
        % Ángulo de incidencia simplificado (similar al original)
        anguloIncidencia = 0; % Asumimos paneles siempre apuntando al Sol
        
        potenciaGenerada(i) = potenciaMaximaActual * variacionAnual * cos(anguloIncidencia * pi/180);
    end

    % Actualizar el estado de la batería
    deltaT = (tiempoSimulacion(i) - tiempoSimulacion(i-1)) * 24; % horas
    capacidadActual = capacidadBateriaTiempo(tiempoDias);

    if estadoEnEclipse(i) == 0
        energiaGenerada = potenciaGenerada(i) * deltaT;
        energiaConsumida = consumoPromedio * deltaT;
        energiaNeta = energiaGenerada * eficienciaCarga - energiaConsumida;
        energiaBateria(i) = min(capacidadActual, energiaBateria(i-1) + energiaNeta);
    else
        energiaConsumida = consumoPromedio * deltaT;
        energiaBateria(i) = max(capacidadActual * profundidadDescargaMaxima, energiaBateria(i-1) - energiaConsumida / eficienciaDescarga);
        energiaBateria(i) = min(capacidadActual, energiaBateria(i));
        % En el bucle de simulación, después de actualizar energiaBateria[i]
if estadoEnEclipse(i) == 0
    energiaGenerada = potenciaGenerada(i) * deltaT;
    energiaConsumida = consumoPromedio * deltaT;
    energiaNeta = energiaGenerada * eficienciaCarga - energiaConsumida;
    energiaBateria(i) = energiaBateria(i-1) + energiaNeta;
    % Aplicar límites min/max
    energiaBateria(i) = min(capacidadActual, energiaBateria(i));
    energiaBateria(i) = max(capacidadActual * profundidadDescargaMaxima, energiaBateria(i));
else
    energiaConsumida = consumoPromedio * deltaT;
    energiaBateria(i) = energiaBateria(i-1) - energiaConsumida / eficienciaDescarga;
    % Aplicar límites min/max
    energiaBateria(i) = min(capacidadActual, energiaBateria(i));
    energiaBateria(i) = max(capacidadActual * profundidadDescargaMaxima, energiaBateria(i));
end

    end
end

%% Cálculo de promedios móviles de 30 días
window_size = 30 * puntosOrbitasPorDia; % 30 días
potenciaPromedioMovil = movmean(potenciaGenerada, window_size);
energiaPromedioMovil = movmean(energiaBateria, window_size);

%% Análisis de capacidad de batería necesaria
energiaMinima = min(energiaBateria);
descargaMaximaObservada = (capacidadBateria - energiaMinima) / capacidadBateria;

% Cálculo de capacidad mínima requerida
if energiaMinima < capacidadBateria * profundidadDescargaMaxima
    capacidadRequerida = capacidadBateria * descargaMaximaObservada / (1 - profundidadDescargaMaxima);
else
    capacidadRequerida = capacidadBateria;
end

% Cálculo de superficie mínima de paneles considerando degradación al final de la vida útil
eficienciaFinalVidaUtil = eficienciaPanelesTiempo(diasSimulacion);
superficieMinimaPaneles = (consumoPromedio * (1 + fraccionEclipse/(1-fraccionEclipse))) / (irradiancia * eficienciaFinalVidaUtil * (1 - fraccionEclipse));

%% Visualización de resultados
figure(1);
plot(tiempoSimulacion, potenciaGenerada, 'LineWidth', 0.5, 'Color', [0.7 0.7 1]);
hold on;
plot(tiempoSimulacion, potenciaPromedioMovil, 'LineWidth', 2, 'Color', [0 0 0.8]);
title('Potencia generada por los paneles solares durante 8 años');
xlabel('Tiempo (días)');
ylabel('Potencia (W)');
legend('Potencia instantánea', 'Promedio móvil 30 días');
grid on;

figure(2);
plot(tiempoSimulacion, energiaBateria, 'LineWidth', 0.5, 'Color', [0.7 1 0.7]);
hold on;
plot(tiempoSimulacion, energiaPromedioMovil, 'LineWidth', 2, 'Color', [0 0.6 0]);
plot(tiempoSimulacion, ones(size(tiempoSimulacion)) * capacidadBateria * profundidadDescargaMaxima, 'r--', 'LineWidth', 1.5);
title('Estado de carga de la batería durante 8 años');
xlabel('Tiempo (días)');
ylabel('Energía (Wh)');
legend('Energía instantánea', 'Promedio móvil 30 días', 'Límite de descarga mínima permitida');
ylim([0, capacidadBateria * 1.1]);
grid on;

figure(3);
indiceOrbita = round(puntosOrbitasPorDia * periodoDias);
subplot(3,1,1);
yyaxis left;
plot(tiempoSimulacion(1:indiceOrbita), potenciaGenerada(1:indiceOrbita), 'b-', 'LineWidth', 1.5);
ylabel('Potencia (W)');
yyaxis right;
plot(tiempoSimulacion(1:indiceOrbita), energiaBateria(1:indiceOrbita), 'g-', 'LineWidth', 1.5);
hold on;
plot(tiempoSimulacion(1:indiceOrbita), estadoEnEclipse(1:indiceOrbita) * capacidadBateria * 0.5, 'r-', 'LineWidth', 1.5);
title('Ciclos orbitales durante una órbita');
xlabel('Tiempo (días)');
ylabel('Energía (Wh) / Estado');
legend('Potencia generada', 'Energía batería', 'Estado de eclipse');
grid on;

indiceDia = puntosOrbitasPorDia;
subplot(3,1,2);
yyaxis left;
plot(tiempoSimulacion(1:indiceDia), potenciaGenerada(1:indiceDia), 'b-', 'LineWidth', 1.5);
ylabel('Potencia (W)');
yyaxis right;
plot(tiempoSimulacion(1:indiceDia), energiaBateria(1:indiceDia), 'g-', 'LineWidth', 1.5);
hold on;
plot(tiempoSimulacion(1:indiceDia), estadoEnEclipse(1:indiceDia) * capacidadBateria * 0.5, 'r-', 'LineWidth', 1.5);
title('Ciclos orbitales durante un día');
xlabel('Tiempo (días)');
ylabel('Energía (Wh) / Estado');
legend('Potencia generada', 'Energía batería', 'Estado de eclipse');
grid on;

indiceSemana = puntosOrbitasPorDia * 7;
subplot(3,1,3);
yyaxis left;
plot(tiempoSimulacion(1:indiceSemana), potenciaGenerada(1:indiceSemana), 'b-', 'LineWidth', 1.5);
ylabel('Potencia (W)');
yyaxis right;
plot(tiempoSimulacion(1:indiceSemana), energiaBateria(1:indiceSemana), 'g-', 'LineWidth', 1.5);
hold on;
plot(tiempoSimulacion(1:indiceSemana), estadoEnEclipse(1:indiceSemana) * capacidadBateria * 0.5, 'r-', 'LineWidth', 1.5);
title('Ciclos orbitales durante una semana');
xlabel('Tiempo (días)');
ylabel('Energía (Wh) / Estado');
legend('Potencia generada', 'Energía batería', 'Estado de eclipse');
grid on;

%% Exportación de resultados a archivos TXT y CSV

% Crear tabla con resultados principales
T = table(tiempoSimulacion', potenciaGenerada, energiaBateria, estadoEnEclipse, ...
          'VariableNames', {'Tiempo_dias', 'Potencia_W', 'Energia_Bateria_Wh', 'Estado_Eclipse'});

% Exportar a CSV
writetable(T, 'resultados_simulacion.csv');

% Exportar a TXT (formato más legible)
fileID = fopen('resultados_simulacion.txt', 'w');
fprintf(fileID, 'RESULTADOS DE SIMULACIÓN DE POTENCIA PARA SATÉLITE\n');
fprintf(fileID, '=================================================\n\n');
fprintf(fileID, 'Parámetros orbitales:\n');
fprintf(fileID, '- Altitud: %.2f km\n', altitud);
fprintf(fileID, '- Inclinación: %.2f grados\n', inclinacion);
fprintf(fileID, '- Período orbital: %.2f minutos\n', periodoOrbital/60);
fprintf(fileID, '- Tiempo en eclipse por órbita: %.2f minutos (%.1f%%)\n', tiempoEclipsePorOrbita/60, fraccionEclipse*100);
fprintf(fileID, '- Órbitas por día: %.2f\n\n', orbitasPorDia);

fprintf(fileID, 'Rendimiento energético:\n');
fprintf(fileID, '- Potencia máxima generada (inicio de vida): %.2f W\n', potenciaMaxima);
fprintf(fileID, '- Potencia máxima generada (fin de vida): %.2f W\n', potenciaMaxima * (1 - anosSimulacion * tasaDegradacionPaneles));
fprintf(fileID, '- Degradación anual estimada de los paneles solares: %.2f%%\n', tasaDegradacionPaneles*100);
fprintf(fileID, '- Descarga máxima de batería observada: %.1f%%\n\n', descargaMaximaObservada*100);

fprintf(fileID, 'Dimensionamiento de la batería:\n');
fprintf(fileID, '- Capacidad mínima de batería requerida: %.2f Wh\n', capacidadRequerida);
fprintf(fileID, '- Peso estimado de la batería: %.2f kg\n\n', capacidadRequerida/densidadEnergeticaBateria);

fprintf(fileID, 'Dimensionamiento de los paneles solares:\n');
fprintf(fileID, '- Superficie mínima de paneles solares requerida: %.2f m^2\n', superficieMinimaPaneles);
fprintf(fileID, '- Peso estimado de los paneles solares: %.2f kg\n', superficieMinimaPaneles * densidadPaneles);
fclose(fileID);

% Exportar datos de resumen a CSV
resumen = table(...
    {'Altitud (km)'; 'Inclinación (grados)'; 'Período orbital (min)'; 'Tiempo en eclipse (min)'; ...
     'Potencia máxima inicio (W)'; 'Potencia máxima fin (W)'; 'Descarga máxima batería (%)'; ...
     'Capacidad batería requerida (Wh)'; 'Peso batería (kg)'; 'Superficie paneles (m^2)'; 'Peso paneles (kg)'}, ...
    [altitud; inclinacion; periodoOrbital/60; tiempoEclipsePorOrbita/60; ...
     potenciaMaxima; potenciaMaxima*(1-anosSimulacion*tasaDegradacionPaneles); descargaMaximaObservada*100; ...
     capacidadRequerida; capacidadRequerida/densidadEnergeticaBateria; superficieMinimaPaneles; superficieMinimaPaneles*densidadPaneles], ...
    'VariableNames', {'Parámetro', 'Valor'});
writetable(resumen, 'resumen_simulacion.csv');

% Mostrar mensaje de confirmación
fprintf('\nResultados del dimensionamiento para 8 años de misión:\n');
fprintf('----------------------------------------------------------\n');
fprintf('Parámetros orbitales:\n');
fprintf('- Período orbital: %.2f minutos\n', periodoOrbital/60);
fprintf('- Tiempo en eclipse por órbita: %.2f minutos (%.1f%%)\n', tiempoEclipsePorOrbita/60, fraccionEclipse*100);
fprintf('- Órbitas por día: %.2f\n\n', orbitasPorDia);

fprintf('Rendimiento energético:\n');
fprintf('- Potencia máxima generada (inicio de vida): %.2f W\n', potenciaMaxima);
fprintf('- Potencia máxima generada (fin de vida): %.2f W\n', potenciaMaxima * (1 - anosSimulacion * tasaDegradacionPaneles));
fprintf('- Degradación anual estimada de los paneles solares: %.2f%%\n', tasaDegradacionPaneles*100);
fprintf('- Descarga máxima de batería observada: %.1f%%\n\n', descargaMaximaObservada*100);

fprintf('Dimensionamiento de la batería:\n');
fprintf('- Capacidad mínima de batería requerida: %.2f Wh\n', capacidadRequerida);
fprintf('- Peso estimado de la batería: %.2f kg\n\n', capacidadRequerida/densidadEnergeticaBateria);

fprintf('Dimensionamiento de los paneles solares:\n');
fprintf('- Superficie mínima de paneles solares requerida: %.2f m^2\n', superficieMinimaPaneles);
fprintf('- Peso estimado de los paneles solares: %.2f kg\n', superficieMinimaPaneles * densidadPaneles);
fprintf('\nLos resultados han sido exportados a archivos CSV y TXT.\n');
