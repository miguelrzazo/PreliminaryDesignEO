% =========================================================================
% SCRIPT PRINCIPAL PARA LA SIMULACION DE POTENCIA DE UN SATELITE
% =========================================================================
clear; clc; close all;

% 1. Cargar la configuracion de la mision
params = configuracion();
fprintf('>> Configuracion cargada: superficie a simular = %.3f m^2, consumo = %.1f W\n', ...
    params.superficiePaneles, params.consumoPromedio);
fprintf('>> Tipo de Orbita: LTAN a las %d:00\n', params.LTAN);

% 2. Realizar calculos orbitales
orbita = calculos_orbitales(params);

% 3. Comprobacion de viabilidad energetica
eficiencia_eol = params.eficienciaPaneles * (1 - params.tasaDegradacionPaneles * params.anosSimulacion);
potencia_generada_media_eol = params.superficiePaneles * eficiencia_eol * params.irradiancia * (1 - orbita.fraccionEclipse);

fprintf('\n--- Comprobacion de Viabilidad Energetica ---\n');
if potencia_generada_media_eol < params.consumoPromedio
    superficie_minima_teorica = params.consumoPromedio / ((1 - orbita.fraccionEclipse) * eficiencia_eol * params.irradiancia);
    warning('Diseno energeticamente NO VIABLE con la superficie configurada.');
    fprintf('La generacion media (EOL) de %.2f W es INFERIOR al consumo de %.2f W.\n', ...
        potencia_generada_media_eol, params.consumoPromedio);
    fprintf('La simulacion mostrara un agotamiento de la bateria.\n');
    fprintf('La superficie minima teorica para este consumo es de %.3f m^2.\n', superficie_minima_teorica);
else
    fprintf('Diseno energeticamente viable: Generacion media (EOL) de %.2f W > Consumo de %.2f W.\n', ...
        potencia_generada_media_eol, params.consumoPromedio);
end
fprintf('-------------------------------------------\n\n');


% 4. Ejecutar la simulacion
disp('Iniciando simulacion de 8 anos...');
resultados = simulacion_energia(params, orbita);
disp('Simulacion completada.');

% 5. Analizar y mostrar los resultados de diseno
analisis_resultados(params, orbita, resultados);

% 6. Generar las visualizaciones graficas
visualizacion(params, orbita, resultados);



function params = configuracion()
    % Devuelve una estructura con todos los parametros de configuracion.
    
    % --- Parametros de la Mision y Orbita ---
    params.altitud = 760;            % km, altura orbital
    params.inclinacion = 98;         % grados, inclinacion (heliosincrona)
    params.anosSimulacion = 8;       % anos
    params.LTAN = 6;                 % Horas (6 o 18 para Dawn/Dusk, 12 o 0 para Noon/Midnight)

    % --- Parametros de los Paneles Solares ---
    params.superficiePaneles = 0.05; % m^2 
    params.eficienciaPaneles = 0.30; % 30%
    params.tasaDegradacionPaneles = 0.045; % 4.5% anual
    params.densidadPaneles = 15;     % kg/m^2

    % --- Parametros de Consumo y Bateria ---
    params.consumoPromedio = 20;      % W 
    params.capacidadBateria = 10;    % Wh
    params.profundidadDescargaMaxima = 0.9; 
    params.densidadEnergeticaBateria = 200; % Wh/kg
    params.eficienciaCarga = 0.9;
    params.eficienciaDescarga = 0.95;
    params.tasaDegradacionAnualBateria = 0.02; % 2% anual

    % --- Constantes Fisicas ---
    params.radioTierra = 6371;       % km
    params.mu = 398600;              % km^3/s^2
    params.irradiancia = 1366;       % W/m^2

    % --- Parametros de Simulacion ---
    params.puntosPorDia = 100; % Puntos de calculo por dia
end

function orbita = calculos_orbitales(params)
    % Calcula y devuelve una estructura con las propiedades orbitales,
    % considerando el LTAN para el calculo del eclipse.
    
    radioOrbita = params.altitud + params.radioTierra;
    
    orbita.periodoOrbital = 2 * pi * sqrt(radioOrbita^3 / params.mu); % segundos
    orbita.periodoDias = orbita.periodoOrbital / 86400;
    orbita.orbitasPorDia = 24 * 3600 / orbita.periodoOrbital;
    
    % --- Calculo del tiempo en eclipse basado en LTAN ---
    if params.LTAN == 6 || params.LTAN == 18
        % Orbita Dawn-Dusk: el plano orbital es perpendicular al vector Sol-Tierra.
        % Idealmente, no hay eclipse.
        orbita.fraccionEclipse = 0;
    elseif params.LTAN == 12 || params.LTAN == 0
        % Orbita Noon-Midnight: el plano orbital contiene el vector Sol-Tierra.
        % Se experimenta el maximo eclipse posible.
        beta_angle = asin(params.radioTierra / radioOrbita);
        orbita.fraccionEclipse = (2 * beta_angle) / (2 * pi);
    else
        % Para otros LTAN, se necesitaria un modelo mas complejo (usando el angulo beta).
        % Por ahora, se asume el peor caso (Noon-Midnight) para LTANs no especificados.
        warning('LTAN no estandar. Asumiendo peor caso de eclipse (Noon-Midnight).');
        beta_angle = asin(params.radioTierra / radioOrbita);
        orbita.fraccionEclipse = (2 * beta_angle) / (2 * pi);
    end
    
    orbita.tiempoEclipsePorOrbita = orbita.periodoOrbital * orbita.fraccionEclipse; % segundos
end

function resultados = simulacion_energia(params, orbita)
    % Realiza la simulación energética de forma vectorizada.

    % --- Preparación del Vector de Tiempo ---
    diasSimulacionTotal = params.anosSimulacion * 365;
    puntosTotales = diasSimulacionTotal * params.puntosPorDia;
    tiempoDias = linspace(0, diasSimulacionTotal, puntosTotales);
    deltaT_horas = (tiempoDias(2) - tiempoDias(1)) * 24;

    % --- Vectorización de Cálculos Orbitales y de Potencia ---
    % Determinar estado de eclipse para todos los puntos
    fraccionOrbita = mod(tiempoDias / orbita.periodoDias, 1);
    anguloOrbita_rad = fraccionOrbita * 2 * pi;
    estadoEnEclipse = abs(sin(anguloOrbita_rad)) < orbita.fraccionEclipse / 2; % Condición simplificada y robusta

    % Calcular degradación y variación de irradiancia
    degradacionPaneles = (1 - params.tasaDegradacionPaneles * tiempoDias / 365);
    eficienciaActual = params.eficienciaPaneles * degradacionPaneles;
    diaDelAno = mod(tiempoDias, 365);
    variacionIrradiancia = 1 + 0.034 * cos(2*pi*(diaDelAno-4)/365);

    % Calcular potencia generada (0 en eclipse)
    potenciaGenerada = params.superficiePaneles * eficienciaActual .* params.irradiancia .* variacionIrradiancia;
    potenciaGenerada(estadoEnEclipse) = 0;

    % --- Bucle Secuencial para el Estado de la Batería ---
    energiaBateria = zeros(1, puntosTotales);
    energiaBateria(1) = params.capacidadBateria * 0.8; % Carga inicial
    degradacionBateria = (1 - params.tasaDegradacionAnualBateria * tiempoDias / 365);
    capacidadBateriaActual = params.capacidadBateria * degradacionBateria;

    for i = 2:puntosTotales
        energiaConsumida = params.consumoPromedio * deltaT_horas;
        if ~estadoEnEclipse(i) % En sol
            energiaNeta = (potenciaGenerada(i) * deltaT_horas) - energiaConsumida;
            energiaBateria(i) = energiaBateria(i-1) + energiaNeta * params.eficienciaCarga;
        else % En eclipse
            energiaBateria(i) = energiaBateria(i-1) - energiaConsumida / params.eficienciaDescarga;
        end
        % Aplicar límites de la batería (0 y capacidad máxima actual)
        energiaBateria(i) = max(0, min(energiaBateria(i), capacidadBateriaActual(i)));
    end

    % --- Guardar resultados en una estructura ---
    resultados.tiempoDias = tiempoDias;
    resultados.potenciaGenerada = potenciaGenerada;
    resultados.energiaBateria = energiaBateria;
    resultados.estadoEnEclipse = estadoEnEclipse;
    resultados.capacidadBateriaActual = capacidadBateriaActual;
end
function analisis_resultados(params, orbita, resultados)
    % Calcula e imprime un resumen, usando (1 - orbita.fraccionEclipse).
    
    fprintf('\n========== ANALISIS DE RESULTADOS ==========\n');
    
    % --- Resultados de la SIMULACION ---
    fprintf('\n1. Resultados de la Simulacion con Superficie Configurada (%.3f m^2):\n', params.superficiePaneles);
    potencia_ini = params.superficiePaneles * params.eficienciaPaneles * params.irradiancia;
    eficiencia_fin_sim = params.eficienciaPaneles * (1 - params.tasaDegradacionPaneles * params.anosSimulacion);
    potencia_fin_sim = params.superficiePaneles * eficiencia_fin_sim * params.irradiancia;
    energiaMinimaObservada = min(resultados.energiaBateria);
    descargaMaximaObservada = (params.capacidadBateria - energiaMinimaObservada) / params.capacidadBateria;

    fprintf('- Potencia maxima generada (BOL): %.2f W\n', potencia_ini);
    fprintf('- Potencia maxima generada (EOL): %.2f W\n', potencia_fin_sim);
    fprintf('- Descarga maxima de bateria observada: %.1f %% (Energia minima: %.2f Wh)\n', descargaMaximaObservada * 100, energiaMinimaObservada);

    % --- Calculos para el DIMENSIONAMIENTO TEORICO ---
    fprintf('\n2. Calculos de Dimensionamiento Teorico:\n');
    
    % ** LA CORRECCION ESTA AQUI **
    % Se usa (1 - orbita.fraccionEclipse) en lugar del inexistente 'fraccionSol'
    fraccionSol = 1 - orbita.fraccionEclipse;
    
    eficiencia_fin_teo = params.eficienciaPaneles * (1 - params.tasaDegradacionPaneles * params.anosSimulacion);
    superficieMinimaRequerida = params.consumoPromedio / (fraccionSol * eficiencia_fin_teo * params.irradiancia);
    pesoPanelesRequerido = superficieMinimaRequerida * params.densidadPaneles;
    
    fprintf('- Superficie MINIMA de paneles requerida (EOL): %.3f m^2\n', superficieMinimaRequerida);
    fprintf('- Peso estimado para esta superficie minima (%.1f kg/m^2): %.2f kg\n', params.densidadPaneles, pesoPanelesRequerido);

    energia_en_eclipse = (params.consumoPromedio * orbita.tiempoEclipsePorOrbita / 3600);
    if energia_en_eclipse > 0
        capacidadRequeridaBateria = energia_en_eclipse / (params.profundidadDescargaMaxima * params.eficienciaDescarga);
        pesoBateriaRequerido = capacidadRequeridaBateria / params.densidadEnergeticaBateria;
        fprintf('- Capacidad MINIMA de bateria requerida (DoD %.0f%%): %.2f Wh\n', params.profundidadDescargaMaxima*100, capacidadRequeridaBateria);
        fprintf('- Peso estimado para esta capacidad minima (%.0f Wh/kg): %.2f kg\n', params.densidadEnergeticaBateria, pesoBateriaRequerido);
    else
        fprintf('- Capacidad MINIMA de bateria requerida: 0 Wh (No hay eclipse)\n');
    end
    fprintf('===========================================\n');
end


function visualizacion(params, orbita, resultados)
    % Genera todas las graficas de resultados con colores suaves y formato LaTeX.

    %% --- Definicion de Colores Suaves ---
    colores.potenciaInst = [0.6 0.75 1.0];  % Azul claro
    colores.potenciaAvg = [0.1 0.3 0.7];     % Azul oscuro
    colores.energiaInst = [0.6 1.0 0.6];    % Verde claro
    colores.energiaAvg = [0.1 0.5 0.1];     % Verde oscuro
    colores.limiteDoD = [1.0 0.5 0.2];      % Naranja
    
    tiempo = resultados.tiempoDias;
    
    % --- Grafica 1: Potencia Generada ---
    figure('Name', 'Rendimiento de Paneles Solares');
    plot(tiempo, resultados.potenciaGenerada, 'LineWidth', 0.5, 'Color', colores.potenciaInst);
    hold on;
    window_size = 30 * params.puntosPorDia; % Promedio de 30 dias
    potenciaPromedioMovil = movmean(resultados.potenciaGenerada, window_size);
    plot(tiempo, potenciaPromedioMovil, 'LineWidth', 2, 'Color', colores.potenciaAvg);
    
    title('Evolucion de la Potencia Generada', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Tiempo (dias)', 'Interpreter', 'latex');
    ylabel('Potencia (W)', 'Interpreter', 'latex');
    legend({'Potencia Instantanea', 'Promedio movil (30 dias)'}, 'Interpreter', 'latex');
    grid on; hold off;

    % --- Grafica 2: Energia de la Bateria ---
    figure('Name', 'Estado de la Bateria');
    plot(tiempo, resultados.energiaBateria, 'LineWidth', 0.5, 'Color', colores.energiaInst);
    hold on;
    energiaPromedioMovil = movmean(resultados.energiaBateria, window_size);
    plot(tiempo, energiaPromedioMovil, 'LineWidth', 2, 'Color', colores.energiaAvg);
    
    limiteDescarga = params.capacidadBateria * (1 - params.profundidadDescargaMaxima);
    line([0 tiempo(end)], [limiteDescarga limiteDescarga], 'Color', colores.limiteDoD, 'LineStyle', '--', 'LineWidth', 1.5);
    
    title('Estado de Carga de la Bateria', 'Interpreter', 'latex', 'FontSize', 14);
    xlabel('Tiempo (dias)', 'Interpreter', 'latex');
    ylabel('Energia (Wh)', 'Interpreter', 'latex');
    legend({'Energia Instantanea', 'Promedio movil (30 dias)', 'Limite de Descarga (DoD$_{max}$)'}, 'Interpreter', 'latex');
    ylim([0, params.capacidadBateria * 1.1]);
    grid on; hold off;
end

