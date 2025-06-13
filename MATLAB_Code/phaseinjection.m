% --- 1. CONSTANTES Y PARÁMETROS ORBITALES ---
mu = 3.986004418e14; % Parámetro gravitacional estándar de la Tierra (m^3/s^2)
R_tierra = 6371e3;   % Radio de la Tierra (m)

% Órbita final (circular)
altitud_final = 520e3;
r_final = R_tierra + altitud_final;
T_final = 2 * pi * sqrt(r_final^3 / mu); % Periodo para el print
n_final = 2 * pi / T_final;             % Movimiento medio

% Órbita de fase (elíptica)
altitud_apogeo_fase = 720e3;
r_apogeo_fase = R_tierra + altitud_apogeo_fase;
r_perigeo_fase = r_final;
a_fase = (r_perigeo_fase + r_apogeo_fase) / 2;
e_fase = (r_apogeo_fase - r_perigeo_fase) / (r_apogeo_fase + r_perigeo_fase);
T_fase = 2 * pi * sqrt(a_fase^3 / mu); % Periodo para el print
n_fase = 2 * pi / T_fase;             % Movimiento medio

% --- 2. ANÁLISIS PREVIO Y RESULTADOS EN CONSOLA ---
delta_n = n_final - n_fase; % Diferencia de movimiento medio
t_desfase_seg = pi / abs(delta_n); % Tiempo teórico para desfase de 180°
num_orbitas_fase = t_desfase_seg / T_fase; % Número de órbitas en faseo


fprintf('--- Análisis de la Maniobra de Faseo ---\n');
fprintf('Periodo de la órbita final: %.2f minutos\n', T_final / 60);
fprintf('Periodo de la órbita de fase: %.2f minutos\n', T_fase / 60);
fprintf('Tiempo requerido para un desfase de 180°: %.2f horas\n', t_desfase_seg / 3600);
fprintf('Número de órbitas en la trayectoria de fase: %.2f\n', num_orbitas_fase);
fprintf('-------------------------------------------\n\n');

% --- 3. CONFIGURACIÓN DE LA SIMULACIÓN ---
figure('Name', 'Simulación de Desfase Orbital - Versión Final', 'NumberTitle', 'off');
hold on;
axis equal;
grid on;

% Dibujar Tierra y Órbitas
[x, y, z] = sphere;
surf(x*R_tierra/1e3, y*R_tierra/1e3, z*R_tierra/1e3, 'FaceColor', '#0077be', 'EdgeColor', 'none', 'DisplayName', 'Tierra');
ang = linspace(0, 2*pi, 300);
plot(r_final/1e3 * cos(ang), r_final/1e3 * sin(ang), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Órbita Final (520 km)');
r_fase_plot = (a_fase*(1-e_fase^2))./(1+e_fase*cos(ang));
plot(r_fase_plot/1e3 .* cos(ang), r_fase_plot/1e3 .* sin(ang), 'b', 'LineWidth', 1.5, 'DisplayName', 'Órbita de Fase');

% Inicializar objetos gráficos
h_sat1 = plot(NaN, NaN, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Satélite 1');
h_etapa = plot(NaN, NaN, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Etapa/Satélite 2');
h_sat2 = plot(NaN, NaN, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'DisplayName', 'Satélite 2 (Desplegado)', 'Visible', 'off');
title_handle = title('Inicio de Misión');
legend('Location', 'northeastoutside', 'AutoUpdate', 'off');
xlim([-1.2*r_apogeo_fase/1e3, 1.2*r_apogeo_fase/1e3]);
ylim([-1.2*r_apogeo_fase/1e3, 1.2*r_apogeo_fase/1e3]);
xlabel('X (km)');
ylabel('Y (km)');

% --- 4. BUCLE DE ANIMACIÓN HASTA 180 GRADOS ---
t = 0;
dt = 120;
desfase_acumulado_rad = 0;

while desfase_acumulado_rad < pi
    % Posición Satélite 1
    theta1 = n_final * t;
    x1 = r_final * cos(theta1);
    y1 = r_final * sin(theta1);
    
    % Posición Satélite 2
    theta_etapa = n_fase * t;
    r_etapa_actual = (a_fase * (1 - e_fase^2)) / (1 + e_fase * cos(theta_etapa));
    x_etapa = r_etapa_actual * cos(theta_etapa);
    y_etapa = r_etapa_actual * sin(theta_etapa);
    
    % Calcular desfase acumulado
    desfase_acumulado_rad = delta_n * t;
    
    % Actualizar título
    title_handle.String = sprintf('Fase de Deriva | Tiempo: %.1f h | Desfase: %.1f°', t / 3600, rad2deg(desfase_acumulado_rad));
    
    % Actualizar posiciones
    set(h_sat1, 'XData', x1/1e3, 'YData', y1/1e3);
    set(h_etapa, 'XData', x_etapa/1e3, 'YData', y_etapa/1e3);
    
    drawnow;
    pause(0.01);
    t = t + dt;
end

% --- 5. FASE FINAL: DESPLIEGUE ---
title_handle.String = sprintf('Desfase de 180° alcanzado en %.1f horas. Desplegando Satélite 2.', t / 3600);
set(h_sat2, 'XData', x_etapa/1e3, 'YData', y_etapa/1e3, 'Visible', 'on');
set(h_etapa, 'Visible', 'off');

hold off;
