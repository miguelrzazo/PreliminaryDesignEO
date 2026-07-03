%% SIMULACION DE CICLO DE VIDA ORBITAL
%
% Descripcion:
% Este script simula el perfil completo de una mision satelital, incluyendo:
% 1. Decaimiento orbital por arrastre atmosferico.
% 2. Impulsos de mantenimiento para restaurar la altitud.
% 3. Impulso final de deorbitacion.
% 4. Curva de reentrada atmosferica hasta la desintegracion.
%
% El script genera dos visualizaciones:
%   - Figura 1: El ciclo de vida completo de la mision.
%   - Figura 2: Una comparacion de la reentrada con y sin impulso,
%               

clear; clc; close all;

%% ========================================================================
%  PARAMETROS INICIALES
%  ========================================================================

% --- Parametros de la Mision y Orbita ---
h0_km = 520;            % Altitud orbital inicial y objetivo (km)
mission_duration_years = 8; % Duracion de la mision (años)
decay_percentage_drop = 2.0;  % Porcentaje de decaimiento para activar impulso (%)
h_perigeo_deorbit_km = 120; % Altitud del perigeo para la maniobra de deorbitacion (km)
h_reentry_km = 100;      % Altitud a la que se considera reentrada/desintegracion (km)

% --- Parametros del Satelite ---
masa_seca_kg = 1.34;    % Masa seca del satelite (kg)
Area_m2 = 0.07;           % Area transversal efectiva para resistencia aerodinamica (m^2)
Cd = 2.5;               % Coeficiente de resistencia aerodinamica (adimensional)
Isp_s = 209;            % Impulso especifico ECAPS HPGP 100 mN (s)


%% ========================================================================
%  CONSTANTES Y CONFIGURACION
%  ========================================================================
mu = 3.986004418e14; R_earth_m = 6378e3; g0 = 9.80665;
h0_m = h0_km * 1e3;
h_umbral_m = h0_km * (1 - decay_percentage_drop / 100) * 1e3;
h_perigeo_deorbit_m = h_perigeo_deorbit_km * 1e3;
mission_duration_s = mission_duration_years * 365.25 * 24 * 3600;
hist_tiempo_s = [0]; hist_altitud_m = [h0_m];
t_actual_s = 0; h_actual_m = h0_m; masa_actual_kg = masa_seca_kg;
total_delta_v_mps = 0; num_impulsos = 0;

%% ========================================================================
%  1. PRE-CALCULO DEL COMBUSTIBLE PARA MANTENIMIENTO
%  ========================================================================
t_sim = 0; h_sim = h0_m; dv_sim = 0;
while t_sim < mission_duration_s
    if h_sim <= h_umbral_m
        r1 = R_earth_m + h_sim; r2 = R_earth_m + h0_m;
        dv_sim = dv_sim + hohmannDeltaV(r1, r2, mu); h_sim = h0_m;
    end
    a_sim = R_earth_m + h_sim;
    h_sim = h_sim + decayODE(0, a_sim, R_earth_m, mu, Cd, Area_m2, masa_seca_kg) * 3600;
    t_sim = t_sim + 3600;
end
fuel_mass_maintenance_kg = masa_seca_kg * (exp(dv_sim / (Isp_s * g0)) - 1) * 1.1;
masa_total_inicial_kg = masa_seca_kg + fuel_mass_maintenance_kg;
masa_actual_kg = masa_total_inicial_kg;
fprintf('Calculo preliminar:\n  - Masa total inicial: %.2f kg\n\n', masa_total_inicial_kg);

%% ========================================================================
%  2. SIMULACION DETALLADA DE LA MISION
%  ========================================================================
fprintf('Iniciando simulacion de la mision de %d a~nos...\n', mission_duration_years);
a_target_decay = R_earth_m + h_umbral_m;
while t_actual_s < mission_duration_s
    a_initial = R_earth_m + h_actual_m;
    options = odeset('Events', @(t,a) decayEvent(t, a, a_target_decay), 'RelTol', 1e-7);
    [t_ode, a_ode] = ode45(@(t,a) decayODE(t, a, R_earth_m, mu, Cd, Area_m2, masa_actual_kg), ...
                           [0, mission_duration_s - t_actual_s], a_initial, options);
    hist_tiempo_s = [hist_tiempo_s; t_actual_s + t_ode(2:end)];
    hist_altitud_m = [hist_altitud_m; a_ode(2:end) - R_earth_m];
    t_actual_s = t_actual_s + t_ode(end); h_actual_m = a_ode(end) - R_earth_m;
    if t_actual_s >= mission_duration_s; break; end
    num_impulsos = num_impulsos + 1;
    fprintf('Impulso #%d en t = %.2f a~nos. Altitud: %.1f km\n', num_impulsos, t_actual_s / (365.25*24*3600), h_actual_m/1000);
    r1 = R_earth_m + h_actual_m; r2 = R_earth_m + h0_m;
    dv_impulso = hohmannDeltaV(r1, r2, mu);
    masa_actual_kg = masa_actual_kg / exp(dv_impulso / (Isp_s * g0));
    h_actual_m = h0_m;
    hist_tiempo_s = [hist_tiempo_s; t_actual_s]; hist_altitud_m = [hist_altitud_m; h_actual_m];
end
fprintf('Simulacion de mision completada.\n\n');

%% ========================================================================
%  3. SIMULACIONES POST-MISION
%  ========================================================================
r_apogeo_deorbit = R_earth_m + hist_altitud_m(end);
r_perigeo_deorbit = R_earth_m + h_perigeo_deorbit_m;
v_circular_final = sqrt(mu / r_apogeo_deorbit);
a_transfer_deorbit = (r_apogeo_deorbit + r_perigeo_deorbit) / 2;
v_apogeo_transfer = sqrt(mu * (2/r_apogeo_deorbit - 1/a_transfer_deorbit));
dv_deorbit = v_circular_final - v_apogeo_transfer;
masa_antes_impulso = masa_actual_kg;
masa_despues_impulso = masa_antes_impulso / exp(dv_deorbit / (Isp_s * g0));
fprintf('Maniobra de Deorbitacion:\n  - Delta-V: %.2f m/s\n', dv_deorbit);

a_reentry_target = R_earth_m + h_reentry_km * 1e3;
options_reentry = odeset('Events', @(t,a) decayEvent(t, a, a_reentry_target), 'RelTol', 1e-7);

fprintf('Simulando reentrada con impulso...\n');
[t_reentry_con_impulso, a_reentry_con_impulso] = ode45(@(t,a) decayODE(t, a, R_earth_m, mu, Cd, Area_m2, masa_despues_impulso), ...
                               [0, inf], a_transfer_deorbit, options_reentry);
tiempo_reentrada_dias = t_reentry_con_impulso(end) / (24 * 3600);
fprintf('  - Reentrada estimada en %.2f dias.\n\n', tiempo_reentrada_dias);

fprintf('Simulando decaimiento natural (sin impulso)...\n');
a_initial_natural = R_earth_m + hist_altitud_m(end);
[t_reentry_sin_impulso, a_reentry_sin_impulso] = ode45(@(t,a) decayODE(t, a, R_earth_m, mu, Cd, Area_m2, masa_antes_impulso), ...
                               [0, inf], a_initial_natural, options_reentry);
tiempo_natural_dias = t_reentry_sin_impulso(end) / (24 * 3600);
fprintf('  - Reentrada natural estimada en %.2f dias (%.2f a~nos).\n\n', tiempo_natural_dias, tiempo_natural_dias/365.25);


%% ========================================================================
%  4. GENERACION DE GRAFICA PRINCIPAL (CICLO DE VIDA)
%  ========================================================================
hist_tiempo_total_s = [hist_tiempo_s; t_actual_s + t_reentry_con_impulso(2:end)];
hist_altitud_total_m = [hist_altitud_m; a_reentry_con_impulso(2:end) - R_earth_m];
hist_tiempo_years = hist_tiempo_total_s / (365.25 * 24 * 3600);
hist_altitud_km = hist_altitud_total_m / 1000;

figure('Name', 'Evolucion de la Altitud con Reentrada');
hold on; grid on;

plot(hist_tiempo_years, hist_altitud_km, 'b-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Trayectoria ($M_{seca}$=%dkg, A=%.1fm$^2$)', masa_seca_kg, Area_m2));
yline(h0_km, '--g', 'LineWidth', 1.5, 'DisplayName', 'Altitud Inicial');
yline(h_umbral_m/1e3, '--r', 'LineWidth', 1.5, 'DisplayName', 'Altitud Umbral');
xline(mission_duration_years, '--m', 'LineWidth', 1.5, 'DisplayName', 'Fin de Mision');
plot(mission_duration_years, hist_altitud_km(find(hist_tiempo_years >= mission_duration_years, 1)), ...
     'r*', 'MarkerSize', 10, 'LineWidth', 2, ...
     'DisplayName', sprintf('Impulso Deorbit ($\\Delta V$=%.1f m/s)', dv_deorbit));

title(sprintf('Evolucion de la Altitud durante %d a~nos y Reentrada', mission_duration_years), 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Tiempo (a~nos)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Altitud (km)', 'FontSize', 14, 'Interpreter', 'latex');
xlim([0, hist_tiempo_years(end) * 1.05]); ylim([0, h0_km * 1.05]);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
hold off;

%% ========================================================================
%  5. GENERACION DE GRAFICA COMPARATIVA DE REENTRADA 
%  ========================================================================
figure('Name', 'Comparativa de Tiempos de Reentrada');
hold on; grid on;

semilogx(t_reentry_con_impulso / (24*3600), (a_reentry_con_impulso - R_earth_m) / 1000, ...
    'r-', 'LineWidth', 2, 'DisplayName', 'Con Impulso de Deorbitacion');
semilogx(t_reentry_sin_impulso / (24*3600), (a_reentry_sin_impulso - R_earth_m) / 1000, ...
    'k--', 'LineWidth', 2, 'DisplayName', 'Decaimiento Natural (sin impulso)');

title('Comparativa de Tiempos de Reentrada', 'FontSize', 16, 'Interpreter', 'latex');
xlabel('Tiempo tras fin de mision (dias)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Altitud (km)', 'FontSize', 14, 'Interpreter', 'latex');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
hold off;


%% ========================================================================
%  FUNCIONES AUXILIARES
%  ========================================================================
function dadt = decayODE(~, a, R_earth, mu, Cd, A, mass)
    h = a - R_earth; if h < 0; h = 0; end
    rho = getAtmosphericDensity(h);
    dadt = -Cd * A * rho * sqrt(mu * a) / mass;
end
function [value, isterminal, direction] = decayEvent(~, a, a_target)
    value = a - a_target; isterminal = 1; direction = -1;
end
function dv = hohmannDeltaV(r1, r2, mu)
    v1 = sqrt(mu / r1); v2 = sqrt(mu / r2); a_transfer = (r1 + r2) / 2;
    v_transfer_1 = sqrt(mu * (2/r1 - 1/a_transfer));
    v_transfer_2 = sqrt(mu * (2/r2 - 1/a_transfer));
    dv = abs(v_transfer_1 - v1) + abs(v2 - v_transfer_2);
end
function rho = getAtmosphericDensity(h)
    h_km = max(h / 1e3, 0);
    h_base_km = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, ...
                 130, 140, 150, 180, 200, 250, 300, 350, 400, 450, ...
                 500, 600, 700, 800, 900, 1000];
    rho_base = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, ...
                3.206e-4, 8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7, ...
                9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9, ...
                5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, ...
                9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13, ...
                1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15];
    H_km = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, ...
            5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149, ...
            22.523, 29.740, 37.105, 45.546, 53.628, 53.298, ...
            58.515, 60.828, 63.822, 71.835, 88.667, 124.64, ...
            181.05, 268.00];

    idx = find(h_km >= h_base_km, 1, 'last');
    if isempty(idx)
        idx = 1;
    elseif idx == numel(h_base_km)
        idx = numel(h_base_km);
    end

    rho = rho_base(idx) * exp(-(h_km - h_base_km(idx)) / H_km(idx));
end
