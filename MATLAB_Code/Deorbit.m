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
decay_percentage_drop = 1.5;  % Porcentaje de decaimiento para activar impulso (%)
h_perigeo_deorbit_km = 120; % Altitud del perigeo para la maniobra de deorbitacion (km)
h_reentry_km = 100;      % Altitud a la que se considera reentrada/desintegracion (km)

% --- Parametros del Satelite ---
masa_seca_kg = 1.34;    % Masa seca del satelite (kg)
Area_m2 = 0.07;           % Area transversal efectiva para el arrastre (m^2)
Cd = 2.5;               % Coeficiente de arrastre (adimensional)
Isp_s = 220;            % Impulso especifico del sistema de propulsion (s)


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
    h_km = h / 1e3;
    if h_km < 100,      rho0=1.225;   H=6.7e3; elseif h_km < 150,  rho0=4.79e-7; H=9.5e3;
    elseif h_km < 200,  rho0=1.81e-9; H=25.5e3;elseif h_km < 250,  rho0=2.53e-10;H=37.5e3;
    elseif h_km < 300,  rho0=6.24e-11;H=44.8e3;elseif h_km < 350,  rho0=7.40e-9; H=50.3e3;
    elseif h_km < 400,  rho0=6.98e-12;H=54.8e3;elseif h_km < 450,  rho0=2.72e-12;H=58.2e3;
    elseif h_km < 500,  rho0=1.13e-12;H=61.3e3;elseif h_km < 600,  rho0=5e-13;   H=70e3;
    elseif h_km < 700,  rho0=1e-13;   H=80e3;  elseif h_km < 800,  rho0=2e-14;   H=90e3;
    elseif h_km < 900,  rho0=5e-15;   H=100e3; elseif h_km < 1000, rho0=1e-15;   H=110e3;
    elseif h_km < 1500, rho0=1e-16;   H=150e3; elseif h_km < 2000, rho0=1e-17;   H=200e3;
    else; rho = 1e-18; return; end
    base_h_km = floor(h_km/50)*50;
    rho = rho0 * exp(-(h - base_h_km*1e3) / H);
end
