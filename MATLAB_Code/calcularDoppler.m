function [t, doppler_Hz, maxShift_Hz, maxRate_Hz_s] = calcularDoppler(h_km, f0_Hz, el_min_deg, exportPath)
% CALCULARDOPPLER  Perfil de desplazamiento Doppler durante un pase de satélite.
%
%   [t, doppler_Hz, maxShift_Hz, maxRate_Hz_s] = calcularDoppler(h_km, f0_Hz, el_min_deg)
%
%   Entradas:
%     h_km        - Altitud orbital [km]  (por defecto: 500)
%     f0_Hz       - Frecuencia portadora [Hz]  (por defecto: 8.1e9, banda X)
%     el_min_deg  - Elevación mínima de contacto [°]  (por defecto: 10)
%     exportPath  - Ruta donde guardar la figura PDF  (vacío = no exportar)
%
%   Salidas:
%     t            - Vector de tiempo del pase [s]
%     doppler_Hz   - Desplazamiento Doppler a lo largo del pase [Hz]
%     maxShift_Hz  - Desplazamiento máximo absoluto [Hz]
%     maxRate_Hz_s - Tasa de variación máxima (en TCA) [Hz/s]
%
%   Modelo: pase cenital (elevación máxima = 90°), Tierra plana local.
%   Válido para estimaciones de orden de magnitud en diseño preliminar.

if nargin < 1 || isempty(h_km),       h_km       = 500;   end
if nargin < 2 || isempty(f0_Hz),      f0_Hz      = 8.1e9; end
if nargin < 3 || isempty(el_min_deg), el_min_deg = 10;     end
if nargin < 4,                         exportPath = '';     end

%% === Constantes y parámetros orbitales ===
mu  = 3.986e14;          % [m³/s²] parámetro gravitacional terrestre
R_e = 6371e3;            % [m]     radio medio terrestre
c   = 3e8;               % [m/s]   velocidad de la luz

h = h_km * 1e3;                         % [m]    altitud orbital
r = R_e + h;                            % [m]    radio orbital
v = sqrt(mu / r);                       % [m/s]  velocidad orbital (circular)

% Semiduración del pase: tiempo desde el inicio (elevación = el_min) hasta TCA
t_max = h / (v * tand(el_min_deg));     % [s]  modelo de Tierra plana local

%% === Cálculo del perfil de Doppler ===
t = linspace(-t_max, t_max, 2000);     % [s]  vector de tiempo centrado en TCA

% Rango satélite-GS (aproximación planar)
d = sqrt(h^2 + (v .* t).^2);          % [m]

% Velocidad radial (+ = alejándose)
d_rate = v^2 .* t ./ d;               % [m/s]

% Desplazamiento Doppler: Δf = -f0 · v_r / c (< 0 cuando se acerca)
doppler_Hz = -f0_Hz .* d_rate / c;

maxShift_Hz  = max(abs(doppler_Hz));
maxRate_Hz_s = max(abs(diff(doppler_Hz) ./ diff(t)));

%% === Resumen en consola ===
fprintf('\n=== Desplazamiento Doppler — Banda X (%.2f GHz, h = %d km) ===\n', ...
        f0_Hz/1e9, h_km);
fprintf('  Velocidad orbital:            %.1f km/s\n',  v/1e3);
fprintf('  Semiduración del pase:        %.0f s  (elevación mín. %d°)\n', t_max, el_min_deg);
fprintf('  Desplazamiento máximo:        ±%.0f Hz  (±%.1f kHz)\n', maxShift_Hz, maxShift_Hz/1e3);
fprintf('  Tasa de variación máx. (TCA): %.1f Hz/s  (%.2f kHz/s)\n', maxRate_Hz_s, maxRate_Hz_s/1e3);

%% === Figura ===
fig = figure('Name', 'Perfil Doppler', 'NumberTitle', 'off');
plot(t, doppler_Hz/1e3, 'b-', 'LineWidth', 1.5);
hold on;
yline(0, 'k--', 'LineWidth', 0.8);
xline(0, 'r--', 'LineWidth', 0.8, 'Label', 'TCA', 'LabelHorizontalAlignment', 'left');
xlabel('Tiempo relativo al TCA [s]');
ylabel('\Deltaf  [kHz]');
title(sprintf('Perfil de desplazamiento Doppler — Banda X (%.1f GHz, h = %d km)', ...
              f0_Hz/1e9, h_km));
legend(sprintf('\\Deltaf_{max} = \\pm%.0f kHz', maxShift_Hz/1e3), 'Location', 'northeast');
grid on;
box on;

if ~isempty(exportPath)
    exportgraphics(fig, exportPath, 'ContentType', 'vector');
    fprintf('  Figura exportada a: %s\n', exportPath);
end
