function CoverageRevisitCalc(GSD, alturas_orbitales, swaths_km, Cov_Requirement, current_Npix, current_fov_limit, N_sat, N_telescopes, detector_id, telescope_name, solapamiento, cobertura_nubes, max_detectores)

% --- Constantes y Parámetros Orbitales ---
R_E = 6378.137e3; % Radio de la Tierra [m]
mu_E = 3.986004418e14; % Constante gravitacional terrestre [m3/s2]
J2 = 1.082629989052e-3; % Achatamiento J2
SSO_nodal_rate = 360/(365.2421897*24*3600) * (pi/180); % Tasa de precesión para SSO [rad/s]

% --- Configuración de la Simulación ---
config_name = sprintf('%d Satelite(s); %d Telescopio(s): %s; Detector %d', N_sat, N_telescopes, strrep(telescope_name, ' ', ''), detector_id);
lat = 40; % Latitud objetivo [grados]
dayLimit = [1, Cov_Requirement]; % Límites para el cálculo de revisita
fprintf('Iniciando análisis para: %s\n', strrep(config_name, '_', ' '));

% --- Cálculo de Límites de Swath para Visualización ---
detector_swath_limits = (1:max_detectores) * GSD * current_Npix / 1000; % [km]

% --- Inicialización de la Matriz de Resultados ---
coverage_days = zeros(length(alturas_orbitales), length(swaths_km));

% --- Bucle Principal de Cálculo ---
for h = 1:length(alturas_orbitales)
    for s = 1:length(swaths_km)
        height = alturas_orbitales(h);
        swath = swaths_km(s);

        % --- Verificación de Restricciones Físicas ---
        effective_fov_limit = current_fov_limit;
        if N_telescopes == 2
            effective_fov_limit = current_fov_limit * 2;
        end
        max_swath_by_fov = 2 * height * tand(effective_fov_limit / 2);
        max_swath_by_detector = max_detectores * GSD * current_Npix / 1000; % [km]
        
        if swath > max_swath_by_fov || swath > max_swath_by_detector
            coverage_days(h,s) = NaN;
            continue;
        end

        % --- Cálculo del Swath Efectivo (con solapamiento) ---
        effective_swath = swath * (1 - solapamiento);

        % --- Parámetros Orbitales para SSO ---
        a = (R_E + height * 1000); % Semieje mayor [m]
        n = sqrt(mu_E / a^3); % Movimiento medio [rad/s]
        cos_i = -SSO_nodal_rate * (a/R_E)^(7/2) / (1.5 * n * J2);
        if abs(cos_i) > 1
            coverage_days(h,s) = NaN; % SSO no posible a esta altura
            continue;
        end
        inc = acos(cos_i);
        coes = [a, 0, inc, 0, 0]; % [sma, ecc, inc, RAAN, AoP]

        % --- Llamada a RevisitCalc ---
        half_fov_angle = atand(effective_swath / (2 * height)); % [grados]
        try
            maxRevisit = RevisitCalc(coes, lat, 1, 0, N_sat, 1, 0, 1, 0, dayLimit, half_fov_angle);
            final_revisit = maxRevisit / (1 - cobertura_nubes);
            if final_revisit > Cov_Requirement
                coverage_days(h,s) = NaN;
            else
                coverage_days(h,s) = final_revisit;
            end
        catch
            coverage_days(h,s) = NaN;
        end
    end
end

% --- Generación de Gráficos y Guardado de Datos ---
output_dir = 'coverage';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
fig = figure('Visible', 'off');
custom_cmap = [1 1 1; parula(256)];
colormap(custom_cmap);
imagesc(swaths_km, alturas_orbitales, coverage_days);
set(gca, 'Color', [1 1 1]); % Fondo blanco para áreas NaN
hold on;

% Calcula y dibuja la línea de límite de FOV
effective_fov_limit_plot = current_fov_limit;
if N_telescopes == 2
    effective_fov_limit_plot = current_fov_limit * 2;
end
max_swath_by_fov_curve = 2 * alturas_orbitales .* tand(effective_fov_limit_plot / 2);
plot(max_swath_by_fov_curve, alturas_orbitales, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Límite FOV');

% Dibujar líneas de referencia para los límites de swath de los detectores
for i = 1:max_detectores
    swath_limit = detector_swath_limits(i);
    if swath_limit <= max(swaths_km) && swath_limit >= min(swaths_km)
        xline(swath_limit, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        text(swath_limit + 2, max(alturas_orbitales) - 15, sprintf('%d Det', i), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'Color', 'black', 'FontWeight', 'bold');
    end
end

legend('show', 'Location', 'southeast');
hold off;

% Formato del gráfico
title_text = sprintf('Tiempo Revisita - %s\n(Blanco: $\\geq$ %d dias)', strrep(config_name, '_', ' '), Cov_Requirement);
title(title_text, 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Swath (km)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Altura orbital (km)', 'Interpreter', 'latex', 'FontSize', 12);
c = colorbar;
ylabel(c, 'Tiempo revisita (dias)', 'Interpreter', 'latex', 'FontSize', 11);
clim([1 Cov_Requirement]);
set(gca, 'YDir', 'normal');
axis tight;

% Guardar archivos
png_filename = fullfile(output_dir, sprintf('heatmap_%s.jpg', config_name));
csv_filename = fullfile(output_dir, sprintf('coverage_%s.csv', config_name));
print(fig, png_filename, '-dpng', '-r300');
writematrix(coverage_days, csv_filename);
close(fig);
fprintf('Análisis completado para: %s. Archivos guardados en ''%s''.\n', strrep(config_name, '_', ' '), output_dir);
end
