% calcularFrozen.m - Frozen orbit condition analysis for SSO orbits
clear; close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'lib'));

mu     = 3.986004418e14;
R_E    = 6378.137e3;
J2     = 1.08263e-3;
J3     = -2.5327e-6;
w_sol  = 1.991e-7;

h_vals       = 200:5:1000;
n_h          = length(h_vals);
i_vals       = zeros(1, n_h);
e_f_vals     = zeros(1, n_h);
dh_vals      = zeros(1, n_h);
deadband_vals = zeros(1, n_h);
rho_vals     = zeros(1, n_h);

% piecewise exponential density model (kg/m^3)
h_rho = [200, 300, 400, 500, 600, 700, 800];
rho_ref = [2.5e-12, 1.0e-12, 4e-13, 1.5e-13, 6e-14, 2e-14, 1e-14];

for k = 1:n_h
    h    = h_vals(k);
    a    = R_E + h*1e3;
    r    = a;

    % SSO inclination
    ci = -(2 * r^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
    ci = max(-1, min(1, ci));
    i_vals(k) = acosd(ci);

    % First-order frozen eccentricity for e << 1 and omega = 90 deg
    e_f_vals(k) = -J3 / (2*J2) * (R_E / a) * sind(i_vals(k));

    dh_vals(k)      = 2 * a * e_f_vals(k) / 1e3;
    deadband_vals(k) = 0.02 * h;

    % interpolate density
    rho_vals(k) = exp(interp1(h_rho, log(rho_ref), h, 'linear', 'extrap'));
end

% --- print selected-altitude results ---
h_selected = 630;
idx_selected = find(h_vals == h_selected);
fprintf('--- Frozen condition at %.0f km ---\n', h_selected);
fprintf('e_frozen  = %.2e\n', e_f_vals(idx_selected));
fprintf('Delta h   = %.1f km\n', dh_vals(idx_selected));
fprintf('Deadband half-width = %.1f km\n', deadband_vals(idx_selected));
fprintf('Deadband full width = %.1f km\n', 2*deadband_vals(idx_selected));
fprintf('i_SSO     = %.2f deg\n', i_vals(idx_selected));
fprintf('Density   = %.1e kg/m^3\n', rho_vals(idx_selected));
fprintf('\n');

% --- altitudes where frozen IS compatible ---
% dh is the full apogee-perigee excursion; deadband is one-sided.
compat = dh_vals <= 2*deadband_vals;
if any(compat)
    h_first = h_vals(find(compat, 1, 'first'));
    fprintf('--- Frozen compatible (dh <= deadband) ---\n');
    fprintf('First compatible altitude: %.0f km\n', h_first);
    fprintf('Compatible altitudes: %.0f to %.0f km\n', ...
        h_first, h_vals(find(compat, 1, 'last')));
else
    fprintf('No altitudes with dh <= deadband in range.\n');
end

% density ratio
fprintf('\n--- Density comparison ---\n');
idx800 = find(h_vals == 800);
if ~isempty(idx800)
    fprintf('Density at %.0f km = %.1e kg/m^3\n', h_selected, rho_vals(idx_selected));
    fprintf('Density at 800 km = %.1e kg/m^3\n', rho_vals(idx800));
    fprintf('Ratio (%.0f/800)  = %.1f\n', h_selected, rho_vals(idx_selected)/rho_vals(idx800));
end

% --- Height variation over the acquisition arc at the critical latitude ---
% The coverage model uses the lower reference latitude of continental USA
% and the effective swath selected by the payload analysis.
lat_critical = 25;             % geodetic latitude [deg]
swath_nominal = 190;           % nominal swath [km]
swath_effective = 0.95 * swath_nominal;
psi = atand(swath_effective / (2 * h_selected));
a_selected = R_E + h_selected * 1e3;
inc_selected = deg2rad(i_vals(idx_selected));
omega_perigee = deg2rad(90);   % frozen-orbit orientation
coes_circular = [a_selected, 0, inc_selected, 0, omega_perigee];
lat_model = atand((1 - 1/298.257223563)^2 * tand(lat_critical));
[theta, ~, ~, ~] = calcTheta(lat_model, coes_circular, 0, 1, 0, [], psi);
lat_arc = linspace(lat_model - theta, lat_model + theta, 100);

height_cases = struct( ...
    'name', {'frozen', 'circular'}, ...
    'eccentricity', {e_f_vals(idx_selected), 0});
height_rows = zeros(2, 9);
height_curves = cell(2, 2);

for case_idx = 1:numel(height_cases)
    e_case = height_cases(case_idx).eccentricity;
    pass_heights = zeros(2, numel(lat_arc));

    for pass_idx = 1:2
        if pass_idx == 1
            lat_pass = lat_arc;
            nu = deg2rad(asind(sind(lat_pass) / sind(i_vals(idx_selected)))) - omega_perigee;
        else
            lat_pass = fliplr(lat_arc);
            nu = pi - deg2rad(asind(sind(lat_pass) / sind(i_vals(idx_selected)))) - omega_perigee;
        end

        r_orbit = a_selected * (1 - e_case^2) ./ ...
            (1 + e_case * cos(nu));
        R_ellipsoid = ellipsoidRadius(R_E, lat_pass);
        pass_heights(pass_idx, :) = (r_orbit - R_ellipsoid) / 1e3;
        height_curves{case_idx, pass_idx} = pass_heights(pass_idx, :);
    end

    h_min_obs = min(pass_heights, [], 'all');
    h_max_obs = max(pass_heights, [], 'all');
    delta_h_obs = h_max_obs - h_min_obs;
    max_deviation = max(abs(pass_heights - h_selected), [], 'all');
    height_rows(case_idx, :) = [h_selected, height_cases(case_idx).eccentricity, ...
        theta, h_min_obs, h_max_obs, delta_h_obs, max_deviation, ...
        delta_h_obs / h_selected, max_deviation / h_selected];

    fprintf('\n--- Height variation over acquisition arc (%s) ---\n', height_cases(case_idx).name);
    fprintf('Critical latitude: %.1f deg N, effective swath: %.1f km\n', ...
        lat_critical, swath_effective);
    fprintf('Acquisition half-angle theta: %.4f deg\n', theta);
    fprintf('Observed height range: %.6f - %.6f km\n', h_min_obs, h_max_obs);
    fprintf('Observed height variation: %.6f km\n', delta_h_obs);
    fprintf('Maximum deviation from nominal: %.6f km\n', max_deviation);
end

% Save reproducible results and a compact figure for the thesis.
analysis_dir = fullfile(script_dir, '..', 'Compiled data', 'OptimumConfigs');
latex_figure_dir = fullfile(script_dir, '..', 'Latex_Code', '5.Mission');
if ~exist(analysis_dir, 'dir'), mkdir(analysis_dir); end
if ~exist(latex_figure_dir, 'dir'), mkdir(latex_figure_dir); end

height_table = array2table(height_rows, 'VariableNames', { ...
    'NominalAltitude_km', 'Eccentricity', 'AcquisitionHalfAngle_deg', ...
    'ObservedMinAltitude_km', 'ObservedMaxAltitude_km', ...
    'ObservedVariation_km', 'MaxDeviationFromNominal_km', ...
    'RelativeVariation', 'RelativeMaxDeviation'});
writetable(height_table, fullfile(analysis_dir, ...
    'frozen_coverage_height_analysis.csv'));

figure('Position', [100 100 760 460]);
hold on;
plot(lat_arc, height_curves{1, 1}, 'b-', 'LineWidth', 1.4, ...
    'DisplayName', 'Frozen, ascendente');
plot(lat_arc, height_curves{1, 2}, 'b--', 'LineWidth', 1.4, ...
    'DisplayName', 'Frozen, descendente');
plot(lat_arc, height_curves{2, 1}, 'r-', 'LineWidth', 1.4, ...
    'DisplayName', 'Circular, ascendente');
plot(lat_arc, height_curves{2, 2}, 'r--', 'LineWidth', 1.4, ...
    'DisplayName', 'Circular, descendente');
yline(h_selected, 'k:', 'LineWidth', 1.1, 'DisplayName', 'Altura nominal');
xline(lat_model, 'k-.', 'LineWidth', 1.0, 'HandleVisibility', 'off');
hold off;
grid on;
xlabel('Latitud del arco de adquisición [grados]');
ylabel('Altura sobre el elipsoide [km]');
title('Variación de altura sobre la zona de adquisición');
legend('Location', 'eastoutside');
set(gca, 'FontSize', 10);
exportgraphics(gcf, fullfile(latex_figure_dir, ...
    'frozen_coverage_height_analysis.pdf'), 'ContentType', 'vector');
close(gcf);

fprintf('\nSaved frozen_coverage_height_analysis.csv and frozen_coverage_height_analysis.pdf\n');

% --- Mean eccentricity-vector evolution over the mission ---
% First-order J2/J3 mean-element model. The semimajor axis decay and its
% maintenance impulses are evaluated independently in calcularMasaTotal.
mission_days = 8 * 365.25;
t_days = (0:mission_days).';
t_seconds = t_days * 86400;
p_selected = a_selected * (1 - e_f_vals(idx_selected)^2);
e_center = -J3 / (2 * J2) * (R_E / p_selected) * sin(inc_selected);

evolution_cases = struct( ...
    'name', {'frozen', 'circular'}, ...
    'initial_e', {e_center, 0}, ...
    'initial_omega', {omega_perigee, omega_perigee});
evolution_rows = zeros(numel(t_days) * numel(evolution_cases), 11);
summary_rows = zeros(numel(evolution_cases), 10);
summary_names = cell(numel(evolution_cases), 1);
row = 0;

for case_idx = 1:numel(evolution_cases)
    [ex, ey] = meanEccentricityVector(t_seconds, ...
        evolution_cases(case_idx).initial_e, ...
        evolution_cases(case_idx).initial_omega, ...
        a_selected, inc_selected, mu, R_E, J2, J3);
    e_history = hypot(ex, ey);
    omega_history = atan2(ey, ex);
    omega_history(e_history < 1e-12) = omega_perigee;

    h_min_history = zeros(size(t_days));
    h_max_history = zeros(size(t_days));
    delta_h_history = zeros(size(t_days));
    max_deviation_history = zeros(size(t_days));

    for time_idx = 1:numel(t_days)
        [h_min_history(time_idx), h_max_history(time_idx)] = ...
            acquisitionHeightEnvelope(e_history(time_idx), omega_history(time_idx), ...
            lat_arc, i_vals(idx_selected), a_selected, R_E);
        delta_h_history(time_idx) = h_max_history(time_idx) - h_min_history(time_idx);
        max_deviation_history(time_idx) = max(abs([h_min_history(time_idx), ...
            h_max_history(time_idx)] - h_selected));

        row = row + 1;
        evolution_rows(row, :) = [t_days(time_idx), case_idx, ex(time_idx), ey(time_idx), ...
            e_history(time_idx), rad2deg(omega_history(time_idx)), ...
            h_min_history(time_idx), h_max_history(time_idx), delta_h_history(time_idx), ...
            max_deviation_history(time_idx), max_deviation_history(time_idx) <= deadband_vals(idx_selected)];
    end

    summary_names{case_idx} = evolution_cases(case_idx).name;
    summary_rows(case_idx, :) = [evolution_cases(case_idx).initial_e, ...
        min(e_history), max(e_history), min(h_min_history), max(h_max_history), ...
        max(delta_h_history), max(max_deviation_history), ...
        sum(max_deviation_history > deadband_vals(idx_selected)), ...
        max(e_history) - min(e_history), ...
        mean(max_deviation_history <= deadband_vals(idx_selected))];

    fprintf('\n--- Mean eccentricity evolution (%s) ---\n', evolution_cases(case_idx).name);
    fprintf('e range: %.8e - %.8e\n', min(e_history), max(e_history));
    fprintf('Maximum observed height variation: %.6f km\n', max(delta_h_history));
    fprintf('Maximum deviation from nominal: %.6f km\n', max(max_deviation_history));
    fprintf('Days outside maintenance half-band: %.0f\n', ...
        sum(max_deviation_history > deadband_vals(idx_selected)));
end

assert(all(isfinite(evolution_rows), 'all'), 'Non-finite secular evolution result.');
assert(summary_rows(1, 7) < summary_rows(2, 7), ...
    'Frozen case should have the smaller maximum height deviation.');
assert(all(summary_rows(:, 8) == 0), ...
    'A propagated case leaves the maintenance half-band.');

evolution_table = array2table(evolution_rows, 'VariableNames', { ...
    'MissionDay', 'Case', 'EccentricityVectorX', 'EccentricityVectorY', ...
    'Eccentricity', 'ArgumentOfPerigee_deg', 'ObservedMinAltitude_km', ...
    'ObservedMaxAltitude_km', 'ObservedVariation_km', ...
    'MaxDeviationFromNominal_km', 'WithinMaintenanceHalfBand'});
writetable(evolution_table, fullfile(analysis_dir, ...
    'frozen_secular_evolution_8yr.csv'));

summary_table = array2table(summary_rows, 'VariableNames', { ...
    'InitialEccentricity', 'MinEccentricity', 'MaxEccentricity', ...
    'GlobalObservedMinAltitude_km', 'GlobalObservedMaxAltitude_km', ...
    'MaximumObservedVariation_km', 'MaximumDeviationFromNominal_km', ...
    'DaysOutsideMaintenanceHalfBand', 'EccentricityRange', ...
    'FractionOfDaysWithinHalfBand'});
summary_table.Case = summary_names;
summary_table = movevars(summary_table, 'Case', 'Before', 1);
writetable(summary_table, fullfile(analysis_dir, ...
    'frozen_secular_evolution_8yr_summary.csv'));

figure('Position', [100 100 820 650]);
tiledlayout(2, 1);
nexttile;
hold on;
for case_idx = 1:numel(evolution_cases)
    case_rows = evolution_rows(evolution_rows(:, 2) == case_idx, :);
    plot(case_rows(:, 1) / 365.25, case_rows(:, 5), 'LineWidth', 1.3, ...
        'DisplayName', evolution_cases(case_idx).name);
end
hold off;
grid on;
xlabel('Tiempo desde el inicio [años]');
ylabel('Excentricidad media');
legend('Location', 'best');
title('Evolución de la excentricidad media');

nexttile;
hold on;
for case_idx = 1:numel(evolution_cases)
    case_rows = evolution_rows(evolution_rows(:, 2) == case_idx, :);
    plot(case_rows(:, 1) / 365.25, case_rows(:, 10), 'LineWidth', 1.3, ...
        'DisplayName', evolution_cases(case_idx).name);
end
yline(deadband_vals(idx_selected), 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'Semibanda de mantenimiento');
hold off;
grid on;
xlabel('Tiempo desde el inicio [años]');
ylabel('Desviación máxima de altura [km]');
legend('Location', 'best');
title('Desviación de altura sobre la zona de adquisición');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 10);
exportgraphics(gcf, fullfile(latex_figure_dir, ...
    'frozen_secular_evolution_8yr.pdf'), 'ContentType', 'vector');
close(gcf);

fprintf('\nSaved frozen_secular_evolution_8yr.csv, summary CSV and figure\n');

% --- plot ---
figure('Position', [100 100 700 450]);
yyaxis left;
plot(h_vals, e_f_vals, 'b-', 'LineWidth', 1.2);
ylabel('e_{frozen}');
hold on;
xline(h_selected, 'r--', 'LineWidth', 1.2);
ylim auto;

yyaxis right;
plot(h_vals, dh_vals, 'k-', 'LineWidth', 1.2);
hold on;
plot(h_vals, deadband_vals, 'k--', 'LineWidth', 1.2);
ylabel('Altura [km]');

xlabel('Altitud [km]');
title('Condición frozen para órbita SSO');
legend('e_f', sprintf('h = %.0f km', h_selected), '\Delta h', 'Deadband half-width (\pm2% h)', 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'Frozen_e_vs_h.pdf', 'ContentType', 'vector');
fprintf('\nSaved Frozen_e_vs_h.pdf\n');

function radius = ellipsoidRadius(equatorial_radius, latitude)
    flattening = 1 / 298.257223563;
    polar_radius = equatorial_radius * (1 - flattening);
    radius = sqrt(((equatorial_radius^2 .* cosd(latitude)).^2 + ...
        (polar_radius^2 .* sind(latitude)).^2) ./ ...
        ((equatorial_radius .* cosd(latitude)).^2 + ...
        (polar_radius .* sind(latitude)).^2));
end

function [ex, ey] = meanEccentricityVector(t_seconds, e_initial, omega_initial, ...
        a, inc, mu, R_E, J2, J3)
    p = a * (1 - e_initial^2);
    e_frozen = -J3 / (2 * J2) * (R_E / p) * sin(inc);
    center_y = e_frozen;
    n = sqrt(mu / a^3);
    alpha_rate = -3 * J2 * (R_E / p)^2 * ...
        (5/4 * sin(inc)^2 - 1) * n;
    alpha = alpha_rate * t_seconds;

    initial_x = e_initial * cos(omega_initial);
    initial_y = e_initial * sin(omega_initial);
    displacement_x = initial_x;
    displacement_y = initial_y - center_y;

    ex = displacement_x .* cos(alpha) - displacement_y .* sin(alpha);
    ey = center_y + displacement_x .* sin(alpha) + displacement_y .* cos(alpha);
end

function [h_min, h_max] = acquisitionHeightEnvelope(eccentricity, omega, ...
        lat_arc, inclination_deg, a, R_E)
    h_passes = zeros(2, numel(lat_arc));
    for pass_idx = 1:2
        if pass_idx == 1
            lat_pass = lat_arc;
            u = deg2rad(asind(sind(lat_pass) / sind(inclination_deg)));
        else
            lat_pass = fliplr(lat_arc);
            u = pi - deg2rad(asind(sind(lat_pass) / sind(inclination_deg)));
        end
        nu = u - omega;
        r_orbit = a * (1 - eccentricity^2) ./ ...
            (1 + eccentricity * cos(nu));
        R_ellipsoid = ellipsoidRadius(R_E, lat_pass);
        h_passes(pass_idx, :) = (r_orbit - R_ellipsoid) / 1e3;
    end
    h_min = min(h_passes, [], 'all');
    h_max = max(h_passes, [], 'all');
end
