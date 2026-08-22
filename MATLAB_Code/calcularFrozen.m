% calcularFrozen.m - Frozen orbit condition analysis for SSO orbits
clear; close all;

mu     = 3.986004418e14;
R_E    = 6378.137e3;
J2     = 1.08263e-3;
J3     = -2.5327e-6;
w_sol  = 1.991e-7;
lat_coverage = [25, 49]; % Limites de la franja de cobertura [deg N]
omega_frozen = pi/2;     % Argumento del perigeo usado en el caso frozen [rad]

h_vals       = 200:5:1000;
n_h          = length(h_vals);
i_vals       = zeros(1, n_h);
e_f_vals     = zeros(1, n_h);
dh_vals      = zeros(1, n_h);
dh_coverage_vals = zeros(1, n_h);
dh_pass_vals = zeros(1, n_h);
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

    % frozen eccentricity for w = 90 deg
    e_f_vals(k) = -J3 / (2*J2) * (R_E / a) * sind(i_vals(k));

    dh_vals(k)      = 2 * a * e_f_vals(k) / 1e3;
    u_coverage = deg2rad(asind(sind(lat_coverage) / sind(i_vals(k))));
    u_crossings = [u_coverage, pi - u_coverage];
    nu_crossings = u_crossings - omega_frozen;
    p_frozen = a * (1 - e_f_vals(k)^2);
    h_crossings = (p_frozen ./ (1 + e_f_vals(k) * cos(nu_crossings)) - R_E) / 1e3;
    dh_coverage_vals(k) = max(h_crossings) - min(h_crossings);
    dh_pass_vals(k) = max(abs(diff(h_crossings(1:2))), abs(diff(h_crossings(3:4))));
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
fprintf('Delta h coverage, e0 = 0       = %.2f km\n', 0);
fprintf('Delta h coverage, frozen       = %.2f km\n', dh_coverage_vals(idx_selected));
fprintf('Delta h within one pass, frozen = %.2f km\n', dh_pass_vals(idx_selected));
fprintf('Deadband  = %.1f km\n', deadband_vals(idx_selected));
fprintf('i_SSO     = %.2f deg\n', i_vals(idx_selected));
fprintf('Density   = %.1e kg/m^3\n', rho_vals(idx_selected));
fprintf('\n');

% --- altitudes where frozen IS compatible ---
compat = dh_coverage_vals <= deadband_vals;
if any(compat)
    h_first = h_vals(find(compat, 1, 'first'));
    fprintf('--- Frozen compatible over coverage (dh <= deadband) ---\n');
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

% --- plot ---
if ~strcmp(getenv('TFG_SKIP_FROZEN_PLOT'), '1')
figure('Position', [100 100 700 450]);
yyaxis left;
plot(h_vals, e_f_vals, 'b-', 'LineWidth', 1.2);
ylabel('e_{frozen}');
hold on;
xline(h_selected, 'r--', 'LineWidth', 1.2);
ylim auto;

yyaxis right;
plot(h_vals, dh_coverage_vals, 'k-', 'LineWidth', 1.2);
hold on;
plot(h_vals, deadband_vals, 'k--', 'LineWidth', 1.2);
ylabel('Altura [km]');

xlabel('Altitud [km]');
title('Condición frozen para órbita SSO');
legend('e_f', sprintf('h = %.0f km', h_selected), '\Delta h coverage', 'Deadband (2% h)', 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'Frozen_e_vs_h.pdf', 'ContentType', 'vector');
fprintf('\nSaved Frozen_e_vs_h.pdf\n');
end
