% calcularFrozen.m - Frozen orbit condition analysis for SSO orbits
clear; close all;

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

    % frozen eccentricity for w = 90 deg
    e_f_vals(k) = -J3 / (2*J2) * (R_E / a) * sind(i_vals(k));

    dh_vals(k)      = 2 * a * e_f_vals(k) / 1e3;
    deadband_vals(k) = 0.02 * h;

    % interpolate density
    rho_vals(k) = exp(interp1(h_rho, log(rho_ref), h, 'linear', 'extrap'));
end

% --- print 520 km results ---
idx520 = find(h_vals == 520);
fprintf('--- Frozen condition at 520 km ---\n');
fprintf('e_frozen  = %.2e\n', e_f_vals(idx520));
fprintf('Delta h   = %.1f km\n', dh_vals(idx520));
fprintf('Deadband  = %.1f km\n', deadband_vals(idx520));
fprintf('i_SSO     = %.2f deg\n', i_vals(idx520));
fprintf('Density   = %.1e kg/m^3\n', rho_vals(idx520));
fprintf('\n');

% --- altitudes where frozen IS compatible ---
compat = dh_vals <= deadband_vals;
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
    fprintf('Density at 520 km = %.1e kg/m^3\n', rho_vals(idx520));
    fprintf('Density at 800 km = %.1e kg/m^3\n', rho_vals(idx800));
    fprintf('Ratio (520/800)  = %.1f\n', rho_vals(idx520)/rho_vals(idx800));
end

% --- plot ---
figure('Position', [100 100 700 450]);
yyaxis left;
plot(h_vals, e_f_vals, 'b-', 'LineWidth', 1.2);
ylabel('e_{frozen}');
hold on;
xline(520, 'r--', 'LineWidth', 1.2);
ylim auto;

yyaxis right;
plot(h_vals, dh_vals, 'k-', 'LineWidth', 1.2);
hold on;
plot(h_vals, deadband_vals, 'k--', 'LineWidth', 1.2);
ylabel('Altura [km]');

xlabel('Altitud [km]');
title('Condición frozen para órbita SSO');
legend('e_f', 'h = 520 km', '\Delta h', 'Deadband (2% h)', 'Location', 'northwest');
grid on;
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'Frozen_e_vs_h.pdf', 'ContentType', 'vector');
fprintf('\nSaved Frozen_e_vs_h.pdf\n');
