% calcularRGT.m - RGT exact altitudes for SSO orbits with J2
clear; close all;

mu        = 3.986004418e14;
R_E       = 6378.137e3;
J2        = 1.08263e-3;
omega_E   = 7.292115373194e-5; % Earth rotation [rad/s]
omega_sol = 1.991e-7;           % Apparent solar motion [rad/s]
T_day     = 86400;              % Reference day for drift [s]

h_vals  = 200:2:1000;  % km
n_h     = length(h_vals);
i_SSO      = zeros(1,n_h);
omega_nodal = zeros(1,n_h);
dotOmega   = zeros(1,n_h);
T_nodal    = zeros(1,n_h);

for k = 1:n_h
    h = h_vals(k);
    [i_SSO(k), omega_nodal(k), dotOmega(k)] = ...
        rgt_rates(h*1e3, mu, R_E, J2, omega_sol);
    T_nodal(k) = 2*pi / omega_nodal(k);
end

fprintf('--- Reference at 630 km ---\n');
h_ref = 630;
[i_ref, omega_nodal_ref, dotOmega_ref] = ...
    rgt_rates(h_ref*1e3, mu, R_E, J2, omega_sol);
Tn_ref = 2*pi / omega_nodal_ref;
fprintf('i_SSO(630 km)      = %.2f deg\n', i_ref);
fprintf('T_nodal(630 km)     = %.2f s\n', Tn_ref);
fprintf('omega_nodal/(omega_E-dotOmega) = %.4f\n', ...
    omega_nodal_ref / (omega_E - dotOmega_ref));
fprintf('\n');

fprintf('--- Exact RGT altitudes (D <= 30) ---\n');
fprintf('  h [km]    i_SSO [deg]    N     D    T_nodal [s]    N/D\n');
fprintf('---------------------------------------------------------\n');

tol   = 1e-10;
maxit = 100;
roots = [];

ratio_rgt = omega_nodal ./ (omega_E - dotOmega);
for D = 1:30
    N_min = ceil(D * min(ratio_rgt));
    N_max = floor(D * max(ratio_rgt));
    for N = N_min:N_max
        ha = 200e3;  % lower bound in m
        hb = 1000e3; % upper bound in m

        fa = f_RGT(ha, N, D, mu, R_E, J2, omega_E, omega_sol);
        fb = f_RGT(hb, N, D, mu, R_E, J2, omega_E, omega_sol);

        if fa * fb > 0
            continue;
        end

        for iter = 1:maxit
            hc = (ha + hb) / 2;
            fc = f_RGT(hc, N, D, mu, R_E, J2, omega_E, omega_sol);
            if abs(fc) < tol
                break;
            end
            if fa * fc < 0
                hb = hc; fb = fc;
            else
                ha = hc; fa = fc;
            end
        end

        h_exact = hc / 1e3;
        [i_exact, omega_nodal_exact] = ...
            rgt_rates(hc, mu, R_E, J2, omega_sol);
        Tn_exact = 2*pi / omega_nodal_exact;

        roots = [roots; h_exact, i_exact, N, D, Tn_exact, N/D];
        fprintf('%9.3f   %8.2f     %4d   %2d   %10.2f    %.4f\n', ...
            h_exact, i_exact, N, D, Tn_exact, N/D);
    end
end

roots = sortrows(roots, 1);

h_nearest = 630;
[~, idx] = min(abs(roots(:,1) - h_nearest));
fprintf('\n--- Nearest to %.0f km ---\n', h_nearest);
fprintf('h = %.3f km, i = %.2f deg, N = %d, D = %d, T_n = %.2f s\n', ...
    roots(idx,1), roots(idx,2), roots(idx,3), roots(idx,4), roots(idx,5));
fprintf('Delta h = %.1f m\n', (roots(idx,1) - h_nearest)*1e3);

% Daily longitudinal drift
residual = zeros(1, n_h);
for k = 1:n_h
    ratio_k = ratio_rgt(k);
    best_err = inf;
    for D = 1:30
        N_cand = round(ratio_k * D);
        err = abs(omega_nodal(k) - (omega_E - dotOmega(k)) * ...
            N_cand / D);
        if err < best_err
            best_err = err;
        end
    end
    residual(k) = best_err * R_E * T_day / 1e3;  % km/day
end

figure('Position', [100 100 700 420]);
semilogy(h_vals, residual, 'b-', 'LineWidth', 1.2);
hold on;
for k = 1:size(roots,1)
    [~, ih] = min(abs(h_vals - roots(k,1)));
    plot(roots(k,1), residual(ih), 'g.', 'MarkerSize', 18);
end
xline(h_nearest, 'r--', 'LineWidth', 1.5);
xlabel('Altura [km]');
ylabel('Deriva longitudinal diaria [km/día]');
title('Condición RGT para órbita SSO con J_2');
grid on;
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'RGT_residual.pdf', 'ContentType', 'vector');
fprintf('\nSaved RGT_residual.pdf\n');

function [inc, omega_nodal, dotOmega] = rgt_rates(h_m, mu, R_E, J2, omega_sol)
    a = R_E + h_m;
    n = sqrt(mu / a^3);
    ci = -(2 * a^(7/2) * omega_sol) / (3 * J2 * R_E^2 * sqrt(mu));
    ci = max(-1, min(1, ci));
    inc = acosd(ci);

    p = a; % Circular mission orbit: e = 0
    dotOmega = -(3 * J2 * n * R_E^2 / (2 * p^2)) * cosd(inc);
    dotomega = (3 * J2 * n * R_E^2 / (4 * p^2)) * ...
        (4 - 5 * sind(inc)^2);
    dotM1 = (3 * J2 * n * R_E^2 / (4 * p^2)) * ...
        (2 - 3 * sind(inc)^2);
    omega_nodal = n + dotM1 + dotomega;
end

function out = f_RGT(h_m, N, D, mu, R_E, J2, omega_E, omega_sol)
    [~, omega_nodal, dotOmega] = ...
        rgt_rates(h_m, mu, R_E, J2, omega_sol);
    out = omega_nodal - (omega_E - dotOmega) * ...
        N / D;
end
