% calcularRGT.m - RGT exact altitudes for SSO orbits with J2
clear; close all;

mu    = 3.986004418e14;
R_E   = 6378.137e3;
J2    = 1.08263e-3;
T_rot = 86164.0905; % Periodo sideral de rotacion terrestre [s]
omega_E = 2*pi/T_rot; % Velocidad angular sideral terrestre [rad/s]
w_sol = 1.991e-7;     % Velocidad angular aparente del Sol [rad/s]
e = 0;                % Excentricidad de la orbita nominal

h_vals  = 200:2:1000;  % km
n_h     = length(h_vals);
i_SSO   = zeros(1,n_h);
T_orb_eff = zeros(1,n_h);
rgt_ratio = zeros(1,n_h);
rate_ground_vals = zeros(1,n_h);

for k = 1:n_h
    h  = h_vals(k);
    a  = (R_E + h*1e3);
    r  = a;
    ci = -(2 * r^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
    ci = max(-1, min(1, ci));
    i_SSO(k) = acosd(ci);

    [rate_orb, rate_ground] = rgt_rates(a, mu, R_E, J2, e, deg2rad(i_SSO(k)), omega_E);
    T_orb_eff(k) = 2*pi / rate_orb;
    rgt_ratio(k) = rate_orb / rate_ground;
    rate_ground_vals(k) = rate_ground;
end

fprintf('--- Reference at 630 km ---\n');
h_ref = 630;
a_ref = R_E + h_ref*1e3;
ci_ref = -(2 * a_ref^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
ci_ref = max(-1, min(1, ci_ref));
i_ref = acosd(ci_ref);
[rate_ref, rate_ground_ref, n_ref, dM_ref, domega_ref, dOmega_ref] = ...
    rgt_rates(a_ref, mu, R_E, J2, e, deg2rad(i_ref), omega_E);
Torb_ref = 2*pi / rate_ref;
Tground_ref = 2*pi / rate_ground_ref;
fprintf('i_SSO(630 km)      = %.2f deg\n', i_ref);
fprintf('n_k(630 km)         = %.10e rad/s\n', n_ref);
fprintf('Mdot_J2(630 km)     = %.10e rad/s\n', dM_ref);
fprintf('omegadot_J2(630 km) = %.10e rad/s\n', domega_ref);
fprintf('Omegadot_J2(630 km) = %.10e rad/s\n', dOmega_ref);
fprintf('T_orb_eff(630 km)   = %.2f s\n', Torb_ref);
fprintf('T_ground_eff        = %.2f s\n', Tground_ref);
fprintf('rate_orb / rate_ground = %.4f\n', rate_ref / rate_ground_ref);
fprintf('\n');

fprintf('--- Exact RGT altitudes (D <= 30) ---\n');
fprintf('  h [km]    i_SSO [deg]    N     D    T_orb_eff [s]    N/D\n');
fprintf('---------------------------------------------------------\n');

tol   = 1e-10;
maxit = 100;
roots = [];

for D = 1:30
    N_min = ceil(D * rgt_ratio(end));   % at h=1000 km
    N_max = floor(D * rgt_ratio(1));    % at h=200 km
    for N = N_min:N_max
        ha = 200e3;  % lower bound in m
        hb = 1000e3; % upper bound in m

        fa = f_RGT(ha, N, D, mu, R_E, J2, e, omega_E, w_sol);
        fb = f_RGT(hb, N, D, mu, R_E, J2, e, omega_E, w_sol);

        if fa * fb > 0
            continue;
        end

        for iter = 1:maxit
            hc = (ha + hb) / 2;
            fc = f_RGT(hc, N, D, mu, R_E, J2, e, omega_E, w_sol);
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
        a_exact = R_E + hc;
        ci = -(2 * a_exact^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
        ci = max(-1, min(1, ci));
        i_exact = acosd(ci);
        [rate_exact, rate_ground_exact] = rgt_rates(a_exact, mu, R_E, J2, e, deg2rad(i_exact), omega_E);
        T_orb_exact = 2*pi / rate_exact;
        ratio_exact = rate_exact / rate_ground_exact;

        roots = [roots; h_exact, i_exact, N, D, T_orb_exact, N/D];
        fprintf('%9.3f   %8.2f     %4d   %2d   %10.2f    %.4f\n', ...
            h_exact, i_exact, N, D, T_orb_exact, N/D);
    end
end

roots = sortrows(roots, 1);

h_nearest = 630;
[~, idx] = min(abs(roots(:,1) - h_nearest));
fprintf('\n--- Nearest to %.0f km ---\n', h_nearest);
fprintf('h = %.3f km, i = %.2f deg, N = %d, D = %d, T_orb_eff = %.2f s\n', ...
    roots(idx,1), roots(idx,2), roots(idx,3), roots(idx,4), roots(idx,5));
fprintf('Delta h = %.1f m\n', (roots(idx,1) - h_nearest)*1e3);

% Daily longitudinal drift
residual = zeros(1, n_h);
for k = 1:n_h
    T_ground_eff_k = 2*pi / rate_ground_vals(k);
    best_err = inf;
    for D = 1:30
        N_cand = round(rgt_ratio(k) * D);
        err = abs(N_cand * T_orb_eff(k) - D * T_ground_eff_k) / D;
        if err < best_err
            best_err = err;
        end
    end
    residual(k) = best_err * (2*pi*R_E / T_rot) / 1e3;  % km/day
end

if ~strcmp(getenv('TFG_SKIP_RGT_PLOT'), '1')
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
end

function out = f_RGT(h_m, N, D, mu, R_E, J2, e, omega_E, w_sol)
    a = R_E + h_m;
    ci = -(2 * a^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
    ci = max(-1, min(1, ci));
    i = acos(ci);
    [rate_orb, rate_ground] = rgt_rates(a, mu, R_E, J2, e, i, omega_E);
    out = rate_orb - rate_ground * N / D;
end

function [rate_orb, rate_ground, n_k, dM, domega, dOmega] = rgt_rates(a, mu, R_E, J2, e, i, omega_E)
    n_k = sqrt(mu / a^3);
    p = a * (1 - e^2);
    s2 = sin(i)^2;
    dOmega = -1.5 * n_k * J2 * (R_E / p)^2 * cos(i);
    domega = 0.75 * n_k * J2 * (R_E / p)^2 * (4 - 5*s2);
    dM = 0.75 * n_k * J2 * sqrt(1 - e^2) * R_E^2 / ...
        (a^2 * (1 - e^2)^2) * (2 - 3*s2);
    rate_orb = n_k + dM + domega;
    rate_ground = omega_E - dOmega;
end
