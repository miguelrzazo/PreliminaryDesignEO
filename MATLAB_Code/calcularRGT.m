% calcularRGT.m - RGT exact altitudes for SSO orbits with J2
clear; close all;

mu    = 3.986004418e14;
R_E   = 6378.137e3;
J2    = 1.08263e-3;
T_sol = 86400;
w_sol = 1.991e-7;

h_vals  = 200:2:1000;  % km
n_h     = length(h_vals);
i_SSO   = zeros(1,n_h);
T_nodal = zeros(1,n_h);

for k = 1:n_h
    h  = h_vals(k);
    a  = (R_E + h*1e3);
    r  = a;
    ci = -(2 * r^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
    ci = max(-1, min(1, ci));
    i_SSO(k) = acosd(ci);

    T_k = 2*pi * sqrt(a^3 / mu);
    s2  = sind(i_SSO(k))^2;
    ratio = (R_E / a)^2;
    T_nodal(k) = T_k / (1 + 0.75 * J2 * ratio * (6 - 8*s2));
end

fprintf('--- Reference at 520 km ---\n');
h_ref = 520;
a_ref = R_E + h_ref*1e3;
ci_ref = -(2 * a_ref^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
ci_ref = max(-1, min(1, ci_ref));
i_ref = acosd(ci_ref);
Tk_ref = 2*pi * sqrt(a_ref^3 / mu);
s2_ref = sind(i_ref)^2;
ratio_ref = (R_E / a_ref)^2;
Tn_ref = Tk_ref / (1 + 0.75 * J2 * ratio_ref * (6 - 8*s2_ref));
fprintf('i_SSO(520 km)      = %.2f deg\n', i_ref);
fprintf('T_k(520 km)         = %.2f s\n', Tk_ref);
fprintf('T_nodal(520 km)     = %.2f s\n', Tn_ref);
fprintf('T_sol / T_nodal      = %.4f\n', T_sol / Tn_ref);
fprintf('\n');

fprintf('--- Exact RGT altitudes (D <= 30) ---\n');
fprintf('  h [km]    i_SSO [deg]    N     D    T_nodal [s]    N/D\n');
fprintf('---------------------------------------------------------\n');

tol   = 1e-10;
maxit = 100;
roots = [];

for D = 1:30
    N_min = ceil(D * T_sol / T_nodal(end));   % at h=1000 km (smallest T_nodal -> largest ratio)
    N_max = floor(D * T_sol / T_nodal(1));     % at h=200 km
    for N = N_min:N_max
        ha = 200e3;  % lower bound in m
        hb = 1000e3; % upper bound in m

        fa = f_RGT(ha, N, D, mu, R_E, J2, T_sol, w_sol);
        fb = f_RGT(hb, N, D, mu, R_E, J2, T_sol, w_sol);

        if fa * fb > 0
            continue;
        end

        for iter = 1:maxit
            hc = (ha + hb) / 2;
            fc = f_RGT(hc, N, D, mu, R_E, J2, T_sol, w_sol);
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
        Tk_exact = 2*pi * sqrt(a_exact^3 / mu);
        s2_exact = sind(i_exact)^2;
        ratio_exact = (R_E / a_exact)^2;
        Tn_exact = Tk_exact / (1 + 0.75 * J2 * ratio_exact * (6 - 8*s2_exact));

        roots = [roots; h_exact, i_exact, N, D, Tn_exact, N/D];
        fprintf('%9.3f   %8.2f     %4d   %2d   %10.2f    %.4f\n', ...
            h_exact, i_exact, N, D, Tn_exact, N/D);
    end
end

roots = sortrows(roots, 1);

h_nearest = 520;
[~, idx] = min(abs(roots(:,1) - h_nearest));
fprintf('\n--- Nearest to %.0f km ---\n', h_nearest);
fprintf('h = %.3f km, i = %.2f deg, N = %d, D = %d, T_n = %.2f s\n', ...
    roots(idx,1), roots(idx,2), roots(idx,3), roots(idx,4), roots(idx,5));
fprintf('Delta h = %.1f m\n', (roots(idx,1) - h_nearest)*1e3);

% Daily longitudinal drift
residual = zeros(1, n_h);
for k = 1:n_h
    ratio_k = T_sol / T_nodal(k);
    best_err = inf;
    for D = 1:30
        N_cand = round(ratio_k * D);
        err = abs(N_cand * T_nodal(k) - D * T_sol) / D;
        if err < best_err
            best_err = err;
        end
    end
    residual(k) = best_err * (2*pi*R_E / T_sol) / 1e3;  % km/day
end

figure('Position', [100 100 700 420]);
semilogy(h_vals, residual, 'b-', 'LineWidth', 1.2);
hold on;
for k = 1:size(roots,1)
    [~, ih] = min(abs(h_vals - roots(k,1)));
    plot(roots(k,1), residual(ih), 'g.', 'MarkerSize', 18);
end
xline(520, 'r--', 'LineWidth', 1.5);
xlabel('Altura [km]');
ylabel('Deriva longitudinal diaria [km/día]');
title('Condición RGT para órbita SSO con J_2');
grid on;
set(gca, 'FontSize', 11);
exportgraphics(gcf, 'RGT_residual.pdf', 'ContentType', 'vector');
fprintf('\nSaved RGT_residual.pdf\n');

function out = f_RGT(h_m, N, D, mu, R_E, J2, T_sol, w_sol)
    a = R_E + h_m;
    ci = -(2 * a^(7/2) * w_sol) / (3 * J2 * R_E^2 * sqrt(mu));
    ci = max(-1, min(1, ci));
    s2 = (sqrt(1 - ci^2))^2;
    T_k = 2*pi * sqrt(a^3 / mu);
    ratio = (R_E / a)^2;
    T_n = T_k / (1 + 0.75 * J2 * ratio * (6 - 8*s2));
    out = N * T_n - D * T_sol;
end
