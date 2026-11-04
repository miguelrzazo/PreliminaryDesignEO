function coverage_days = referenceCoverage(GSD, alturas_orbitales, swaths_km, Cov_Requirement, current_Npix, current_fov_limit, N_sat, N_telescopes, solapamiento, cobertura_nubes, max_detectores, revisitSolver)
%REFERENCECOVERAGE Frozen scalar equations from the pre-refactor implementation.
R_E = 6378.137e3; % Radio de la Tierra [m]
mu_E = 3.986004418e14; % Constante gravitacional terrestre [m3/s2]
J2 = 1.082629989052e-3; % Achatamiento J2
SSO_nodal_rate = 360/(365.2421897*24*3600) * (pi/180); % Tasa de precesión para SSO [rad/s]

lat = 25;
dayLimit = [1, Cov_Requirement];
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

        % --- Cálculo del Swath Efectivo (W_ef, con solapamiento) ---
        effective_swath = swath * (1 - solapamiento);

        % --- Parámetros Orbitales para SSO ---
        a = (R_E + height * 1000); % Semieje mayor [m]
        n = sqrt(mu_E / a^3); % Movimiento medio [rad/s]
        cos_i = -SSO_nodal_rate / (1.5 * n * J2 * (R_E/a)^2);
        if abs(cos_i) > 1
            coverage_days(h,s) = NaN; % SSO no posible a esta altura
            continue;
        end
        inc = acos(cos_i);
        coes = [a, 0, inc, 0, 0]; % [sma, ecc, inc, RAAN, AoP]

        % --- Llamada a RevisitCalc ---
        % psi: semiángulo geométrico desde el nadir hasta el borde del swath.
        % Se pasa a RevisitCalc como el parámetro de FoV equivalente.
        half_fov_angle = atand(effective_swath / (2 * height)); % psi [grados]
        try
            maxRevisit = revisitSolver(coes, lat, 1, 0, N_sat, 1, 0, 1, 0, dayLimit, half_fov_angle);
            % Penalización media del enunciado: cobertura_nubes = 1/6.
            final_revisit = maxRevisit / (1 - cobertura_nubes);
            if final_revisit > Cov_Requirement
                coverage_days(h,s) = NaN;
            else
                coverage_days(h,s) = final_revisit;
            end
        catch ME
            fprintf('    RevisitCalc error h=%g km swath=%g km: %s\n', height, swath, ME.message);
            coverage_days(h,s) = NaN;
        end
    end
end


end
