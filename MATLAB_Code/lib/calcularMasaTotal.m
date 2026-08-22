function [total_mass, fuel_mass, num_impulses, total_delta_v] = calcularMasaTotal(h0_array, masa_seca, Am)
    %% Parametros iniciales
    mu = 3.986004418e14; % Parametro gravitacional
    R_earth = 6378e3; % Radio terrestre (m)
    Cd = 2.5; % Coeficiente de resistencia aerodinamica
    Isp = 220; % Impulso especifico, valor tipico para propulsion quimica (s)
    g0 = 9.80665; % Gravedad 
    mission_duration = 8 * 365.25 * 24 * 3600; % duracion mision, 8 años en seg
    
    %% Prealocar resultados
    num_impulses = zeros(size(h0_array));
    total_delta_v = zeros(size(h0_array));
    fuel_mass = zeros(size(h0_array));
    total_mass = zeros(size(h0_array));
    
    for i = 1:length(h0_array)
        h0 = h0_array(i);
        h_target = h0 * 1e3; % Altura orbital nominal en metros
        h_threshold = 0.98 * h_target; % Umbral inferior de mantenimiento
        h_recovery = 1.02 * h_target; % Altura de recuperacion tras el impulso
        dry_mass = getIndexedValue(masa_seca, i);
        drag_area_margin = 1.2;
        A = getIndexedValue(Am, i) * drag_area_margin;
        
        cycle_time = decayCycleTime(h_threshold, h_recovery, R_earth, mu, Cd, A, dry_mass);
        if isfinite(cycle_time) && cycle_time > 0
            impulses = floor((mission_duration - eps(mission_duration)) / cycle_time);
        else
            impulses = 0;
        end

        % Calcular delta-V para cada impulso (Hohmann transfer)
        r1 = R_earth + h_threshold; % Radio despues del decaimiento
        r2 = R_earth + h_recovery;  % Radio de la orbita de recuperacion
        dv = hohmannDeltaV(r1, r2, mu);
        total_dv = impulses * dv;
        
        % Calculo masa
        num_impulses(i) = impulses;
        total_delta_v(i) = total_dv;
        fuel = dry_mass * (exp(total_dv / (Isp * g0)) - 1);
        fuel_mass(i) = fuel*1.1;
        total_mass(i) = fuel_mass(i) + dry_mass;
        
        % Aplicar filtro para valores extremos
        if total_mass(i) > 100000
            total_mass(i) = NaN;
            fuel_mass(i) = NaN;
        end
    end
end
    
    %% Funciones internas
    function value = getIndexedValue(values, idx)
        if isscalar(values)
            value = values;
        else
            value = values(idx);
        end
    end

    function t_decay = decayCycleTime(h_lower, h_upper, R_earth, mu, Cd, A, mass)
        if h_upper <= h_lower || Cd <= 0 || A <= 0 || mass <= 0
            t_decay = Inf;
            return;
        end

        [h_base_km, ~, ~] = atmosphereTable();
        h_lower_km = h_lower / 1e3;
        h_upper_km = h_upper / 1e3;
        split_points_km = h_base_km(h_base_km > h_lower_km & h_base_km < h_upper_km);
        bounds_m = [h_lower, split_points_km * 1e3, h_upper];
        t_decay = 0;

        for k = 1:(numel(bounds_m) - 1)
            h_a = bounds_m(k);
            h_b = bounds_m(k + 1);
            if h_b <= h_a
                continue;
            end
            integrand = @(h) mass ./ (Cd .* A .* getAtmosphericDensity(h) .* sqrt(mu .* (R_earth + h)));
            t_decay = t_decay + integral(integrand, h_a, h_b, ...
                'RelTol', 1e-8, 'AbsTol', 1e-3, 'ArrayValued', true);
        end
    end

    function dadt = decayODE(~, a, R_earth, mu, Cd, A, mass)
        h = a - R_earth;
        rho = getAtmosphericDensity(h);
        dadt = -Cd * A * rho * sqrt(mu * a) / mass;
    end
    
    function [value, isterminal, direction] = decayEvent(~, a, a_target)
        value = a - a_target;
        isterminal = 1; % Detener integracion
        direction = -1; % Detectar decrecimiento
    end
    
    function dv = hohmannDeltaV(r1, r2, mu)
        % Calcular delta-V para transferencia Hohmann (m/s)
        a_transfer = (r1 + r2) / 2;
        v1 = sqrt(mu / r1);
        v_p = sqrt(mu * (2/r1 - 1/a_transfer));
        dv1 = v_p - v1;
        v2 = sqrt(mu / r2);
        v_a = sqrt(mu * (2/r2 - 1/a_transfer));
        dv2 = v2 - v_a;
        dv = dv1 + dv2;
    end
    
    function rho = getAtmosphericDensity(h)
        % Exponential atmosphere table from Vallado/Curtis, consistent with USSA76.
        h_km = max(h / 1e3, 0);
        [h_base_km, rho_base, H_km] = atmosphereTable();

        rho = zeros(size(h_km));
        for j = 1:numel(h_km)
            idx = find(h_km(j) >= h_base_km, 1, 'last');
            if isempty(idx)
                idx = 1;
            elseif idx == numel(h_base_km)
                idx = numel(h_base_km);
            end
            rho(j) = rho_base(idx) * exp(-(h_km(j) - h_base_km(idx)) / H_km(idx));
        end
    end

    function [h_base_km, rho_base, H_km] = atmosphereTable()
        h_base_km = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, ...
                     130, 140, 150, 180, 200, 250, 300, 350, 400, 450, ...
                     500, 600, 700, 800, 900, 1000];
        rho_base = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, ...
                    3.206e-4, 8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7, ...
                    9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9, ...
                    5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, ...
                    9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13, ...
                    1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15];
        H_km = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, ...
                5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149, ...
                22.523, 29.740, 37.105, 45.546, 53.628, 53.298, ...
                58.515, 60.828, 63.822, 71.835, 88.667, 124.64, ...
                181.05, 268.00];
    end
