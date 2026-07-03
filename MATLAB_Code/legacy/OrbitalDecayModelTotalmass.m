function OrbitalDecayModelTotalmass()
    %% Parámetros iniciales
    mu = 3.986004418e14;      % Parámetro gravitacional (m³/s²)
    R_earth = 6378e3;         % Radio terrestre (m)
    Cd = 2.5;                 % Coeficiente de arrastre

    Isp = 220;                % Impulso específico (s)
    g0 = 9.80665;             % Gravedad estándar (m/s²)
    mission_duration = 8 * 365.25 * 24 * 3600;  % Duración de misión (8 años en segundos)
    datos = readtable('valores_mediossat2.csv');
    h0_array = 350:10:1400;
    Lm = datos.Longitud_media;
    Am= 0.0267*1.2;
    masa_seca = datos.Masa_seca;
    %% Prealocar resultados
    num_impulses = zeros(size(h0_array));
    total_delta_v = zeros(size(h0_array));
    fuel_mass = zeros(size(h0_array));

    for i = 1:length(h0_array)
        h0 = h0_array(i);
        h_current = h0 * 1e3;  % Convertir a metros
        t = 0;
        total_dv = 0;
        impulses = 0;
        dry_mass=masa_seca(i);
        A= Am(i)*1.2;
        while t < mission_duration
            % Calcular el tiempo de decaimiento hasta el 98% de la altitud actual
            a_initial = R_earth + h_current;
            h_target = 0.98 * h_current / 1e3;  % 98% de la altitud en km
            a_target = R_earth + h_target * 1e3;

            options = odeset('Events', @(t,a) decayEvent(t,a,a_target), 'RelTol', 1e-6);
            [t_ode, a_ode] = ode45(@(t,a) decayODE(t,a,R_earth,mu,Cd,A,dry_mass), [0, mission_duration], a_initial, options);

         

            % Actualizar tiempo y altitud
            t = t + t_ode(end);
            if t > mission_duration
                break;
            end

            % Calcular delta-V para el impulso (Hohmann transfer)
            r1 = R_earth + 0.98 * h_current;
            r2 = R_earth + 1.00 * h_current;
            dv = hohmannDeltaV(r1, r2, mu);

            % Acumular delta-V y combustible
            total_dv = total_dv + dv;
            impulses = impulses + 1;

            % Actualizar altitud
            h_current = 1.0 * h_current;
        end

        % Calcular masa de combustible
        fuel = dry_mass * (exp(total_dv / (Isp * g0)) - 1);

        % Almacenar resultados
        num_impulses(i) = impulses;
        total_delta_v(i) = total_dv;
        fuel_mass(i) = fuel*10;
        total_mass(i) = fuel_mass(i)+ masa_seca(i);
    end

    %% Graficar resultados
    figure;
    subplot(3,1,1);
    plot(h0_array, num_impulses, 'b-o');
    xlabel('Altura Inicial (km)');
    ylabel('Número de Impulsos');
    grid on;

    subplot(3,1,2);
    plot(h0_array, total_delta_v, 'r-o');
    xlabel('Altura Inicial (km)');
    ylabel('\DeltaV Total (m/s)');
    grid on;

    subplot(3,1,3);
    plot(h0_array, fuel_mass, 'g-o');
    xlabel('Altura Inicial (km)');
    ylabel('Masa de Combustible (kg)');
    grid on;

    figure;
    hold on
    plot(h0_array, fuel_mass, 'b-o');
    plot(h0_array, masa_seca, 'g-o');
    plot(h0_array, total_mass, 'r-o');
    
    xlabel('Altura Inicial (km)');
    ylabel('Masa Total (kg)');
    grid on;
    hold off
end

function dadt = decayODE(~, a, R_earth, mu, Cd, A, mass)
    h = a - R_earth;
    rho = getAtmosphericDensity(h);
    dadt = -Cd * A * rho * sqrt(mu * a) / mass;
end

function [value, isterminal, direction] = decayEvent(~, a, a_target)
    value = a - a_target;
    isterminal = 1;  % Detener integración
    direction = -1;  % Detectar decrecimiento
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
    % Modelo exponencial por tramos extendido hasta 2000 km
    h_km = h / 1e3;  % Altura en kilómetros

    if h_km < 100
        H = 6.7e3;  % Escala de altura (m)
        rho0 = 1.225;
    elseif h_km < 150
        H = 9.5e3;
        rho0 = 4.79e-7;
    elseif h_km < 200
        H = 25.5e3;
        rho0 = 1.81e-9;
    elseif h_km < 250
        H = 37.5e3;
        rho0 = 2.53e-10;
    elseif h_km < 300
        H = 44.8e3;
        rho0 = 6.24e-11;
    elseif h_km < 350
        H = 50.3e3;
        rho0 = 7.40e-9;
    elseif h_km < 400
        H = 54.8e3;
        rho0 = 6.98e-12;
    elseif h_km < 450
        H = 58.2e3;
        rho0 = 2.72e-12;
    elseif h_km < 500
        H = 61.3e3;
        rho0 = 1.13e-12;
    elseif h_km < 600
        H = 70e3;
        rho0 = 5e-13;
    elseif h_km < 700
        H = 80e3;
        rho0 = 1e-13;
    elseif h_km < 800
        H = 90e3;
        rho0 = 2e-14;
    elseif h_km < 900
        H = 100e3;
        rho0 = 5e-15;
    elseif h_km < 1000
        H = 110e3;
        rho0 = 1e-15;
    elseif h_km < 1500
        H = 150e3;
        rho0 = 1e-16;
    elseif h_km < 2000
        H = 200e3;
        rho0 = 1e-17;
    else
        rho = 1e-18;  % Densidad insignificante para alturas >2000 km
        return;
    end

    % Ajuste fino usando la fórmula exponencial por tramos
    rho = rho0 * exp(-(h_km - floor(h_km/50)*50) / (H / 1e3));
end
