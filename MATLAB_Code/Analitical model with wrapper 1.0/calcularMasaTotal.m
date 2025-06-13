function [total_mass, fuel_mass, num_impulses, total_delta_v] = calcularMasaTotal(h0_array, masa_seca, Am)    %% ParÃ¡metros iniciales
    %% Paraimetros iniciales
    mu = 3.986004418e14; % ParÃ¡metro gravitacional (mÂ³/sÂ²)
    R_earth = 6378e3; % Radio terrestre (m)
    Cd = 2.5; % Coeficiente de arrastre
    Isp = 220; % Impulso especÃ­fico (s)
    g0 = 9.80665; % Gravedad estÃ¡ndar (m/sÂ²)
    mission_duration = 8 * 365.25 * 24 * 3600; % DuraciÃ³n de misiÃ³n (8 aÃ±os en segundos)
    
    %% Prealocar resultados
    num_impulses = zeros(size(h0_array));
    total_delta_v = zeros(size(h0_array));
    fuel_mass = zeros(size(h0_array));
    total_mass = zeros(size(h0_array));
    
    for i = 1:length(h0_array)
        h0 = h0_array(i);
        h_target = h0 * 1e3; % Altura objetivo original en metros
        h_current = h_target; % Iniciar en la altura objetivo
        
        t = 0;
        total_dv = 0;
        impulses = 0;
        dry_mass = masa_seca(i);
        A = Am(i) * 1.2;
        
        while t < mission_duration
            % Calcular el tiempo de decaimiento hasta el 98% de la altura OBJETIVO ORIGINAL
            a_initial = R_earth + h_current;
            h_threshold = 0.98 * h_target; % 98% de la altura objetivo ORIGINAL
            a_target = R_earth + h_threshold;
            
            % Configurar opciones con tolerancias mÃ¡s estrictas y MaxStep
            options = odeset('Events', @(t,a) decayEvent(t,a,a_target), ...
                           'RelTol', 1e-8, 'AbsTol', 1e-10, ...
                           'MaxStep', 86400); % 1 dÃ­a como paso mÃ¡ximo
            
            try
                [t_ode, a_ode] = ode45(@(t,a) decayODE(t,a,R_earth,mu,Cd,A,dry_mass), ...
                                    [0, mission_duration-t], a_initial, options);
                
                % Actualizar tiempo y altitud
                t = t + t_ode(end);
                h_current = a_ode(end) - R_earth; % Altura despuÃ©s del decaimiento
            catch ME
                % Si hay un error, usar integraciÃ³n paso a paso
                fprintf('Error en ode45 para altura %d km: %s\n', h0, ME.message);
                
                % Enfoque alternativo con pasos pequeÃ±os
                dt_manual = 3600; % Paso de 1 hora
                t_elapsed = 0;
                a_current = a_initial;
                
                while a_current > a_target && t + t_elapsed < mission_duration
                    % Calcular derivada
                    dadt = decayODE(0, a_current, R_earth, mu, Cd, A, dry_mass);
                    
                    % Actualizar con paso de Euler
                    a_current = a_current + dadt * dt_manual;
                    t_elapsed = t_elapsed + dt_manual;
                end
                
                % Actualizar tiempo y altitud 
                t = t + t_elapsed;
                h_current = a_current - R_earth;
            end
            
            if t >= mission_duration
                break;
            end
            
            % Calcular delta-V para el impulso (Hohmann transfer)
            r1 = R_earth + h_current; % Radio despuÃ©s del decaimiento
            r2 = R_earth + h_target;  % Radio de la altura objetivo ORIGINAL
            dv = hohmannDeltaV(r1, r2, mu);
            
            % Acumular delta-V y combustible
            total_dv = total_dv + dv;
            impulses = impulses + 1;
            
            % Restaurar altura a la altura objetivo original
            h_current = h_target;
        end
        
        % CÃ¡lculo normal
        num_impulses(i) = impulses;
        total_delta_v(i) = total_dv;
        fuel = dry_mass * (exp(total_dv / (Isp * g0)) - 1);
        fuel_mass(i) = fuel*1.1;
        total_mass(i) = fuel_mass(i) + masa_seca(i);
        
        % Aplicar filtro para valores extremos
        if total_mass(i) > 100000
            total_mass(i) = NaN;
            fuel_mass(i) = NaN;
        end
    end
end
    
    %% Funciones internas (igual que antes)
    function dadt = decayODE(~, a, R_earth, mu, Cd, A, mass)
        h = a - R_earth;
        rho = getAtmosphericDensity(h);
        dadt = -Cd * A * rho * sqrt(mu * a) / mass;
    end
    
    function [value, isterminal, direction] = decayEvent(~, a, a_target)
        value = a - a_target;
        isterminal = 1; % Detener integraciÃ³n
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
        % Modelo exponencial
        h_km = h / 1e3; % Altura en kilÃ³metros
        
        if h_km < 100
            H = 6.7e3; % Escala de altura (m)
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
            rho = 1e-18; % Densidad insignificante para alturas >2000 km
            return;
        end
        
        % Ajuste fino usando la fÃ³rmula exponencial por tramos
        rho = rho0 * exp(-(h_km - floor(h_km/50)*50) / (H / 1e3));
    end