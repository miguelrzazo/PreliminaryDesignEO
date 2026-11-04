function MTF_total = referenceMTFfunction(lambda, pixel_size, MTF_Detector, MTF_alineamiento, GSD, R, alturas_orbitales, diametros_pupila, filename_prefix, telescope_name, detector_idx,MTF_req)

%% Parámetros MTF adicionales
MTF_aberraciones = 0.95;
MTF_fabricacion = 0.98;
MTF_vibraciones = 0.99;
MTF_Termoelastico = 0.95;
MTF_Margen = 0.9;
MTF_resto = MTF_Margen * MTF_vibraciones * MTF_fabricacion * MTF_Termoelastico * MTF_aberraciones * MTF_Detector;

%% Inicializar matriz MTF
MTF_total = zeros(length(alturas_orbitales), length(diametros_pupila));

%% Cálculo MTF para cada combinación altura-diámetro
for i = 1:length(alturas_orbitales)
    for j = 1:length(diametros_pupila)
        altura_m = alturas_orbitales(i) * 1e3;
        diametro_m = diametros_pupila(j) * 1e-3;
        
        % Distancia focal
        distancia_focal = (altura_m * pixel_size) / GSD;
        
        % Frecuencias de corte
        f_co = diametro_m / lambda;
        f_Nyquist = 1 / (2 * pixel_size);
        f_x = f_Nyquist * distancia_focal;
        
        if f_x > f_co
            f_x = f_co;
        end
        
        % Cálculo para sistema con obscuración central (ecuaciones 6-10 a 6-15)
        X = f_x / f_co;
        Y = X / R;
        
        % Cálculo de alpha según la ecuación 6-12
        if (1 + R^2 - 4*X^2) / (2*R) >= -1 && (1 + R^2 - 4*X^2) / (2*R) <= 1
            alpha = acos((1 + R^2 - 4*X^2) / (2*R));
        else
            alpha = 0;
        end
        
        % Cálculo de A (ecuación 6-13)
        if 0 <= X && X <= 1
            A = (2/pi) * (acos(X) - X*sqrt(1 - X^2));
        else
            A = 0;
        end
        
        % Cálculo de B (ecuación 6-14)
        if 0 <= Y && Y <= 1
            B = (2*R^2/pi) * (acos(Y) - Y*sqrt(1 - Y^2));
        else
            B = 0;
        end
        
        % Cálculo de C (ecuaciones 6-15)
        if 0 < X && X <= (1-R)/2
            C = -2*R^2;
        elseif (1-R)/2 < X && X < (1+R)/2
            C = (2*R/pi)*sin(alpha) + ((1+R^2)/pi)*alpha - ((2*(1-R^2))/pi)*atan((1+R)/(1-R)*tan(alpha/2)) - 2*R^2;
        elseif X >= (1+R)/2
            C = 0;
        else
            C = 0;
        end
        
        % Cálculo del OTF difracción según ecuación 6-10
        MTF_difraccion = (A + B + C) / (1 - R^2);
        
        % Si el valor es negativo o NaN, ajustarlo a 0
        if isnan(MTF_difraccion) || MTF_difraccion < 0
            MTF_difraccion = 0;
        end
        
        % MTF total
        MTF_value = MTF_difraccion * MTF_alineamiento * MTF_resto;
        
        % Aplicar requirement: valores < 0.25 se convierten a NaN
        if MTF_value < MTF_req
            MTF_total(i, j) = NaN;
        else
            MTF_total(i, j) = MTF_value;
        end
    end
end


end
