function MTFfunction(lambda, pixel_size, MTF_Detector, MTF_alineamiento, GSD, R, alturas_orbitales, diametros_pupila, filename_prefix, telescope_name, detector_idx)
MTF_aberraciones = 0.95;
MTF_fabricacion = 0.98;
MTF_vibraciones = 0.99;
MTF_Termoelastico = 0.95;
MTF_Margen = 0.9;
MTF_resto = MTF_Margen*MTF_vibraciones*MTF_fabricacion*MTF_Termoelastico*MTF_aberraciones*MTF_Detector;

%% Inicializar matriz MTF
MTF_total = zeros(length(alturas_orbitales), length(diametros_pupila));

% CÃ¡lculo MTF
for i = 1:length(alturas_orbitales)
    for j = 1:length(diametros_pupila)
        altura_m = alturas_orbitales(i) * 1e3;
        diametro_m = diametros_pupila(j) * 1e-3;
        
        distancia_focal = (altura_m * pixel_size) / GSD;
        
        f_co = diametro_m / lambda;
        f_Nyquist = 1 / (2 * pixel_size);
        f_x = f_Nyquist * distancia_focal;
        
        if f_x > f_co
            f_x = f_co;
        end
        
        % CÃ¡lculo para sistema con obscuraciÃ³n central
        X = f_x / f_co;
        Y = X / R;
        
        % CÃ¡lculo de alpha segÃºn la ecuaciÃ³n 6-12
        if (1 + R^2 - 4*X^2) / (2*R) >= -1 && (1 + R^2 - 4*X^2) / (2*R) <= 1
            alpha = acos((1 + R^2 - 4*X^2) / (2*R));
        else
            alpha = 0;
        end
        
        % CÃ¡lculo de A (ecuaciÃ³n 6-13)
        if 0 <= X && X <= 1
            A = (2/pi) * (acos(X) - X*sqrt(1 - X^2));
        else
            A = 0;
        end
        
        % CÃ¡lculo de B (ecuaciÃ³n 6-14)
        if 0 <= Y && Y <= 1
            B = (2*R^2/pi) * (acos(Y) - Y*sqrt(1 - Y^2));
        else
            B = 0;
        end
        
        % CÃ¡lculo de C (ecuaciones 6-15)
        if 0 < X && X <= (1-R)/2
            C = -2*R^2;
        elseif (1-R)/2 < X && X < (1+R)/2
            C = (2*R/pi)*sin(alpha) + ((1+R^2)/pi)*alpha - ((2*(1-R^2))/pi)*atan((1+R)/(1-R)*tan(alpha/2)) - 2*R^2;
        elseif X >= (1+R)/2
            C = 0;
        else
            C = 0;
        end
        
        % CÃ¡lculo del OTF difracciÃ³n segÃºn ecuaciÃ³n 6-10
        MTF_difraccion = (A + B + C) / (1 - R^2);
        
        % Si el valor es negativo o NaN, ajustarlo a 0
        if isnan(MTF_difraccion) || MTF_difraccion < 0
            MTF_difraccion = 0;
        end
        
        MTF_total(i, j) = MTF_difraccion * MTF_alineamiento * MTF_resto;
    end
end


end