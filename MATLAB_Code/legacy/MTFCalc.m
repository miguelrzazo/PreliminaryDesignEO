function MTF_total = MTFCalc(lambda, pixel_size, MTF_Detector, MTF_alineamiento, GSD, R, alturas_orbitales, diametros_pupila, detector_num, telescope_name, MTF_Requirement, lambda_name)

MTF_aberraciones = 0.95;
MTF_fabricacion = 0.98;
MTF_vibraciones = 0.99;
MTF_Termoelastico = 0.95;
MTF_Margen = 0.9;

MTF_resto = MTF_Margen*MTF_vibraciones*MTF_fabricacion*MTF_Termoelastico*MTF_aberraciones*MTF_Detector;

%% Inicializar matriz MTF
MTF_total = zeros(length(alturas_orbitales), length(diametros_pupila));

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
        
        X = f_x / f_co;
        
        % CORREGIDO: Manejar caso cuando R = 0 (telescopios no obscurados)
        if R == 0
            % Para telescopios no obscurados (Refractivo y TMA)
            if 0 <= X && X <= 1
                MTF_difraccion = (2/pi) * (acos(X) - X*sqrt(1 - X^2));
            else
                MTF_difraccion = 0;
            end
        else
            % Para telescopios obscurados (Korsch y Cassegrain)
            Y = X / R;
            
            if (1 + R^2 - 4*X^2) / (2*R) >= -1 && (1 + R^2 - 4*X^2) / (2*R) <= 1
                alpha = acos((1 + R^2 - 4*X^2) / (2*R));
            else
                alpha = 0;
            end
            
            if 0 <= X && X <= 1
                A = (2/pi) * (acos(X) - X*sqrt(1 - X^2));
            else
                A = 0;
            end
            
            if 0 <= Y && Y <= 1
                B = (2*R^2/pi) * (acos(Y) - Y*sqrt(1 - Y^2));
            else
                B = 0;
            end
            
            if 0 < X && X <= (1-R)/2
                C = -2*R^2;
            elseif (1-R)/2 < X && X < (1+R)/2
                C = (2*R/pi)*sin(alpha) + ((1+R^2)/pi)*alpha - ((2*(1-R^2))/pi)*atan((1+R)/(1-R)*tan(alpha/2)) - 2*R^2;
            elseif X >= (1+R)/2
                C = 0;
            else
                C = 0;
            end
            
            MTF_difraccion = (A + B + C) / (1 - R^2);
        end
        
        if isnan(MTF_difraccion) || MTF_difraccion < 0
            MTF_difraccion = 0;
        end
        
        MTF_total(i, j) = MTF_difraccion * MTF_alineamiento * MTF_resto;
    end
end

%% Crear matriz para heatmap (NaN para valores que no cumplen requisito)
MTF_heatmap = MTF_total;
MTF_heatmap(MTF_total < MTF_Requirement) = NaN;

%% Generar heatmap
fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600]);
h = heatmap(MTF_heatmap);

% Configurar etiquetas (mostrar cada 50 valores)
h.XDisplayLabels = repmat({''}, 1, length(diametros_pupila));
h.YDisplayLabels = repmat({''}, 1, length(alturas_orbitales));

for i = 1:50:length(diametros_pupila)
    h.XDisplayLabels{i} = sprintf('%d', diametros_pupila(i));
end

for i = 1:100:length(alturas_orbitales)
    h.YDisplayLabels{i} = sprintf('%d', alturas_orbitales(i));
end

h.Title = sprintf('MTF - Det%d, %s, GSD=%dm, R=%.1f, %s', ...
    detector_num, telescope_name, GSD, R, lambda_name);
h.XLabel = 'Diámetro Pupila (mm)';
h.YLabel = 'Altura Orbital (km)';

h.Colormap = parula(256);
h.ColorLimits = [0, 1.0];
h.MissingDataColor = [1, 1, 1]; % Blanco para NaN
h.GridVisible = 'off';

%% Guardar archivos
filename_base = sprintf('MTF_Det%d_%s_GSD%d_R%.1f_%s', ...
    detector_num, strrep(telescope_name, ' ', ''), GSD, R*10, lambda_name);

% Guardar heatmap
saveas(fig, fullfile('MTF', [filename_base '.
