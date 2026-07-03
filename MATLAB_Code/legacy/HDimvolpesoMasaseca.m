%% Parámetros base
TM_aperture = 400;     % mm
TM_length = 1;         % valor normalizado
TM_weight = 240;       % kg
TM_power = 280;        % W

SEOSAT_aperture = 250; % mm
SEOSAT_length = 1;     % valor normalizado
SEOSAT_weight = 100;   % kg
SEOSAT_power = 100;    % W

%% Leer datos desde CSV
datos = readtable('hvsDminsat.csv');
alturas_km = datos.Altura_km;
diametros_mm = datos.D_min_mm;
n = length(diametros_mm);

%% Inicialización
Li_mean = zeros(n,1); Wi_mean = zeros(n,1); Vi_mean = zeros(n,1); Pi_mean = zeros(n,1); Si_mean = zeros(n,1);

for i = 1:n
    % Escalado TM
    [Li_TM, Si_TM, Vi_TM, Wi_TM, Pi_TM] = scaleInstrument(diametros_mm(i), TM_aperture, TM_length, TM_weight, TM_power);
    
    % Escalado SEOSAT
    [Li_SEOSAT, Si_SEOSAT, Vi_SEOSAT, Wi_SEOSAT, Pi_SEOSAT] = scaleInstrument(diametros_mm(i), SEOSAT_aperture, SEOSAT_length, SEOSAT_weight, SEOSAT_power);
    
    % Promedios
    Li_mean(i) = (Li_TM + Li_SEOSAT) / 2;
    Wi_mean(i) = (Wi_TM + Wi_SEOSAT) / 2;
    Vi_mean(i) = (Vi_TM + Vi_SEOSAT) / 2;
    Pi_mean(i) = (Pi_TM + Pi_SEOSAT) / 2;
    Si_mean(i) = (Si_TM + Si_SEOSAT) / 2;
end

%% Masa seca total (4 instrumentos)
masa_seca = 4 * Wi_mean;

%% Guardar CSV final
tabla_final = table(alturas_km, diametros_mm, Li_mean,Si_mean, Wi_mean, Vi_mean, Pi_mean, masa_seca, ...
    'VariableNames', {'Altura_km', 'Diametro_pupila', 'Longitud_media', 'Sup_media',...
                      'Peso_medio', 'Volumen_medio', 'Potencia_media', 'Masa_seca'});
writetable(tabla_final, 'instrumentos_escalados.csv');

%% Gráficas
figure;
subplot(3,2,1);
plot(alturas_km, Li_mean, 'k-o', 'LineWidth', 1.5);
title('Longitud Media'); xlabel('Altura orbital (km)'); ylabel('Longitud'); grid on;

subplot(3,2,2);
plot(alturas_km, Wi_mean, 'k-o', 'LineWidth', 1.5);
title('Peso Medio'); xlabel('Altura orbital (km)'); ylabel('Peso (kg)'); grid on;

subplot(3,2,3);
plot(alturas_km, Vi_mean, 'k-o', 'LineWidth', 1.5);
title('Volumen Medio'); xlabel('Altura orbital (km)'); ylabel('Volumen (m³)'); grid on;

subplot(3,2,4);
plot(alturas_km, Pi_mean, 'k-o', 'LineWidth', 1.5);
title('Potencia Media'); xlabel('Altura orbital (km)'); ylabel('Potencia (W)'); grid on;

subplot(3,2,[5 6]);
plot(alturas_km, masa_seca, 'm-s', 'LineWidth', 1.5);
title('Masa Seca Total (4×Peso Medio)'); xlabel('Altura orbital (km)'); ylabel('Masa Seca (kg)'); grid on;

set(gcf, 'Position', [100, 100, 1000, 900]);
sgtitle('Escalamiento Promedio de Instrumentos vs. Altura Orbital');

%% --------- Función de escalado -------------
function [Li, Si, Vi, Wi, Pi] = scaleInstrument(Ai, Ao, Lo, Wo, Po, K)
    if nargin < 6
        K = 1;
    end
    R = Ai / Ao;
    Li = R * Lo;
    Si = Li^2;
    Vi = Li^3;
    Wi = K * R^3 * Wo;
    Pi = K * R^3 * Po;
end