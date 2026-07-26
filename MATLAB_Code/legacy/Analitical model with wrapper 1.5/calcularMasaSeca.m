function [masa_seca, masa_seca_TM, masa_seca_SEOSAT, Li_instrumento, Si_instrumento, Vi_instrumento, Wi_instrumento, Pi_instrumento, V_sat, L_sat, S_sat, U_sat] = calcularMasaSeca(diametros_mm, N_telescopes)
% calcularMasaSeca: Calcula las propiedades de la masa seca de forma vectorizada.

% Entradas:
%   diametros_mm - Vector de diámetros de pupila [mm].
%   N_telescopes - Número de telescopios en la configuración (1 o 2).

% Salidas (Instrumento):
%   Li_instrumento, Si_instrumento, Vi_instrumento, Wi_instrumento, Pi_instrumento
% Salidas (Satélite):
%   masa_seca - Masa seca total del satélite [kg].
%   V_sat, L_sat, S_sat - Volumen, longitud y área del satélite.
%   U_sat - Estándar CubeSat 'U'.

%% Parámetros base de los instrumentos de referencia
TM_aperture = 400; % Diámetro Thematic Mapper [mm]
TM_length = 1; % Longitud normalizada
TM_weight = 240; % Masa Thematic Mapper [kg]
TM_power = 280; % Potencia Thematic Mapper [W]

SEOSAT_aperture = 250; % Diámetro SEOSAT [mm]
SEOSAT_length = 1; % Longitud normalizada
SEOSAT_weight = 100; % Masa SEOSAT [kg]
SEOSAT_power = 100; % Potencia SEOSAT [W]

%% Escalado vectorizado de propiedades del instrumento
[Li_TM, Si_TM, Vi_TM, Wi_TM, Pi_TM] = scaleInstrument(diametros_mm, TM_aperture, TM_length, TM_weight, TM_power);
[Li_SEOSAT, Si_SEOSAT, Vi_SEOSAT, Wi_SEOSAT, Pi_SEOSAT] = scaleInstrument(diametros_mm, SEOSAT_aperture, SEOSAT_length, SEOSAT_weight, SEOSAT_power);

%% Promedio de las propiedades de los instrumentos escalados
Li_instrumento = (Li_TM + Li_SEOSAT) / 2;
Si_instrumento = (Si_TM + Si_SEOSAT) / 2;
Vi_instrumento = (Vi_TM + Vi_SEOSAT) / 2;
Wi_instrumento = (Wi_TM + Wi_SEOSAT) / 2;
Pi_instrumento = (Pi_TM + Pi_SEOSAT) / 2;

%% Ajuste por configuración de múltiples telescopios
if N_telescopes == 2
    Wi_instrumento = Wi_instrumento * 1.5;
end

%% Cálculo de masa seca del satélite
% Se estima que la masa de la carga útil (instrumento) es un 25% de la masa seca total.
masa_seca_TM = 4 * Wi_TM;
masa_seca_SEOSAT = 4 * Wi_SEOSAT;
masa_seca = 4 * Wi_instrumento;

%% ---- Propiedades del Satélite ----
% Se asume una densidad para el satélite y se calculan sus dimensiones como un CubeSat.
densidad_sat = 79; % Densidad del satélite [kg/m^3]

V_sat = masa_seca ./ densidad_sat; % Volumen del satélite [m^3]
L_sat = V_sat.^(1/3); % Longitud del satélite (lado del cubo) [m]
S_sat = L_sat.^2; % Área de la sección transversal del satélite [m^2]
U_sat = V_sat / (0.1^3); % Estándar 'U' del satélite (1U = 10x10x10 cm)

end

function [Li, Si, Vi, Wi, Pi] = scaleInstrument(Ai, Ao, Lo, Wo, Po)
    R = Ai / Ao;       % Ratio de aperturas

    % Definir K según la condición de R
    K = ones(size(R)); % Inicializa K en 1
    K(R <= 0.5) = 2;   % Asigna K=2 donde R<=0.5

    Li = R .* Lo; 
    Si = Li.^2;
    Vi = Li.^3;
    Wi = (R.^3) .* Wo .* K; % Multiplica por K la masa
    Pi = (R.^3) .* Po .* K; % Multiplica por K la potencia
end
