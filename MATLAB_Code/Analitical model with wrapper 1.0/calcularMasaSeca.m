function [masa_seca, Li_mean, Si_mean, Wi_mean, Vi_mean, Pi_mean] = calcularMasaSeca(diametros_mm, N_telescopes)
% calcularMasaSeca: Calcula las propiedades de la masa seca de forma vectorizada.
%
%   Esta versión optimizada elimina los bucles para un cálculo más eficiente.
%   Aplica un factor de escala al peso medio (Wi_mean) si hay 2 telescopios.
%
%   Inputs:
%     diametros_mm - Vector de diámetros de pupila [mm].
%     N_telescopes - Número de telescopios en la configuración.
%
%   Outputs:
%     masa_seca  - Vector de masa seca total [kg].
%     Li_mean    - Vector de longitudes medias [m normalizado].
%     Si_mean    - Vector de superficies medias [m^2 normalizado].
%     Wi_mean    - Vector de pesos medios del instrumento [kg].
%     Vi_mean    - Vector de volúmenes medios [m^3 normalizado].
%     Pi_mean    - Vector de potencias medias [W].

%% Parámetros base
TM_aperture = 400; % mm
TM_length = 1;     % valor normalizado
TM_weight = 240;   % kg
TM_power = 280;    % W

SEOSAT_aperture = 250; % mm
SEOSAT_length = 1;     % valor normalizado
SEOSAT_weight = 100;   % kg
SEOSAT_power = 100;    % W

%% Escalado vectorizado
% El bucle 'for' se elimina y los cálculos se aplican a todo el vector de diámetros.
[Li_TM, Si_TM, Vi_TM, Wi_TM, Pi_TM] = scaleInstrument(diametros_mm, TM_aperture, TM_length, TM_weight, TM_power);
[Li_SEOSAT, Si_SEOSAT, Vi_SEOSAT, Wi_SEOSAT, Pi_SEOSAT] = scaleInstrument(diametros_mm, SEOSAT_aperture, SEOSAT_length, SEOSAT_weight, SEOSAT_power);

%% Promedios vectorizados
Li_mean = (Li_TM + Li_SEOSAT) / 2;
Wi_mean = (Wi_TM + Wi_SEOSAT) / 2;
Vi_mean = (Vi_TM + Vi_SEOSAT) / 2;
Pi_mean = (Pi_TM + Pi_SEOSAT) / 2;
Si_mean = (Si_TM + Si_SEOSAT) / 2;

%% Aplicar factor si hay 2 telescopios
% Se multiplica el peso medio del instrumento (Wi_mean), que luego se propaga a la masa seca.
if N_telescopes == 2
    Wi_mean = Wi_mean * 1.5;
end

%% Masa seca total
% La masa seca del satélite se estima como 4 veces el peso medio del instrumento.
masa_seca = 4 * Wi_mean;

end

function [Li, Si, Vi, Wi, Pi] = scaleInstrument(Ai, Ao, Lo, Wo, Po)
    % Las operaciones (/, *, .^2, .^3) son vectoriales por defecto en MATLAB.
    % Esta función interna funciona con escalares o vectores sin cambios.
    R = Ai / Ao;
    Li = R .* Lo;
    Si = Li.^2;
    Vi = Li.^3;
    Wi = (R.^3) .* Wo;
    Pi = (R.^3) .* Po;
end
