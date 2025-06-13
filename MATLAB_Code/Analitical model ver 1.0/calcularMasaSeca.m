

function [masa_seca, tabla_final] = calcularMasaSeca(alturas_km, diametros_mm, N_telescopes)
    %% ParÃ¡metros base
    TM_aperture = 400; % mm
    TM_length = 1; % valor normalizado
    TM_weight = 240; % kg
    TM_power = 280; % W
    
    SEOSAT_aperture = 250; % mm
    SEOSAT_length = 1; % valor normalizado
    SEOSAT_weight = 100; % kg
    SEOSAT_power = 100; % W
    
    %% InicializaciÃ³n
    n = length(diametros_mm);
    Li_mean = zeros(n,1); Wi_mean = zeros(n,1); 
    Vi_mean = zeros(n,1); Pi_mean = zeros(n,1); Si_mean = zeros(n,1);
    
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
    
    %% Masa seca total
    masa_seca = 4 * Wi_mean;
    
    %% Aplicar factor adicional si hay 2 telescopios
    if nargin >= 3 && N_telescopes == 2
        masa_seca = masa_seca * 1.5;
    end
    
    %% Crear tabla final
    tabla_final = table(alturas_km, diametros_mm, Li_mean, Si_mean, Wi_mean, Vi_mean, Pi_mean, masa_seca, ...
        'VariableNames', {'Altura_km', 'Diametro_pupila', 'Longitud_media', 'Sup_media',...
        'Peso_medio', 'Volumen_medio', 'Potencia_media', 'Masa_seca'});
    
    %% FunciÃ³n interna de escalado
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
end