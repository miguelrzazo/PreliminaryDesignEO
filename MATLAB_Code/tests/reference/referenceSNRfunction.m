function SNR_table = referenceSNRfunction(lambda_c, pixel_size, eta, tau, GSD, r_obs, altitudes, diameters, filename_prefix, telescope_name, detector_idx,SNR_req)

%% Constantes físicas
h = 6.626e-34; % Constante de Planck (J·s)
c = 3e8; % Velocidad de la luz (m/s)
bandwidth_m = 20e-3; % Ancho de banda espectral Δλ (20 nm)
TDI = 1; % Número de etapas TDI
rad_ref = 100; % Radiancia de referencia [W/(m²·sr·μm)]

%% Ruidos del sistema (en e⁻ RMS)
Ndark = 50; % Ruido térmico alto por mal enfriamiento
Nread = 100; % Sistema de lectura básico o rápido
Npreamp = 5; % Electrónica con ganancia moderada
Nvideo = 10; % Línea de vídeo sin filtrado
Njitter = 5; % Plataforma con jitter notable
Nemc = 5; % Alta interferencia electromagnética
Nquant = 2; % ADC de baja resolución
Nnonlin = 2; % Mala calibración o no linealidad no compensada

%% Inicialización
SNR_table = zeros(length(altitudes), length(diameters));

%% Cálculo SNR para cada combinación altura-diámetro
for i = 1:length(altitudes)
    h_orb = altitudes(i) * 1e3; % convertir a metros

    % Focal length basado en GSD, altura y pixel size
    focal_length = (pixel_size * h_orb) / GSD;

    % Velocidad sobre el suelo
    v_orb = sqrt(3.986e14 / (6371e3 + h_orb));

    % Tiempo de integración
    integration_time = GSD / v_orb;

    for j = 1:length(diameters)
        D = diameters(j) / 1000; % Diámetro de la pupila en metros

        % Irradiancia en el detector (W/m²)
        irradiance = (pi * tau * bandwidth_m * rad_ref) / (1 + 4 * (focal_length / sqrt(D^2*(1-r_obs^2)))^2);

        % Área del píxel
        pixel_area = pixel_size^2;

        % Cálculo de Ne (electrones generados)
        Ne = (irradiance * pixel_area * eta * lambda_c * TDI * integration_time) / (h * c);

        % Ruido total
        N_total = sqrt(Ndark^2 + Nread^2 + Npreamp^2 + Nvideo^2 + ...
            Njitter^2 + Nemc^2 + Nquant^2 + Nnonlin^2 + Ne);

        % SNR
        SNR_value = 0.9 * Ne / N_total;

        % Aplicar requirement: valores < 400 se convierten a NaN
        if SNR_value < SNR_req
            SNR_table(i, j) = NaN;
        else
            SNR_table(i, j) = SNR_value;
        end
    end
end


end
