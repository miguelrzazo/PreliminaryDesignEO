function SNRfunction(lambda_c, pixel_size, eta, tau, GSD, r_obs, altitudes, diameters, filename_prefix, telescope_name, detector_idx)

% === Constantes fÃ­sicas ===
h = 6.626e-34; % Constante de Planck (JÂ·s)
c = 3e8; % Velocidad de la luz (m/s)
bandwidth_m = 20e-3; % Ancho de banda espectral Î”Î» (20 nm)
TDI = 1; % NÃºmero de etapas TDI
rad_ref = 100; % Radiancia de referencia [W/(mÂ²Â·srÂ·Âµm)]

% === Ruidos del sistema (en eâ» RMS) ===
Ndark = 50; % Ruido tÃ©rmico alto por mal enfriamiento
Nread = 100; % Sistema de lectura bÃ¡sico o rÃ¡pido
Npreamp = 5; % ElectrÃ³nica con ganancia moderada
Nvideo = 10; % LÃ­nea de vÃ­deo sin filtrado
Njitter = 5; % Plataforma con jitter notable
Nemc = 5; % Alta interferencia electromagnÃ©tica
Nquant = 2; % ADC de baja resoluciÃ³n
Nnonlin = 2; % Mala calibraciÃ³n o no linealidad no compensada

%% === InicializaciÃ³n ===
SNR_table = zeros(length(altitudes), length(diameters));

for i = 1:length(altitudes)
    h_orb = altitudes(i) * 1e3; % convertir a metros
    
    % Focal length basado en GSD, altura y pixel size
    focal_length = (pixel_size * h_orb) / GSD;
    
    % Velocidad sobre el suelo
    v_orb = sqrt(3.986e14 / (6371e3 + h_orb));
    
    % Tiempo de integraciÃ³n
    integration_time = GSD / v_orb;
    
    for j = 1:length(diameters)
        D = diameters(j) / 1000; % DiÃ¡metro de la pupila en metros
        
        % Irradiancia en el detector (W/mÂ²)
        irradiance = (pi * tau * bandwidth_m * rad_ref) / (1 + 4 * (focal_length /sqrt(D^2*(1-r_obs^2)))^2) ;
        
        % Ãrea del pÃ­xel
        pixel_area = pixel_size^2;
        
        % CÃ¡lculo de Ne (electrones generados)
        Ne = (irradiance * pixel_area * eta * lambda_c * TDI * integration_time) / (h * c);
        
        % Ruido total
        N_total = sqrt(Ndark^2 + Nread^2 + Npreamp^2 + Nvideo^2 + ...
            Njitter^2 + Nemc^2 + Nquant^2 + Nnonlin^2 + Ne);
        
        % SNR
        SNR = Ne / N_total;
        
        SNR_table(i,j) = SNR;
    end
end


end