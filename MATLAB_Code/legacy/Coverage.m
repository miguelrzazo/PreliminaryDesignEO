% ---------------- CONFIGURACIÓN ----------------
clear; clc;

% Nuevos parámetros
Npix = 1000;           % Número de píxeles por detector
GSD = 30;              % Ground Sample Distance (m)
max_detectores = 3;    % Límite máximo de detectores

% Parámetros del área de interés (EE.UU. continental)
lat_min_US = 24.5; lat_max_US = 49.5;
lon_min_US = -125; lon_max_US = -66;
ancho_US_km = 4500; alto_US_km = 2500; % aprox.

% Parámetros orbitales
alturas_km = 200:50:2000;
swaths_km = 20:10:200;
inclinaciones_deg = 50:5:140;
fov_limit_deg = 10;      % Límite de FOV (grados)

% Números de satélites a evaluar
sat_counts = [1];

% Constantes
R_earth_km = 6378.137;   % Radio ecuatorial terrestre
overlap_factor = 0.95;   % 5% de solapamiento
effective_pass_fraction = 5/6;
daytime_only = false;     % Switch para solo pasadas diurnas

% Inicialización de resultados
numAlt = numel(alturas_km);
numSwath = numel(swaths_km);
numInc = numel(inclinaciones_deg);
numSat = numel(sat_counts);

% Matriz para almacenar resultados: [h, sw, inc, n_sat, days, fov, n_detectores]
results = [];

% ---------------- CÁLCULO PRINCIPAL ----------------
fprintf('Iniciando simulación...\n');
for iSat = 1:numSat
    n_sat = sat_counts(iSat);
    RAAN_shift = 360 / n_sat;
    for iAlt = 1:numAlt
        h = alturas_km(iAlt);
        for iSw = 1:numSwath
            swath_nominal = swaths_km(iSw);
            swath_eff = swath_nominal * overlap_factor;
            
            % Calcula número de detectores
            n_detectores = swath_nominal * 1000 / (Npix * GSD);
            
            % Salta si excede el límite de detectores
            if n_detectores > max_detectores
                continue;
            end
            
            % Calcula FOV correspondiente al swath y altura
            fov_calc = 2 * atand(swath_nominal/(2*h));
            if fov_calc > fov_limit_deg
                continue;  % Skip si excede límite de FOV
            end
            for iInc = 1:numInc
                inc = inclinaciones_deg(iInc);
                if inc < lat_max_US
                    continue;
                end
                % Días para cobertura
                dias = calcular_dias_cobertura(h, swath_eff, inc, alto_US_km, effective_pass_fraction, n_sat, daytime_only);
                % Solo guardar si <= 7 días
                if dias <= 7
                    results = [results; h, swath_nominal, inc, n_sat, dias, fov_calc, n_detectores];
                end
            end
        end
    end
end

% ---------------- EXPORTAR CSV FILTRADO ----------------
fprintf('Guardando CSV filtrado...\n');
T = array2table(results, 'VariableNames', ...
    {'Altitude_km','Swath_km','Inclination_deg','NumSat','CoverageDays','FOV_deg','NumDetectores'});
% Guardar un único CSV con todos los casos válidos
writetable(T, 'optimal_configs_filtered.csv');

% ---------------- HEATMAPS (multiplo de 5°) ----------------
fprintf('Generando heatmaps...\n');
for i = 1:numInc
    inc = inclinaciones_deg(i);
    if mod(inc,5)~=0, continue; end
    mask = (results(:,3)==inc);
    if ~any(mask), continue; end
    % Construir matriz de días
    M = nan(numAlt, numSwath);
    for j = find(mask)'
        [~, sw_idx] = ismember(results(j,2), swaths_km);
        [~, alt_idx] = ismember(results(j,1), alturas_km);
        M(alt_idx, sw_idx) = results(j,5);
    end
    figure;
    imagesc(swaths_km, alturas_km, M);
    xlabel('Swath (km)'); ylabel('Altura (km)');
    title(sprintf('CoverageDays <=7 (Inc %d°)',inc));
    colorbar; caxis([0 7]); set(gca,'YDir','normal');
end

% ---------------- GROUND TRACK + SWATH ----------------
fprintf('Generando ground track...\n');
% Usar configuración de visualización: primer elemento de results
default = results(1,:);
plot_groundtrack_con_swath(default(1), default(3), default(2), default(4), lat_min_US, lat_max_US, lon_min_US, lon_max_US);

% ---------------- ESCENA 3D ----------------
%fprintf('Generando visualización 3D...\n');
%visualizar_sat3D(default(1), default(3));

fprintf('Simulación finalizada.\n');


%% ---------------- FUNCIONES ----------------
function dias = calcular_dias_cobertura(h, swath, inc, alto, nublado_factor, n_sat, solo_dia)
    mu = 398600.4418; R = 6378.137;
    % Periodo orbital min
    T_orb_min = 2*pi*sqrt((R + h)^3 / mu) / 60;
    n_orbitas_dia = floor(1440 / T_orb_min);
    % Fracción de órbitas que cruzan banda latitudinal
    lat_min_US = 24.5; lat_max_US = 49.5;
    frac_inc = (sin(lat_max_US*pi/180)-sin(lat_min_US*pi/180))/sin(inc*pi/180);
    frac_inc = min(max(frac_inc,0),1);
    pasadas = n_orbitas_dia * n_sat * frac_inc;
    if solo_dia, pasadas = pasadas * 0.5; end
    pasadas = pasadas * nublado_factor;
    n_strips = ceil(alto/swath);
    dias = n_strips / pasadas;
end

function plot_groundtrack_con_swath(h, inc, swath, RAAN_shift, latmin, latmax, lonmin, lonmax)
    mu = 398600.4418; R = 6378.137; e=0;
    t = 0:60:86400*7;
    lat=zeros(size(t)); lon=zeros(size(t));
    for k=1:length(t)
        [r_eci,~]=keplerorbit(mu,R+h,e,deg2rad(inc),deg2rad(RAAN_shift),0,0,t(k));
        theta=deg2rad(360.985647*t(k)/86400);
        R_e=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
        r_rot=R_e*r_eci;
        lla=eci2lla(r_rot'*1000,datevec(datetime(2020,1,1)+seconds(t(k))));
        lat(k)=lla(1); lon(k)=lla(2);
    end
    delta=(swath/2)/111.32;
    u=lat+delta; l=lat-delta;
    figure; worldmap([latmin latmax],[lonmin lonmax]); load coastlines; plotm(coastlat,coastlon,'k'); hold on;
    plotm(lat,lon,'r-','LineWidth',1.5); plotm(u,lon,'b-'); plotm(l,lon,'b-');
    title(sprintf('GT & Swath Inc=%.1f° Alt=%.0fkm',inc,h)); legend('Coast','Track','Upper','Lower'); hold off;
end

function visualizar_sat3D(h, inc)
    % Visualización 3D usando satelliteScenario con órbita circular
    % Crear escenario con valores por defecto
    ss = satelliteScenario;

    % Crear satélite en órbita circular
    semiMajor = 6378.137 + h;
    ecc = 0;
    incl = deg2rad(inc);
    RAAN = 0;
    argPer = 0;
    trueAnom = 0;
    sat = satellite(ss, semiMajor, ecc, incl, RAAN, argPer, trueAnom, 'Name', 'DemoSat');

    % Crear estación de tierra en EE.UU.
    groundStation(ss, 38, -97, 'Name', 'USCenter');

    % Iniciar visor y reproducir
    viewer = satelliteScenarioViewer(ss);
    % Reproducir 1 día de simulación
    play(ss, 86400);

    % Dibujar Tierra
    [sx, sy, sz] = sphere(50);
    figure; hold on;
    surf(R*sx, R*sy, R*sz, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    colormap([0.6 0.8 1]);
    % Dibujar órbita
    plot3(X, Y, Z, 'r-', 'LineWidth', 2);
    axis equal;
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    title(sprintf('Órbita Circular (Inc=%.1f°, Alt=%.0f km)', inc, h));
    camlight; lighting gouraud;
    hold off;
end

function [r_eci,v_eci]=keplerorbit(mu,a,e,i,RAAN,argp,nu0,t)
    n=sqrt(mu/a^3); M=n*t; nu=mod(nu0+M,2*pi);
    r_pf=[a*cos(nu);a*sin(nu);zeros(1,numel(t))]; v_pf=[-a*n*sin(nu);a*n*cos(nu);zeros(1,numel(t))];
    R1=[cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0;0 0 1];
    R2=[1 0 0;0 cos(i) -sin(i);0 sin(i) cos(i)]; R3=[cos(argp) -sin(argp) 0;sin(argp) cos(argp) 0;0 0 1];
    Q=R1*R2*R3;
    r_eci=Q*r_pf; v_eci=Q*v_pf;
end
