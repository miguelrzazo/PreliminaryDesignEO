function budget = calcularLinkBudget(Pt_dBW, Gt_dBi, Gr_dBi, freq_GHz, range_km, T_sys_K, rate_bps, L_misc_dB)
% CALCULARLINKBUDGET  Balance de enlace RF (ecuación de Friis).
%
%   budget = calcularLinkBudget(Pt_dBW, Gt_dBi, Gr_dBi, freq_GHz, ...
%                                range_km, T_sys_K, rate_bps, L_misc_dB)
%
%   Entradas:
%     Pt_dBW    - Potencia de transmisión [dBW]
%     Gt_dBi    - Ganancia antena transmisora (satélite) [dBi]
%     Gr_dBi    - Ganancia antena receptora (estación) [dBi]
%     freq_GHz  - Frecuencia portadora [GHz]
%     range_km  - Distancia satélite-estación [km]
%     T_sys_K   - Temperatura de ruido del sistema receptor [K]
%     rate_bps  - Velocidad de datos [bps]
%     L_misc_dB - Pérdidas de sistema adicionales (atmosf., cables, etc.) [dB, positivo]
%
%   Salida:
%     budget - struct con todos los términos del balance
%
%   Eb/N0 requerida: 10.5 dB (QPSK, BER = 1e-6, sin FEC).
%   Modificar EbN0_req_dB para otra modulación o código.

%% === Valores por defecto (enlace de descarga banda X, Fairbanks) ===
if nargin < 1,  Pt_dBW   = 5;        end   % 3 W ≈ 5 dBW
if nargin < 2,  Gt_dBi   = 8;        end   % parche banda X, satélite pequeño
if nargin < 3,  Gr_dBi   = 58;       end   % antena 13 m Fairbanks, banda X
if nargin < 4,  freq_GHz = 8.1;      end
if nargin < 5,  range_km = 700;      end   % elevación ~60°, 500 km de altitud
if nargin < 6,  T_sys_K  = 50;       end
if nargin < 7,  rate_bps = 150e6;    end   % 150 Mbps
if nargin < 8,  L_misc_dB = 0.5;    end   % atenuación atmosférica nominal

EbN0_req_dB = 10.5;  % QPSK, BER = 1e-6

%% === Cálculo ===
c   = 3e8;                              % [m/s]
k_dB = -228.6;                          % [dBW/Hz·K]  constante de Boltzmann

lambda_m  = c / (freq_GHz * 1e9);
FSPL_dB   = 20 * log10(4 * pi * range_km * 1e3 / lambda_m);
EIRP_dBW  = Pt_dBW + Gt_dBi;
GT_dBK    = Gr_dBi - 10 * log10(T_sys_K);

CN0_dBHz  = EIRP_dBW - FSPL_dB - L_misc_dB + GT_dBK - k_dB;
Rb_dBHz   = 10 * log10(rate_bps);
EbN0_dB   = CN0_dBHz - Rb_dBHz;
margin_dB = EbN0_dB - EbN0_req_dB;

%% === Almacenamiento en struct ===
budget.Pt_dBW        = Pt_dBW;
budget.Gt_dBi        = Gt_dBi;
budget.EIRP_dBW      = EIRP_dBW;
budget.FSPL_dB       = FSPL_dB;
budget.L_misc_dB     = L_misc_dB;
budget.Gr_dBi        = Gr_dBi;
budget.T_sys_K       = T_sys_K;
budget.GT_dBK        = GT_dBK;
budget.CN0_dBHz      = CN0_dBHz;
budget.rate_bps      = rate_bps;
budget.Rb_dBHz       = Rb_dBHz;
budget.EbN0_dB       = EbN0_dB;
budget.EbN0_req_dB   = EbN0_req_dB;
budget.link_margin_dB = margin_dB;

%% === Tabla en consola ===
fprintf('\n========== BALANCE DE ENLACE (Link Budget) ==========\n');
fprintf('  Frecuencia:   %.2f GHz    Rango: %.0f km\n', freq_GHz, range_km);
fprintf('------------------------------------------------------\n');
fprintf('  %-40s %+8.1f dBW\n',  'Potencia de transmisión (Pt)',        Pt_dBW);
fprintf('  %-40s %+8.1f dBi\n',  'Ganancia antena satélite (Gt)',       Gt_dBi);
fprintf('  %-40s %+8.1f dBW\n',  '>>> EIRP',                            EIRP_dBW);
fprintf('------------------------------------------------------\n');
fprintf('  %-40s %+8.1f dB\n',   'Pérdida de espacio libre (FSPL)',     -FSPL_dB);
fprintf('  %-40s %+8.1f dB\n',   'Pérdidas de sistema (atm., etc.)',    -L_misc_dB);
fprintf('------------------------------------------------------\n');
fprintf('  %-40s %+8.1f dBi\n',  'Ganancia antena receptora (Gr)',      Gr_dBi);
fprintf('  %-40s %+8.0f K\n',    'Temperatura de ruido (Tsys)',         T_sys_K);
fprintf('  %-40s %+8.1f dB/K\n', 'Figura de mérito (G/T)',              GT_dBK);
fprintf('------------------------------------------------------\n');
fprintf('  %-40s %+8.1f dBW/Hz\n','Cte. Boltzmann (k)',                k_dB);
fprintf('  %-40s %+8.1f dB·Hz\n', '>>> C/N0',                          CN0_dBHz);
fprintf('------------------------------------------------------\n');
fprintf('  %-40s %+8.1f dB·Hz\n', sprintf('Velocidad de datos (%.0f Mbps)', rate_bps/1e6), Rb_dBHz);
fprintf('  %-40s %+8.1f dB\n',   'Eb/N0 obtenida',                     EbN0_dB);
fprintf('  %-40s %+8.1f dB\n',   'Eb/N0 requerida (QPSK, BER 1e-6)',   EbN0_req_dB);
fprintf('======================================================\n');
if margin_dB >= 3
    fprintf('  MARGEN DE ENLACE:  %+.1f dB  ✓ (≥ 3 dB)\n', margin_dB);
else
    fprintf('  MARGEN DE ENLACE:  %+.1f dB  ✗ (< 3 dB — revisar parámetros)\n', margin_dB);
end
fprintf('======================================================\n\n');
