# Cálculos tras la revisión de la tutora

MATLAB R2025b Update 7 con Aerospace Toolbox ejecutó `MATLAB_Code/review_calculations.m` sobre el commit `0fa3c4db02dc232c3c642403efa4a974eb9a81c1`. El [workflow](https://github.com/miguelrzazo/PreliminaryDesignEO/actions/runs/41822472658) terminó con éxito. El [artefacto completo](https://github.com/miguelrzazo/PreliminaryDesignEO/actions/runs/41822472658/artifacts/13317933494) contiene los CSV y el diario de ejecución. `summary.json` conserva las cifras que imprimió MATLAB.

| Parámetro | Resultado |
|---|---:|
| Satélites / telescopios por satélite | 2 / 2 |
| Altura nominal | 690 km |
| Inclinación SSO | 98,147265° |
| Apertura / focal | 38 mm / 155,25 mm |
| Barrido adquirido / disponible por FoV | 240 km / 243,331 km |
| Peor revisita geométrica muestreada | 5,534958 días, a 41° N |
| Peor revisita efectiva muestreada | 6,641950 días, a 41° N |
| Masa seca por satélite | 3,341705 kg |
| Combustible de mantenimiento por satélite | 0,050304 kg |
| Mantenimiento | 2 maniobras, 29,324714 m/s |
| Reserva de deorbitado por satélite | 0,146061 kg |
| Masa por satélite / constelación con reservas | 3,538070 kg / 7,076140 kg |
| Mayor tasa mínima de descarga con SGP4 | 5,014169 Mbps |
| Mayor ocupación de memoria con SGP4 | 1,389472 GB |

La selección compara la masa seca y el combustible de mantenimiento de la malla original; añade la reserva de deorbitado después. La comprobación de cobertura usa 13 paralelos de 25 a 49° N, con paso de 2°, y la corrección media por nubosidad del enunciado. El factor por nubosidad es `T_efectiva = T_geometrica / (1 - 1/6) = 1,2 * T_geometrica`; las cifras de `revisit_days` y `Dias_Cobertura` incluyen ese factor. Esta comprobación no equivale a una validación continua de todo CONUS. La comparación frozen/circular usa el modelo de excentricidad media J2/J3 durante ocho años y arcos completos de latitud, sin recorte longitudinal. El cálculo de deorbitado conserva la aproximación circular equivalente del modelo original: propaga el semieje mayor de la transferencia, sin integrar el arrastre de una órbita excéntrica. Los contactos y la memoria corresponden a la semana que comienza el 28 de abril de 2024, con paso de 10 s.

El TeX y las figuras publicadas conservan la versión anterior. Las sustituciones de texto y cifras quedan pendientes de revisión con el autor. Para el enlace S, la guía NASA 453-NENUG, revisión 4, tabla 7-2, documenta 10 Mbps con FM/PCM para ASF AS1. Esa referencia identifica ASF, una instalación distinta de NOAA FCDAS.
