# Guión de Defensa — TFG: Diseño Preliminar de un Satélite de Observación de la Tierra

**Autor:** Miguel Rosa Zazo  
**Fecha:** Febrero 2026  
**Duración estimada:** 20 minutos + preguntas  
**Tiempo total de trabajo:** ~287 horas  

---

## Herramientas de IA y distribución del tiempo

### Herramientas de IA utilizadas

| Herramienta | Uso principal |
|---|---|
| Claude Sonnet 3.7 / 4.0 | Redacción de memoria, corrección de LaTeX, análisis |
| Claude Thinking (Extended) 3.7 / 4.0 | Resolución de problemas complejos, diseño de subsistemas |
| Claude Opus (via Claude Code) | Programación MATLAB, Python, LaTeX |
| Perplexity Deep Research | Revisión bibliográfica, búsqueda de misiones comparables |
| Google Gemini Pro 2 | Síntesis de documentación técnica de componentes COTS |

### Distribución de las ~287 horas

| Fase | Descripción | Horas |
|---|---|---|
| 1 | Revisión bibliográfica: misiones comparables (OCO-2, GOSAT, GHGSat, CO2M) | 28 |
| 2 | Desarrollo del modelo MTF (MATLAB) | 18 |
| 3 | Desarrollo del modelo SNR (MATLAB) | 18 |
| 4 | Desarrollo del modelo de cobertura: wrapper RevisitTime | 22 |
| 5 | Desarrollo del modelo de masa seca y total | 20 |
| 6 | Integración del pipeline analítico (CrossData, ConfigAnalysis) | 15 |
| 7 | Modelo de propagación numérica (Propagation ver 2.0) | 20 |
| 8 | Análisis y validación de resultados (tablas, heatmaps) | 25 |
| 9 | Redacción de la memoria en LaTeX | 55 |
| 10 | Diseño de subsistemas (ADCS, propulsión, comunicaciones) | 22 |
| 11 | Análisis de lanzadores y segmento tierra | 15 |
| 12 | Puerto Python (`eo-mission` package) | 18 |
| 13 | Preparación de la presentación Beamer | 22 |
| 14 | Revisión final, correcciones y GitHub | 9 |
| **Total** | | **287** |

---

## Guión por diapositivas

> **Notación de tiempo:** cada bloque indica el tiempo acumulado aproximado desde el inicio.
> Total presentación: ~20 min.

---

### [0:00] Diapositiva 1 — Portada

Buenos días a todos los miembros del tribunal.

Me llamo Miguel Rosa Zazo y presento mi Trabajo de Fin de Grado: *Diseño Preliminar de un Satélite de Observación de la Tierra*.

El trabajo aborda el reto de monitorizar gases de efecto invernadero desde el espacio usando nanosatélites de muy baja masa. La presentación dura aproximadamente 20 minutos seguidos de un turno de preguntas.

---

### [0:30] Diapositiva 2 — Introducción y Objetivos

El objetivo principal es demostrar que es posible monitorizar CO₂ y O₂ atmosférico sobre EE.UU. con un GSD de 80 m —20 veces mejor que OCO-2 de la NASA— usando satélites con una masa total de tan solo **1,71 kg**.

Importante: el proyecto no fabrica el satélite. Establece un diseño preliminar analítico verificable, con código abierto, que puede servir de punto de partida para una misión real. La restricción de masa es la que guía todas las decisiones de diseño.

Los tres requisitos técnicos fundamentales son:
- MTF ≥ 0,25 (calidad óptica mínima)
- SNR ≥ 400 (relación señal-ruido)
- Cobertura completa de EE.UU. en ≤ 7 días con nubosidad 1/6

---

### [1:30] Diapositiva 3 — Misiones Semejantes

Antes de justificar las elecciones de diseño, es útil ver dónde se sitúa este trabajo respecto a las misiones existentes.

OCO-2, GOSAT, TanSat y CO₂M tienen masas entre 200 y 1750 kg con resoluciones entre 2 y 10 km. GHGSat llega a 25 m de GSD pero pesa 15 kg. Este diseño ocupa un nicho que no existe todavía: alta resolución espacial (80 m) con masa ultra-baja (1,71 kg).

**Mensaje clave:** ninguna misión operativa combina resolución sub-100 m con masa menor de 2 kg. Este TFG sitúa esa capacidad al alcance de un microsatélite comercial.

---

### [2:30] Diapositiva 4 — Flujo de trabajo del modelo analítico

El corazón del trabajo es el modelo analítico implementado en MATLAB, disponible en GitHub.

El pipeline sigue esta secuencia: a partir de los parámetros de entrada (GSD, bandas, detectores, altitudes), se calculan MTF y SNR, luego se calcula la cobertura con el toolkit RevisitTime, y finalmente se estiman la masa seca y la masa total. La función CrossData intersecta los tres mapas para extraer el espacio de diseño válido.

Como validación existe un modelo numérico de propagación de trazas terrestres, más lento pero más preciso para la geometría de cobertura real.

---

### [3:30] Diapositiva 5 — GSD

El GSD es el parámetro de entrada más importante del diseño, porque determina directamente la distancia focal del telescopio: f = p·h/GSD. Con GSD=80 m, h=520 km y p=15 µm, la focal resulta ser 97,5 mm —valor razonable para un telescopio miniaturizado.

**La elección de 80 m es deliberada.** Se probaron valores menores, como GSD=30 m, pero no dieron ninguna solución válida en todo el espacio de diseño: los telescopios necesarios eran demasiado grandes o la SNR era insuficiente para cumplir el requisito de 400 con los detectores disponibles.

80 m es 20 veces mejor que OCO-2 y permite atribuir emisiones a escala de ciudad o instalación industrial.

---

### [4:30] Diapositiva 6 — Bandas espectrales

Se monitorizan tres bandas: CO₂ a 1,61 µm y 2,01 µm, y O₂-A a 0,76 µm. Estas son exactamente las mismas bandas que OCO-2 y CO₂M, lo que permite la comparación directa y el intercalibrado en vuelo.

El ancho de banda de 20 nm por banda es el valor **óptimo**: suficientemente estrecho para la precisión espectral requerida, y no tan estrecho como para perder señal. Se probaron valores de 10 nm, pero la SNR caía por debajo del umbral de 400 en buena parte del espacio de diseño.

La banda O₂-A actúa como referencia de presión de columna para corregir el scattering atmosférico —la misma técnica que OCO-2 y GOSAT.

---

### [5:30] Diapositiva 7 — Detectores evaluados

Se evaluaron 6 tipos de detectores distintos. El criterio de selección principal fue la eficiencia cuántica en las tres bandas de interés —especialmente a 1,61 y 2,01 µm donde solo los detectores de InGaAs o HgCdTe son eficientes.

El detector H2RG (Hawaii-2RG) de Teledyne es el resultado de la optimización: ofrece la mejor combinación de eficiencia cuántica, ruido de lectura bajo y tamaño de píxel de 15 µm. Se usa en instrumentos espaciales como NIRSpec del JWST.

Los detectores CCD estándar se descartaron porque no son sensibles en el infrarrojo de onda corta (SWIR) donde están las bandas de CO₂.

---

### [6:00] Diapositiva 8 — Modo de escaneo Pushbroom

El modo pushbroom es el estándar en todos los instrumentos EO modernos: Sentinel-2, DESIS, OCO-2. No hay partes móviles, lo que aumenta la fiabilidad en el espacio.

El swath total es 2288 píxeles × 80 m = 183 km. Este valor cubre el requisito de cobertura semanal con la configuración de 2 satélites.

---

### [6:30] Diapositiva 9 — Tipo de órbita y LTAN

La órbita seleccionada es una SSO de amanecer-anochecer con LTAN de 6h. Esta elección maximiza la exposición solar: el panel solar ve el sol casi continuamente, minimizando los requisitos de batería.

El ángulo beta elevado también minimiza los eclipses. A 520 km la deorbitación natural ocurre en 257 días. Se exploró el rango completo de altitudes de 200 a 1000 km para encontrar el óptimo de masa.

---

### [7:15] Diapositiva 10 — MTF: Modelo

La MTF cuantifica la capacidad del sistema óptico de reproducir detalles espaciales. El requisito es MTF ≥ 0,25 a la frecuencia de Nyquist. Un valor menor de 0,25 implica imágenes con contraste insuficiente para distinguir emisores adyacentes.

El modelo factoriza la MTF total en cuatro componentes: difracción óptica, detector, alineamiento mecánico y margen de degradación. La componente dominante es la difracción cuando la pupila es pequeña.

El parámetro crítico es la relación entre la frecuencia de Nyquist del detector y la frecuencia de corte óptica: si la pupila es muy pequeña, el sistema no puede resolver el píxel y la MTF cae.

---

### [8:00] Diapositiva 11 — MTF: Mapas de calor

Estos mapas de calor muestran el valor de MTF en función de la altura orbital y el diámetro de pupila. Las celdas en blanco son combinaciones que no cumplen MTF ≥ 0,25 —son diseños descartados.

A altitudes bajas (200–400 km) se puede usar una pupila más pequeña y aun así cumplir el requisito, porque la frecuencia de Nyquist es menor relativa a la frecuencia de corte óptica. Los dos mapas corresponden a bandas espectrales distintas.

---

### [8:30] Diapositiva 12 — SNR: Modelo

La SNR es el factor limitante más importante: con tiempos de integración muy cortos —solo 10,5 ms a 520 km— la señal recogida por cada píxel es pequeña.

La SNR se calcula como la raíz cuadrada de los electrones de señal dividido por la suma cuadrática de todos los ruidos: shot noise, dark current, ruido de lectura, jitter, EMC y cuantización.

La irradiancia de referencia es 100 W·m⁻²·sr⁻¹·µm⁻¹, que corresponde a una escena terrestre típica bajo iluminación solar normalizada.

---

### [9:00] Diapositiva 13 — SNR: Mapas de calor

Estos mapas de SNR son análogos a los de MTF pero muestran la relación señal-ruido. La zona válida (SNR ≥ 400) es más restrictiva que la de MTF: requiere pupilas mayores, especialmente a altitudes altas.

La intersección de las zonas válidas de MTF y SNR define el subespacio de diseño. La banda CO₂-2 (2,01 µm) es la banda limitante porque la luminancia de referencia es menor en el SWIR.

---

### [9:30] Diapositiva 14 — Tiempo de revisita: Metodología

El tercer requisito es la cobertura: cubrir el 100% de EE.UU. continental en 7 días o menos, con factor de nubosidad de 1/6. Este factor significa que solo 1 de cada 6 pasadas se asume libre de nubes, multiplicando por 6 el tiempo efectivo de revisita.

Se usa el RevisitTime Toolkit de la Universidad de Manchester (Crisp & Livadiotti, 2017), que calcula analíticamente el tiempo de revisita para una constelación SSO. Los resultados del modelo numérico de propagación son consistentes con el analítico.

---

### [10:00] Diapositiva 15 — Tiempo de revisita: Resultados

Los mapas de revisita muestran el número de días para cubrir EE.UU. completo. La configuración óptima es 2 satélites con 2 telescopios cada uno.

Con 2 satélites separados 180° en la órbita, el tiempo de revisita se reduce a la mitad. Con 2 telescopios, el swath total es 183 km —suficiente para cubrir EE.UU. en 7 días desde 520 km.

A altitudes bajas, el swath es insuficiente. A altitudes altas, el telescopio es demasiado grande para MTF y SNR.

---

### [10:30] Diapositiva 16 — Cruce de datos: h vs D_min

La función CrossData intersecta los tres mapas de MTF, SNR y cobertura para extraer, para cada altura h, el diámetro mínimo de pupila D_min(h) que satisface **simultáneamente** los tres criterios.

Existe una altitud óptima que minimiza la masa total: la combinación de masa del instrumento (que decrece con h) y masa de combustible (que crece con altitudes bajas). Este equilibrio aparece a h = 520 km.

---

### [11:00] Diapositiva 17 — Comparativa de detectores

Esta tabla compara el rendimiento de todos los detectores evaluados con el telescopio refractivo, mostrando SNR, MTF y masa para cada combinación.

El H2RG resulta el detector óptimo para todas las bandas, especialmente en el SWIR. Los CCD se descartan por su insensibilidad a 1,61 y 2,01 µm.

---

### [11:20] Diapositiva 18 — Masa seca: Método de escalado

El método de escalado geométrico extrae la masa del instrumento interpolando entre dos misiones de referencia: el Thematic Mapper (0,4 m, 240 kg) y SEOSAT (0,25 m, 100 kg). Es el método estándar del SMAD para diseño preliminar.

El coeficiente K=2 cuando R<0,5 refleja que miniaturizar mucho penaliza la masa específica. Para nuestros diámetros (~28 mm), K=2.

La masa total del satélite se estima como 4× la masa del instrumento —relación estándar SMAD para pequeños satélites.

---

### [11:50] Diapositiva 19 — Masa seca: Resultado

El mínimo de masa seca aparece a 520 km con D_min = 28 mm: la pupila es suficientemente pequeña para que el instrumento sea ligero. Con la configuración de 2 telescopios, la masa seca total es 1,349 kg.

---

### [12:10] Diapositiva 20 — Masa total: Modelo de decaimiento

El modelo calcula el ΔV acumulado necesario para mantener la órbita durante 8 años usando el modelo atmosférico estándar USSA76 integrado con ode45. Los impulsos de Hohmann se aplican cuando la altitud cae un 2%.

La masa de combustible usa Tsiolkovsky con Isp=209 s del motor ECAPS HPGP y un margen del 10%.

---

### [12:40] Diapositiva 21 — Masa total: Resultado

El mínimo absoluto es a 520 km: **1,71 kg por satélite**. La constelación completa pesa 3,42 kg en órbita.

El mínimo es robusto: el rango 480–560 km da masas dentro del 5% del valor óptimo.

---

### [13:00] Diapositiva 22 — Otras soluciones: verificación del mínimo absoluto

Esta diapositiva verifica que 1,71 kg es el mínimo absoluto global, no un mínimo local. Se muestran todas las combinaciones evaluadas: distintos detectores, tipos de telescopio y configuraciones de constelación.

La configuración H2RG + refractivo + 2 sat/2 tel es la ganadora en todas las comparaciones. Las reflectivas son ~15% más pesadas por el oscurecimiento central que penaliza la MTF.

---

### [13:20] Diapositiva 23 — Disposición final del instrumento

La disposición física del instrumento es compatible con un formato 1,5U CubeSat (10×10×15 cm). Los dos telescopios están orientados en paralelo para cubrir el swath total de 183 km.

Esta disposición es conceptual —el diseño mecánico detallado queda como trabajo futuro.

---

### [13:40] Diapositiva 24 — Punto de diseño final

Esta tabla resume el punto de diseño final. Todos los requisitos se cumplen simultáneamente: MTF ≥ 0,25, SNR ≥ 400, revisita ≤ 7 días. La masa de 1,71 kg incluye el propelente para 8 años y la maniobra de phasing.

---

### [14:00] Diapositiva 25 — Deorbitación natural

A 520 km, la simulación de decaimiento natural muestra que el satélite reingresará en **257 días** sin ninguna maniobra activa. Esto cumple la normativa NASA/ESA de menos de 1 año.

No se necesita un kit de deorbitación ni propelente adicional —ventaja directa de operar a 520 km.

---

### [14:20] Diapositiva 26 — Subsistemas: Componentes seleccionados

Todos los componentes son COTS, disponibles en catálogo. El diseño de subsistemas es **preliminar**: la gestión térmica, el presupuesto de potencia detallado y el ADCS se han dimensionado con criterios de primer orden, no con un diseño de ingeniería completo. Los números son coherentes y demuestran la viabilidad del concepto.

El ADCS usa gradiente gravitatorio + magnetopares GMAT-1 en lugar de ruedas de reacción: apuntamiento nadir ≤1°, solo 50 g, sin consumo activo. Para análisis de tendencias atmosféricas con corrección geométrica en tierra, 1° de error de apuntamiento (9 km de geolocalización a 520 km) es aceptable.

---

### [15:10] Diapositiva 27 — Desglose de masa por subsistema

El propelente (21%) y el instrumento (20%) son los dos subsistemas dominantes. El sistema de comunicaciones (15,8%) ocupa el tercer lugar —la antena X-Band es el componente más masivo.

La batería pesa solo 13 g porque se dimensiona únicamente para el déficit durante las transmisiones (10,38 W × 0,25 h = 2,6 Wh). La órbita amanecer-anochecer elimina la necesidad de grandes baterías.

---

### [15:40] Diapositiva 28 — Lanzadores evaluados

Se evaluaron 6 lanzadores por cuatro criterios: fiabilidad, precio/kg en SSO, compatibilidad CubeSat/PPOD y disponibilidad rideshare para cargas < 5 kg.

---

### [16:00] Diapositiva 29 — Bases de lanzamiento evaluadas

Vandenberg SFB (California) es la única base con azimut compatible con SSO (198°) sin sobrevolar tierra habitada. Kennedy y Kourou no pueden alcanzar la inclinación polar de 97,6°.

---

### [16:20] Diapositiva 30 — Selección final del lanzador: Falcon 9

El Falcon 9 Transporter (SpaceX) se selecciona por su máxima fiabilidad (99,5%, >300 vuelos), precio rideshare más competitivo (~6.500 $/kg) y compatibilidad total con PPOD.

El coste estimado para los dos satélites (3,42 kg) es de ~22.000 $. Cada gramo adicional tiene un coste real —de ahí la importancia de minimizar la masa.

---

### [16:40] Diapositiva 31 — Vandenberg SFB

El complejo SLC-4E está dedicado al Falcon 9 con historial de misiones SSO: Planet Labs, Sentinel-6, SMAP. El azimut de 198° garantiza la inclinación de 97,6° pasando sobre el Pacífico.

---

### [17:00] Diapositiva 32 — Maniobra de phasing

Los dos satélites se lanzan juntos y se separan con una maniobra de phasing: Sat-2 baja a una órbita de drift, avanza 180° en anomalía media y vuelve a 520 km. La separación de 180° maximiza la cobertura temporal.

El ΔV total es ~3 m/s, incluido en el presupuesto de combustible. El tiempo de separación es ~3 semanas.

---

### [17:20] Diapositiva 33 — Estaciones terrestres evaluadas

Criterios: latitud alta (más pases/día), banda X disponible, infraestructura existente (sin CAPEX). Una estación propia se descartó por el coste de instalación.

---

### [17:40] Diapositiva 34 — Selección de Fairbanks

Fairbanks (64,8°N, red NOAA) con 22,1 h/semana de contacto es más que suficiente para descargar los 6,9 GB semanales. La infraestructura banda X existe y opera actualmente.

---

### [18:00] Diapositiva 35 — Velocidad de descarga de datos

150 Mbps × 22,1 h/semana = >14 TB de capacidad, vs 6,9 GB generados. Factor de margen ×2000. La descarga no es el cuello de botella de la misión.

---

### [18:15] Diapositiva 36 — Gestión de datos a bordo

Pipeline: instrumento → OBC (compresión ×3) → almacenamiento 512 GB → descarga en Fairbanks → procesamiento en tierra (corrección geométrica, calibración radiométrica, corrección atmosférica) → productos de columna CO₂/O₂.

---

### [18:30] Diapositiva 37 — Conclusiones

En resumen, el trabajo demuestra que es posible diseñar una constelación de 2 nanosatélites capaz de monitorizar CO₂ y O₂ sobre EE.UU. con GSD=80 m —20 veces mejor que OCO-2— con una masa total de solo **1,71 kg/satélite**. Todos los requisitos técnicos se cumplen simultáneamente.

El código está completamente disponible en GitHub bajo licencia abierta.

Quiero agradecer explícitamente las herramientas de inteligencia artificial utilizadas en este TFG: **Claude Sonnet/Thinking 3.7 y 4.0** y **Claude Opus** para redacción y código; **Perplexity Deep Research** para revisión bibliográfica; y **Google Gemini Pro 2** para síntesis de documentación técnica. Estas herramientas permitieron completar el trabajo en aproximadamente **287 horas**.

---

### [19:30] Diapositiva 38 — Gracias

Muchas gracias por su atención. Quedo a disposición del tribunal para responder preguntas.

El repositorio contiene todos los ficheros: código MATLAB, paquete Python `eo-mission`, la memoria completa en LaTeX, datos de salida en CSV y estas diapositivas.

---

## Líneas de mejora / Trabajo futuro

Posibles mejoras identificadas durante el trabajo:

### Modelo analítico
- Extender el análisis de MTF para incluir degradación por jitter real medido (en vez de estimado)
- Incorporar el efecto del ángulo de iluminación solar en la SNR (variación estacional del ángulo beta)
- Añadir modelo de transmitancia atmosférica dependiente de la banda (actualmente constante)

### Modelo numérico de cobertura
- Completar la integración del modelo de propagación numérica (Propagation ver 2.0) con el pipeline principal
- Validar los resultados del toolkit RevisitTime con el modelo numérico para todo el rango de altitudes
- Extender la malla de EE.UU. a resolución de 0,1° para mayor precisión geográfica

### Instrumento óptico
- Diseño detallado del sistema óptico refractivo con tolerancias de fabricación
- Análisis de calibración radiométrica en vuelo usando vicarious calibration sobre sites conocidos
- Estudio de la función de respuesta espectral del filtro de banda y su impacto en la recuperación de columnas

### Subsistemas
- Diseño térmico detallado: modelo nodal del satélite, identificación de zonas críticas
- Presupuesto de potencia completo incluyendo todos los modos operativos (safe, nominal, transmisión)
- Análisis de controlabilidad del ADCS magnetorquer en pases polares (geometría del campo magnético)
- Evaluación de un actuador de rueda de reacción de bajo perfil como alternativa para misiones que requieran < 0,1° de apuntamiento

### Sistema
- Análisis de fiabilidad y FMEA a nivel sistema
- Estudio de autonomía de vuelo (detección y respuesta a fallos on-board)
- Validación del presupuesto de enlace completo con modelo de propagación de señal (link margin)
- Extensión de la misión a otras regiones geográficas (Europa, Asia)

---

## Posibles preguntas del tribunal

### P1: ¿Por qué GSD = 80 m y no mayor o menor?

**Respuesta:** El GSD es la variable de diseño más importante porque determina la distancia focal. Con GSD=80 m, la focal resulta 97,5 mm a 520 km —manejable para un telescopio miniaturizado. Se probaron valores menores como GSD=30 m, pero no generaron ninguna solución válida en todo el espacio de diseño: el telescopio resultante era demasiado grande o la SNR insuficiente. Valores mayores (120 m) sí daban soluciones pero con menor resolución. El 80 m es además un valor significativo: es 20× mejor que OCO-2 y permite la atribución de emisiones a escala de ciudad, que es el caso de uso más valioso para inventarios de emisiones.

---

### P2: ¿Cómo se valida que SNR ≥ 400 es suficiente para recuperar concentraciones de CO₂?

**Respuesta:** El umbral de SNR ≥ 400 proviene de los estudios de sensibilidad de misiones como OCO-2 y CO₂M (Taylor et al., 2016; Butz et al., 2011). Para recuperar concentraciones de CO₂ con una precisión de < 1 ppm usando métodos de inversión óptima, se requiere una SNR mínima de 200–400 en las bandas SWIR. Nosotros usamos 400 como umbral conservador. La validación rigurosa requeriría simulaciones end-to-end con modelos de transferencia radiativa (como VLIDORT o DISORT), lo que queda como trabajo futuro.

---

### P3: ¿El modelo de escalado de masa es suficientemente preciso para un diseño preliminar?

**Respuesta:** El escalado geométrico del SMAD es el método estándar en la fase de diseño preliminar (Fase A/B), donde las incertidumbres son del orden del ±30%. Las dos misiones de referencia (SEOSAT y Thematic Mapper) tienen parámetros similares al diseño objetivo, lo que reduce la extrapolación. Para la fase de diseño detallado, el siguiente paso sería un análisis de masa de componentes bottom-up con datos reales de los fabricantes de los componentes COTS seleccionados.

---

### P4: ¿El ADCS con gradiente gravitatorio + magnetopares puede mantener el apuntamiento durante toda la misión?

**Respuesta:** Para esta misión de monitorización de tendencias atmosféricas con procesamiento geométrico en tierra, un apuntamiento nadir de ≤1° es suficiente. El gradiente gravitatorio proporciona apuntamiento pasivo de 3–5° y los magnetopares GMAT-1 lo refinan a ≤1°. La geometría de la órbita amanecer-anochecer es favorable para los magnetopares porque el campo geomagnético es casi perpendicular al plano orbital durante la mayor parte de la órbita. En pases polares hay una zona de reducida controlabilidad de un eje, pero la duración es breve. Si la misión evolucionara hacia productos con geolocalización de sub-píxel sin corrección post-procesado, habría que considerar un sistema con rueda de reacción (BCT XACT o similar, ~200 g).

---

### P5: ¿Por qué Isp = 209 s y no el valor habitual de 220 s para HPGP?

**Respuesta:** El valor nominal del catálogo de ECAPS para el motor HPGP 100 mN con propelente LMP-103S es de 209 s en condiciones de referencia (300 K, 20 bar de presión de entrada). El valor de 220 s que aparece en algunas referencias corresponde a condiciones óptimas de precalentamiento o a versiones anteriores del propelente. Se optó por el valor conservador de 209 s para tener un presupuesto de combustible más realista, lo que resulta en un ~5% más de propelente respecto a usar 220 s.

---

### P6: ¿Por qué no se usa un lanzador dedicado para mayor flexibilidad orbital?

**Respuesta:** Para una carga de 3,42 kg, un lanzador dedicado de pequeña capacidad (como Electron) costaría del orden de 7–10 M$ —dos órdenes de magnitud más que el rideshare. El programa Transporter del Falcon 9 ofrece lanzamientos frecuentes (~cada 6 meses) a SSO, con lo que la ventana de lanzamiento es flexible y el coste es mínimo. El único compromiso es la incertidumbre en la fecha exacta de lanzamiento, que es aceptable para una misión de monitorización que no tiene dependencia temporal crítica.

---

### P7: ¿Cómo se asegura la coherencia de los datos entre los dos satélites de la constelación?

**Respuesta:** Los dos satélites tienen la misma configuración hardware (H2RG, refractivo 28 mm, misma óptica). El intercalibrado en vuelo se realizaría usando vicarious calibration sobre sitios terrestres de reflectancia conocida (por ejemplo, RVUS de RadCalNet) y comparando con medidas cruzadas cuando ambos satélites observan la misma zona en días consecutivos. Esta es la técnica estándar usada en constelaciones como Planet Labs y Sentinel-2.

---

### P8: ¿El modelo de cobertura considera la variación estacional de la nubosidad?

**Respuesta:** El modelo usa un factor de nubosidad constante de 1/6, que es el valor promedio anual sobre EE.UU. No considera la variación estacional ni geográfica. En la realidad, la nubosidad sobre los Grandes Llanos es mucho menor que sobre las costas, y varía con la estación. Un modelo más detallado utilizaría datos de nubosidad de MODIS o CERES para calcular el porcentaje de cielo despejado por celda de la malla. Esto queda como una de las líneas de mejora identificadas.

---

### P9: ¿Por qué 8 años de vida útil de la misión?

**Respuesta:** Los 8 años de vida útil son una elección de diseño representativa para una misión de monitorización atmosférica. Misiones como OCO-2 (>10 años operativa) y GOSAT (>15 años) demuestran que los instrumentos de este tipo son muy longevos. Para el diseño de la masa de combustible, 8 años es un valor conservador: si la misión resulta ser más larga, el propelente se agota antes del fin de misión y el satélite simplemente entra en la fase de deorbitación natural pasiva.

---

### P10: ¿Se ha considerado la degradación del detector con la dosis de radiación en 8 años?

**Respuesta:** No en detalle —este es uno de los aspectos que quedan como trabajo futuro. A 520 km en SSO, la dosis de radiación total en 8 años es del orden de 5–10 krad dependiendo del blindaje. El detector H2RG tiene una tolerancia a radiación de >50 krad en las versiones diseñadas para espacio (como las del JWST), por lo que la degradación no debería ser crítica con un blindaje adecuado. Sin embargo, el análisis SHIELDOSE detallado no está incluido en este diseño preliminar.

---

### P11: ¿Por qué una sola estación terrena? ¿Qué pasa si hay una falla de comunicaciones?

**Respuesta:** Una sola estación es suficiente para esta misión porque el margen de descarga de datos es enorme (factor 2000×). La memoria a bordo de 512 GB puede almacenar semanas de datos en caso de pérdida de contacto. Para una misión real, se añadiría Wallops (Virginia) como estación de respaldo dentro de la misma red NOAA, sin coste adicional de infraestructura. Esta es una línea de mejora de robustez del sistema.

---

### P12: ¿Cómo se compara la relación coste-resolución de esta misión con las alternativas comerciales?

**Respuesta:** GHGSat-D ofrece 25 m de GSD (mejor que nosotros) con una masa de ~15 kg, un coste de misión estimado en varios millones de dólares y datos de pago. Este diseño, con GSD=80 m y masa de 1,71 kg, costaría del orden de 200.000–500.000 $ en total (incluyendo hardware COTS, integración y lanzamiento), lo que lo sitúa en un nicho de monitorización de bajo coste accesible a agencias regionales o universidades. El compromiso es la resolución de 80 m vs los 25 m de GHGSat, que es aceptable para inventarios de emisiones nacionales.

---

*Fin del guión*
