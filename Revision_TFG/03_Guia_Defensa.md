# Guía de defensa — Diseño Preliminar de un Satélite de Observación de la Tierra

> Guion hablado por diapositiva (2–4 frases) + banco de preguntas del tribunal.
> **Todos los números usan los valores YA corregidos** (ver `05_Cambios_para_Defensa.md` para actualizar el .pptx).
> Tiempo objetivo típico: ~15–20 min. Marca con ⏱️ las diapositivas donde puedes acelerar si vas justo.

---

## Cifras maestras (ten estas en la cabeza)
- Constelación: **2 satélites**, SSO *dawn-dusk*, **520 km**, LTAN 06:00, i = **97,48°**, periodo **95 min**.
- Carga útil: GSD **80 m**; bandas **0,76 / 1,61 / 2,01 µm**, ancho 20 nm; detector **Teledyne H2RG**; telescopio **refractivo** D = **28 mm**, f = 12 cm, **F/4,2**, FoV combinado 20°.
- Calidad: **MTF ≥ 0,25** (banda restrictiva 2,01 µm); **SNR ≥ 400** (se obtiene 717 en 0,76 µm).
- Revisita: efectiva **6,6 días**, física **5,47 días**.
- Masa/satélite: seca **1,349 kg** + combustible **0,361 kg** = **1,710 kg** → constelación **3,42 kg**.
- Propulsión: ECAPS HPGP, **Isp 209 s**, **Δv 436 m/s**, 76 impulsos, deorbitado natural **257 días**.
- Comms: **IQ Spacecom XLink-X** banda X (200 g, 16 W, 200 Mbps) + UHF EnduroSat (85 g).
- Lanzador: **Falcon 9** desde **Vandenberg**, rideshare **6 500 USD/kg**; fase 180° en 37 h (órbita 520×720 km).
- Segmento tierra: **Fairbanks**, **79 620 s/semana** (22,1 h) de contacto, **90 Mbps** reales; datos **6,9 GB/cobertura → 8,8 GB/semana**, descarga 786 s, memoria 512 MB.

---

## Guion por diapositiva

**1. Portada.** Preséntate y enmarca: «Diseño preliminar de un satélite de observación de la Tierra para monitorizar CO₂ antropogénico sobre EE. UU.». Nombra tutores.

**2. Objetivos.** «El cliente —U.S. Observatory Against Global Change— pide un mapa semanal de CO₂ sobre EE. UU. continental, con MTF ≥ 0,25 y SNR ≥ 400, minimizando la masa puesta en órbita.» Menciona el factor de nubes (1 día cubierto por cada 5).

**3. Misiones semejantes (portadilla).** «Antes de diseñar, reviso el estado del arte para inspirar la solución.»

**4. Comparativa de misiones.** «Ninguna misión existente resuelve directamente el problema: van de satélites grandes (GOSAT-2, ~1 800 kg) a constelaciones de nanosatélites (GHGSat). Todas usan SSO. Esto orienta hacia una constelación ligera.»

**5. Conceptos previos (portadilla).** ⏱️ «Defino los conceptos que sustentan los requisitos: calidad de imagen, detectores, óptica y órbitas.»

### Carga de pago (Payload) — slides 6–27 (núcleo técnico)

**6. Carga de pago (portadilla).** «Diseño el sistema óptico y elijo el detector partiendo de los requisitos del cliente.»

**7. GSD.** «El GSD lo fija implícitamente el cliente: hay que distinguir fuentes puntuales. Pixelando una imagen de Houston a distintos GSD, **80 m** distingue zonas industriales de residenciales sin gastar masa de más.» Es la justificación visual clave.

**8. Bandas espectrales.** «Elijo tres bandas: **O₂ a 0,76 µm** para calibrar la columna de aire (el O₂ tiene concentración conocida, 20,95 %), y **CO₂ a 1,61 y 2,01 µm**. Ancho de banda 20 nm: equilibrio entre resolución espectral y SNR.»

**9. Detectores y telescopios.** «Evalúo tres detectores SWIR (CAPYORK, H2RG, Saturn VISIR) y cuatro telescopios (refractivo, Cassegrain, Korsch, TMA), cada uno con su MTF de alineamiento y FoV.»

**10. Filtros y escaneo.** «Escaneo **pushbroom**: sin partes móviles, mayor tiempo de integración → mejor SNR. Filtros **microstrip** integrados en el detector: sin mecanismos, robustos para los 8 años.»

**11. Órbita y hora de paso.** «SSO por iluminación constante (los detectores dependen de la luz solar reflejada). **Dawn-dusk (LTAN 06:00)**: el satélite bordea el terminador, evita eclipses y maximiza potencia.»

**12. Flujo de scripts.** «Desarrollé un paquete propio en MATLAB: a partir de los requisitos del cliente y el espacio de diseño (alturas y diámetros), calcula MTF, SNR, revisita, masa seca y total, y cruza los datos para dar las soluciones óptimas.» (open-source en GitHub).

**13–14. Cálculo de MTF.** «La MTF total es el producto de difracción, aberraciones, fabricación, alineamiento, vibraciones, termoelástico y detector, con un margen del 10 %. La banda más restrictiva es la de mayor λ, **2,01 µm**. Los mapas de calor muestran qué combinaciones altura–diámetro cumplen ≥ 0,25.»

**15–16. Cálculo de SNR.** «Modelo la señal (electrones por píxel vía radiancia, transmisión óptica y eficiencia cuántica) y todas las fuentes de ruido. La banda restrictiva aquí es la de menor λ, **0,76 µm**. Es el requisito menos exigente: se cumple con holgura (717 > 400), así que descarto el TDI.»

**17–18. Revisita.** «La revisita de 7 días es lo más limitante. Uso un código semianalítico (con *wrapper* propio) que, para cada altura y *swath*, da la revisita sobre EE. UU., descontando un 1/6 de pasadas por nubes y un 5 % de solapamiento.»

**19. Cruce de datos (h vs Dmin).** «Para cada altura, busco el **diámetro mínimo** que cumple MTF, SNR y revisita simultáneamente.»

**20. Comparativa de detectores (refractivo).** «Con telescopio refractivo, el **H2RG** da el menor diámetro para todas las alturas → es el elegido. El TMA necesitaría diámetros casi 3× mayores.»

**21. Optimización de masa.** «Con la relación h–Dmin, estimo la masa de cada solución para elegir la altura que minimiza la masa total.»

**22. Masa seca.** «Estimo la masa seca por escalado a partir de dos misiones de referencia (Thematic Mapper y SEOSAT), asumiendo instrumento = 25 % del satélite y +50 % por el segundo telescopio.»

**23. Masa de combustible.** «Modelo el arrastre atmosférico (densidad exponencial USSA76, Cd = 2,5), integro con ode45, y cuando la altura cae un 2 % aplico un impulso tipo Hohmann. Con Tsiolkovski (**Isp 209 s**) obtengo el combustible para 8 años.» ⚠️ *La diapositiva pone Isp 220 s en un sitio: corrige a 209.*

**24. Masa total.** «El mínimo de masa total está en **520 km**: seca 1,349 kg + combustible 0,361 kg = **1,710 kg/satélite**, diámetro de pupila 28 mm.»

**25. Disposición del instrumento.** «Dos telescopios refractivos embarcados, dispuestos para duplicar el FoV; *swath* final 183 km sobre 170 km requeridos.»

**26. Punto de diseño.** «Aquí se ve el punto final cumpliendo las tres restricciones: MTF en 2,01 µm, SNR en 0,76 µm y cobertura.»

**27. Deorbitación.** «Tras fin de misión, el satélite reentra de forma natural en **257 días**, muy por debajo del año que me impuse → no hace falta combustible extra de deorbitado.»

### Subsistemas — slides 28–29

**28. Subsistemas (portadilla).** «Con la masa y volumen definidos, dimensiono los subsistemas con componentes comerciales (*off-the-shelf*).»

**29. Componentes.** «Desglose para 1,710 kg: instrumento, eléctrico, comunicaciones, estructura (15 %), ADCS, propulsión, OBC y térmica.» ⚠️ *La diapositiva tiene datos antiguos (comms en banda S, 40 g; combustible 0,334): actualiza a **XLink-X banda X 200 g / 16 W**, UHF 85 g, combustible 0,361 kg.* (ver `05_Cambios_para_Defensa.md`).

### Lanzadores — slides 30–35

**30. Lanzadores (portadilla).** «Selecciono el vehículo y la base de lanzamiento.»

**31. Lanzadores evaluados.** «Comparo seis: Falcon 9, Firefly Alpha, Electron, Vulcan, Vega-C, Ariane 6, por fiabilidad, coste/kg y rideshare.»

**32. Bases de lanzamiento.** «Tres bases: Kennedy (azimuts no aptos para SSO), Vandenberg y Kourou.»

**33. Selección del lanzador.** «**Falcon 9**: fiabilidad 99,35 %, programa rideshare a ~6 500 USD/kg, reignición de 6 h para afinar la órbita, y autonomía estratégica de EE. UU.»

**34. Base de lanzamiento.** «**Vandenberg**: azimut compatible con SSO y operaciones previas con F9.»

**35. Puesta en órbita.** «Despliego los dos satélites con una maniobra de fase: la etapa superior sube a una órbita 520×720 km y, tras alcanzar 180° de desfase (**37 h**), recircula y suelta el segundo satélite.» *Azimut **189°** (corregido; antes figuraba 198° por error de transcripción). Compatible con el rango WTR de Vandenberg. Regenera la figura corriendo `Azimuth.m`.*

### Segmento tierra — slides 36–41

**36. Segmento tierra (portadilla).** «Genero y descargo los datos.»

**37. Estaciones evaluadas.** «Fairbanks, Wallops, red EUMETSAT o estación propia.»

**38. Selección de Fairbanks.** «**Fairbanks**: latitud alta → muchos contactos para órbita polar; banda X; territorio EE. UU. Simulación con Aerospace Toolbox: **22,1 h/semana** de contacto (elevación mínima 10°).»

**39. Generación de datos.** «3 bandas × 12 bits sobre 9,83 M km² a 80 m → **6,9 GB por cobertura**; con revisita física 5,47 días, **8,8 GB/semana**.»

**40. Velocidad de descarga.** «Banda X ~150 Mbps teóricos; con el modelo meteorológico (ITU-R P.838-3) la media anual es **137 Mbps**.» ⚠️ *La tabla de la diapositiva tiene cifras antiguas (lluvia 75–100, tormenta 22–38): actualiza a las de la memoria.*

**41. Descarga.** «Con factores correctores de hardware (12,5 %) y baja elevación (25 %) → **90 Mbps reales**. Descargar 8,8 GB lleva 786 s frente a 79 620 s de contacto: **sobra margen**. Una memoria de 512 MB basta.»

### Cierre — slides 42–43

**42. Conclusiones.** → *Contenido propuesto abajo (la diapositiva está vacía).*

**43. ¡Gracias!** «Gracias. El código está disponible en GitHub. Quedo a su disposición para preguntas.»

---

## Contenido para la diapositiva 42 (Conclusiones) — actualmente vacía

Puntos para rellenar la diapositiva y para el cierre hablado:

**Cumplimiento de requisitos del cliente:**
- ✅ Cobertura semanal de EE. UU. continental (constelación de 2 satélites, revisita efectiva 6,6 días con nubes).
- ✅ **MTF ≥ 0,25** y **SNR ≥ 400** (se alcanza 717), con márgenes del 10 %.
- ✅ Masa minimizada: **1,71 kg/satélite**, 3,42 kg la constelación — un orden de magnitud por debajo de las misiones comparables.

**Lecciones / valor del trabajo:**
- Metodología **iterativa** y herramienta MATLAB propia que recorre el espacio de diseño y da la solución óptima.
- Síntesis **multidisciplinar**: óptica, mecánica orbital, propulsión, potencia, comunicaciones.

**Limitaciones (diseño preliminar):**
- Apuntamiento ≤ 2° (geolocalización refinada en postproceso); aberración cromática del refractivo; modelo de arrastre simplificado; sin análisis estructural/térmico detallado.

**Líneas futuras:**
- Rueda de reacción para cartografía de precisión (<0,1°); validación óptica con software de trazado de rayos; análisis de costes y de ciclo de vida.

---

## Banco de preguntas del tribunal (con respuesta)

**P1. ¿Por qué 80 m de GSD?** Lo fija el cliente implícitamente: hay que distinguir fuentes puntuales. El análisis visual sobre Houston muestra que 80 m separa instalaciones industriales de zonas residenciales; más resolución gastaría masa sin aportar al objetivo.

**P2. ¿Por qué el detector H2RG y no un CCD?** El H2RG (HgCdTe) cubre las tres bandas (0,4–2,5 µm), tiene la mejor MTF (0,5) y, en el cruce de datos, da el **menor diámetro de pupila** para todas las alturas. Un CCD de silicio no llega al SWIR (1,61 y 2,01 µm).

**P3. ¿Por qué telescopio refractivo y no TMA/Cassegrain?** El refractivo no tiene obscuración central → mayor MTF de difracción, y su penalización de alineamiento es baja (0,90). El TMA exigiría diámetros casi 3× mayores. La aberración cromática no penaliza porque las bandas son estrechas (20 nm) y se filtran por separado.

**P4. Justifique SNR ≥ 400 y MTF ≥ 0,25.** Son requisitos del cliente. La MTF es el requisito **limitante** (banda 2,01 µm, mayor λ → menor frecuencia de corte). La SNR se cumple con holgura (717 en 0,76 µm), por eso descarto el TDI.

**P5. ¿Por qué SSO dawn-dusk a 520 km?** SSO para iluminación solar constante (imprescindible para detectores pasivos). Dawn-dusk (LTAN 06:00) evita eclipses y maximiza generación. 520 km es el **mínimo de la curva masa–altitud**: más bajo dispara el combustible por arrastre, más alto exige mayor pupila.

**P6. ¿Por qué dos satélites?** La revisita de 7 días es inalcanzable con uno (las configuraciones de 1 satélite no bajan de 7 días en simulación). La configuración mínima viable es 2 satélites con 2 telescopios cada uno, desfasados 180° en anomalía verdadera.

**P7. ¿Por qué Falcon 9 desde Vandenberg?** Fiabilidad 99,35 %, rideshare a ~6 500 USD/kg, etapa superior con reignición de 6 h para afinar la órbita, y autonomía estratégica de EE. UU. Vandenberg da azimut compatible con SSO e infraestructura existente para F9.

**P8. ¿Por qué Fairbanks?** Latitud alta (64,8°N) → múltiples contactos diarios con órbita cuasi-polar; antenas de banda X; en territorio EE. UU. La simulación da 22,1 h/semana de contacto, frente a los 786 s necesarios para descargar.

**P9. Supuestos del modelo de arrastre.** Densidad exponencial por capas basada en USSA76, Cd = 2,5, integración con ode45; impulso (Hohmann) cuando la altura cae un 2 %; combustible por Tsiolkovski con Isp = 209 s y margen del 10 %.

**P10. ¿Qué márgenes tiene el diseño?** 10 % en MTF, SNR y combustible; componentes comerciales con datos de fabricante; la suma bottom-up de subsistemas coincide con la estimación por escalado (1,71 kg).

**P11. Banda 2,01 vs 2,06 µm.** El diseño se computó a 2,01 µm (próxima a la de MicroCarb). OCO-2/GOSAT usan 2,06 µm; ambas están en la región de absorción del CO₂ en SWIR. *(Si el tutor prefiere 2,06, habría que recomputar la MTF.)*

**P12. ¿Por qué no hay análisis de enlace (link budget) ni Doppler?** Se consideró fuera del alcance de un diseño **preliminar**; las pérdidas relevantes (meteorología, hardware, baja elevación) se cubren con los factores correctores aplicados a la velocidad de descarga.

**P13. ¿Es realista una masa de 1,71 kg?** Es coherente con nanosatélites de 1–10 kg; cada subsistema usa componentes off-the-shelf reales (XLink-X, EnduroSat UHF, ECAPS, iOBC, GMAT-1). El instrumento (28 mm de pupila) es pequeño porque el GSD es modesto (80 m).

**P14. ¿Cómo se controla la actitud con tan poco ADCS?** Estabilización pasiva por gradiente gravitatorio + magnetopares triaxiales (B-dot para detumbling), magnetómetro y sensor solar. Logra ~1–2°, suficiente porque la geolocalización se refina en postproceso.
