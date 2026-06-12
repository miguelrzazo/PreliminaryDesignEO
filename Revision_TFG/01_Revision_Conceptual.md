# Pasada 1 — Revisión conceptual y fundamental

> **Naturaleza de esta pasada:** solo reporte. **No** he editado ningún `.tex`.
> Revisa los puntos, marca los que apruebas y los aplico en el paso de correcciones.
> Clasificación: **🔴 Crítico** (afecta a un resultado o afirmación central) · **🟠 Importante** (incoherencia o imprecisión real) · **🟡 Menor** (matiz / aproximación).
>
> Verifiqué la aritmética con Python; los valores recalculados aparecen entre paréntesis.

> ### ESTADO DE APLICACIÓN (actualizado)
> **✅ Aplicados en los `.tex`:** C2 (banda→2,01 µm), C3 (fórmula→T#²), I2, I3, I4, I5, I6, I7, I8, I9, I10, y menores M1, M2, M3, M4, M5, M6, M8, M9, M10.
> **⏸️ Pendientes de tu decisión:** **C1** (masa: "escalado oficial" no es aplicable sin forzar componentes — ver nota abajo), **I1** (azimut 198° hardcodeado en `Azimuth.m`, no se deriva de la fórmula), parte de **I11** (celdas de coste Vulcan/Ariane en la tabla — no las cambio sin valor verificado).
> **⏭️ Omitido a propósito:** M7 (FoV), M11 (Vulcan ya matizado en prosa).

---

## 🔴 Críticos

### C1 — Dos masas secas/totales distintas en el documento
- **Dónde:** `5.Mission/00AnalisisMision.tex:136` y Resumen Ejecutivo (`1.Introduccion/Introduccion.tex:66-67`) dicen **masa seca 1,349 kg / total 1,710 kg**. La tabla de subsistemas `5.Mission/00AnalisisMision.tex:226-228` (suma componente a componente) da **masa seca 1,455 kg / total 1,816 kg**.
- **Problema:** la estimación *top-down* por escalado (1,349) y el presupuesto *bottom-up* (1,455) difieren en ~0,106 kg (~8 %). El Resumen Ejecutivo y el capítulo de optimización usan 1,710; el desglose real de subsistemas suma 1,816. La constelación sería 3,42 kg o 3,63 kg según cuál se tome.
- **Sugerencia:** decidir cuál es la cifra oficial. Lo más honesto: presentar la de escalado como estimación preliminar y la suma de subsistemas como presupuesto detallado, y **reconciliar el Resumen Ejecutivo con la tabla** (o añadir una frase explicando el margen entre ambas). Hay que elegir antes de la defensa porque el tribunal puede cruzar las dos cifras.

### C2 — Banda de CO₂: 2,01 µm vs 2,06 µm (inconsistente en toda la memoria)
- **Dónde:** la tabla de bandas `4.Payload/00Payload.tex:108` y el texto `4.Payload/00Payload.tex:92` dicen **2,06 µm (2060 nm)**. En cambio, el cálculo de MTF (`4.Payload/00Payload.tex:228,247,...`, mapas "Lambda2"), el Resumen Ejecutivo (`1.Introduccion/Introduccion.tex:42,46`) y Conceptos Previos (`3.Conceptos_Previos/00Conceptos_Previos.tex:282`) usan **2,01 µm**.
- **CONFIRMADO en el código:** `MATLAB_Code/master.m:28` → `lambda = [1.61e-6, 2.01e-6, 0.76e-6]`. **Toda la computación de MTF/SNR/diseño se hizo a 2,01 µm.** Por tanto el "2,06" de la tabla de bandas (y de `4.Payload:92`) es el valor atípico, no al revés.
- **Problema adicional:** en `3.Conceptos_Previos:282` se atribuye 2,01 µm a OCO-2 y GOSAT, que en realidad usan **2,06 µm** (así lo dice tu propio Cap. 2, líneas 20 y 67).
- **Decisión necesaria (ver pregunta al final):**
  - **Opción A (barata):** unificar todo a **2,01 µm** (coincide con el código, sin recomputar; justificar citando MicroCarb, que usa ~2,01 µm). Solo hay que tocar la tabla de bandas y `4.Payload:92`, y corregir la atribución a OCO-2/GOSAT en Cap. 3.
  - **Opción B (rigurosa):** unificar a **2,06 µm** (coincide con la física de OCO-2/GOSAT/TanSat) y **recomputar** la banda restrictiva de MTF a 2,06 µm (la MTF empeora a mayor λ, así que cambian los mapas y posiblemente algún margen).

### C3 — Fórmula de irradiancia mal transcrita ✅ RESUELTO (solo errata de fórmula)
- **Dónde:** `3.Conceptos_Previos/00Conceptos_Previos.tex:159`: `I = π·τ·Δλ·L_ref / (1 + 4·F#)`.
- **Verificado en el código:** `MATLAB_Code/SNRCalc.m:41`, `BSignalNoiseRatio.m:49`, `SNRfunction.m:39` usan `(1 + 4 * (focal/sqrt(D²(1-r_obs²)))²)`, es decir **(1 + 4·T#²)** con el término **al cuadrado** y usando la relación focal equivalente T# (no F#).
- **Conclusión:** el **cálculo es correcto**; los SNR son válidos. Lo único erróneo es la fórmula *escrita* en la memoria, a la que le falta el cuadrado y el uso de T#.
- **Corrección (segura, aplicable):** escribir `I = π·τ·Δλ·L_ref / (1 + 4·T#²)` (coherente con la ec. \ref{relfoc} de T# que ya está definida en el Cap. 3).

---

## 🟠 Importantes

### I1 — Azimut de lanzamiento 198°
- **Dónde:** `6.Lanzadores/00Lanzador.tex:106-110` y Resumen Ejecutivo (`1.Introduccion:106`).
- **Problema:** con `cos(i)=sin(Az)·cos(lat)`, i=97,48°, lat=34,7°, sale **Az≈189°** (rama retrógrada), no 198° (recalculado: arcsin → −9,1° → 189,1°).
- **Sugerencia:** si aplicaste corrección por velocidad de rotación terrestre (azimut real ≠ inercial), indícalo en el texto; si no, revisar el cálculo. En cualquier caso 189–198° cae dentro del rango de Vandenberg (147°–201°), así que la conclusión (viable) se mantiene, pero el número debe cuadrar con la fórmula mostrada.

### I2 — Periodo orbital 97 min
- **Dónde:** Resumen Ejecutivo `1.Introduccion:90` dice **97 min**. La simulación de fase (`6.Lanzadores:137`) da **94,88 min** y el valor teórico a 520 km es **95,0 min**.
- **Sugerencia:** unificar a ~95 min.

### I3 — "Datos generados 6,9 GB/semana"
- **Dónde:** Resumen Ejecutivo `1.Introduccion:110`.
- **Problema:** 6,909 GB es el volumen **por cobertura** (`7.Segmento_Tierra:164`). El volumen **semanal** real es **8,843 GB** (`7.Segmento_Tierra:178`). La etiqueta "/semana" es incorrecta.
- **Sugerencia:** poner "6,9 GB/cobertura" o "8,8 GB/semana".

### I4 — Telemetría: "Banda S" vs UHF
- **Dónde:** Resumen Ejecutivo `1.Introduccion:109` ("Banda S (telemetría)") y diapositiva de subsistemas.
- **Problema:** el componente elegido es el **EnduroSat UHF II a 437 MHz** (`5.Mission:300`, `7.Segmento_Tierra:29`). El TM/TC es UHF, no banda S.
- **Sugerencia:** corregir el Resumen a "Banda UHF (telemetría)". (Ver también nota para la defensa: la diapositiva 29 lista masas de comunicaciones distintas — banda S, 40 g/17 g — que no coinciden con la memoria.)

### I5 — Falcon 9 "segunda etapa criogénica"
- **Dónde:** `6.Lanzadores:12` y `6.Lanzadores:95`.
- **Problema:** la 2.ª etapa del Falcon 9 es **LOX/RP-1 (queroseno)**, no criogénica en el sentido habitual (que implica hidrógeno líquido, como Centaur o la etapa superior de Ariane). El LOX es criogénico, pero llamar "etapa criogénica" al conjunto es impreciso.
- **Sugerencia:** "etapa superior LOX/queroseno con capacidad de reignición".

### I6 — OCO-2 "fotodiodos de silicio HAWAII-1RG"
- **Dónde:** `2.Misiones_Semejantes:20` y tabla `2.Misiones_Semejantes:122`.
- **Problema:** los arrays **HAWAII-1RG (Teledyne) son HgCdTe**, no silicio; OCO-2 usa silicio solo en la banda O₂ (0,76 µm) y HgCdTe en las bandas SWIR. Describir todos como "fotodiodos de silicio HAWAII-1RG" es incorrecto.
- **Sugerencia:** "matrices HgCdTe (HAWAII-1RG) refrigeradas, con silicio en la banda O₂".

### I7 — 14–16 µm presentada como "ventana atmosférica"
- **Dónde:** `3.Conceptos_Previos:271` lista 14–16 µm como ventana; pero `3.Conceptos_Previos:285` dice (correctamente) que a 15 µm hay **fuerte absorción de CO₂**.
- **Problema:** contradicción interna. 14–16 µm **no** es ventana atmosférica (es banda de absorción del CO₂). Las ventanas reales son ~8–12 µm.
- **Sugerencia:** eliminar 14–16 µm de la lista de ventanas (o sustituir por la ventana de microondas/otra), dejando 8–12 µm.

### I8 — Tabla comparativa Cap. 2 incoherente con el texto
- **Dónde:** tabla `2.Misiones_Semejantes:118-136` vs texto de cada misión.
- **Problemas:** GOSAT-2 masa **2900 kg** (tabla) vs **1700–2000 kg** (texto, línea 35; real ≈1800 kg); altura **666 km** (tabla) vs **613 km** (texto); revisita **3 días** (tabla) vs **6 días** (texto); detector **InSb/MCT** (tabla) vs **Si/InGaAs/MCT** (texto). También CO2M altura 735 (tabla) sin valor en texto.
- **Sugerencia:** revisar fila por fila y unificar tabla↔texto (priorizar los valores del texto, que citan fuente).

### I9 — Masa de combustible 0,361 vs 0,334 kg
- **Dónde:** tabla subsistemas y Resumen usan **0,361 kg**; el dimensionado del depósito (`5.Mission:356`) usa **0,334 kg**.
- **Sugerencia:** unificar (decidir si 0,334 es pre-margen y 0,361 con margen del 10 %, o si es un valor de otra altura), y usar la misma en todo.

### I10 — Error aritmético en la suma de masa de propulsión
- **Dónde:** `5.Mission:362-364`: `0,014 + 0,040 = 0,0567 kg`.
- **Problema:** 0,014+0,040 = 0,054, no 0,0567. El primer sumando debería ser el depósito **0,0167** (calculado dos líneas antes), y entonces 0,0167+0,040 = 0,0567 ✓.
- **Sugerencia:** cambiar `0,014` → `0,0167`.

### I11 — Costes de lanzadores incoherentes
- **Dónde:** `6.Lanzadores`, tabla `:47-52` vs texto.
- **Problemas:** Vulcan **4.044 USD/kg** (tabla) vs texto "110 M$ / 14.400 kg" ≈ **7.640 USD/kg**; Ariane 6 **5.500 USD/kg** (tabla) vs **11.164 USD/kg** (texto línea 28). Además, "bajar el coste de lanzamiento en torno a los **26.000 USD**" (`:95`) no cuadra con 6.500 USD/kg × 3,42 kg ≈ **22.000 USD** (ni con el coste/kg de Falcon 9; 26.984 es el coste/kg de Firefly en la tabla — posible confusión).
- **Sugerencia:** revisar las cifras y aclarar si "26.000 USD" es coste total de lanzamiento de la constelación o un error.

---

## 🟡 Menores

| # | Dónde | Observación | Sugerencia |
|---|-------|-------------|------------|
| M1 | `3.Conceptos_Previos:191` vs `:716` | Radio terrestre **6371 km** (señal) vs **6378,1 km** (SSO) | Unificar (o aclarar medio vs ecuatorial) |
| M2 | `3.Conceptos_Previos:174,182` | `t_int = GSD/v_orb` usa la velocidad **orbital**, no la de la traza en tierra (~7 % menor a 520 km) | Es conservador; opcional notar la aproximación |
| M3 | `4.Payload:331` | "la longitud de onda es directamente proporcional a la irradiancia" | Es proporcional al **nº de electrones/fotones** (E=hc/λ), no a la irradiancia; reformular |
| M4 | `5.Mission:174` | Solapamiento "10 %" con swath 183 vs 170 km → en realidad ~7,6 % | Revisar el porcentaje |
| M5 | `7.Segmento_Tierra:32` vs `:25` | Latitud Fairbanks **64,84°** vs **64,9°** | Unificar |
| M6 | `3.Conceptos_Previos:159` vs `:167` | Símbolo irradiancia **I** vs **E** | Usar el mismo símbolo |
| M7 | `1.Introduccion:51` | "Field of View 20°" se presenta como FoV del telescopio, pero es el **FoV combinado** de 2 telescopios (cada uno ≈10° a 520 km) | Aclarar "FoV combinado" |
| M8 | `5.Mission:341` | El epígrafe **"Dimensionado de la planta de potencia"** trata en realidad de **propulsión**, justo tras la sección eléctrica → confunde | Renombrar a "Propulsión" / "Sistema propulsivo" |
| M9 | `4.Payload:348-353` | Las gráficas SNR son de Detector 4/5/6 (variantes banda O₂) pero los pies dicen "Detector 1/2/3" | Aclarar la correspondencia |
| M10 | `4.Payload:27` | GSD Houston 8,516 (14251/1673 = 8,518) | Trivial, redondeo |
| M11 | `6.Lanzadores:20,49` | Vulcan "fiabilidad 100 %" sin historial operativo | Matizar ("proyectada, sin historial") |

---

## ✔ Comprobaciones que dan correctas (para tu tranquilidad)
- Inclinación SSO **97,48°** a 520 km ✓ (fórmula y valor verificados).
- Volumen por cobertura **6,909 GB** y semanal **8,843 GB** ✓.
- Tiempo de descarga **786 s** ✓; velocidad real 137×0,875×0,75 = **90 Mbps** ✓.
- Potencia específica panel 217,8 W/m², generada **12,47 W** ✓; energía batería 2,60 Wh ✓.
- F# = 120/28 = **4,29 (≈F/4,2)** ✓; θ_pix = 1,54·10⁻⁴ rad ✓; error geoloc. ≈18 km ✓.
- Maniobra de fase: 180° en 36,98 h ✓ (coherente con periodos 96,95/94,88 min).
- Revisita física 6,56×5/6 = **5,47 días** ✓ (dirección correcta: física < efectiva).
- Suma de la tabla de subsistemas = **1,455 kg** ✓ (consistente consigo misma, ver C1).
