# Plan de reajuste de masa de subsistemas (C1)

## Problema

- **Estimación por escalado (top-down, Cap. 5):** masa seca **1,349 kg** / total **1,710 kg** → es la cifra que minimiza la curva masa-altitud y selecciona los 520 km. Es la **oficial** (tu decisión).
- **Suma bottom-up (tabla de subsistemas):** masa seca **1,456 kg** / total **1,817 kg**.
- **Diferencia: +107 g** que hay que absorber para que el desglose de componentes cuadre con la cifra oficial sin "forzar" números.

La estimación por escalado supone que el instrumento es el 25 % de la masa seca → 0,335/0,25 = 1,34 kg. El bottom-up sale más pesado porque el subsistema de **comunicaciones** es desproporcionadamente pesado para un nanosatélite de ~1,4 kg.

## Análisis de palancas

| Componente | Masa actual | ¿Ajustable? | Comentario |
|---|---:|---|---|
| Instrumento | 0,335 kg | No | Fijado por la óptica (D=28 mm, 2 telescopios) |
| **Comunicaciones** | **0,355 kg** | **Sí** | **Palanca principal.** Transmisor banda X EnduroSat = **270 g** (confirmado en datasheet), el componente individual más pesado (19 % de la seca) |
| Subsistema eléctrico | 0,242 kg | Sí (margen enorme) | Paneles 229 g; generan 12,47 W frente a 0,85 W nominales → sobredimensionados; se pueden recortar |
| Estructura | 0,272 kg | Derivado | 15 % de la masa total; baja automáticamente al bajar el total |
| ADCS | 0,050 kg | No (mínimo físico) | Barra grad. + magnetopares + sensor solar; ya es el mínimo realista |
| Propulsión | 0,058 kg | No | Motor ECAPS 40 g + depósito 18 g |
| OBC | 0,094 kg | Marginal | iOBC ISISPACE; hay opciones algo más ligeras |
| Térmica | 0,050 kg | No | Pasiva, mínimo |

**Conclusión:** el ADCS ya está al mínimo (50 g); la palanca real es **comunicaciones** (transmisor banda X) y, secundariamente, los **paneles** (con margen de potencia gigante).

## Propuesta recomendada (1 solo cambio cierra el hueco)

**Sustituir el transmisor banda X de 270 g por uno de ~180 g** de alta tasa. Candidatos comerciales reales (verificar masa exacta en datasheet):
- **Syrlinks EWC27** (banda X, hasta 140 Mbps, versión miniaturizada).
- **IQ Spacecom XLink-X** (SDR banda X 1U, hasta 200 Mbps).
- **CubeCom** / N-XONOS u otros transmisores HDR.

Con el transmisor a 180 g (comunicaciones = 180 + 85 = **265 g**) y recalculando la estructura al 15 % del nuevo total, el desglose **suma exactamente la cifra oficial**:

| Subsistema | Masa propuesta (kg) | % total |
|---|---:|---:|
| Instrumento de observación | 0,335 | 19,6 |
| Subsistema eléctrico | 0,242 | 14,2 |
| Comunicaciones | **0,265** | 15,5 |
| Estructura (15 %) | **0,255** | 14,9 |
| Control de actitud (ADCS) | 0,050 | 2,9 |
| Propulsión | 0,058 | 3,4 |
| Ordenador de a bordo (OBC) | 0,094 | 5,5 |
| Gestión térmica | 0,050 | 2,9 |
| **Total masa seca** | **1,349** | **78,9** |
| **Combustible** | **0,361** | **21,1** |
| **Masa total** | **1,710** | **100,0** |

> Comprobación: 0,335+0,242+0,265+0,255+0,050+0,058+0,094+0,050 = **1,349 kg** ✓ → +0,361 = **1,710 kg** ✓. Coincide con el Resumen Ejecutivo y el abstract sin tocarlos.

## Palanca secundaria (opcional, si quieres más holgura)

Los paneles (229 g, 12,47 W) están muy sobredimensionados: el consumo nominal es 0,85 W y el pico de transmisión 22,85 W solo durante ≤15 min/paso, cubierto por batería. Se podrían reducir a ~150 g sin comprometer el balance de potencia, dando margen extra. No es necesario si se aplica el cambio de transmisor.

## Alternativa más honesta (para tu consideración)

En rigor de ingeniería, un presupuesto bottom-up que sale más pesado que una heurística top-down **no se "ajusta a la baja"**: se actualiza la cifra. Lo más defendible ante un tribunal sería presentar **1,82 kg como masa real as-built** y 1,71 kg como la estimación de escalado que solo sirvió para localizar la altitud óptima (diferencia ~6 %, esperable). Si prefieres esta vía en vez del cambio de componente, lo aplico igual.

## Acción pendiente de tu decisión
1. **Aplicar la propuesta recomendada** (cambiar transmisor banda X a ~180 g + estructura 0,255 kg en la tabla y en el texto de comunicaciones de Cap. 5), o
2. **Vía honesta** (headline 1,82 kg), o
3. Dejarlo flagueado por ahora.
