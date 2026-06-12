# Informe de cambios — para sincronizar la defensa (Defensa.pptx)

> Resumen de TODO lo que cambió en la memoria durante esta sesión.
> **NOTA:** las diapositivas 4, 23, 29, 40 y 42 de `Defensa.pptx` ya se han actualizado automáticamente (backup en `Defensa_backup_20260612.pptx`). Esta sección queda como registro de qué se cambió.

---

## 1. Diapositivas que DEBES actualizar

### ⚠️ Slide 29 (Subsistemas. Componentes) — cambios grandes
| Antes (slide) | Ahora (memoria) |
|---|---|
| Banda X "Transmitter" 0,04 kg | **IQ Spacecom XLink-X, banda X, 200 g, 16 W, 200 Mbps** |
| "Transceiver Banda S" 0,017 kg | **Transceptor UHF EnduroSat, 437 MHz, 85 g** |
| Masa de combustible 0,334 kg | **0,361 kg** |
| (Comms total implícito) | **Comunicaciones = 285 g** |
| Estructura 0,272 kg | **0,258 kg** (15 %) |
| — | Paneles 212 g, baterías 7 g (eléctrico 219 g) |
| Masa seca / total | **1,349 / 1,710 kg** (igual titular, ahora el desglose cuadra) |

> Motivo del cambio de transmisor: el EnduroSat banda X pesaba 270 g (19 % de la masa seca) y hacía que la suma de subsistemas (1,82 kg) no cuadrara con el óptimo por escalado (1,71 kg). El XLink-X (200 g, 16 W, 200 Mbps) es más ligero, más rápido y consume menos → la batería y los paneles bajan, y el desglose **suma exactamente 1,710 kg**.

### ⚠️ Slide 23 (Masa de combustible)
- Isp **220 s → 209 s** (la memoria usa 209 s, el del propulsor ECAPS).

### ⚠️ Slide 40 (Velocidad de descarga) — tabla meteorológica desactualizada
- Sustituir por la tabla actual de la memoria: favorables 130–150, lluvia moderada 120–140, lluvia intensa 90–120, nubes densas 125–145 Mbps; **media anual 137 Mbps**. (El "Tormenta 22–38" antiguo ya no aplica.)

### ⚠️ Slide 4 (Misiones semejantes)
- GOSAT: masa **2000 → ~1800 kg**; (en la tabla de la memoria también altura 613 km y revisita 6 días).

### ⚠️ Slide 39 / cualquier "6,9 GB/semana"
- 6,9 GB es **por cobertura**; el volumen **semanal** es **8,8 GB**.

### Azimut — corregido a 189°
- Era un **error de transcripción** (198 ↔ 189). Corregido en memoria y en `Azimuth.m`. **Regenera la figura** `azimuth.jpg` corriendo `Azimuth.m`.

---

## 2. Cambios conceptuales aplicados en la memoria (coherencia)
- **Banda CO₂** unificada a **2,01 µm** en todo el documento (era 2,06 en la tabla de bandas; el código calculó a 2,01).
- **Telemetría = UHF** (437 MHz), no banda S, en el Resumen Ejecutivo.
- **Periodo orbital 97 → 95 min**.
- **Falcon 9**: "segunda etapa criogénica" → "etapa superior LOX/queroseno".
- **OCO-2**: "fotodiodos de silicio H1RG" → **HgCdTe** (con Si en la banda O₂).
- **Ventana atmosférica 14–16 µm** eliminada (es absorción de CO₂, no ventana).
- **Fórmula de irradiancia** corregida a `(1 + 4·T#²)` (el código ya lo hacía bien; era errata de transcripción).
- **Suma de propulsión** corregida (0,0181 + 0,040 = 0,0581) y combustible unificado a **0,361 kg**.
- **Coste lanzamiento** "26.000 USD" confuso → **6.500 USD/kg**.
- Menores: radio terrestre etiquetado como ecuatorial, solapamiento swath ~8 %, latitud Fairbanks 64,84°, etiquetas de detectores en gráficas SNR, etc.

## 3. Eliminado (a petición tuya)
- **CSP** (CubeSat Space Protocol): quitado de Cap. 5 y Cap. 7 y de la bibliografía.
- **Link budget / Doppler / polarización**: ya estaban fuera del texto (commit previo); limpié 2 referencias colgantes a "22,8 dB" y borré los scripts huérfanos `calcularDoppler.m` y `calcularLinkBudget.m`.
- **No había nada de IoT ni MQTT**.

## 4. Git
- Borradas en **local** las 4 ramas fusionadas (`latex-fixes`, `latex-frontmatter`, `matlab-refactor`, `python-port`).
- **Pendiente que ejecutes tú** (falló por auth):
  ```
  git push origin --delete feature/latex-fixes feature/latex-frontmatter feature/matlab-refactor feature/python-port
  ```
- `feature/beamer-defence` (notas de presentador) se conserva — útil para la defensa.

## 5. Pendiente de tu decisión / verificación
- **Azimut**: corregido a 189°; solo queda **regenerar `azimuth.jpg`** corriendo `Azimuth.m`.
- **Costes de lanzadores** (tabla Cap. 6): Vulcan 4.044 vs ~7.640 USD/kg; Ariane 5.500 vs 11.164 USD/kg → revisa con datasheet.
- **Banda 2,01 vs 2,06 µm**: quedó en 2,01 (coincide con el código). Si el tutor prefiere 2,06, hay que **recomputar la MTF**.
- **XLink-X**: verifica en datasheet la masa/consumo exactos antes de imprimir (usé 200 g / 16 W / 200 Mbps de CubeSatShop).
- **Barrido completo de tildes**: aquí no hay corrector español; pasa el spellchecker de Overleaf/TeXstudio para cazar tildes sueltas restantes (mas→más, orbita→órbita, calculo→cálculo, linea→línea, etc.).
- **Recompila** la memoria en tu entorno (aquí no hay LaTeX): la verificación estática no encontró citas ni labels colgantes, pero conviene un build limpio.

## 6. Entregables de esta revisión (carpeta `Revision_TFG/`)
- `01_Revision_Conceptual.md` — hallazgos conceptuales (críticos/importantes/menores) y estado.
- `02_Revision_Escritura.md` — propuestas de reescritura (stop-slop), **sin aplicar**, para que elijas.
- `03_Guia_Defensa.md` — guion por diapositiva + banco de Q&A del tribunal.
- `04_Plan_Masa_Subsistemas.md` — justificación del reajuste de masa.
- `05_Cambios_para_Defensa.md` — este documento.
