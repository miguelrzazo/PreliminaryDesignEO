# Guía rápida de defensa

## Chuleta de una hoja

### Qué problema resuelvo
- Diseñé una misión preliminar para medir CO2 antropogénico sobre Estados Unidos continental.
- El cliente pide cobertura semanal, MTF mayor o igual que 0,25 y SNR mayor o igual que 400.
- La variable que más aprieta el diseño no es la SNR. Es la cobertura.

### Solución final
- Constelación: 2 satélites.
- Órbita: SSO dawn-dusk.
- Altura: 520 km.
- LTAN: 06:00.
- Inclinación: 97,48 grados.
- GSD: 80 m.
- Bandas: 0,76 / 1,61 / 2,01 micras.
- Ancho de banda: 20 nm.
- Detector: Teledyne H2RG.
- Telescopio: refractivo.
- Pupila: 28 mm.
- Masa por satélite: 1,710 kg.
- Masa total en órbita: 3,42 kg.

### Prestaciones que tengo que decir de memoria
- MTF de diseño: 0,25 en 2,01 micras.
- SNR de diseño: 717 en 0,76 micras.
- Revisita física: 5,47 días.
- Revisita efectiva con nubes: 6,6 días.
- Swath requerido: 170 km.
- Swath logrado: 183 km.
- Deorbitación natural: 257 días.

### Segmento tierra y datos
- Estación elegida: Fairbanks.
- Tiempo de contacto semanal: 79 620 s, unas 22,1 h.
- Velocidad teórica en banda X: 150 Mbps.
- Velocidad media anual: 137 Mbps.
- Velocidad final con pérdidas: 90 Mbps.
- Datos por cobertura: 6,9 GB.
- Datos por semana: 8,8 GB.
- Tiempo necesario de descarga: 786 s.

### Decisiones que tengo que defender
- 80 m porque separa zona industrial y residencial sin disparar masa.
- Dos satélites porque con uno no llego a la cobertura semanal.
- 520 km porque minimiza la masa total.
- H2RG porque cubre las tres bandas y da el menor diámetro de pupila.
- Refractivo porque evita obscuración central y reduce el diámetro necesario.
- Fairbanks porque maximiza contactos en órbita polar y ya tiene banda X.
- Falcon 9 porque combina fiabilidad, rideshare y compatibilidad con SSO.

### Frases de cierre
- Este trabajo no entrega un satélite listo para fabricar. Entrega una solución preliminar coherente y defendible.
- La aportación fuerte no es una cifra aislada. Es el método para cruzar carga útil, órbita, cobertura y masa.
- El resultado demuestra viabilidad a nivel preliminar. El siguiente paso sería subir detalle en propagación, subsistemas y costes.

## Guion por diapositiva

### 1. Portada
Presento el diseño preliminar de una misión de observación de la Tierra para monitorizar CO2 antropogénico sobre Estados Unidos continental. El objetivo del trabajo no es cerrar un diseño de detalle, sino demostrar viabilidad técnica con una solución ligera.

### 2. Objetivos del proyecto
El cliente pide tres cosas al mismo tiempo: calidad óptica, sensibilidad radiométrica y cobertura. En números, eso significa MTF mayor o igual que 0,25, SNR mayor o igual que 400 y cobertura semanal de Estados Unidos continental.

### 3. Misiones semejantes
Antes de diseñar, revisé qué hacen las misiones existentes. Me interesaba ver qué compromisos aceptan entre masa, resolución y cobertura.

### 4. Comparativa de misiones
La conclusión es que no existe una solución única. Hay misiones grandes con mucha capacidad y misiones pequeñas con mucha agilidad, pero ninguna resuelve exactamente este problema con esta combinación de resolución y cobertura. Eso justifica hacer un proceso iterativo de diseño.

### 5. Conceptos previos
Aquí solo fijo la base teórica. No me paro mucho porque lo importante viene después, cuando esos conceptos entran en el espacio de diseño real.

### 6. Carga de pago
Empiezo por la carga de pago porque condiciona casi todo lo demás. Si aquí cambias GSD, bandas o detector, cambias también cobertura, masa y subsistemas.

### 7. GSD
Probé varios GSD sobre una imagen de Houston. El punto de diseño es 80 m porque sigue permitiendo distinguir zonas industriales de residenciales y, al mismo tiempo, evita penalizar el swath y la masa.

### 8. Bandas espectrales y ancho de banda
Trabajo con tres bandas. La de 0,76 micras sirve para calibrar la columna de aire con O2, y las de 1,61 y 2,01 micras atacan la absorción de CO2 en SWIR. Elegí 20 nm porque da un equilibrio razonable entre selectividad espectral y SNR.

### 9. Detectores y telescopios
Comparé varios detectores y varias arquitecturas ópticas. No elegí una solución por intuición, sino por el cruce posterior entre MTF, SNR, revisita y masa.

### 10. Filtros y escaneo
Elegí pushbroom porque evita mecanismos móviles y mejora el tiempo de integración. Elegí filtros microstrip porque reducen masa y riesgo mecánico, y permiten adquirir bandas de forma simultánea.

### 11. Tipo de órbita y LTAN
La misión pide observación con iluminación estable, así que la opción natural es una órbita heliosíncrona. La seleccioné en configuración dawn-dusk, LTAN 06:00, para mantener condiciones de iluminación consistentes y reducir eclipses.

### 12. Flujo de trabajo de scripts
Desarrollé una cadena de scripts en MATLAB para recorrer el espacio de diseño. Ese bloque calcula MTF, SNR, cobertura, masa seca, combustible y masa total, y luego cruza todo para encontrar soluciones viables.

### 13. Cálculo de MTF
Aquí calculo la MTF total del sistema. La banda más restrictiva es la de 2,01 micras porque la difracción penaliza más al crecer la longitud de onda.

### 14. Gráficas de MTF
Estas gráficas muestran, para cada combinación, qué pares altura-diámetro cumplen el requisito. La idea clave es que no todas las arquitecturas pagan el mismo coste en pupila.

### 15. Cálculo de SNR
La SNR la calculo con señal y ruido a nivel de detector. En este caso el requisito de SNR no es el cuello de botella del diseño. Se cumple con margen.

### 16. Gráficas de SNR
La lectura buena aquí es que la SNR no me gobierna la solución final. Me sirve para descartar configuraciones pobres, pero el filtro fuerte aparece cuando meto cobertura.

### 17. Tiempo de revisita
La cobertura semanal es el requisito más duro. Por eso necesitaba un modelo específico de revisita y no solo una intuición geométrica.

### 18. Gráficas de revisita
Estas curvas dejan claro que con un solo satélite no llego. La constelación de dos satélites aparece porque la cobertura obliga, no porque yo la quisiera desde el principio.

### 19. Cruce de los datos
Aquí junto MTF, SNR y cobertura. Para cada altura busco el diámetro mínimo que cumple todo a la vez.

### 20. Comparativa de detectores con refractivo
Con telescopio refractivo, el H2RG sale mejor parado que el resto. Me permite cumplir antes y con menos pupila.

### 21. Optimización de la masa
Una vez tengo la relación entre altura y diámetro mínimo, traduzco eso a masa. Ese paso convierte una solución óptica viable en una solución de misión defendible.

### 22. Masa seca
La masa seca la estimo por escalado con dos referencias de misión. No vendo esto como un cálculo de detalle. Lo uso como una aproximación preliminar consistente.

### 23. Masa de combustible
Modelo el arrastre atmosférico durante ocho años y aplico impulsos de reposición cuando la altura cae un 2 por ciento. A partir del delta-v acumulado calculo la masa de combustible con Tsiolkovski.

### 24. Masa total
El mínimo de masa total aparece a 520 km. Ese es el punto donde el compromiso entre arrastre y apertura óptica sale mejor.

### 25. Disposición final del instrumento
Aquí cierro la disposición física del instrumento. Dos telescopios me permiten alcanzar el swath necesario sin disparar el resto del sistema.

### 26. Punto de diseño
Esta es la slide que resume el proyecto. En el punto final cumplo MTF, SNR y cobertura con una constelación de dos satélites de 1,71 kg cada uno.

### 27. Deorbitación natural
No necesito reservar combustible adicional para fin de vida porque el satélite reentra de forma natural en 257 días. Eso me deja dentro del criterio de sostenibilidad que fijé.

### 28. Subsistemas
Con la carga útil ya cerrada, paso a validar si el resto del satélite cabe en masa y volumen. La idea aquí es demostrar factibilidad, no cerrar ingeniería de detalle.

### 29. Componentes
Seleccioné componentes comerciales para comprobar que la solución se puede materializar. El mensaje importante es que la suma cuadra con la masa estimada y no aparece ninguna incoherencia gruesa.

### 30. Lanzadores
Con una masa total tan baja, el factor dominante no es capacidad de carga. Es acceso fiable y barato a una SSO.

### 31. Lanzadores evaluados
Comparé varias opciones de rideshare. No basta con mirar coste por kilo. También miré fiabilidad y compatibilidad operativa.

### 32. Bases de lanzamiento
La base importa porque necesito un azimut compatible con SSO. No todas las bases sirven aunque el lanzador sí sirva.

### 33. Selección final del lanzador
El Falcon 9 sale mejor por fiabilidad, precio, experiencia en rideshare y compatibilidad con la órbita objetivo. Para esta misión, es la opción más limpia.

### 34. Selección de la base
Vandenberg encaja porque permite el azimut requerido y ya opera este tipo de misiones. No necesito forzar una solución rara.

### 35. Puesta en órbita
Los dos satélites se despliegan con una maniobra de fase hasta quedar separados 180 grados. Eso me da la geometría que necesito para la cobertura.

### 36. Segmento tierra
No basta con captar datos. También tengo que demostrar que puedo descargarlos.

### 37. Estaciones evaluadas
Comparé varias opciones y me quedé con las que tenían sentido por latitud e infraestructura.

### 38. Selección de Fairbanks
Fairbanks gana porque multiplica el número de contactos en órbita polar y ya cuenta con infraestructura de banda X. La simulación me da 22,1 horas de contacto a la semana.

### 39. Generación de datos
La constelación genera 6,9 GB por cobertura y unos 8,8 GB por semana. La cifra importante es la semanal, porque es la que después comparo con la capacidad de descarga.

### 40. Velocidad de descarga
Parto de 150 Mbps teóricos en banda X y los bajo con un modelo meteorológico anual. La media queda en 137 Mbps antes de pérdidas adicionales.

### 41. Descarga de datos
Aplicando pérdidas de hardware y baja elevación, la velocidad final útil queda en 90 Mbps. Con ese valor, descargar la producción semanal lleva solo 786 segundos, así que hay margen de sobra.

### 42. Conclusiones
La solución final cumple los requisitos del cliente con una constelación muy ligera. El resultado fuerte del trabajo no es solo el número final, sino el método para cruzar óptica, órbita, cobertura, masa y operaciones. El paso siguiente sería profundizar en propagación, nubes, potencia, térmico, ADCS y coste.

### 43. Cierre
Gracias por la atención. El código está disponible en GitHub y ahora quedo a disposición del tribunal.

## Preguntas probables y respuestas

### Objetivo y alcance

**¿Qué aporta este TFG si no es un diseño final?**  
Aporta una solución preliminar coherente y un método reproducible para llegar a ella. En esta fase, eso es más valioso que cerrar detalle falso.

**¿Por qué centrarte en Estados Unidos continental?**  
Porque el cliente fija esa región como zona de interés y porque concentra emisiones antropogénicas relevantes. Eso define cobertura, órbita y segmento tierra.

**¿Por qué cobertura semanal?**  
Porque es el requisito operativo del cliente. Además, obliga a introducir constelación y hace el problema interesante.

**¿Qué partes consideras dentro y fuera de alcance?**  
Dentro: carga útil, órbita, cobertura, masa, lanzador y descarga. Fuera: diseño estructural detallado, térmico detallado, link budget completo y coste de programa cerrado.

### Carga útil

**¿Por qué 80 m y no 50 m?**  
Porque 50 m mejora el detalle, pero castiga swath, revisita y masa. A 80 m sigo distinguiendo fuentes industriales con una solución mucho más ligera.

**¿Por qué 80 m y no 100 m?**  
Porque a 100 m empiezo a perder capacidad para separar zonas industriales y residenciales. 80 m es el punto donde el compromiso sigue siendo útil.

**¿Por qué tres bandas y no una sola de CO2?**  
Porque necesito calibrar la columna de aire con O2 y además comparar dos bandas de CO2 para robustez y recuperación de concentración.

**¿Por qué 2,01 micras?**  
Porque es una banda útil de absorción de CO2 en SWIR y es la que usé de forma coherente en el cálculo de MTF y en la memoria final.

**¿Por qué 20 nm?**  
Porque un ancho más estrecho mejora selectividad, pero castiga señal. Un ancho mayor mete radiación no deseada. Con 20 nm consigo un equilibrio razonable.

**¿Por qué H2RG?**  
Porque cubre visible y SWIR con buen rendimiento y, en el cruce de datos, minimiza el diámetro de pupila necesario.

**¿Por qué telescopio refractivo?**  
Porque evita obscuración central y mejora el comportamiento de MTF para este caso. El precio es la aberración cromática, pero aquí las bandas son estrechas y separadas.

**¿Por qué pushbroom?**  
Porque simplifica mecánica, mejora SNR y es una arquitectura natural para esta misión.

**¿Por qué filtros microstrip?**  
Porque evitan mecanismos móviles, reducen masa y permiten adquirir bandas a la vez.

**¿Por qué dos telescopios?**  
Porque con uno no llego al swath necesario de forma limpia. Dos telescopios cierran el requisito de cobertura.

### Prestaciones

**¿Qué requisito manda de verdad?**  
La cobertura. La SNR se cumple con margen y la MTF condiciona el diámetro, pero la revisita es lo que fuerza constelación y swath.

**¿Por qué MTF 0,25?**  
Porque es el requisito del cliente. No lo inventé yo.

**¿Por qué la banda restrictiva para MTF es la de 2,01 micras?**  
Porque al aumentar la longitud de onda baja la frecuencia de corte por difracción. Esa banda castiga más la transmisión de detalle espacial.

**¿Por qué la SNR no gobierna el diseño final?**  
Porque sale holgada en las configuraciones viables. El filtro fuerte llega cuando cruzo cobertura y masa.

### Órbita y cobertura

**¿Por qué SSO?**  
Porque necesito observación con iluminación consistente. Es la elección natural para sensores pasivos de este tipo.

**¿Por qué dawn-dusk?**  
Porque mantiene mejor estabilidad de iluminación y reduce eclipses, lo que también simplifica la parte eléctrica.

**¿Por qué 520 km?**  
Porque ahí aparece el mínimo de masa total. Más abajo crece el coste en combustible. Más arriba crece la pupila óptica.

**¿Por qué dos satélites y no tres?**  
Porque dos ya cumplen el requisito. Tres añadirían coste sin necesidad para este objetivo.

**¿Por qué no vale un satélite?**  
Porque no alcanza la revisita semanal exigida sobre la región de interés.

**¿Qué limita tu modelo de revisita?**  
Es un modelo válido para diseño preliminar, pero no sustituye una propagación numérica completa a largo plazo con más perturbaciones.

**¿Cómo metiste las nubes?**  
Como una corrección simple sobre la revisita efectiva. Es útil para primera aproximación, pero se puede mejorar con un modelo probabilístico espacial y temporal.

**¿Por qué no aplicaste repeating ground track?**  
Porque quise resolver primero la factibilidad general. Sería una mejora clara en trabajo futuro.

### Masa y subsistemas

**¿Cómo sacas la masa seca?**  
Por escalado con dos referencias de misión y con hipótesis explícitas. No es un presupuesto de detalle. Es una estimación preliminar.

**¿Es creíble una masa de 1,71 kg?**  
Sí, porque la apertura óptica es pequeña y el resto del sistema usa componentes comerciales ligeros. Aun así, esa cifra debe validarse con más detalle en fases posteriores.

**¿No falta análisis térmico?**  
Sí. Lo reconozco como limitación. Aquí solo cierro factibilidad global.

**¿No falta análisis de potencia?**  
Sí, aunque sí comprobé que el orden de magnitud encaja con la masa asignada al subsistema eléctrico.

**¿No falta ADCS detallado?**  
Sí. Hice una validación preliminar de requerimiento y arquitectura, no una simulación fina de control.

**¿Cómo justificas el volumen del instrumento?**  
Con la disposición preliminar y el escalado geométrico. En la siguiente fase tocaría pasar a CAD y empaquetado real.

### Lanzamiento y operaciones

**¿Por qué Falcon 9?**  
Porque combina fiabilidad alta, coste bajo en rideshare y experiencia operativa en SSO.

**¿Por qué Vandenberg?**  
Porque permite el azimut requerido para la órbita heliosíncrona y ya tiene experiencia con el lanzador elegido.

**¿Cómo separas 180 grados los satélites?**  
Con una maniobra de fase: subes temporalmente el apogeo, dejas derivar el segundo y recircularizas al alcanzar el desfase.

**¿Cuál es el coste aproximado de lanzamiento?**  
Con el dato de rideshare de Falcon 9, la masa total de la constelación sale muy baja. En esta fase me interesa más la viabilidad que cerrar una cifra contractual.

### Segmento tierra y datos

**¿Por qué Fairbanks y no Wallops?**  
Porque la latitud alta da más pases útiles para una órbita casi polar. Eso mejora mucho el margen de descarga.

**¿Cómo justificas 90 Mbps finales?**  
Parto de la capacidad teórica, aplico un modelo meteorológico anual y luego dos correctores prácticos: hardware y pases de baja elevación.

**¿Qué margen real tienes en descarga?**  
Muchísimo. Necesito 786 segundos por semana y dispongo de 79 620 segundos de contacto.

**¿No te preocupa la memoria a bordo?**  
No en esta configuración. Con 512 MB ya cierro el caso y con 1 GB gano robustez.

### Metodología

**¿Qué validaste y qué estimaste?**  
Validé coherencia física del resultado y la parte de cobertura y contactos con herramientas numéricas. Estimé lo que corresponde a una fase preliminar: masa seca, detalle de subsistemas y parte de operaciones.

**¿Qué software usaste?**  
MATLAB para la cadena de diseño, Aerospace Toolbox para escenarios orbitales y simulación de contactos, y herramientas de apoyo para la memoria y la defensa.

**¿Qué harías si tuvieras tres meses más?**  
Primero mejoraría propagación y cobertura. Después cerraría potencia, térmico, ADCS y coste.

### Críticas difíciles

**Esto parece una misión demasiado pequeña para medir CO2. ¿No es optimista?**  
La masa es muy baja, sí. Por eso presento el resultado como diseño preliminar y no como misión lista para fabricar. El valor del trabajo está en demostrar que el orden de magnitud no es absurdo y en señalar dónde hace falta profundizar.

**Si el modelo de nubes es simple, ¿cómo confías en la cobertura?**  
La cobertura geométrica sí está bien capturada. Lo que simplifico es la degradación meteorológica. Para una fase preliminar me sirve para ver si el concepto vive o muere.

**¿No estás mezclando precisión climática con una arquitectura demasiado pequeña?**  
Lo que demuestro aquí es capacidad de observación preliminar y cierre de requisitos de misión. La precisión final de producto geofísico exigiría una cadena de calibración y validación más profunda.

**¿Tu resultado depende demasiado del escalado de masa?**  
Sí, depende de esa hipótesis, y lo digo de forma explícita. Por eso acompaño el valor final con un desglose de subsistemas y con referencias comerciales.

### IA y declaración responsable

**¿Usaste IA?**  
Sí. La usé como apoyo de redacción, revisión y organización.

**¿Para qué no la usaste?**  
No la usé para sustituir el criterio técnico ni para aceptar números sin comprobación.

**¿Qué validaste tú?**  
Las decisiones de diseño, las cifras finales, las ecuaciones, las comparaciones y la coherencia global del trabajo.

**Si te preguntan por ética o autoría, qué dices?**  
Que la IA fue una herramienta de apoyo, igual que un corrector, un editor o un asistente de búsqueda. La responsabilidad técnica y académica del contenido final es mía.
