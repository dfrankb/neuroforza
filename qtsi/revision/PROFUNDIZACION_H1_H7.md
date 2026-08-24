# Profundización: hallazgos 1 y 7

Continuación de `REVISION_CRITICA.md`. Todo lo afirmado aquí está verificado en
`deep_h1_theta_and_unitarity.py`, `deep_h7_free_sector_and_hydrogen.py` y
`deep_h7_born_rule_trilemma.py`.

**Resumen de lo que cambió respecto de la primera pasada.** El hallazgo 1 resultó
*reparable*, y en el camino apareció un resultado que mejora el paper: la Θ-covarianza
sola, sin ansatz, es más fuerte de lo que el paper reclama. El hallazgo 7 resultó *peor*
de lo que había dicho: no es solo que falte una predicción, es que falta la regla de Born,
y sin ella la teoría no dice qué predice.

---

# Hallazgo 1 — el pivote de §5.3

## 1.1 Un resultado que el paper podría reclamar y no reclama

El paper supone el ansatz `V = I + w σ₁` y deriva "w real" (5.9). **La Θ-covarianza sola,
sin ningún ansatz, da más.**

Con `V = [[a,b],[c,d]]`:

```
Θ V Θ⁻¹ = σ₁ V* σ₁ = [[d*, c*], [b*, a*]]
```

Igualar a `V` da `d = a*` y `c = b*`, es decir

```
V = [[a,  b ],
     [b*, a*]]
```

que es **exactamente la forma U(1,1)**, con `det V` real automáticamente. Verificado sobre
vértices aleatorios de esa forma: `|Θ(V) − V| = 0` y determinante real.

Esto es una mejora directa a §5.3: se puede pasar de (5.7) a (5.10) en un solo paso, sin
postular linealidad en w y sin la expansión a primer orden. El vértice como boost no es una
elección de parametrización, es lo único que la simetría admite.

## 1.2 La incompatibilidad con unitariedad es estructural, no perturbativa

Mi objeción original (`w = 0`) era el caso particular del ansatz lineal. El enunciado
general es más fuerte. Con `V = [[a,b],[b*,a*]]`:

```
V†V = [[|a|²+|b|², 2a*b], [2ab*, |a|²+|b|²]]
```

La unitariedad exige `|a|²+|b|² = 1` **y** `a*b = 0`, o sea `a = 0` o `b = 0`. Verificado
qué son esos dos casos:

| caso | qué es | verificación |
|---|---|---|
| `b = 0` | V diagonal, las ramas nunca se mezclan | acoplamiento entre ramas = 0.0000 |
| `a = 0` | swap puro, `T² ∝ I`, período dos | tr T = 0.0000, nada propaga |

Así que **Θ-covarianza + unitariedad admite exactamente dos vértices, ambos degenerados**:
uno que nunca pliega y otro que pliega con certeza en cada paso y no propaga nada. No hay
teoría en la intersección. La incompatibilidad no es un artefacto de la expansión a primer
orden en w.

## 1.3 La justificación que el paper necesita — y que casi tiene

Mi objeción era que `ΘVΘ⁻¹ = V` equivale a `ΘHΘ⁻¹ = −H`, una *anti*simetría, y no a la
condición T estándar `+H`. **Hay una respuesta buena, y el paper no la da.**

Θ no es inversión temporal. Es una **redundancia de etiquetado**. Tres hechos verificados
lo sostienen:

- `Θ² = I`. Es una involución, como debe serlo un reetiquetado.
- `D` es Θ-invariante para **todo** E (verificado en E = 0.1, 1.0, 7.3, residuo 0). Es decir:
  intercambiar las etiquetas de rama *y* conjugar la convención de fase amarrada a ellas
  devuelve el mismo D. Θ es la afirmación de que el corte en "adelante" y "atrás" es una
  etiqueta, no un dato físico.
- El swap puro σ₁, sin conjugación, **no** es simetría de D: `|σ₁Dσ₁ − D| = 0.59`. O sea el
  K no es opcional — la etiqueta de rama y la convención de fase están amarradas y hay que
  girarlas juntas.

Y la razón de que sea `ΘTΘ⁻¹ = T` y no `T⁻¹`: **Θ no toca τ**. El contador de pasos n crece
a lo largo del camino tanto si éste avanza como si retrocede en t. La condición estándar
`ΘU(t)Θ⁻¹ = U(−t)` aplica a operaciones que invierten el parámetro de evolución; Θ no lo
invierte, intercambia ramas. Es correcto — pero hay que escribirlo, porque es donde vive
todo el contenido de la derivación.

## 1.4 El puente que falta: §4.2 hace trabajo positivo

Aquí está la reparación de fondo. §4.2 está presentada como resultado negativo ("what the
counting does not settle"), pero **hace trabajo positivo que el paper no le atribuye**.

§4.2 establece que el conteo total de cruces no es invariante de la línea de universo. La
consecuencia que el paper no extrae: entonces `|ψ₊|² + |ψ₋|²` no puede ser la densidad de
nada que el plegado conserve, y por tanto **no hay ninguna condición de normalización que
imponer**. El warrant de la unitariedad se retira, no se viola.

El contraste con el caso espacial es exactamente la asimetría que justifica usar dos
criterios distintos, y el paper ya lo tiene escrito en §5.3: a un instante fijo un camino
espacial ocupa una celda en una rama, así que `|ψ₊|²+|ψ₋|²` **es** todo el contenido de
probabilidad de esa celda. Para pliegues en el tiempo eso falla, y falla por §4.2.

**Cadena reparada:**

| paso | contenido | estado en el paper |
|---|---|---|
| (a) | los caminos plegados cruzan un instante un número impar de veces; el conteo total no es invariante | ✅ §4.2 |
| (b) | ⟹ no hay normalización que imponer ⟹ la unitariedad en τ queda sin motivación | ❌ **falta** |
| (c) | Θ = σ₁K es redundancia de etiquetado; D es Θ-invariante idénticamente | ✅ (5.8), sin la interpretación |
| (d) | Θ-covarianza ⟹ V en forma U(1,1) | ⚠️ (5.9), en versión débil con ansatz |
| (e) | normalización det = 1 ⟹ SU(1,1), σ₃ conservado, y todo lo demás | ✅ (5.10)–(5.12) |

Con (b) explícito, §5.3 deja de ser un cambio de postulados sin anunciar y pasa a ser un
argumento cerrado. El paper apunta al puente sin cruzarlo: §5.3 dice solo *"For folds in
time no such quantity is available in advance, and Section 4.2 gives the reason the counting
cannot supply one."* Esa frase señala el argumento; hay que hacerlo.

## 1.5 La objeción que sobrevive: Stueckelberg–Horwitz–Piron

Hay que decirlo porque debilita la reparación anterior. La teoría estándar parametrizada
—Stueckelberg–Horwitz–Piron, y el paper cita a Stueckelberg [14] precisamente para el
parámetro τ— **sí toma la evolución en τ como unitaria**, con generador hermítico K y
`∫d⁴x |ψ|²` conservada en τ.

O sea existe un marco establecido, para exactamente este parámetro, donde la unitariedad se
retiene. El paso (b) retira una motivación para la unitariedad (la normalización por celda)
pero no la otra (que la τ-evolución sea evolución física con probabilidad conservada).

El paper tiene que **elegir explícitamente entre dos lecturas de T**:

- **T como matriz de transferencia** de una suma sobre caminos. Las matrices de
  transferencia no son unitarias a priori: en mecánica estadística son positivas, en QFT
  euclídea hermíticas positivas. Lo que se les exige es invertibilidad y preservar el
  bilineal que la teoría declare físico. Esta lectura sostiene la construcción.
- **T como evolución física en τ**, à la SHP. Esta lectura exige unitariedad y mata la
  construcción.

El paper opera con la primera sin decirlo, y cita [14], que trabaja en la segunda. Declarar
y defender esa elección es, en mi lectura, la tarea principal que le queda a §5.

---

# Hallazgo 7 — el sector libre

## 7.1 La equivalencia es exacta pero no uniforme

Encontré la forma cerrada del condicionamiento de la métrica:

```
cond(η₊) = (E + Δ) / (E − Δ)      exacto
```

Verificado a 5 dígitos en ocho valores de E/Δ entre 50 y 1.0001. Para un modo libre con
`E = √(Δ² + c²p²)` se tiene `E − Δ ≈ p²/2m`, de donde

```
cond(η₊) ≈ 4 (mc/p)²
```

verificado: p/mc = 0.1 → 402 (predicho 400); p/mc = 0.001 → 4.0000×10⁶ (predicho 4×10⁶).

Así que **la similitud a una teoría hermítica existe para todo p ≠ 0 pero se degrada como
p⁻²**. Es puntual, no uniforme. Un paquete con peso cerca de p = 0 explora la región donde
el mapa está arbitrariamente mal condicionado — y eso es exactamente lo que produce el
crecimiento que el paper reporta en §8.

## 7.2 El crecimiento del §8 es real, y su origen es ese mal condicionamiento

Verifiqué la pendiente `2π|φ(0)|²` que el paper reporta para el paquete gaussiano de ancho
0.3 mc: **predicho 8.3554, medido 8.3554, razón 1.0000**. La derivación analítica: cerca de
p = 0 el corchete va como `(2Δ²/c²p²) sin²(cpτ/ħ)`, y con `∫ sin²(ap)/p² dp = π|a|` sale la
linealidad en τ con esa pendiente.

Pero **lo que crece es el bilineal positivo, que no es el conservado**. Q se queda quieto.
De modo que el crecimiento no es una predicción: es un síntoma de que `⟨ψ|ψ⟩` no es la
cantidad correcta. Lo cual lleva directo al hueco central.

## 7.3 El trilema de la regla de Born

**El paper nunca enuncia una regla de Born.** Los tres bilineales disponibles fallan cada uno
de manera distinta, y verifiqué los tres:

| bilineal | positivo | conservado | distingue de MC ordinaria |
|---|---|---|---|
| `⟨ψ\|ψ⟩` | sí | **no** — crece | — |
| `Q = ⟨ψ\|σ₃\|ψ⟩` | **no** | sí | — |
| `⟨ψ\|η₊\|ψ⟩` | sí | sí | **no** |

Detalle de cada fallo, verificado:

- `⟨ψ|ψ⟩` para el estado (1,0) en E = 2Δ: 1.0000 → 1.6495 → 1.0670 → 1.5224 en τ = 0,1,2,3.
  No conservado, luego no es probabilidad.
- `Q`: el estado `(1,1)/√2`, perfectamente normalizable, da **Q = 0**. El estado
  `(1, 1.3)` normalizado da **Q = −0.2565**. Ambos conservados a 10⁻¹⁰. Q es una carga, tal
  como el paper mismo dice en §1.4 y §11 — pero entonces no puede ser la probabilidad.
- `⟨ψ|η₊|ψ⟩` sí funciona: positivo y conservado. Pero solo existe para E > Δ, se degrada
  como (E+Δ)/(E−Δ), y **hace que el sector libre sea MC ordinaria por construcción**.
  Verifiqué las probabilidades de transición contra la imagen hermítica:

  | τ | lectura η₊ | imagen hermítica | ¿iguales? |
  |---|---|---|---|
  | 0.5 | 0.8549304327 | 0.8549304327 | sí |
  | 1.0 | 0.7564446005 | 0.7564446005 | sí |
  | 2.0 | 0.9748861239 | 0.9748861239 | sí |

  Coinciden a diez decimales.

Esto es más grave que "falta una predicción". Sin regla de Born la teoría **no dice qué
predice**. §11 concede que Q no es probabilidad, pero no pone nada en su lugar, y el hueco
queda abierto en todo el manuscrito.

## 7.4 De dónde sale realmente la predicción de hidrógeno

Reconstruí la fórmula citada en §1.9. El espectro de `H = Eₙσ₃ + iΔ₀σ₁` es
`±√(Eₙ² − Δ₀²) ≈ |Eₙ|(1 − Δ₀²/2Eₙ²)`, y con `Eₙ = −Ry/n²`:

```
ΔEₙ/|Eₙ| = Δ₀² n⁴ / (2 Ry²)
```

que reproduce exactamente lo citado. La reconstrucción funciona, así que la fórmula es
consistente con la estructura.

**El punto clave:** `√(E² − Δ₀²)` es **un invariante de similitud** — es el espectro de H.
Por eso la predicción sí es observable. Pero por eso mismo **no proviene de la métrica ni de
la estructura de dos ramas**, que son covariantes bajo similitud y por tanto inobservables
en el sector libre (7.1–7.3). Proviene enteramente de postular `Δ₀ ≠ 0`.

Y §11 concede exactamente eso: *"The lattice supplies the ceiling Δ < Ep and the invariant
E² − Δ², and does not supply the value of Δ for a given system."*

De modo que la única predicción falsable del programa descansa sobre un parámetro que la
derivación no fija. Eso pertenece al abstract, no a §11.

## 7.5 La identificación natural ya está excluida

Si se toma `Δ = mc²` para el electrón ligado, como sugiere la analogía con la checkerboard
espacial que el paper usa en todo el §7:

```
ΔEₙ/|Eₙ| = (mc²/Ry)² n⁴/2 = 7.05×10⁸  en n = 1
```

Un corrimiento relativo de ~10⁹. Absurdo por cualquier medida.

Cotas sobre Δ₀ a partir de espectroscopía:

| n | precisión relativa | cota sobre Δ₀ | Δ₀/mc² |
|---|---|---|---|
| 2 | 10⁻¹⁴ | 4.8×10⁻⁷ eV | 9.4×10⁻¹³ |
| 10 | 10⁻¹² | 1.9×10⁻⁷ eV | 3.8×10⁻¹³ |
| 50 | 10⁻¹¹ | 2.4×10⁻⁸ eV | 4.8×10⁻¹⁴ |
| 100 | 10⁻¹⁰ | 1.9×10⁻⁸ eV | 3.8×10⁻¹⁴ |

Es decir, Δ₀ para un electrón ligado está más de **13 órdenes de magnitud** por debajo de
mc². La teoría necesita una regla para qué es Δ en un sistema ligado, y no la tiene.

> **Advertencia de honestidad.** Las precisiones relativas de esa tabla son estimaciones mías
> de orden de magnitud sobre espectroscopía de Rydberg, no cifras tomadas de mediciones
> específicas. La conclusión cualitativa (Δ = mc² excluido por márgenes enormes) es robusta
> frente a varios órdenes de error en esas cifras; **las cotas numéricas exactas no lo son**
> y habría que rehacerlas contra mediciones reales antes de citarlas.

## 7.6 Dónde sí podría haber física — lo constructivo

Dado 7.1–7.3, el sector libre no distingue. Los únicos lugares donde puede haber contenido
invariante son tres:

1. **La fase rota `E < Δ`**, donde no existe η₊ y el crecimiento es invariante. Pero §1.8
   establece que la evolución libre nunca llega ahí. La conclusión que sigue es favorable al
   paper y vale la pena adoptarla: **QTSI no es una modificación de la MC libre, es una
   teoría de la medición.** Ese reencuadre es consistente con todo lo derivado y elimina de
   un golpe la objeción de "¿qué predice distinto?" para el sector libre.

2. **El escalamiento n⁴.** La sensibilidad a Δ₀ va como `√(precisión)/n²`, así que subir n
   un factor 10 gana un factor 100 en Δ₀. Eso es lo que hace de la predicción de hidrógeno
   *"the most accessible falsifiable prediction"*, y conviene decirlo cuantitativamente en
   vez de solo citarlo de [9,10].

3. **Interacciones locales en las variables originales**, que dejan de serlo tras la
   similitud dependiente de p. Ahí sí puede haber contenido no covariante, pero hay que
   exhibirlo.

---

# Recomendaciones revisadas

Reemplazan a las #1 y #6 de `REVISION_CRITICA.md`.

**Para §5 (hallazgo 1) — reparable, y con ganancia:**

1. Sustituir (5.9) por el resultado general: la Θ-covarianza sola fuerza la forma U(1,1) sin
   ansatz. Es más limpio y más fuerte que lo que el paper reclama.
2. Añadir el paso (b): §4.2 retira el warrant de la unitariedad, no la viola. Con eso §5.3
   pasa de ser un cambio de postulados sin anunciar a ser un argumento cerrado.
3. Escribir la interpretación de Θ como redundancia de etiquetado, y decir por qué la
   condición es `ΘTΘ⁻¹ = T` y no `T⁻¹` (porque Θ no toca τ).
4. Declarar y defender la lectura de T como matriz de transferencia frente a la lectura SHP,
   citando que [14] trabaja en la segunda.

**Para el hallazgo 7 — más serio de lo que parecía:**

5. **Enunciar la regla de Born.** Es prerequisito de todo lo demás y ahora no está en ningún
   lado. Si es `⟨ψ|η₊|ψ⟩`, decir explícitamente que el sector libre es MC ordinaria por
   construcción y que el contenido físico vive en la fase rota.
6. Mover al abstract la concesión de §11 sobre Δ: la única predicción falsable depende de un
   parámetro que la derivación no fija.
7. Decir que `Δ = mc²` está excluido para sistemas ligados, y dar la regla que lo reemplaza.
8. Adoptar el reencuadre de 7.6.1: QTSI como teoría de la medición, no como modificación de
   la MC libre.

---

## Reproducir

```
python3 deep_h1_theta_and_unitarity.py        # forma U(1,1), degeneración, Theta como reetiquetado
python3 deep_h7_free_sector_and_hydrogen.py   # cond(eta+), pendiente 2pi|phi(0)|^2, hidrógeno, cotas
python3 deep_h7_born_rule_trilemma.py         # los tres bilineales y sus tres fallos
```
