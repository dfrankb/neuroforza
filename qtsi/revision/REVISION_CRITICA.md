# Revisión crítica

**Paper:** *Forward Feynman Gives Schrödinger, Reversible Feynman Gives QTSI — The postulates of
quantum mechanics as theorems of a sum over histories on an invariant lattice*
**Autor:** Alejandro Frank (ICN / C3, UNAM; El Colegio Nacional). Agosto 2026.
**Base:** `Forward_Feynman_QTSI_revised-6.docx` (~13,265 palabras, 12 secciones, 11 figuras).

Toda afirmación matemática citada abajo fue verificada numérica o simbólicamente. Los scripts
están en este mismo directorio y son reproducibles (`python3 check_*.py`).

---

## 0. Veredicto

El álgebra del paper es correcta. Verifiqué (1.4), (1.5), (1.6), (5.10)–(5.12), (6.1)–(6.6),
(8.3), (8.4), la expansión de la traza de §10, el exponente de Lyapunov y todas las constantes
físicas de §8: todo reproduce, en varios casos a precisión de máquina. La expansión de cuarto
orden (8.4), que es el resultado técnico más delicado, es correcta — el residuo escala como
O(s⁶) bajo reescalamiento conjunto.

El problema no está en el álgebra sino en el **balance de postulados**, que es la tesis central
del trabajo. La afirmación del abstract — *"the list of postulates is shorter at the end of the
derivation than at its beginning"* — no queda establecida. En los dos lugares donde se supone
que ocurre la reducción (§3.4 y §5.3) lo que hay es una **sustitución** de supuestos, no una
eliminación. Los hallazgos 1 y 2 de abajo son, en mi lectura, los que el autor tiene que
responder antes que nada; son también los más fáciles de arreglar, porque en gran medida son un
problema de cómo está enunciado el resultado, no de qué se calculó.

Aparte de eso hay cuatro inconsistencias internas (hallazgos 4–7) que sí requieren decisiones de
fondo, y una lista de erratas menores.

Un mérito que hay que decir de entrada: **§4.2 es intelectualmente honesta de una manera
infrecuente.** El autor construye un argumento de conteo, y acto seguido lo desmonta él mismo
explicando en dos razones distintas por qué no puede fijar el vértice, en vez de dejarlo pasar.
La sección §11 ("Not established") concede limitaciones reales. Eso da confianza en el resto.

---

## 1. §5.3 — el pivote intercambia postulados, no los reduce

Es el punto de quiebre del argumento.

El vértice es `V = I + w σ₁`. El paper aplica dos criterios distintos en los dos casos:

- **Caso espacial (5.6):** se impone unitariedad `V†V = I` ⟹ `w + w* = 0` ⟹ *w imaginaria*, vértice en SU(2).
- **Caso temporal (5.9):** se impone Θ-covarianza `ΘVΘ⁻¹ = V` con `Θ = σ₁K` ⟹ `w = w*` ⟹ *w real*, vértice en U(1,1).

Verificado (`check_algebra_and_pivot.py`): las dos condiciones son **mutuamente excluyentes**.
Imponer ambas fuerza `w = 0`, es decir, ningún pliegue, ningún acoplamiento y `H = Eσ₃` — la
teoría trivial. No hay un régimen donde convivan.

Por tanto QTSI no se obtiene *añadiendo* reversibilidad a la estructura usual. Se obtiene
**retirando la unitariedad y poniendo Θ-covarianza en su lugar.** El conteo de postulados no
baja: cambia de nombre. Esto es distinto de lo que el abstract, §1.10 y §12 afirman.

**Sub-punto, igual de importante.** A nivel del generador, la condición `ΘVΘ⁻¹ = V` con Θ
antiunitaria equivale a `ΘHΘ⁻¹ = −H`. Verificado: `σ₁(Eσ₃ + iΔσ₁)*σ₁ = −Eσ₃ − iΔσ₁ = −H`. Pero
la condición estándar de invariancia bajo inversión temporal en MC es `ΘHΘ⁻¹ = +H`, no `−H`.
Lo que se impone es una **anti**simetría, no la simetría T habitual.

Hay una defensa razonable — τ crece monótonamente a lo largo del camino aunque t no, así que Θ
no es inversión temporal en el parámetro sino intercambio de ramas — pero el paper no la da. Y
como *todo* el contenido de la derivación vive en ese signo, hay que argumentarlo explícitamente
en vez de presentarlo como lectura obvia de "reversibilidad".

Nótese también que el vértice espacial (unitario, SU(2)) **no** es Θ-covariante: verificado,
`|Θ(V) − V| = 0.78`. Como la checkerboard espacial de Feynman describe la ecuación de Dirac, que
sí es invariante bajo inversión temporal en el sentido ordinario, la Θ elegida aquí no es el
operador de inversión temporal de esa teoría. Conviene decir qué es, entonces.

---

## 2. §3.4 — "cuatro nombres para un supuesto" es demasiado fuerte

§3.4 afirma que tiempo monótono, primer orden, unitariedad y positividad de Born son cuatro
nombres de un solo hecho geométrico. El abstract lo repite: *"one geometric fact read three
times"*.

De la univaluación de los caminos se sigue **primer orden** y la propiedad de semigrupo. No se
sigue unitariedad ni positividad.

**Contraejemplo decisivo:** la fórmula de Feynman–Kac. Suma sobre caminos *igualmente
univaluados*, con peso `e^{−S_E/ħ}` (medida de Wiener), y produce la ecuación de calor. Primer
orden en t ✓. Semigroup ✓. Caminos univaluados ✓. Unitariedad ✗ — lo conservado es ∫ψ, no
∫|ψ|², y ∫|ψ|² decae. La geometría de los caminos es idéntica; lo que cambia es el peso.

Lo que produce la unitariedad es entonces el **peso de fase pura** supuesto en (3.1), no la
geometría de §3.3.

Y la justificación que se da para la fase pura está mal argumentada. Cita textual de §3.1:

> *"The weight must be a pure phase: were its magnitude to vary from path to path, the classical
> limit would select the route of largest weight, and no reason exists for that route to coincide
> with the route of extremal action."*

Con peso euclídeo `e^{−S_E/ħ}` el peso máximo es exactamente el de **acción mínima** — sí
coincide con la ruta extremal, por descenso más pronunciado. El motivo dado es falso. El motivo
real de exigir fase pura es la conservación de la probabilidad, es decir, la unitariedad misma.
Con lo cual §3.4 deriva la unitariedad de un supuesto que ya la contenía.

En descargo: §3.4 dice literalmente *"Because it is first order **with a self-adjoint
generator**"*, reconociendo el segundo insumo, y §3.2 exige explícitamente autoadjunción para
Trotter–Kato. El problema es que el abstract, §1.10 y §12 omiten esa cláusula y presentan el
resultado como reducción a un solo hecho geométrico. Es un problema de enunciado, y por eso es
arreglable sin tocar el cálculo.

---

## 3. §7 no es una ruta independiente

§7 dice: *"obtained a second time by an independent route... The two derivations agree, which is
the check that matters."*

No es independiente. (6.6) se dedujo de H en §6; §7 recupera H de (6.6). Es una inversión, no una
corroboración. En concreto: el único insumo que fuerza `Γ₂² = −I` — y con ello la firma
indefinida, la métrica η y toda la estructura su(1,1) — es el signo menos en `E² − Δ²`, que entró
en §5 con w real. El acuerdo entre §5 y §7 está garantizado por construcción.

Lo que §7 sí establece, y es valioso, es que **dada** la firma, la linealización es esencialmente
única en 2×2. Eso vale la pena decirlo así.

---

## 4. Contradicción interna sobre la longitud mínima

El paper usa y repudia una longitud mínima en secciones distintas.

**Usa:**
- §2.3: *"A particle cannot be a point, since nothing is localisable below ℓp. An electron occupies at least one cell."*
- §8.1: *"A path moves exactly one spatial cell per temporal cell, so the maximum speed is c = δ/ε = ℓp/tp"*, con δ = ℓp.

**Repudia:**
- §9.2: una cota inferior de longitud es violada por cualquier observador con γ > ℓ₀/ℓp, y restaurarla exige relatividad doblemente especial.
- §9.3: *"No proper length is attached to a worldline"* y *"The Planck length enters Section 8 in exactly two places, and neither is a coordinate bound: through the ratio c = ℓp/tp, which is the statement that a path advances one spatial cell per tick."*

La defensa de §9.3 es circular: *"un paso espacial por tick"* **es** una discretización del
espacio con celda ℓp. No se puede tener una celda espacial sin discretizar el espacio. Y §2.3
afirma directamente una extensión mínima del electrón, que es exactamente la clase de enunciado
que §9.2 declara insostenible bajo boost.

O bien se retiran §2.3 y la lectura espacial de §8.1, o bien §9.2/§9.3 no aplican a esta
construcción. Ahora mismo el paper sostiene las dos cosas.

---

## 5. §9.4 se refuta a sí mismo

§9.4 cita [22] (Bombelli–Lee–Meyer–Sorkin) y [23] (Bombelli–Henson–Sorkin) para argumentar que la
discreción no rompe Lorentz, y concede explícitamente:

> *"What breaks Lorentz invariance in a lattice is regularity, not discreteness, and the two are
> routinely conflated."*

Correcto — y es precisamente el problema. §5.1 construye *"a chain of cells of size ε"*: una red
**regular**. El teorema citado cubre el sprinkling de Poisson a densidad fija, que es justamente
el caso que **no** es el de esta construcción. La sección enuncia el criterio que la condena y no
lo aplica a sí misma.

Si el argumento de covariancia va a apoyarse en [23], la construcción tiene que ser un sprinkling,
no una cadena regular — y entonces el conteo binomial de (5.3), que depende de N vértices
equiespaciados, hay que rehacerlo.

---

## 6. Dos criterios de medición incompatibles

- **§1.8:** la medición es lo que empuja el modo a través del borde **inferior**, hacia `E < Δ`
  (fase rota, crecimiento invariante). *"Reaching Δ > E requires something to drive the mode
  across the band edge, and in the measurement context that something is the coupling to the
  apparatus."*
- **§8:** *"A collapse criterion phrased as the driving of τRC down to its minimum therefore
  selects the **upper** edge, and identifies measurement with the exhaustion of the available
  temporal resolution."*

Son bordes opuestos de la banda y el paper no los reconcilia.

Además verifiqué dónde está el borde superior. La condición es `cos θ = −sech χ`, de donde para
Δ ≪ Ep resulta **E ≈ π·Ep**:

| Δ/Ep | E borde inferior | E borde superior |
|---|---|---|
| 4.19×10⁻²³ (electrón) | ≈ Δ | 3.1416 Ep ≈ 3.8×10¹⁹ GeV |
| 0.15 (Figura 9) | 0.1506 Ep | 2.9910 Ep ≈ 3.7×10¹⁹ GeV |
| 0.9 | 1.1198 Ep | 2.0218 Ep |

Identificar la medición con el borde superior exige excitar el modo a energía **transplanckiana**
(~3 Ep). El LHC llega a ~1.4×10⁴ GeV, quince órdenes de magnitud por debajo de Ep. Como criterio
de colapso de laboratorio no es viable, y conviene decirlo o abandonar esa lectura.

---

## 7. El sector libre es indistinguible de la MC ordinaria — por el propio argumento del paper

Dos afirmaciones del paper, ambas correctas, se combinan en algo que el paper no extrae:

- §1.8: todo modo libre con Δ = mc² cumple `E ≥ Δ`, con igualdad solo en reposo.
- §1.7: para `E > Δ` existe una métrica positiva η₊ en la que la evolución es genuinamente unitaria.

Verificado (`check_free_sector_and_constants.py`): `h = η₊^{1/2} H η₊^{−1/2}` es hermítica a
precisión de máquina (residuo ~10⁻¹⁵) para E/Δ ∈ {1.001, 1.5, 2, 5, 50}. Es decir, **todo modo
libre es una transformación de similitud de un sistema hermítico de dos niveles.**

Consecuencia: en todo el régimen libre, QTSI y la MC ordinaria coinciden en espectros y
probabilidades de transición. El paper dice la mitad de esto (*"growth there is a choice of inner
product"*) pero no saca la conclusión.

Queda entonces sin explicar de dónde sale la desviación de hidrógeno
`ΔEn/|En| = (Δ₀²/2Ry²)·n⁴` que §1.9 cita de [9,10] y describe como *"the most accessible
falsifiable prediction of the scheme"*. Si el régimen libre es equivalente por similitud, la
predicción tiene que venir del parámetro nuevo Δ₀ y no de la estructura derivada aquí. Vale
aclararlo, porque de ello depende que el programa sea falsable.

**Matiz justo, para no sobrepasarme:** la similitud η₊ depende del momento, así que es no-local en
posición. Lo que cambia entre las dos descripciones no son los espectros sino **qué operadores
cuentan como observables**. Esa es la salida disponible para el autor — pero entonces hay que
decir cuáles son, explícitamente.

---

## 8. Detalles técnicos verificados y erratas

**(7.4) no fuerza Γ₂ anti-hermítica.** §7 afirma: *"Γ₂ squares to −I and must therefore be
anti-Hermitian"*, y de ahí concluye que la métrica η = Γ₁ *"comes with the algebra rather than
being imposed"*. Falso como está enunciado. `Γ₂ = [[0, b], [−1/b, 0]]` cumple `Γ₂² = −I` y
`{σ₃, Γ₂} = 0` para todo b ≠ 0, y es anti-hermítica **solo si |b| = 1**. Verifiqué b = 2 y
b = 0.5: ambos satisfacen (7.4) completa y `H² = (E²−Δ²)I`, sin ser anti-hermíticos. Hace falta
una elección adicional de hermiticidad, lo cual afecta directamente la frase *"neither the metric,
nor the indefiniteness of the conserved form, nor the algebra ... has been assumed at any point."*

**El borde inferior de banda no está en E = Δ.** En la red es `E = Ep·arcsin(Δ/Ep)`, no `E = Δ`.
Coinciden a primer orden. Para Δ = 0.5 Ep el corrimiento es 4.7%; para Δ = 0.9 Ep es 24%. §8 y
§11 dicen *"the lower edge, E = Δ"* sin la salvedad — es un enunciado del continuo trasladado a la
red sin marcar.

**Notación de w.** (5.6) concluye *"w purely imaginary"*, pero (5.10) escribe *"w = tan φ in the
spatial case"* con `V = cos φ·I + i sin φ·σ₁`, lo que corresponde a `w = i·tan φ`. Inconsistencia
de notación entre las dos ecuaciones.

**§8.2, la cota |w| < 1.** Se sigue de exigir que el vértice sea un boost de parámetro **real**,
no de la invertibilidad. `det(I + wσ₁) = 1 − w²` se anula solo en |w| = 1; para |w| > 1 la matriz
sigue siendo invertible (determinante negativo). Así que *"a maximum energy emerges from the
counting rather than from an additional assumption"* está sobredicho: emerge de la forma de boost.

**§10, dos factores de crecimiento no reproducen.**

| Caso | Paper | Verificado |
|---|---|---|
| E = 2Δ, elíptico (rango de oscilación) | 1.66 | 1.667 ✓ |
| E = Δ, parabólico | 197 | **201** |
| E = 0.5Δ, hiperbólico | 1.87×10⁷ | **2.219×10⁷** |
| Lyapunov en E = 0.5Δ | 0.866026 | 0.866026 ✓ |
| Q conservada, elíptico | 2.6×10⁻¹² | 2.62×10⁻¹² ✓ |
| Q conservada, parabólico | 3.2×10⁻¹² | 3.31×10⁻¹² ✓ |
| Q conservada, hiperbólico | 1.6×10⁻⁸ | 2.56×10⁻⁸ ✓ (mismo orden, en el piso aritmético) |

Probé los tres ordenamientos del paso (VD, DV, simétrico) y además la forma cerrada del continuo
`exp(−iHτ)`: los cuatro dan 201 y 2.219×10⁷. Los valores del paper corresponden a τΔ/ħ ≈ 9.9 y
≈ 9.67 respectivamente, no a 10. Parece un desajuste de bookkeeping en el punto final del
intervalo, no un error de física — el exponente de Lyapunov, que es la cantidad físicamente
significativa, reproduce a seis cifras. Conviene recomputarlos.

**El piso τRC ≥ tp/π está por debajo de tp.** Merece un comentario en el texto: una teoría cuyo
punto de partida es que tp es el intervalo elemental entrega un piso para el tiempo de
interferencia que es ~3.14 veces **menor** que tp. No es necesariamente inconsistente (τRC no es
un intervalo de la red), pero como está queda raro sin explicación.

**(Δ/Ep)²/3 para el electrón** = 5.84×10⁻⁴⁶, que el paper redondea a *"about 10⁻⁴⁵"*. Aceptable
como orden de magnitud; podría decir 10⁻⁴⁶.

---

## 9. Lo que sí se sostiene (verificado)

Para que quede claro qué parte del edificio está firme:

- **(1.4)** σ₃-pseudo-hermiticidad `σ₃Hσ₃ = H†` — residuo exactamente 0.
- **(1.5)** El álgebra su(1,1): las tres relaciones de conmutación, con el signo de `[K₁,K₂] = −iK₃` — residuo 0. Verifiqué además que `{σ₃, iσ₁, iσ₂}` es efectivamente el span real completo de las matrices sin traza que cumplen (1.4).
- **(1.6)** `H² = (E² − Δ²)I` — residuo 0.
- **(5.10)/(5.11)** `V_time†σ₃V_time = σ₃` exacto (2×10⁻¹⁶) y `V_time†V_time ≠ I`; `V_space†V_space = I` exacto y `V_space†σ₃V_space ≠ σ₃`. La Θ-covarianza es exacta para D y para el vértice real, y falla para el imaginario. Todo exactamente como se afirma en §10.
- **(5.12)** `H = Eσ₃ + iΔσ₁` sale correctamente de `V_time·D` a primer orden.
- **(6.1)–(6.6)** La eliminación de la componente es correcta y sin aproximación; los términos en ∂τψ₊ cancelan como se dice.
- **(8.3)** `cos(Ωε) = cosh χ · cos(Eε/ħ)`, vía `tr T = 2 cosh χ cos θ` y `det T = 1`. Correcto.
- **(8.4)** Verificado por escalamiento de residuos: fijando x = 2s, y = s y reduciendo s, el residuo contra la fórmula completa cae como **O(s⁶)** (razón 10⁶ por década) mientras que contra el término principal cae como O(s⁴). Los términos de cuarto orden son correctos. El cambio de signo en E² = 2Δ² también: la corrección es `(Δ²/3Ep²)(E² − 2Δ²)`.
- **§10, expansión de la traza:** `|tr T| − 2 = −3.0000×10⁻⁶` en E = 2Δ y `+7.5000×10⁻⁷` en E = 0.5Δ. Exactos como se reportan.
- **Constantes de §8:** ℓp = 1.6163×10⁻³⁵ m, tp = 5.3912×10⁻⁴⁴ s, mp = 2.1764×10⁻⁸ kg, Ep = 1.22×10¹⁹ GeV, τC(e⁻) = 1.2881×10⁻²¹ s, w(e⁻) = 4.1855×10⁻²³ = m_e/m_p. La identidad `λC·r_s = 2ℓp²` se cumple a 6 cifras. Todas correctas.
- **§1.7, η₊ en E = 2Δ:** autovalores de H en ±√3 ✓; intertwining `η₊H = H†η₊` a 2×10⁻¹⁶; η₊ con dos autovalores positivos en razón exactamente 3, consistente con el "2/3 y 2" reportado salvo normalización.

---

## 10. Recomendaciones

Ordenadas por impacto sobre la tesis.

1. **Rehacer el balance de postulados.** Decir explícitamente que §5.3 sustituye unitariedad por Θ-covarianza y que las dos son mutuamente excluyentes (`w = 0` es la única intersección). Y justificar por qué `ΘHΘ⁻¹ = −H` es la lectura correcta de "reversibilidad" en lugar de la condición T estándar `+H`. El argumento de que τ crece monótonamente aunque t no lo haga es defendible y basta — pero hay que escribirlo.

2. **Corregir §3.1 y matizar §3.4, el abstract y §12.** El argumento para la fase pura es incorrecto tal como está (el peso euclídeo sí selecciona la acción extremal). Reemplazarlo por la razón real: conservación de probabilidad. Y separar lo que se sigue de la univaluación (primer orden, semigrupo) de lo que requiere además el peso de fase pura (unitariedad, positividad). Citar Feynman–Kac como el contraejemplo que delimita el alcance — fortalece el paper en vez de debilitarlo, porque hace preciso un enunciado que ahora es vulnerable.

3. **Presentar §7 como consistencia de la inversión**, no como *"independent route"*. El resultado real y defendible es la unicidad esencial de la linealización dada la firma.

4. **Resolver la longitud mínima.** §2.3/§8.1 y §9.2/§9.3 no pueden coexistir. Mi sugerencia: retirar la afirmación de que el electrón ocupa al menos una celda, y reformular (8.1) de modo que c salga de la razón de escalas sin postular una celda espacial.

5. **Elegir un solo borde de banda para la medición.** Si es el superior, decir que exige ~π Ep y qué significa eso operacionalmente.

6. **Decir qué observable separa QTSI de la MC ordinaria en el régimen libre**, dado §1.7 + §1.8. Es lo que decide si el programa es falsable, y ahora mismo el paper deja abierta la posibilidad de que no lo sea en todo el sector libre.

7. **Ajustar §9.4** o cambiar la red a un sprinkling: el criterio que la propia sección enuncia (regularidad, no discreción, es lo que rompe Lorentz) condena la cadena regular de §5.1.

8. **Menores:** corregir la afirmación sobre Γ₂ anti-hermítica en (7.6); marcar que el borde inferior es `arcsin(Δ/Ep)`; unificar la notación de w entre (5.6) y (5.10); recomputar los factores 197 y 1.87×10⁷; comentar el piso τRC = tp/π < tp.

---

## Apéndice: reproducir la verificación

```
pip install numpy scipy sympy mpmath
python3 check_eq84_expansion.py            # (8.4) por escalamiento de residuos
python3 check_sec10_transfer_matrix.py     # §10 completo + bordes de banda
python3 check_algebra_and_pivot.py         # (1.4)-(1.6), vértices, incompatibilidad unitariedad/Θ, (7.4)
python3 check_free_sector_and_constants.py # equivalencia del sector libre, constantes de §8
python3 check_growth_factors.py            # ordenamientos del paso vs factores reportados
```
