# The Per-Node Rank Theorem: Uniform Spectral Counts for Fresnel-Weighted Stars

**Status**: 🔄 IN PROGRESS | Last updated: 2026-07-31

Replaces the refuted derivation in
[NDF_ANGULAR_SATURATION_DERIVATION](NDF_ANGULAR_SATURATION_DERIVATION.md)
(kept unedited as a record). That derivation claimed per-node star rank
`O(√(ka))` via a Slepian/Landau–Pollak bandlimit argument whose Step C
restricts the star to a circle `r = a` — a rank-*lowering* linear map used
backwards — and thereby discards the radial dimension entirely.

This document supplies (i) a fully proved counting lemma for Gaussian kernel
matrices, uniform in the number of points (Stage 1); (ii) a proved angular
bandwidth lemma and a conditional product theorem for a node's star, together
with a direct empirical decomposition of the star into its angular and radial
factors (Stage 2); (iii) an explanation and out-of-sample test of the
perimeter-versus-interior saturation asymmetry, including the measured
exponent 0.696 (Stage 3); and (iv) what survives at whole-mesh level
(Stage 4). Claims are labelled **proved**, **imported**, **measured**, or
**conjecture** throughout.

Numerical tests: `scripts/verify_gaussian_count_lemma.py` (Stage 1),
`scripts/ndf_star_family_decomposition.py` (Stage 2),
`scripts/ndf_perimeter_two_term.py` (Stage 3). Raw outputs in
`docs/4-research/results/gaussian_count_lemma_verification.json`,
`ndf_star_family_decomposition.json`, `ndf_perimeter_two_term.json`.

---

## 1. Setup and notation

| symbol | meaning |
|---|---|
| `Ω ⊂ ℝ²` | room, width `W`, height `H`, area `A`, radius `a = √(A/π)`, perimeter `P` |
| `λ, k = 2π/λ` | wavelength, wavenumber; `ka` = electrical size |
| `p₁ … p_N` | node positions; `L = N(N−1)/2` links |
| `d_ij` | length of link `(i,j)` |
| `σ_F(d) = ¼√(λd)` | Gaussian cross-section width (σ_F² = λd/16) |
| `w_ij(x)` | row weight: `exp(−d_⊥²/(2σ_F²)) · 4t(1−t)`, clipped to `t∈[0,1]`, truncated below `e^{−4.5}` |
| `W ∈ ℝ^{L×V}` | sensing matrix (`opencsi.geometry.ndf.build_sensing_matrix`) |
| `NDF(τ)` | `#{m : σ_m > τ·σ₁}`, default `τ = 0.01` |
| `S_i` | star of node `i`: the `N−1` rows of links incident on `i` |
| `N_ε(K)` | `#{m : λ_m(K) > ε·λ₁(K)}` for a PSD matrix `K` |

Note `NDF(τ)` of `W` equals `N_{τ²}` of the Gram `G = W Wᵀ`, since
`σ_m(W)² = λ_m(G)`. All experiments below use the reference implementation of
the weight model; the test room is 5×4 m (`a = 2.52` m, `P = 17.2` m after the
0.1 m wall inset).

---

## 2. Stage 1 — the counting lemma (proved)

### 2.1 Statement

> **Lemma 1 (uniform spectral count for Gaussian kernel matrices).**
> Let `s > 0`, let `u₁, …, u_M` be **any** `M` points in an interval of length
> `T`, and let `K_M = [exp(−(u_j − u_{j'})²/(4s²))]`. Then for every
> `ε ∈ (0, ½]`, with `q = ⌈T/s⌉`,
>
> ```
> N_ε(K_M)  ≤  (T/(πs)) · √(ln(2q/ε))  +  (2/π) · ln(8q/ε)  +  4 ,
> ```
>
> **uniformly in `M`** (and trivially `N_ε ≤ M`).

This is the target form `C·(T/s)·√(log(1/ε)) + C′·log(1/ε)` except that
`log(1/ε)` is replaced by `log(q/ε)`: the argument leaves an extra `log(T/s)`
inside the logarithms. That slack comes from the Frobenius step in the proof
(§2.5); Widom's asymptotics for the continuum operator [Widom 1964] indicate
the sharp count is `(T/(πs))·√(ln(1/ε))·(1+o(1))`, but we do not claim it.
For every downstream use in this document the weaker form suffices, since only
polylogarithmic factors differ.

The mechanism to note: uniformity in `M` is **not** free — as points
accumulate, `λ₁(K_M)` grows proportionally (Step 1 below), which raises the
absolute threshold `ε·λ₁` exactly fast enough to absorb the extra points.

### 2.2 Step 1: clustering lower bound on λ₁

`λ₁(K_M) ≥ max(1, e^{−1/4}·M/q)`.

*Proof.* Diagonal entries equal 1, so `λ₁ ≥ 1`. Partition the interval into
`q` subintervals of length `≤ s`. Some subinterval holds `m₀ ≥ M/q` points;
these have pairwise distances `≤ s`, hence pairwise kernel values
`≥ e^{−s²/(4s²)} = e^{−1/4}`. The Rayleigh quotient of the indicator vector of
this subset gives `λ₁ ≥ m₀²e^{−1/4}/m₀ = e^{−1/4}m₀`. ∎

### 2.3 Step 2: an explicit rank-n approximant with entrywise error η

Let `g(t) = e^{−t²/(4s²)}`, with Fourier transform
`ĝ(ω) = 2s√π·e^{−s²ω²} ≥ 0`. For step `h > 0` and cutoff `J ∈ ℕ` define

```
A(t) = (h/2π) Σ_{|j|≤J} ĝ(jh) e^{ijht} .
```

On the points, `A_M = Φ D Φ*` with `Φ_{ℓj} = e^{ijh·u_ℓ}` and
`D = diag((h/2π)ĝ(jh)) ⪰ 0`: a real symmetric PSD matrix of rank `≤ 2J+1`.

Fix `η ∈ (0, 2/e]` and choose

```
Δ = 2s√(ln(8/η)) ,   h = 2π/(T+Δ) ,   Ω = s^{-1}√(ln(2/η)) ,   J = ⌈Ω/h⌉ .
```

Then `|K(u_ℓ,u_{ℓ'}) − A_{ℓℓ'}| ≤ η` for all pairs.

*Proof.* Split via the full lattice sum `S(t) = (h/2π)Σ_{j∈ℤ} ĝ(jh)e^{ijht}`,
which by Poisson summation equals `Σ_{m∈ℤ} g(t + 2πm/h)`.

(i) *Aliasing.* For `|t| ≤ T` and `m ≥ 1`,
`|t + 2πm/h| ≥ m(T+Δ) − T = (m−1)T + mΔ ≥ mΔ`, so
`|g − S|(t) ≤ 2Σ_{m≥1} e^{−m²Δ²/(4s²)} ≤ 2·(η/8)/(1 − η/8) ≤ η/2`
using `e^{−Δ²/(4s²)} = η/8 ≤ 1/8`.

(ii) *Spectral truncation.* With `sJh ≥ sΩ = √(ln(2/η)) ≥ 1`,

```
|S − A|(t) ≤ (2sh/√π) Σ_{j>J} e^{−(sjh)²} ≤ (2/√π) ∫_{sΩ}^∞ e^{−y²} dy
           ≤ e^{−(sΩ)²}/(√π·sΩ) ≤ η/2 . ∎
```

The rank satisfies

```
n = 2J+1 ≤ Ω(T+Δ)/π + 3 = (T/(πs))√(ln(2/η)) + (2/π)√(ln(2/η)·ln(8/η)) + 3
  ≤ (T/(πs))√(ln(2/η)) + (2/π)ln(8/η) + 3 .
```

### 2.4 Step 3: eigenvalue counting

Set `η = ε/q` and `E = K_M − A_M`, so `‖E‖_F ≤ Mη`. Since `A_M ⪰ 0` with
rank `≤ n`, Weyl's inequality gives `λ_{n+j}(K_M) ≤ λ_j(E)` for `j ≥ 1`.
With threshold `θ = ε·λ₁(K_M) ≥ ε·max(1, e^{−1/4}M/q)` (Step 1):

```
#{ m > n : λ_m(K_M) > θ }  ≤  #{ j : λ_j(E) > θ }  ≤  ‖E‖_F²/θ²
  ≤  M²(ε/q)² / (ε²·max(1, e^{−1/2}M²/q²))  ≤  e^{1/2} < 1.65 ,
```

checking both branches of the max. The count is an integer, so it is `≤ 1`.
Adding `n` yields Lemma 1. ∎

### 2.5 Where the log(T/s) slack enters, and tightness

The Frobenius bound `‖E‖_F ≤ Mη` forces `η ≃ ε/q`; a trace bound on a PSD
remainder, or a Landau–Widom two-term analysis of the plunge region
[Landau & Widom 1980], would remove the `q` inside the logs. Conversely,
equispaced points with spacing `2s√(ln(2/ε))` give (Gershgorin) all `M`
eigenvalues within `1 ± 1.1ε` of 1, hence `N_ε ≥ M ≈ T/(2s√(ln(2/ε)))`: the
factor `T/s` is sharp and only the position of the `√log` factor separates the
upper and lower bounds.

### 2.6 Numerical verification (0 violations / 240 configurations)

`scripts/verify_gaussian_count_lemma.py`: `T/s ∈ {5, 20, 80, 320}`,
`M ∈ {30, …, 4000}`, `ε ∈ {10⁻², 10⁻⁴, 10⁻⁶}`, placements uniform-random,
equispaced, clustered (7 clusters), and near-duplicate pairs.

- **0 violations in 240 configurations**; worst tightness `count/bound = 0.80`.
- **Uniformity in M confirmed**: counts freeze while `M` grows 133× (e.g.
  `T/s = 320`, equispaced, `ε = 10⁻²`: counts 30, 100, 220, 220, 220 against
  bound 350).
- Clustered configurations sit far below the bound (count 14 vs bound 350 at
  `T/s = 320`): clustering inflates `λ₁` and demotes eigenvalues, exactly the
  Step 1 mechanism.

---

## 3. Stage 2 — the star: angular × radial structure

### 3.1 Exact polar form (proved, restating the model)

Place node `i` at the origin; partner `j` at bearing `θ_j`, distance `d_j`,
`σ_j = σ_F(d_j)`, `α_j(r) = r²/(2σ_j²)`. In polar coordinates `x = (r, φ)`:

```
w_j(r, φ) = exp(−α_j(r)·sin²(φ−θ_j)) · 4t(1−t) · 1[0 ≤ t ≤ 1] · 1[w ≥ e^{−4.5}],
t = (r/d_j)·cos(φ−θ_j) .
```

Note `t ≤ 1` confines the row essentially to `r ≲ d_j`: on the support,
`r² = r²sin² + r²cos² ≤ 2σ_j²·ln(1/η) + d_j²`, so
`α_j ≤ 8d_j/λ + ln(1/η)` (using `σ_j² = λd_j/16`).

### 3.2 Angular bandwidth lemma (proved)

> **Lemma 3.** For `z > 0` and integer `1 ≤ k`, `e^{−z}I_k(z) ≤ e^{−k²/(4z)}`
> whenever `k ≤ z`. Consequently the angular Fourier coefficients of
> `exp(−α sin²(φ−θ))` satisfy `|c_m| ≤ e^{−m²/(8α)}` (only even `m` occur),
> and the relative `L²` angular tail of `w_j` beyond harmonics
>
> ```
> |m| > M_j = √(8·α_max·[ln(1/η²) + ln(4·α_max)]) + 2 ,   α_max = 8d_j/λ + ln(1/η),
> ```
>
> is at most `η` plus the clip/truncation terms of §3.3.

*Proof sketch.* `e^{−α sin²u} = e^{−α/2}·e^{(α/2)cos 2u}
= e^{−α/2} Σ_k I_k(α/2) e^{2iku}` (Jacobi–Anger for `I_k`). For the Bessel
bound, shift the contour in `I_k(z) = (1/2π)∫ e^{z cos v} e^{−ikv} dv` by
`v → v + iβ`: `|I_k(z)| ≤ e^{−kβ + z(cosh β − 1)}`; choosing
`sinh β = k/z` gives exponent `−k·asinh(k/z) + z(√(1+k²/z²) − 1) ≤ −k²/(4z)`
for `k ≤ z`. (Numerically the bound is conservative by ≥ 4× for
`z ∈ [2, 1000]`; checked in `verify_gaussian_count_lemma.py`.) Summing the
squared tail and dividing by `‖e^{−αsin²}‖² ≥ c/√α` yields the `ln(4α)`
correction. The envelope `4t(1−t) = 4(r/d)cosΔφ − 4(r/d)²cos²Δφ` is exactly a
degree-2 trigonometric polynomial in `φ`, shifting bandwidth by 2. ∎

With `d_j ≤ D_i` (largest in-room distance from node `i`),
`√(8·8d/λ) = 8√(d/λ)`, so the retained harmonic count is

```
n_ang(i) = 2M_i + 1 ≈ 8·√(D_i/λ)·√(ln(1/η²)+ln(32D_i/λ)) = Θ(√(ka·log)) .
```

This replaces the refuted Step B–C: no restriction to a circle, no Slepian
theorem (which concerns strictly bandlimited functions and does not apply),
and the radial dimension is **retained**, as follows.

### 3.3 Clip and truncation remainders (estimated, then measured)

Three non-smooth features fall outside Lemma 3; all are small at `τ = 0.01`:

- *`t = 0` clip* (`|φ−θ_j| = π/2`): jump height `≤ (4r/d)e^{−α(r)}`; total
  relative tail energy `O(σ³/(d³M))` — negligible.
- *`t = 1` clip* (`r > d` wedge): weight vanishes continuously (kink only);
  relative tail `O(M⁻³)` — negligible for `M ≥ 10`.
- *hard truncation at `e^{−4.5}`*: removing the truncation changes each row by
  relative `L²` norm `≈ (e^{−4.5})/√(2·ln(4e^{4.5}))·…` ≈ **0.6 %** — below
  `τ` per row, consistent with the implementation note in `ndf.py` that NDF
  changes by ≤ 2 between weight-model variants. At `τ = 0.01` this is a
  *floor*: model features of amplitude `e^{−4.5} = 1.1τ` are borderline
  significant, and no counting theorem at this `τ` can ignore them entirely.

### 3.4 The product theorem (proved conditional on a radial count)

> **Theorem 4 (star rank).** Fix `η ∈ (0,1)` and let `M_i` be as in Lemma 3.
> For `|m| ≤ M_i` let `ρ_m(η)` be the `η`-rank (relative-threshold eigenvalue
> count) of the radial family `{c_m(·; d) : d ∈ [d_min, D_i]}` of harmonic-`m`
> radial profiles. Then, by the Stage-1 counting mechanism applied to the
> angular-truncation remainder,
>
> ```
> rank_τ(S_i)  ≤  Σ_{|m| ≤ M_i} ρ_m(η)  +  C_rem
>              ≤  (2M_i + 1) · ρ̄(η)  +  C_rem ,     ρ̄ = max_m ρ_m ,
> ```
>
> with `C_rem = O(1)` once `η ≲ τ/√(N·coherence)` (Frobenius step, same slack
> as §2.5). **The angular factor is proved** (Lemma 3). **No analytic bound on
> `ρ̄` is proved here**; it is measured below. Claiming the product without a
> radial bound would repeat the refuted document's error in mirror image, so
> the radial factor is stated as measured fact + conjecture.

### 3.5 What the radial factor actually is (measured, decisive)

`scripts/ndf_star_family_decomposition.py`, room 5×4 m, central node,
`τ = 0.01`, four wavelengths (`ka = 31.7 … 253.7`, an 8× range),
grid-convergence checked (0.0 % change at 1.5× finer grid):

**E1 — collinear fan** (one bearing, distances `d ∈ [0.5, 3.06]` m,
log-ladder, bit-reversed prefix order): rank plateaus at

```
ka = 31.7 : 13     ka = 63.4 : 13     ka = 126.8 : 13     ka = 253.7 : 13
```

**exactly wavelength-independent over 8× in `ka`** (log-log slope 0.000).
This was predicted, not fitted: coaxial beams of widths `σ, σ'` overlap
transversally as `√(2σσ'/(σ²+σ'²))`, a function of `(d'/d)^{1/4}` only —
`λ` cancels — and the axial `4t(1−t)` overlap is `λ`-free as well. The
radial family of this weight model carries **no electrical length scale**.

**E2 — ring fan** (one distance `R = 1.8` m, dense bearings): plateaus

```
ka = 31.7 : 31     ka = 63.4 : 45     ka = 126.8 : 63     ka = 253.7 : 89
n_ang / √(ka) = 5.51, 5.65, 5.60, 5.59   (log-log slope 0.505)
```

The angular mode count is `(5.6 ± 0.1)·√(ka)` — a `√(ka)` law clean to 1 %,
sitting a factor ≈ 4 below the Lemma 3 upper bound (`2M+1 ≈ 250` at
`ka = 126.8`, `d = 1.8` m). The slack splits as ≈ 2× from the coefficient
bound itself and ≈ 2× from the conservative `‖G‖ ≥ |a₀|` norm estimate that
produces the `ln(4α)` term; the *scaling* is exact.

**E3 — random star** (partners uniform in the room, up to 3200 partners,
2 seeds): see §3.6 for the plateaus and the scaling verdict.

### 3.6 Verdict on O(ka) versus O(√ka) for a single star

*(measured; this corrects the prior belief stated at task-issue time)*

The working hypothesis was that the per-node count is a product of an angular
`√(ka)` and a radial `√(ka)` factor, i.e. `O(ka)`. The measurements refute
the second factor **for this weight model**: E1 shows the radial family is
`Θ(1)` in `ka` (13 modes at `τ = 0.01`, exactly `λ`-independent). The physical
intuition behind a radial `√(ka)` — Fresnel-zone radial structure at scale
`λ` — belongs to a *coherent* (phase-carrying) kernel; the implemented weight
model is **phaseless**, and its only `λ`-dependent length is the transverse
width `σ_F ∝ √λ`. Hence:

> **Per-node law (measured, 8× range in ka).**
> `rank_τ(S_i) = C_star(τ) · √(ka)` with `C_star(0.01) ≈ 24 ± 1` for a
> central node — the product of the angular count `5.6·√(ka)` and an
> *effective* radial multiplicity `≈ 4.3 ≤ ρ̄ = 13`.
> E3 plateaus (3200 partners): §3.6a.

Two consequences:

1. The refuted document's *exponent* `√(ka)` happens to be right for a single
   star, but its argument was invalid, its constant (`c₁ ≈ 4.84`, i.e. bound
   109 at `ka = 126.8`) is exceeded by the measurement itself (≥ 244), and its
   restriction step discarded the radial factor that supplies the missing ≈ 4×.
   Finding 2's "no saturation at `NDF = N−1 = 199`, `N = 200`" is fully
   consistent: the central-star ceiling at `ka = 126.8` sits near
   `24·√(126.8) ≈ 270`, above 199 — the mesh simply had not reached it.
2. The middle term of any three-regime bound must read
   `c₁·N·√(ka·log)` with `c₁` including the radial multiplicity — and it is a
   *bound*, vacuous as an additive union (§5).

### 3.6a E3 plateaus (3200 partners, 2 seeds)

| ka | plateau (seed 0 / 1) | last step | NDF/√(ka) | NDF/ka |
|---|---|---|---|---|
| 31.7 | 140 / 141 | +1.4 % / +0.7 % | 25.0 | 4.43 |
| 63.4 | 194 / 196 | +2.6 % / +1.0 % | 24.5 | 3.08 |
| 126.8 | 267 / 270 | +2.3 % / +1.9 % | 23.9 | 2.12 |
| 253.7 | 366 / 372 | +3.4 % / +3.3 % | 23.2 | 1.45 |

`NDF/√(ka)` is constant to ±4 % over 8× in `ka` (and the small residual
growth at large `ka`, where the curves are least converged, closes the gap
further), while `NDF/ka` falls 3.1×. Raw log-log slope 0.464, biased low by
the residual growth at the top end. A star scales as `√(ka)`, not `ka`.

---

## 4. Stage 3 — perimeter versus interior saturation

### 4.1 The mechanism: family exhaustion (proved in outline)

All rows of any mesh whose nodes lie in a set `S` belong to the fixed
continuum family `F_S = {w_{p,q} : p, q ∈ S}`. If every member of `F_S` lies
within relative `L²` distance `η` of a fixed `n`-dimensional space (an
`η`-spanning space for the family), the Stage-1 counting mechanism bounds
`NDF ≤ n + O(1)` **for every `N`**: NDF must saturate in `N`, at a level set
by the family, not the node count.

- `S = ∂Ω` (perimeter): `F_S` is a **2-parameter** family. The `L = N(N−1)/2`
  links sample it at density `~N` per parameter; an `η`-net at parameter scale
  `δ` needs `N ~ P/δ`. Saturation is reached at moderate `N`.
- `S = Ω` (interior): `F_S` is a **4-parameter** family; `L ~ N²` samples must
  net a 4-manifold, so the `N` needed for a net at scale `δ` scales like
  `(1/δ)²` — quadratically further out — **and** the family itself is larger
  (interior pairs reach shorter `d`, hence narrower `σ_F` and finer transverse
  structure). Both effects push interior saturation beyond any tested `N`.

This is the content of Finding 13's observation ("interior nodes are why
random placement keeps growing"), now with a mechanism: it is not that
interior stars are individually richer, but that the interior *pair family*
has twice the parameter dimension and a finer minimal scale.

### 4.2 The ceiling: phase-space count (imported + measured constant)

Each row is `η`-concentrated in spatial support `Ω` and in the 2-D frequency
disc of radius `Ω_c = √(2 ln(1/η))/σ_F(d)` (transverse Gaussian; the axial
spectrum is narrower by `σ_F/d`). The span of any such family is governed by
the space–frequency localisation operator `P_Ω Q_{Ω_c} P_Ω`, whose trace is
`A·Ω_c²/(4π)`; with `σ_F² = λd/16`:

```
NDF_ceiling(η) ≲ A·Ω_c²/(4π) = 8·A·ln(1/η)/(π·λ·d̃) = (4/π)·(a/d̃)·ln(1/η)·ka ,
```

`d̃` a harmonic-type mean link length. This is **linear in `ka`** for any
placement. Two honesty notes: (i) converting "trace" into "count above
`δ ~ τ²`" needs the two-sided plunge estimates for concentration operators
(1-D: Landau–Widom; 2-D: recent concentration-operator spectral estimates,
e.g. Marceca–Romero–Speckbacher); we **import** that step, we do not prove it.
(ii) With `η = τ`, `d̃ ≈ 2.9` m the prefactor is `(4/π)(0.87)(4.6) ≈ 5.1`,
versus the measured bulk coefficient `α ≈ 1.2` (§4.3): the bound holds with
factor ≈ 4 slack.

### 4.3 The two-term law and the 0.696 exponent (measured + out-of-sample)

A phase-space count over a bounded region carries a boundary correction — a
Weyl-type second term. For perimeter placement the natural two-term law is

```
NDF_perim(ka) = α·ka + β·√(ka) ,
```

with the `√(ka)` term counting boundary-layer modes: the beams' transverse
scale `σ_F ∝ √λ` against the fixed perimeter length gives
`P/σ̄_F = 4P/√(λd̄) = (4P/√(2π·d̄·a))·√(ka) ≈ 10·√(ka)` for this room —
the right magnitude and the right `λ`-scaling for the fitted `β` below.

Fitted to the three original plateau points (`ka = 31.7/63.4/126.8` →
`118.9/187.7/311.9`):

```
power law : NDF = 10.64 · ka^0.696            (the fit under test)
two-term  : NDF = 1.17·ka + 14.5·√(ka)        (crossover at ka = (β/α)² ≈ 154)
```

Both are 2-parameter fits and both reproduce the 3 training points to ≤ 1.7 %,
so the training range **cannot** distinguish them — but they diverge outside
it. All three measured points lie *below* the two-term crossover
`ka ≈ 154`, where the boundary term still dominates
(at `ka = 31.7` it contributes 69 % of the total): a genuine
intermediate-exponent regime that mimics `ka^0.7`.

**Out-of-sample test** (new plateau measurements at `ka = 15.9` and
`ka = 253.7`, fits frozen on the original 3 points; both new curves plateaued
cleanly — at `ka = 253.7` the last four values are 545, 543, 543, 543):

| ka | measured | power (3-pt fit) | err | two-term (3-pt fit) | err |
|---|---|---|---|---|---|
| 15.9 | 79.8 | 72.8 | −8.7 % | 75.5 | −5.3 % |
| 253.7 | 543.5 | 500.9 | −7.8 % | 529.3 | −2.6 % |

Joint 5-point fits: power `ka^0.693`, rel-RMS 4.3 % with a systematic
inverted-U residual pattern (−5.1, +2.9, +5.4, +2.6, −4.8 %); two-term
`1.28·ka + 13.7·√(ka)`, rel-RMS 3.0 % with residuals
(−6.5, −1.3, +1.1, +1.3, −0.4 %) — only the `ka = 15.9` point is off, where
the model itself degenerates (`σ_F ≈ 0.43` m is no longer ≪ room size).

The sharpest diagnostic is the **local log-log slope between consecutive
plateau points**, which a power law requires to be constant:

```
ka 15.9→31.7 : 0.575     31.7→63.4 : 0.659     63.4→126.8 : 0.733     126.8→253.7 : 0.801
```

Monotonically rising over a 16× range, heading toward 1 — incompatible with
any fixed exponent, and exactly the crossover behaviour of
`α·ka + β·√(ka)` (predicted local slope `1 − β√(ka)/(2(α·ka+β√(ka)))`,
which rises through these values around the crossover `ka ≈ 140`).

**Verdict: 0.696 is a finite-size crossover artefact, not a genuine
exponent.** All three original points lie below the crossover, where the
`√(ka)` boundary term still carries 50–70 % of the total. The perimeter
plateau is asymptotically `Θ(ka)` with bulk coefficient `α ≈ 1.3` at
`τ = 0.01`, approached from below; the fitted boundary coefficient
`β ≈ 13.7` agrees with the `P/σ̄_F ≈ 10·√(ka)` boundary-layer count to
within 1.4×.

### 4.4 What changes for interior placement

Random interior placement at `N = 90` shows `NDF/ka` still rising
(+0.9 … +2.0 % per step). Under §4.1–4.2 this is *late saturation under a
higher ceiling*, not unbounded growth: the prediction is a plateau at
`≲ (4/π)(a/d̃_int)·ln(1/τ)·ka` with `d̃_int < d̃_perim` (interior pairs are
shorter), reached only at `N` large enough for `N²` links to net a
4-parameter family. The E3 star curves (§3.5) show the same phenomenon in
miniature: slow approach to a family ceiling from below, with the shortfall
largest at large `ka`. **Conjecture** (testable): interior placement plateaus
at `NDF/ka ≈ 6–8` for this room at `τ = 0.01`, at `N` of order several
hundred.

---

## 5. Stage 4 — the whole mesh

The union bound `rank(W) ≤ Σ_i rank(S_i)` is vacuous because stars overlap
almost completely: a **single** central star reaches 268–369 (E3,
`ka = 126.8–253.7`) while the *entire* perimeter mesh plateaus at 312–543 at
the same `ka`. Summing `N` per-star bounds counts that shared span `N` times
(~7× above measurement, as the union-bound arithmetic predicts). The
star-to-mesh plateau ratio falls monotonically with `ka`
(1.18, 1.04, 0.86, 0.68 at the four sizes) — the direct signature of
`√(ka)` per-star against `ka` per-mesh scaling.

What replaces it is the phase-space picture of §4.2, which bounds the span of
*all* rows at once and does not double count:

```
NDF ≤ min( L ,  C(τ)·ka·(a/d̃)·log-terms )        [ceiling: imported + measured]
```

with the transition between the regimes governed by family sampling (§4.1),
not by per-node saturation. A geometric reading consistent with all
measurements: a star's rows occupy the phase-space slice "positions × the one
direction through `p_i`", of volume `∝ √(ka)`; the mesh fills the direction
circle at each point once `N ≳ √(ka)` and caps at `∝ ka`. This also explains
why mesh NDF grows like `N·√(ka) = √2·√(L·ka)` *before* the cap — the
empirical `0.905·√(L·ka)` law — as per-star `√(ka)` contributions that have
not yet saturated the direction circle. **Caution**: the K-sweep found the
NDF/L knee wavelength-independent, which this reading does not yet explain;
the crossover location remains underived (as in the refuted document).

---

## 6. Honest assessment

| stage | status |
|---|---|
| 1. Counting lemma | **Proved** (self-contained, explicit constants), verified 0/240; slack: `log(T/s)` inside logs; tightness factor ~1.3–2.5 |
| 2. Angular lemma | **Proved** (contour shift + Bessel; no Slepian); factor ≈ 2 above measured `n_ang` |
| 2. Product theorem | **Proved conditional** on radial count `ρ̄`; `ρ̄` **measured** = 13, `λ`-independent (predicted); analytic radial bound **open** |
| 2. Star scaling | **Measured**: `≈ C_star(τ)·√(ka)`, `C_star(0.01) ≈ 24`; the prior `O(ka)` belief is **refuted for this phaseless model** — coherent-kernel intuition does not transfer |
| 3. Perimeter law | Two-term `1.28·ka + 13.7·√(ka)`, both constants at predicted magnitudes; out-of-sample at `ka = 15.9, 253.7` beats the power law at both ends; local slope rises 0.575→0.801 ⇒ **0.696 is a crossover artefact** |
| 3. Interior | Mechanism (4-parameter family + finer scales) argued; plateau level **conjecture** |
| 4. Mesh | Union bound replaced by phase-space ceiling (imported plunge step); knee location still underived |

Known gaps: (i) `log(T/s)` slack in Lemma 1; (ii) no analytic radial-family
bound (`ρ̄ = 13` is measured); (iii) the 2-D plunge estimate is imported;
(iv) the `e^{−4.5}` truncation sits at `1.1τ` and caps the precision any
`τ = 0.01` statement can have; (v) everything is 2-D.

---

## 7. Related

- [NDF_ANGULAR_SATURATION_DERIVATION](NDF_ANGULAR_SATURATION_DERIVATION.md) — the refuted derivation (record)
- [MESH_SENSING_NDF_FINDINGS](../4-research/in-progress/MESH_SENSING_NDF_FINDINGS.md) — Findings 12–13
- `docs/4-research/results/spatial_temporal_equilibrium_findings.md` — Finding 2 (star `NDF = N−1` at `N = 200`)
- H. Widom, *Asymptotic behavior of the eigenvalues of certain integral equations II* (1964) — Gaussian-kernel eigenvalue asymptotics on an interval
- H. Landau & H. Widom, *Eigenvalue distribution of time and frequency limiting* (1980)
- Santin & Schaback (2016); Belkin (2018) — RBF kernel eigenvalue decay
- Marceca, Romero, Speckbacher — spectral estimates for concentration operators (2-D plunge)
- Bucci & Franceschetti (1989) — `O(ka)` NDF for *coherent* scattered fields (contrast with the phaseless model here)
