# zeta-explorer

[![tests](https://github.com/louisbelot/zeta-explorer/actions/workflows/tests.yml/badge.svg)](https://github.com/louisbelot/zeta-explorer/actions/workflows/tests.yml)

![Phase portrait of zeta over the critical strip](docs/domain.png)

The Riemann Hypothesis says every non-trivial zero of ζ(s) sits on the line
Re(s) = 1/2. On that line ζ is complex, but multiply it by the phase e^(iθ(t)) and
you get the real-valued Riemann-Siegel Z function, which changes sign at every zero.
This repo computes θ, Z, ζ(s) itself, the zeros and the Gram points straight from
the formulas, checks the results against mpmath, and uses Turing's method to prove
the zero count is complete: all 138,069 zeros of ζ with 0 < Im s ≤ 100000 are on
the critical line and simple. It also draws ζ over the complex plane, compares
100,000 zeros near t = 5·10⁶ with the eigenvalue statistics of large random
Hermitian matrices, rebuilds π(x) from the first 1000 zeros, hunts for Lehmer
pairs (the near misses that bound the de Bruijn-Newman constant), and checks the
value distribution of log|ζ(1/2 + it)| against random unitary matrices. A
simplified Odlyzko-Schönhage evaluator (one nonuniform FFT per window, then
band-limited interpolation) makes Z cost 2.6 μs per point at any height and
reproduces Odlyzko's published zeros near the 10¹²th zero to 3e-9.

The same machinery extends to Dirichlet L-functions: characters with LMFDB's
Conrey labels, L(s, χ) anywhere in the plane, zeros checked against LMFDB and
verified complete with the argument principle, the class number formula, a
smoothed approximate functional equation that needs O(√(q|t|)) terms instead of
O(q|t|), and the low-lying zeros of quadratic characters up to conductor 10⁵.
Those zeros are compared with the symplectic prediction of Katz and Sarnak, with
the lower-order terms of the ratios conjecture (which fit the histograms with
χ² of 16 to 21 on 30 bins, against 298 to 719 for the limit law), and, for the lowest
zero, with gap probabilities of random symplectic matrices computed as Fredholm
determinants.

Elliptic curves come next: L-functions built from point counts, central zeros of
order 0 to 3 for the curves 11a, 37a, 389a and 5077a, and 1,674 quadratic twists
of 11a. Their central values sit exactly on Waldspurger's curves κc²/√d, and that
discreteness explains why their lowest zeros stay away from the center: SO(16)
random matrices with the matching cut on det(I − A) fit the even twists with χ² of
34.6 on 15 bins, against 159 to 182 without the cut. The odd twists repel their
first zero too, but the analogous cut on the central derivative, sized by Gross
and Zagier's height formula and the CFKRS arithmetic factor, removes 0.13% of the
matrices and changes nothing. Fitting the cut instead lands at 50 times that size
and still misses.

Heat flow came next. H_t(z), the heat-flow deformation of ξ behind the
de Bruijn-Newman constant, is computed from its defining integral on a shifted
contour, which keeps double precision at heights where H is 10⁻¹⁹¹. Running the
flow backwards, the close pair of zeros at γ = 564.161 and 564.506 collides at
t = −0.0614041 and leaves the real line, 2·10⁻⁴ above the bound Csordas, Smith and
Varga's criterion gives for that pair. The zeros move like charged particles: the
velocity of each one, predicted from the positions of 4,000 others, matches the
measured motion to 7·10⁻⁷.

The zeros also settle prime races. Under GRH and LI the share of log x where
primes ≡ 3 (mod 4) lead primes ≡ 1 is 0.9959279, computed from 868 zeros of
L(s, χ₋₄) and an exact tail, and the same code reproduces Fiorilli and Martin's
densities for twelve races mod 3 to 24 to six digits. With 3,130 zeros the
explicit formula rebuilds π(x; 4, 3) − π(x; 4, 1) closely enough to see the first
lead change at x = 26861.

Moments and positivity close the loop on the value distribution. The arithmetic
factor of the moment conjecture, read as a correction to the cumulants of
log|ζ|, accounts for the whole gap between the data and random unitary matrices
(variance within 0.015 at five heights, nothing fitted). Integrals of |ζ|^(2k) up
to t = 10¹⁰ follow the full CFKRS polynomials, while their leading term is still
off by a factor 11 for the sixth moment, and recomputing the paper's own tables
turned up a swapped pair of numbers in one row. Li's coefficients λ_n, summed over
the first million zeros and checked against Taylor coefficients of log ξ to
1e−12, are positive for every n tried up to 10⁵.

The extremes of ζ behave like those of a log-correlated field. Over 50,000
windows of length 2π at each height up to 10¹², the maximum of log|ζ| sits far
from the Fyodorov-Hiary-Keating limit law, and so do random unitary matrices of
the same size, even at N = 400. What separates ζ from them is arithmetic: the
window average of the sum over small primes has variance 0.043, and adding it to
the random-matrix variance matches the data to 0.004 from 10⁶ up, with nothing
fitted. Every Jensen polynomial of ξ with degree up to 12 and shift up to 1000
(and degree 40 at a few shifts) is decided hyperbolic in exact rational
arithmetic, as RH requires, and their zeros approach the Hermite limit of
Griffin, Ono, Rolen and Zagier at rate exactly δ(n).

The same numerics also run in the browser at
[louisbelot.github.io/zeta-explorer](https://louisbelot.github.io/zeta-explorer/):
a phase portrait you can pan and zoom, a Z(t) explorer that marks zeros and Gram
points up to t = 10⁹, the explicit formula building ψ(x) one zero at a time, and
Z(t, χ) with its zeros for every primitive Dirichlet character up to modulus 200.

## Install

```bash
pip install "zeta-explorer[plot] @ git+https://github.com/louisbelot/zeta-explorer"
```

## Quick start

```python
from zeta_explorer.riemann_siegel import theta, z_rs, z_exact
from zeta_explorer.zeros import find_zeros, gram_point

find_zeros(10, 30)   # [14.134725141734695, 21.022039638771552, 25.01085758014569]
len(find_zeros(10, 100))   # 29
gram_point(0)        # 17.845599540410863, where theta(t) = 0
z_rs(1000.0)         # 0.9977946375478818 (mpmath: 0.9977946375215866)
z_rs(1000.0, corrections=0)   # 0.9976696568809473, C0 only

from zeta_explorer.zeta import zeta
zeta(-1)             # -1/12
zeta(2 + 3j)         # (0.7980219851462759-0.1137443080529385j)

from zeta_explorer.turing import verify_zeros
verify_zeros(10_000)
# Verification(n=10143, g_n=10000.881777216977, zeros=10144, blocks_after=2)

import numpy as np
zeta(np.array([-3, 0.5 + 14.134725141734693j]))   # vectorized: [1/120, ~1e-16]

from zeta_explorer.spacing import zeros_between, normalized_spacings
z = zeros_between(10**6, 10**6 + 5)   # zeros between Gram points g_1000000 and g_1000005
# [600270.74787059, 600271.48637367, 600271.76148043, 600272.00103217, 600273.02653325]

from zeta_explorer.multieval import zeros_fast
win, u = zeros_fast(10**12 - 2, 10**12 + 4)   # zeros as offsets u from win.t0
win.t0 - 267653395647 + u
# [1.62594824, 1.84752313, 2.36236697, 2.68163092, 2.86198994]
# Odlyzko's table: zero number 10^12 + 1 is at 267653395647 + 1.8475231278

from zeta_explorer.spacing import first_zeros
from zeta_explorer.explicit import pi_from_zeros, riemann_r
gamma = first_zeros(1000)  # the first 1000 zeros, 14.1347 to 1419.4225
pi_from_zeros(np.array([100.5, 10_000.5]), gamma)   # [25.026, 1229.54]; π = 25, 1229
riemann_r(1e15)      # 29844570495886.9

from zeta_explorer.lehmer import lehmer_pairs
pairs = lehmer_pairs(gamma[:700])   # pairs with Δ²g < 4/5
pairs[0]   # LehmerPair(n=34, gamma=111.0295..., delta=0.8451..., gap=0.3868..., score=0.5699..., bound=-0.3954...)

from zeta_explorer.plot import plot_z, plot_domain, plot_spiral, plot_spacing
plot_z(0, 60).savefig("z.png")
plot_domain().savefig("domain.png")

from zeta_explorer.characters import conrey, quadratic, class_number
from zeta_explorer.lfunction import dirichlet_l, l_zeros, verify_zeros
chi = conrey(5, 2)   # LMFDB label 5.2: order 4, odd, conductor 5
chi.root_number()    # (0.85065080835204+0.5257311121191335j)
dirichlet_l(0.5 + 10j, chi)   # (2.1249968234507937+2.163859185370417j)
l_zeros(chi, 0, 13)  # [6.1835782, 8.45722917, 12.67494642]
dirichlet_l(1, quadratic(-23)), class_number(-23)   # 1.9652 = 3π/√23, 3
verify_zeros(conrey(101, 2), 50).ok   # 90 zeros on the line, 90 in the box

from zeta_explorer.afe import l_afe
big = conrey(100003, 2)
l_afe(0.5 + 10j, big)   # (3.9151561336603864-8.036122558438148j), 1557 terms per sum
l_zeros(big, 0, 1.5)    # [0.60517956, 1.49288041]

from zeta_explorer.ratios import quadratic_density
quadratic_density([-100003, 100001], [0.1, 0.5, 1.0])   # [0.1216, 1.4173, 0.8240]

from zeta_explorer.rmt import gap, first_density, symplectic, symplectic_n
gap(symplectic, 0.5), gap(symplectic_n(2), 0.5)   # 0.79479, 0.65935: P(no angle in (0, 1/2))
```

## The math

### Theta

From Stirling's series for log Γ:

```
θ(t) = t/2·log(t/2π) − t/2 − π/8 + 1/(48t) + 7/(5760t³) + …
```

accurate to about 1e-13 for t ≥ 10.

### Riemann-Siegel formula

With a = sqrt(t/2π), N = ⌊a⌋ and p = a − N:

```
Z(t) ≈ 2·Σ_{n=1..N} cos(θ(t) − t·log n)/sqrt(n)  +  (−1)^(N−1)·a^(−1/2)·Σ_{k=0..4} C_k(p)·a^(−k)
```

Only about sqrt(t/2π) terms are needed, so at t = 10⁶ that's 398 terms. The
corrections come from Ψ(p) = cos(2π(p² − p − 1/16)) / cos(2πp) and its
derivatives:

```
C0 = Ψ
C1 = −Ψ⁽³⁾/(96π²)
C2 = Ψ⁽²⁾/(64π²) + Ψ⁽⁶⁾/(18432π⁴)
C3 = −Ψ⁽¹⁾/(64π²) − Ψ⁽⁵⁾/(3840π⁴) − Ψ⁽⁹⁾/(5308416π⁶)
C4 = Ψ/(128π²) + 19Ψ⁽⁴⁾/(24576π⁴) + 11Ψ⁽⁸⁾/(5898240π⁶) + Ψ⁽¹²⁾/(2038431744π⁸)
```

Twelve derivatives of a ratio of cosines is not something to do by hand. Ψ looks
singular at p = 1/4 and 3/4, but the numerator vanishes there too and Ψ is in fact
entire. So the code expands Ψ(1/2 + z) = −cos(2πz² − 5π/8)/cos(2πz) as a Taylor
series in z by dividing the two cosine series. That division runs at 110 digits
in mpmath, because 1/cos has poles at |z| = 1/4 and rounding errors grow like 4^j.
Derivatives of a power series are shifted coefficients, so each C_k becomes a
degree-90 polynomial in z, evaluated with Horner in doubles, and the removable
singularities need no special case. The tests rebuild C1 to C4 from mpmath's
numerical derivatives of Ψ and agree to 1e-15.

![Riemann-Siegel error against t for 0 to 4 correction terms](docs/rs_error.png)

Each extra term should buy another factor of t^(−1/2). Over t in [200, 2000], the
largest value of |error|·t^((2k+3)/4) is:

| terms | error order | measured constant |
|---|---|---|
| C0 | t^(−3/4) | 0.119 |
| C0 to C1 | t^(−5/4) | 0.052 |
| C0 to C2 | t^(−7/4) | 0.0090 |
| C0 to C3 | t^(−9/4) | 0.029 |
| C0 to C4 | t^(−11/4) | 0.015 |

The constants come out about the same over [20, 200], which is the real check that
each C_k is right: a wrong coefficient would leave the previous order behind.
The bottom right of the chart is a different story. Past t ≈ 3000 the C3 and C4
curves hit a floor near 3e-11 set by doubles: the phases θ(t) − t·log n are around
10⁵ radians there, so each cosine is only good to about 1e-11. Gabcke (1979)
proved explicit bounds of this form for t ≥ 200, which is why the code trusts
z_rs from t = 200 up and calls mpmath below.

### ζ(s) everywhere

`zeta.py` has no mpmath in the loop. For Re s ≥ 0 it uses Euler-Maclaurin:

```
ζ(s) = Σ_{n<N} n^(−s) + N^(1−s)/(s−1) + N^(−s)/2
       + Σ_{k=1..M} B_2k/(2k)! · s(s+1)…(s+2k−2) · N^(−s−2k+1) + R
```

with Backlund's bound: |R| is at most |s+2M+1|/(σ+2M+1) times the first omitted
term, for σ > −(2M+1). N = 10 + ⌊|t|/π⌋ makes the Bernoulli terms shrink
geometrically, and M grows until the bound drops below 1e-16 relative. The
Bernoulli numbers are exact `Fraction`s. For Re s < 0 it switches to the
functional equation ζ(s) = χ(s)ζ(1 − s), with
χ(s) = 2^s·π^(s−1)·sin(πs/2)·Γ(1 − s) assembled in log space, since at large |t|
the sine overflows and Γ underflows. log Γ is Stirling's series after shifting
|z| past 15, summing the logs of the shifts so the imaginary part stays on the
continuous branch.

Against mpmath on 300 random points with |Re s| ≤ 30 and |Im s| ≤ 200 the
relative error stays under 1e-12. On the critical line at t = 10⁴ the absolute
error is 1.4e-11 (the same phase floor as above). The tests also check
ζ(1 − 2k) = −B_2k/(2k) for k up to 11, that Backlund's bound really bounds the
error and is within 50x of it, that χ(1/2 + it) = e^(−2iθ(t)), and that
e^(iθ)ζ(1/2 + it) is real and equals the Riemann-Siegel Z.

### Phase portrait

The picture at the top is `plot_domain()`: ζ(σ + it) for −5 ≤ σ ≤ 6 and
0 ≤ t ≤ 60, about 600,000 points. Hue is the argument of ζ, and the brightness
ramps up and resets every time |ζ| doubles, so the shading edges are level curves
of |ζ| (Wegert's enhanced phase portrait). A zero is a point where every color
meets, going around once counterclockwise; the pole at s = 1 goes the other way.
All the zeros in the strip sit on the dashed line σ = 1/2, the trivial zeros at
−2 and −4 sit on the left edge, and to the right of σ = 2 the picture is almost
flat red because ζ(s) is close to 1 there.

This needed ζ on arrays. `zeta`, `chi` and `loggamma` now take any array shape:
every point keeps its own N = 10 + ⌊|t|/π⌋, the main sums are one masked
(points × n) matrix built in chunks, and the Euler-Maclaurin tail drops each point
out of the loop as soon as its Backlund bound is small enough. The whole grid
takes 1.05 s, against 6.5 s for a Python loop over the old scalar version.

### Landscape and a sweep across the strip

![log |zeta| as a surface over the critical strip](docs/landscape.png)

`plot_landscape()` draws log₁₀|ζ(σ + it)| over −0.5 ≤ σ ≤ 2 and 0 ≤ t ≤ 40,
clipped to [−1.5, 1]. The pole at s = 1 is the plateau at the back left, and each
zero is a narrow pit on σ = 1/2. The t grid includes the zeros themselves, so each
pit reaches the floor.

![zeta(sigma + it) as sigma sweeps from 0.1 to 0.9 and back](docs/sweep.gif)

`sweep_sigma()` moves σ from 0.1 to 0.9 and back as σ = 0.5 + 0.4·sin³(2πk/60),
so it lingers near the critical line. Off the line the curve winds around the
origin without touching it and |ζ| stays above 0. At σ = 1/2 both panels touch
zero at the same six heights.

### Analytic continuation

![Partial sums of three series for zeta at the first zero](docs/continuation.png)

Σ n^(−s) only converges for σ > 1, and at the first zero s = 1/2 + 14.1347i the
partial sums spiral outward like N^(1/2)/|s − 1|. Two cheap fixes both land on
ζ(s) = 0. The alternating series η(s) = Σ (−1)^(n−1) n^(−s) converges for σ > 0,
and ζ = η/(1 − 2^(1−s)), but its error only falls like N^(−1/2). Adding the
first two Euler-Maclaurin terms N^(1−s)/(s−1) + N^(−s)/2 to the raw partial sum
gives an error of order |s|·N^(−σ−1), N^(−3/2) on the critical line, and it works
anywhere with σ > −1. The tests check both rates.

The functional equation is easiest to see through
ξ(s) = ½s(s − 1)π^(−s/2)Γ(s/2)ζ(s), which cancels the pole and the trivial zeros.
`xi(s)` evaluates it from the right half plane and reflects with ξ(s) = ξ(1 − s);
the tests check it against mpmath, check ξ(0) = ξ(1) = 1/2 and ξ(2) = π/6, and
check that it is real on the critical line.

![zeta(1/2 + it) traced in the complex plane for t up to 40](docs/spiral.png)

On the critical line itself, ζ(1/2 + it) traces loops that pass through the
origin once per zero. It starts at ζ(1/2) = −1.4604, and between t = 0 and 40 it
goes through 0 six times.

### Zeros

![Z(t) on the critical line with its first 13 zeros](docs/z_function.png)

Sign changes of Z on a grid, refined with the Illinois method
(regula falsi that halves a stuck endpoint's weight). All 29 zeros below t = 100
match mpmath's `zetazero` to 1e-9, and z_rs finds the same 29 to within 4e-6.

### Gram points and Gram's law

These solve θ(g_n) = nπ. θ is increasing and convex there, so Newton's method
started to the right of the root slides monotonically onto it. It's vectorized:
the first 100,001 Gram points take 0.12 s.

![Z(t) around t = 280 with Gram points; g_126 and g_134 break Gram's law](docs/gram_points.png)

At g_n, e^(iθ) = (−1)^n, so if Z behaved like a cosine with one zero per Gram
interval we would always have (−1)^n Z(g_n) > 0. That's Gram's law. It first
fails at g_126 (Hutchinson, 1925). Up to n = 1000, 95.8% of Gram points are good;
up to 10,000, 91.8%; up to 100,000 (t ≈ 74,922), 88.9%.

A Gram block is a run [g_j, g_k) between consecutive good points. Rosser's rule
says a block of length k − j contains at least k − j zeros. The checker samples Z
inside each block, doubling the density until it sees enough sign changes. Up to
n = 100,000 every block passes:

| block length | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| count | 78,694 | 9,445 | 778 | 19 | 1 |

Rosser's rule does fail eventually; Brent found the first failure near
g_13999525. The signs come from z_rs with all four corrections. Any |Z| below 1e-6
is recomputed with mpmath before it's trusted, and below t = 200 mpmath is used
throughout.

### Turing's method

Sign changes of Z give a lower bound on N(T), the number of zeros with
0 < Im s ≤ T. A zero off the critical line, or a double zero, wouldn't show up as
a sign change. The upper bound comes from S(t) = N(t) − θ(t)/π − 1: Lehman proved
|∫ S(t) dt| ≤ 1.7 + 0.114·log(t₂/2π) over any [t₁, t₂] above 168π, so S can't stay
positive for long. Brent turned that into a block count. If K consecutive Gram
blocks [g_n, g_p) obey Rosser's rule and

```
K ≥ (0.114/6π)·log² g_p + ((1.7 − 0.114·log 2π)/6π)·log g_p   (≈ 0.0061 log² + 0.08 log)
```

then N(g_n) ≤ n + 1. `verify_zeros(T)` finds a good Gram point g_n ≥ T, counts
n + 1 sign changes of Z below it, and finds the K blocks past it. The two bounds
meet, so every zero up to g_n is a sign change, which means it is on the line and
simple. For T = 100,000 it runs in 9 seconds, needs K = 2 blocks, and stops at
g_138068 ≈ 100000.29 with N = 138,069 (mpmath's `nzeros` agrees). The tests
cross-check smaller heights the same way.

### Many values of Z at once

![Accuracy and cost of Z from the direct sum and from a NUFFT window](docs/multieval.png)

Every value of Z costs about √(t/2π) terms: 399 at t = 10⁶, 206,393 near the
10¹²th zero. Odlyzko and Schönhage (1988) noticed that the main sum
F(t) = Σ n^(−1/2) e^(−it log n) is the same Dirichlet polynomial for every t in a
window, so its values on an evenly spaced grid can come out of one FFT and
everything in between can be interpolated. `ZWindow(t_lo, t_hi)` is a simplified
version of that.

The grid comes first. With a center t₀ and a spacing δ,
F(t₀ + jδ) = Σ b_n e^(−ij·δ log n) where b_n = n^(−1/2) e^(−it₀ log n). That is a
type-1 nonuniform FFT: arbitrary frequencies δ log n, integer outputs j. `nufft`
does it by Gaussian gridding (Dutt-Rokhlin, Greengard-Lee): spread each term onto
a 2x oversampled grid with a periodic Gaussian, run one FFT, and divide out the
Gaussian's Fourier coefficients e^(−τk²). With w grid points on each side and
τ = πw/(3m²), aliasing costs e^(−2πw/3) and cutting the Gaussian off costs
e^(−3πw/4), both relative to Σ|b_n|, so w = 15 is enough for 1e-13. The tests
compare it with the direct sum at three sizes and three tolerances.

Then interpolation. Multiplied by e^(iβt) with β = ½ log N, F has its spectrum in
[−β, β]. The grid samples at δ = π/α with α = 2β, twice the Nyquist rate. A sinc
of bandwidth α times a Gaussian e^(−ε²u²/2) rebuilds F from the K nearest samples
on each side. Two errors compete: the Gaussian smears the sinc's spectrum, which
lets in e^(−(α − β)²/2ε²) from the aliased copies, and stopping the sum at
distance L = Kδ loses e^(−ε²L²/2). Choosing ε² = (α − β)/L makes both
e^(−(α − β)L/2) = e^(−πK/4), so K = 36 gives 1e-12, and every Z(t) costs 73
multiply-adds whatever t is.

Precision turned out to matter as much as speed. At t = 10¹², t·log n is about
10¹³ radians, and in doubles that has an absolute error near 10⁻³. That is why the
direct sum is off by 3.5e-3 up there (left panel). The window only needs
t₀ log n mod 2π and θ(t₀) once, so those come from mpmath at about 36 digits, and
everything else works with the offset u = t − t₀:

```
θ(t₀ + u) − θ(t₀) = (u/2)·log(t₀/2π) + (t/2)·log1p(u/t₀) − u/2 + (Stirling tail difference)
```

All of those terms are moderate in size. Zeros and Gram points come back as
offsets too, because a double near 2.7·10¹¹ only resolves 6e-5. N goes up by one
each time t passes 2πk², so the grid uses the largest N in the window and the
extra term is subtracted directly for points below the jump.

On this machine, once the window is built, each Z costs 2.6 μs at any height. The
direct sum costs 7.4 μs at t = 10⁶, 795 μs at 10¹⁰ and 7.7 ms at 10¹². Building the
window takes 0.25 s at 10¹⁰ and about 1.3 s at 10¹², mostly for the mpmath phases.
Against mpmath's `siegelz`, the window's error stays between 8e-14 and 3e-13 from
t = 10⁴ to 10¹², while the direct sum's grows from 1e-11 to 3.5e-3.

The check that convinced me was Odlyzko's published table of zeros number
10¹² + 1 to 10¹² + 10⁴. `zeros_fast(10**12 - 2, 10**12 + 10**4 + 2)` finds all of
them in 2.5 s with the same numbering, and the largest difference from his values
is 2.8e-9 (RMS 7e-10; his file promises 1e-8).

Up there Gram blocks start breaking Rosser's rule. The block
[g_(10¹²+13452), g_(10¹²+13454)) has no zeros at all: Z reaches 98 inside it, and
its two zeros sit in the previous block, which is a single Gram interval holding
three. The block finder now deals with that. When a block still comes up short
after sampling it 1024 times more finely, it is merged with the block before it,
then with the blocks after it, until the union has as many sign changes as Gram
intervals. In 100,000 Gram intervals from n = 10¹² that happened 4 times; from
n = 10⁹, never. 20.1% of the Gram points near 10¹² are bad, against 17.5% near
10⁹.

### Zero spacings and random matrices

![Spacing histogram and pair correlation of 100,000 zeros near t = 5 million](docs/spacing.png)

Montgomery conjectured in 1973 that the pair correlation of the zeros, measured in
units of the mean spacing, is 1 − (sin πu/πu)². Dyson recognized it as the pair
correlation of eigenvalues of large random Hermitian matrices (GUE), and Odlyzko's
computations from 1987 on confirmed it numerically.

`zeros_between(n_lo, n_hi)` finds every zero between two good Gram points. It
evaluates Z at the Gram points, samples each Gram block more finely until it has as
many sign changes as the block has intervals, then refines every bracket at once
with a vectorized Illinois method on z_rs. Between n = 10⁷ and 10⁷ + 100,000 it
finds exactly 100,000 zeros in 18 s. That count matters for the statistics: a
missed pair of close zeros would hide exactly the small spacings the plot is
about. The zeros are unfolded with x = θ(γ)/π, which makes the mean spacing 1.

The GUE curve is the exact sine-kernel answer, not the Wigner surmise. The
probability of no eigenvalue in an interval of length s is the Fredholm
determinant E(s) = det(I − K) with K(x, y) = sin π(x − y)/(π(x − y)) on [0, s].
Following Bornemann (2010), Gauss-Legendre quadrature with 40 nodes turns it into
an ordinary 40×40 determinant, and the error falls exponentially with the node
count (20 and 60 nodes agree to 1e-12). The spacing density is p(s) = E''(s); E
is entire, so the code interpolates it with a degree-80 Chebyshev series on
[0, 6] and differentiates the series twice. Checks: p(s) matches
π²s²/3 − 2π⁴s⁴/45 + π⁶s⁶/315 for small s, the derivatives agree with 30-digit
mpmath to 3e-10, p integrates to 1 with mean 1, and the variance comes out
0.17999 (the Wigner surmise gives 3π/8 − 1 = 0.1781). An independent check
samples 300 random 120×120 GUE matrices: the KS distance between their bulk
spacings and the sine-kernel CDF is 0.0044, against 0.28 for Poisson.

The zeros fit GUE well, but not perfectly, and the gap shrinks with height. For
100,000 consecutive zeros starting at five heights (the last two from
`zeros_fast`, which takes 11 s for the n = 10⁷ row against 18 s for
`zeros_between`, and gives the same numbers; at 10⁹ and 10¹² the range between
the outermost good Gram points holds 99,997 and 99,998 zeros):

| first zero | t range | spacing variance | KS distance to GUE | R₂ − Montgomery, 0.3 < u < 0.9 |
|---|---|---|---|---|
| n = 10⁵ | 74,924 to 139,503 | 0.1636 | 0.0174 | −0.012 |
| n = 10⁶ | 600,271 to 654,852 | 0.1668 | 0.0138 | −0.013 |
| n = 10⁷ | 4,992,382 to 5,038,615 | 0.1693 | 0.0110 | −0.013 |
| n = 10⁹ | 371,870,206 to 371,905,313 | 0.1725 | 0.0081 | −0.010 |
| n = 10¹² | 267,653,395,649 to 267,653,421,320 | 0.1763 | 0.0045 | −0.007 |
| GUE | | 0.1800 | 0 | 0 |

The zeros are more regular than GUE at these heights: smaller variance, a taller
peak in the histogram, and slightly fewer close pairs. The variance and the KS
distance move toward GUE as t grows. The dip in the pair correlation barely moves
up to 10⁷, then halves by 10¹². This is the finite-height
effect Bogomolny, Keating and coauthors explained by matching the zeros at height
t to unitary matrices of effective size log(t/2π)/√(12Λ), where Λ ≈ 1.573, which is
only about 3 at t = 5·10⁶. It decays like 1/log² t, so it is still visible in
Odlyzko's data near zero number 10²³.

### The explicit formula

![pi(x) and psi(x) rebuilt from 10, 100 and 1000 zeros](docs/explicit.png)

Riemann's 1859 formula writes the prime counting function as a smooth main term
minus a sum over the zeros. The cleanest version is for Chebyshev's
ψ(x) = Σ_{p^k ≤ x} log p:

```
ψ(x) = x − Σ_ρ x^ρ/ρ − log 2π − ½·log(1 − x⁻²)
```

For π(x) it goes through J(x) = Σ_k π(x^(1/k))/k, which jumps by 1/k at every
prime power p^k:

```
J(x) = li(x) − Σ_ρ li(x^ρ) − log 2 + ∫_x^∞ dt / (t(t² − 1) log t)
π(x) = Σ_n μ(n)/n · J(x^(1/n))
```

The sum over n stops once x^(1/n) < 2, since J is 0 there. Both sums over ρ pair
each zero with its conjugate and run in order of height. At a jump both formulas
converge to the midpoint, so the code gives ψ(97) − ½·log 97 for ψ(97) and 24.5
for π(97).

The awkward part is li(x^ρ). Taking li of the complex number x^ρ lands on the
wrong branch as soon as γ·log x passes π, so it has to be Ei(ρ·log x). `ei(w)` is
written from scratch with three regimes: the power series for |w| < 4, the
continued fraction for E₁ (backward recurrence, depth 20 to 140 depending on |w|)
out to |w| = 40, and the asymptotic series beyond that. It matches mpmath to
1.5e-15 relative on 2,851 random points. The one region it refuses is
4 ≤ |w| < 40 within π/4 of the positive real axis, where neither the continued
fraction nor the series is accurate in doubles. The explicit formula never goes
there, because ρ·log x sits next to the imaginary axis. li(x) for real x, Riemann's
R(x) from the Gram series Σ (log x)^k/(k·k!·ζ(k+1)), and the tail integral
(40-node Gauss-Laguerre after substituting t = e^v) are all checked against
mpmath.

With the first 1000 zeros, on the half-integers 2.5 ≤ x ≤ 100.5:

| zeros | max error ψ | RMS error ψ | max error π | RMS error π |
|---|---|---|---|---|
| 100 | 0.96 | 0.28 | 0.21 | 0.069 |
| 1000 | 0.148 | 0.050 | 0.033 | 0.012 |

π(x) from 1000 zeros also rounds to the right integer at every half-integer from
1000.5 to 1099.5 (worst error 0.35), and π(10000.5) comes out as 1229.54 against
the true 1229. The left panel, 3000 values of x and 1000 zeros, takes 3.3 s. The
zeros come from `first_zeros(n)`, which sign-scans mpmath's Z below the good Gram
point g₇₉ ≈ 201.5 and uses `zeros_between` above it.

### Lehmer pairs

![Z(t) near t = 7005 and the two zeros 0.0377 apart](docs/lehmer.png)

In 1956 Lehmer found two zeros near t = 7005 only 0.0377 apart, with Z rising to
just 0.004 between them. Close pairs like this used to give lower bounds on the
de Bruijn-Newman constant Λ (Rodgers and Tao proved Λ ≥ 0 in 2018). In the
convention of Csordas, Smith and Varga the zeros are x = 2γ, the zeros of
H₀(x) = Ξ(x/2)/8. A pair x_k < x_(k+1) with Δ = x_(k+1) − x_k and

```
g = Σ_{j ≠ k, k+1} [ 1/(x_j − x_k)² + 1/(x_j − x_(k+1))² ]     (j over positive and negative zeros)
```

is a Lehmer pair if Δ²g < 4/5, and then

```
Λ ≥ λ = [ (1 − 5Δ²g/4)^(4/5) − 1 ] / (8g)
```

Δ²g doesn't depend on the scale, so the code works with γ directly and divides
by 2g instead of 8g. `lehmer_sums` adds the 500 nearest zeros on each side
exactly, the mirrored negative zeros exactly over the first 500, and everything
further out through the smooth density log(t/2π)/2π, integrated in closed form.
That density integrates to the main term of θ(t)/π + 1 when it starts from
T₀ = 9.6769, where the smooth count is zero.

The real check doesn't use any list of zeros. Ξ(t) is even with real zeros ±γ,
so −(log Ξ)''(t) = Σ_j [1/(t − γ_j)² + 1/(t + γ_j)²]. Removing the two singular
terms leaves g in terms of Z', Z'' and Z''' at the two zeros plus a trigamma term
from the Gamma factor, all of which mpmath computes directly. The two agree to
3e-6 relative, both on the first pair and on Lehmer's.

Results:

- Below t = 1000 there are exactly 29 Lehmer pairs among the 649 zeros, the count
  Stopple reports in "Lehmer pairs revisited" (2017).
- Lehmer's pair (zeros 6709 and 6710): Δ = 0.0376985 in γ, Δ²g = 0.006933 and
  λ = −7.112e-4. The maximum of Z between them is 0.00397 at t = 7005.0818554,
  next to the zero of ζ' that Stopple lists at 7005.08185555.
- Up to t = 10⁵, 8,262 of the 138,068 consecutive pairs (6.0%) are Lehmer pairs.
  The tightest is zeros 95,248 and 95,249 near t = 71,732.90: 0.0147 apart (0.022
  mean spacings), with Δ²g = 0.0020 and λ = −1.08e-4.
- Starting at t ≈ 10⁶, 7,384 of 114,660 consecutive pairs qualify. Stopple counts
  7,398 of 114,661 near the same height, presumably over a slightly different
  window.
- Csordas, Odlyzko, Smith and Varga (1993) give Δ = 2.171392e-4 in x units for the
  pair van de Lune, te Riele and Winter found near t = 3.9·10⁸. `newman_bound`
  turns that into −5.894e-9; their published bound is −5.895e-9. That pair is out
  of reach here: Z peaks at about 1e-8 between the two zeros, below the rounding
  floor of z_rs at that height.

### Large values and the value distribution

`z_extrema(gamma)` finds the extremum of Z between each pair of consecutive zeros
with a vectorized golden-section search (Z is unimodal between consecutive zeros
in practice). Over all 138,068 gaps below t = 10⁵ (26 s), the running record of
|Z| ends with

| t | Z(t) |
|---|---|
| 50,626.478 | 23.747 |
| 55,104.583 | −24.830 |
| 63,751.863 | −26.073 |
| 74,956.025 | −27.694 |
| 77,403.722 | 28.216 |

and each of those values matches mpmath's `siegelz` to 1e-10.

![Histogram of log|zeta| near t = 5.6e7 with the CUE density, the arithmetic correction and Selberg's normal](docs/log_z.png)

Selberg's central limit theorem says log|ζ(1/2 + it)| divided by
sqrt(½·log log T) tends to a standard normal. Keating and Snaith (2000) proposed a
much sharper model: the characteristic polynomial of a random N×N unitary matrix
(CUE) with N = log(T/2π). Its moment generating function is exact,

```
E|det(I − U)|^s = Π_{j=1..N} Γ(j)·Γ(j + s) / Γ(j + s/2)²
```

so the cumulants are κ_m = (1 − 2^(1−m))·Σ_{j≤N} ψ^(m−1)(j), and the density
comes from Fourier inversion of the product along s = iu (trapezoid rule, step
0.01). The density integrates to 1 within 1e-9, its moments match the cumulants,
and a histogram of 20,000 Haar-random 12×12 unitaries sits within 0.015 of it.

Sampling 200,000 random t in [T, 1.01·T] with T = 2π·e^N, so that N is an
integer:

| N | T | variance | κ₃ | CUE variance | CUE κ₃ | CUE + arithmetic variance | CUE + arithmetic κ₃ | Selberg ½ log log T |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.87·10⁴ | 1.730 | −2.10 | 1.829 | −2.374 | 1.741 | −2.140 | 1.143 |
| 10 | 1.38·10⁵ | 1.864 | −2.20 | 1.940 | −2.393 | 1.852 | −2.159 | 1.236 |
| 12 | 1.02·10⁶ | 1.957 | −2.14 | 2.031 | −2.405 | 1.943 | −2.171 | 1.314 |
| 14 | 7.56·10⁶ | 2.019 | −2.13 | 2.108 | −2.414 | 2.020 | −2.180 | 1.381 |
| 16 | 5.58·10⁷ | 2.093 | −2.10 | 2.175 | −2.421 | 2.087 | −2.187 | 1.441 |

The means are all within 0.011 of zero. Selberg's normal is far too narrow at
these heights: the variance does grow like ½ log log T, but with 0.59 to 0.65 on
top. CUE gets the variance within 0.1 and the strong negative skew roughly right
(log|ζ| has a long left tail from the zeros), but it is too wide by 0.074 to 0.099
at every height, and that gap doesn't shrink with N.

The gap is the primes. In the moment conjecture of Keating and Snaith the CUE
moments get multiplied by an arithmetic factor,

```
a(k) = Π_p (1 − 1/p)^(k²) · E|1 − e^(iθ)/√p|^(−2k)      (θ uniform)
```

which is what you get from the hybrid Euler-Hadamard product of Gonek, Hughes and
Keating: ζ is a product over small primes times a product over nearby zeros, and
only the second part looks like a random matrix. Reading log a(s/2) as an extra
cumulant generating function gives the correction to every cumulant of log|ζ|.
Per prime, `values.arithmetic_cumulants` expands
log Σ_m ((s/2)_m / m!)² p^(−m) + (s²/4)·log(1 − 1/p) as a power series in s,
sums over the 78,498 primes below 10⁶ and adds the tail beyond them in closed form
(the p^(−2) term, −(s/2)²(s/2 − 1)²/4 · Σ_(p>P) p^(−2), with the prime sum from
E1(log P)). The result:

| | κ₂ | κ₃ | κ₄ |
|---|---|---|---|
| arithmetic factor | −0.0881 | +0.2337 | −0.0989 |

The variance correction is ½·Σ_p (Li₂(1/p) + log(1 − 1/p)), which the test suite
checks against mpmath's polylog, and κ₃, κ₄ are checked against the cumulants of
−log|1 − e^(iθ)/√p| by quadrature in θ. With it the predicted variance is within
0.015 of the data at all five heights (the sampling error of a variance here is
0.008), and the skew moves from −2.4 to −2.14..−2.19 against −2.10..−2.20 measured.
Nothing was fitted. The right panel of the chart is the histogram minus the CUE
density next to the same difference for CUE with the four arithmetic cumulants
added (`hybrid_density`): the bump at 0 and the dip near 1.7 are the primes.

The full arithmetic factor can't be Fourier-inverted on its own. Its variance is
negative, so its characteristic function grows like e^(0.044u²) and the product
with the CUE one blows up for large u. Truncating it at the fourth cumulant gives a
transform that decays (κ₄ < 0). At N = 12 the resulting density integrates to 1
within 1e−8, has the predicted variance and skew, and its most negative value is
−3e−8, in the far tail.

### Moments

![Moments of zeta over windows of length 1e5 against the CFKRS polynomials and the leading term](docs/moments.png)

Conrey, Farmer, Keating, Rubinstein and Snaith (2005) conjectured the full
asymptotic for the 2k-th moment,

```
∫_0^T |ζ(1/2 + it)|^(2k) dt = ∫_0^T P_k(log(t/2π)) dt + O(T^(1/2 + ε))
```

with P_k a polynomial of degree k². For k = 1 it is a theorem (Ingham):
P_1(x) = x + 2γ. For k = 2 it matches Heath-Brown's fourth moment, and for k = 3 it
is a conjecture. `moments.py` has P_1 in closed form and P_2, P_3 as printed in the
paper. Their leading coefficient is a_k·G(k + 1)²/G(2k + 1), with G the Barnes
function and a_k the arithmetic factor above at integer k; computing it from the
Euler product (primes to 10⁶ plus the E1 tail) gives 1/(2π²) for k = 2 and
5.708527034652e−6 for k = 3, both within 3e−11 of the paper.

`z_power_integrals(C, D)` integrates Z(t)^(2k) with the trapezoid rule on a grid of
step 1/(2k·log(D/2π)), evaluating Z through the same Odlyzko-Schönhage windows as
the zero finder (own ζ below t = 1000). Z^(2k) is close to band-limited, with
frequencies up to k·log(t/2π), so the interior of the trapezoid rule is
spectrally accurate and all of the error sits at the two ends. An 8-point Gregory
correction at each end removes it: the weights come from solving, in exact
rationals, for the corrections that reproduce the Euler-Maclaurin endpoint terms
on polynomials of degree below 8. Against Gauss-Legendre on mpmath's `siegelz`
over 3-unit blocks at step 0.005, including the blocks around the records
Z(50626.48) = 23.7 and Z(77403.72) = 28.2 where Z⁶ passes 10⁸, the agreement is
1e−11 or better. At the default step the endpoint error is about 0.01 in absolute
terms, which is 2e−8 of a 3-unit block and nothing on a long one: halving the step
moves the integral over [50000, 60000] by 5e−12.

The paper has two tables to compare with. For the smoothed moments
∫_0^∞ |ζ|^(2k) e^(−t/10⁴) dt (`smoothed_moments(1e4)`, integrated to t = 6·10⁵):

| k | measured here | CFKRS, measured | conjecture here | CFKRS, conjecture |
|---|---|---|---|---|
| 1 | 79499.9312639 | 79499.9312635 | 79496.7897047 | 79496.7897047 |
| 2 | 5088336.51209 | 5088336.43654 | 5088332.55512 | 5088332.55512 |
| 3 | 708967539.4 | 708967359.4 | 708965694.5 | 708965694.5 |

The conjecture column uses ∫_0^∞ log(t/2π)^j e^(−t/T) dt = T·j!·[s^j] e^(as)Γ(1 + s)
with a = log(T/2π), a power series exponentiated in closed form. It agrees with
the paper to every printed digit, as long as the two numbers in the paper's k = 2
row are read the other way round: printed as is, they say the conjecture is
5088336.43654 and the measured value 5088332.55512, and my conjecture integral
matches the second one to 12 digits. For the sixth moment over blocks of 50,000:

| block | measured here | CFKRS, measured | relative gap |
|---|---|---|---|
| [0, 50000] | 7231006986.4 | 7231005642.3 | 1.9e−7 |
| [50000, 100000] | 15723925543.3 | 15723919113.6 | 4.1e−7 |
| [100000, 150000] | 21536848443.9 | 21536840937.9 | 3.5e−7 |
| [1000000, 1050000] | 70233469614.7 | 70233393177.0 | 1.1e−6 |

The paper's measured values were computed in Mathematica. The second moment agrees
to 12 digits, so the gap grows with k, which is what a fixed relative error in the
integration on their side would do. My two evaluators of Z agree with each other and
with mpmath to 1e−11 or better on the same ranges.

The chart repeats the comparison over windows [T, T + 10⁵] from T = 10³ to 10¹⁰
(470 s for all fifteen windows). The second moment stays within 8e−4 of P_1 and
the fourth within 0.5% up to 10⁹ (1.1% and 2.3% at 3·10⁹ and 10¹⁰). The sixth moment is within 0.7%
up to 3·10⁶ and then wanders by up to 11%: with |ζ|⁶ a window of 10⁵ is
dominated by a few large peaks, and CFKRS's own blocks at 10⁶ move by ±1%. The
right panel is the point of the full polynomial. Keating and Snaith's leading term
a_k g_k (log T)^(k²) / (k²)! is off by a factor 11 for k = 3 even at 10¹⁰, and the
ratio shrinks only like a polynomial in 1/log T.

### The maximum in a short interval

![Distribution of the maximum of log|zeta| over windows of length 2pi at T = 1e12, with CUE and the two-Gumbel limit, and its mean and spread against N](docs/maxima.png)

Fyodorov, Hiary and Keating (2012) predicted the largest value of |ζ| on a
window of fixed length near height T. With N = log(T/2π), the number of zeros in a
window of length 2π,

```
max_(0 ≤ h ≤ 2π) log|ζ(1/2 + i(T + h))| = log N − ¾ log log N + y
```

where y has a limit law with density 4e^(−2y) K₀(2e^(−y)), the law of half the sum
of two independent standard Gumbel variables. The −¾ log log N (rather than the
−¼ that independent values would give) and the y·e^(−2y) right tail come from the
logarithmic correlations of log|ζ|, the same structure as branching random walks.
The leading order is a theorem (Arguin, Belius, Bourgade, Radziwiłł and
Soundararajan 2019), and Arguin, Bourgade and Radziwiłł proved the −¾ term with
tight tails (2020, 2023). The limit law itself is open.

`maxima.py` computes K₀ and K₁ from ∫_0^∞ e^(−x cosh s) cosh(νs) ds with the
trapezoid rule, which converges double-exponentially (1e−12 against mpmath), so
the limit density and its CDF 2e^(−y) K₁(2e^(−y)) need nothing outside numpy. The
law has mean γ and variance π²/12, and 200,000 simulated pairs of Gumbels sit within
the KS noise of the CDF.

`zeta_maxima(T, count)` takes consecutive windows [T + 2πk, T + 2π(k + 1)],
evaluates Z on 16 points per mean zero spacing through the Odlyzko-Schönhage
windows, with the grid aligned to the window edges, and then runs golden-section
search on every grid cell whose higher end is within 0.02 of its window's best
grid value. The grid misses the top of a typical peak by about 0.005, a quarter
of that margin, so the cell holding the true maximum qualifies. Against scans with 100,001 points per
window the maxima agree to 1e−8. 50,000 windows take 52 s at 10⁴ and 124 s at 10¹².

The random-matrix side is the maximum of log|p_N(θ)| = Σ log|2 sin((θ − φ_j)/2)|
over the circle for Haar unitaries (QR of a complex Ginibre matrix with the phases
of R's diagonal divided out). Each term is concave between eigenangles, so the
maximum on each gap is the unique root of ½Σ cot((θ − φ_j)/2), found by Newton with
a bisection safeguard. My first version evaluated p_N by FFT of its coefficients,
which is fine for N = 26 and useless at N = 100: the coefficients sum to 2^N in
absolute value and the rounding error swamps |p_N| (it reported a mean of 17).

Over 50,000 windows at each height:

| T | N | ζ mean | CUE mean | ζ variance | CUE variance | CUE + primes | KS: CUE, + primes, limit |
|---|---|---|---|---|---|---|---|
| 10⁴ | 7.4 | 0.419 | 0.361 | 0.213 | 0.146 | 0.189 | 0.093, 0.068, 0.168 |
| 10⁵ | 9.7 | 0.419 | 0.375 | 0.220 | 0.172 | 0.215 | 0.065, 0.049, 0.165 |
| 10⁶ | 12.0 | 0.420 | 0.381 | 0.234 | 0.193 | 0.235 | 0.063, 0.048, 0.159 |
| 10⁷ | 14.3 | 0.419 | 0.389 | 0.253 | 0.210 | 0.253 | 0.052, 0.035, 0.153 |
| 10⁸ | 16.6 | 0.419 | 0.390 | 0.266 | 0.225 | 0.268 | 0.049, 0.034, 0.147 |
| 10⁹ | 18.9 | 0.417 | 0.393 | 0.280 | 0.237 | 0.280 | 0.042, 0.029, 0.143 |
| 10¹⁰ | 21.2 | 0.415 | 0.394 | 0.291 | 0.247 | 0.290 | 0.039, 0.028, 0.140 |
| 10¹¹ | 23.5 | 0.414 | 0.395 | 0.302 | 0.255 | 0.298 | 0.037, 0.024, 0.137 |
| 10¹² | 25.8 | 0.412 | 0.396 | 0.309 | 0.264 | 0.307 | 0.034, 0.022, 0.135 |

The CUE columns are interpolated between integer N from 100,000 matrices per
size. The data are nowhere near the limit law: the mean is
0.41 against γ = 0.577 and the variance 0.31 against 0.82 at 10¹², and the KS
distance to the limit only falls from 0.168 to 0.135 across eight decades. Random
matrices show the same lag well past any height ζ can reach: at N = 100, 200 and
400 (20,000, 5,000 and 1,000 matrices) the CUE mean is 0.403, 0.385 and 0.418 and
the variance 0.376, 0.404 and 0.474. The slow convergence belongs to the
log-correlated field, not to ζ. FHK met the same
thing at 10²⁸ and rescaled their data to the predicted variance before comparing.

ζ's maxima also spread more than CUE's at every height, by 0.045 in variance at
10¹². That excess is arithmetic. A random matrix has no structure on scales longer
than the circle, but ζ has primes, and p^(−1/2−it) with small p changes slowly
enough to lift or sink a whole window. Averaged over a window of length L the prime
sum Re Σ_p p^(−1/2−it) has amplitude p^(−1/2)·|sin(½L log p)/(½L log p)| for each
p, and since the log p are linearly independent the phases behave like independent
uniforms over many windows, so the shift has variance
Σ_p sin²(½L log p) / (2p·(½L log p)²) = 0.0430 for L = 2π (0.0420 from p ≤ 30). It
correlates with the measured maximum at 0.42 to 0.44 at every height. Adding it to
the CUE variance, with nothing fitted, lands within 0.0041 of ζ's variance from 10⁶
up (0.005 at 10⁵ and 0.024 at 10⁴), and adding an independent sample of the shift
to the CUE maxima cuts the KS distance to the data by 27 to 37%. What is left is a mean offset that shrinks from
0.058 at 10⁴ to 0.017 at 10¹²; with both means removed, the KS distance between
ζ and CUE plus primes is 0.010 to 0.017 from 10⁵ up.

### Li's criterion

![Li's coefficients from a million zeros and their residual after the n log n trend](docs/li.png)

Li (1997) showed that RH holds if and only if

```
λ_n = Σ_ρ [1 − (1 − 1/ρ)^n]
```

is positive for every n ≥ 1. On the critical line 1 − 1/ρ has modulus 1, so the
pair ρ, ρ̄ contributes 4·sin²(n·arctan(1/(2γ))) ≥ 0, which is why the criterion
works. The same numbers are Taylor coefficients at the point s = 1:
Σ_n λ_n z^n / n = log ξ(1/(1 − z)).

`li.li_exact` takes the second route with no zeros at all. A zero ρ of ξ sits at
z = 1 − 1/ρ, which is on the unit circle when ρ is on the line, and an off-line
zero at height γ would land within about 1/γ² of it. So log ξ(1/(1 − z)) is
analytic well past |z| = 1/2. The code samples it on that circle, unwraps the
phase, and reads the coefficients off a discrete Fourier transform at 40 digits;
the aliasing error falls like 2^(−nodes). λ_1 = 1 + γ/2 − ½ log 4π = 0.0230957089661 to 1e−15, and moving the
circle to radius 0.4 with 220 nodes changes nothing at 1e−12.

`li_coefficients(n, gamma)` takes the first route: the sum over the given zeros
plus the rest of the zeros as the smooth density θ'(t)/π beyond the last one,
∫_T^∞ 4 sin²(n arctan(1/2t)) θ'(t)/π dt, and a boundary term −f(T)·S(T) for the
difference between the true and the smooth count at T (integration by parts of
∫ f dS). The integral is done in φ = arctan(1/2t) with φ = φ_T y⁴ to soften the
logarithm of θ' at φ = 0, 200 Gauss-Legendre nodes. With the 1,000,001 zeros below
the Gram point g_(10⁶) (found by the multi-evaluation zero finder in 148 s, count
checked against the Gram point), the sum matches the Taylor route to 1.3e−12 for
every n ≤ 40. Dropping the top half of the zeros moves λ_n by at most 1.6e−7 for n
up to 10⁵, so the tail does its job.

Every λ_n for n ≤ 20,000 (160 s) is positive and increasing, and so are the 120
log-spaced values I computed up to 10⁵. The trend (n/2)(log n − 1 + γ − log 2π) is
Voros's asymptotic under RH; what is left has standard deviation 0.43
for n ≤ 100, 2.8 for 100 to 1000, 4.4 for 1000 to 5000 and 7.5 for 5000 to 20,000,
so it grows, but slower than √n. If RH failed, a zero off the line at β + iγ would
add a term growing like |1 − 1/ρ|^n, exponentially, but only once n reaches about
γ²/(2β − 1), which for any zero above the verified height is far beyond reach. This
is a check of the computation more than a test of RH.

### Jensen polynomials

![Zeros of renormalized Jensen polynomials of xi moving to Hermite zeros, and the error over delta(n)](docs/jensen.png)

Write the completed zeta function around the center as

```
(4z² − 1)·Λ(1/2 + z) = Σ_n γ(n) z^(2n) / n!,      Λ(s) = π^(−s/2) Γ(s/2) ζ(s)
```

and attach to the sequence γ the Jensen polynomials

```
J^(d,n)(X) = Σ_(j=0..d) C(d, j) γ(n + j) X^j.
```

Pólya showed that RH holds if and only if every J^(d,n) is hyperbolic (has only
real zeros). Griffin, Ono, Rolen and Zagier (2019) proved it for d ≤ 8 and every
n, and for each d once n is large enough, by showing that after the change of
variable X → (δ(n)X − 1)/e^(A(n)) and the scaling δ(n)^(−d)/γ(n), J^(d,n) tends to
the Hermite polynomial H_d (in their normalization Σ H_d t^d/d! = e^(−t² + Xt)).

`jensen.taylor_gamma` gets γ(n) from Riemann's theta integral,
Λ(1/2 + z) = −4/(1 − 4z²) + 2∫_1^∞ t^(−3/4) cosh(z log(t)/2) θ₀(t) dt with
θ₀(t) = Σ e^(−πk²t), which turns into
γ(n) = n!/(2n)! · (32·C(2n, 2)·F(2n − 2) − F(2n)) / 2^(2n−1) with
F(m) = ∫_1^∞ (log t)^m t^(−3/4) θ₀(t) dt. All F(m) share one set of Gauss-Legendre
panels in u = log t, so a thousand of them cost one pass over the nodes; at 130
digits γ(0) to γ(1040) take 35 s. They match the paper's γ(10), γ(100) and γ(1000)
(1.6323380490e−17, 6.5777263785e−205, 3.8760340890e−2567) to every printed digit,
mpmath's Taylor series of the left side for n ≤ 3 to 1e−20, and an independent
50-digit quadrature of F(m) to 1e−35.

Hyperbolicity is decided exactly. The γ(n) are converted to rationals, J^(d,n)
is formed in `Fraction` arithmetic, and Hermite's criterion does the rest: a real
polynomial has d distinct real zeros if and only if the Hankel matrix of the
power sums of its zeros, s_(i+j) for 0 ≤ i, j < d, is positive definite. The
power sums come from Newton's identities and the pivots of Gaussian elimination
decide definiteness, so there is no rounding after the γ(n) themselves. Their
precision still matters: with 15-digit γ, J^(20,100) appears to have a pair of
complex zeros, while 60, 160 and 300 digits all give the same positive pivots.
The sweep below runs at 130 digits.

Every J^(d,n) with d ≤ 12 and n ≤ 1000 is hyperbolic (12,012 polynomials, 335 s),
and so are d = 20, 30 and 40 at n = 0, 10, 100 and 1000 (1.5, 11 and 45 s each).
That is a check, not news: Chasse showed hyperbolicity at n = 0 up to
d = 2·10¹⁷. The more interesting output is how the limit is approached. The
paper's closed form for A(n), δ(n), with L(n) the root of n = L(πe^L + 3/4)
evaluated at 2n − 2, agrees with the values read off the data
(δ² = −½·Δ² log γ(n), A = Δ log γ(n) + δ²) to 9% in δ at n = 20, 0.5% at 170 and
0.14% at 492. With the data values, the largest distance between a zero of the
renormalized J^(d,n) and the matching zero of H_d, divided by δ(n), levels off
near 5.7, 14.6 and 24.3 for d = 3, 6 and 9. So the convergence is first order in
δ(n), which falls like 1/√n up to logs, and it is slow: at n = 1000 every zero of
the d = 6 polynomial is still 0.15 to 0.25 below its Hermite limit (left panel).

### Dirichlet L-functions

![Hardy Z functions for three characters](docs/l_z.png)

A Dirichlet character mod q is a homomorphism (Z/qZ)* → C*, extended by 0 to the
integers that share a factor with q. `conrey(q, n)` builds χ_q(n, ·) with the
Conrey labels LMFDB uses: for an odd prime p the generator g_p is the least
primitive root mod p², and with n ≡ g^a, m ≡ g^b mod p^e,
χ(m) = exp(2πi·ab/φ(p^e)). Powers of 2 use ±5^a instead, and a general q
multiplies the prime-power pieces via CRT. Values are stored as exact numerators
over the order, so χ² = 1 is an integer statement and real characters come out
exactly real. The conductor is the smallest d | q for which χ is trivial on
m ≡ 1 mod d, and `primitive()` builds the character mod d that induces χ.

For σ ≥ 0 the series is split by residue class, L(s, χ) = q^(−s) Σ_a χ(a)·ζ(s, a/q),
and each Hurwitz zeta gets the same Euler-Maclaurin tail as ζ itself, starting at
the term N = 10 + |t|/π so that about q·N terms are summed directly. The tail has
an integral term (N + a/q)^(1−s)/(s − 1), and since Σχ(a) = 0 the 1/(s − 1) parts
cancel. Writing that term as ((N + a/q)^(1−s) − 1)/(s − 1) with `expm1` removes
the cancellation before it happens, so L(1, χ) and its neighborhood come out to full
precision with no special case.

For σ < 0 the code uses the functional equation of the primitive character,

```
L(s, χ) = ε(χ)·2^s·π^(s−1)·q^(1/2−s)·Γ(1 − s)·sin(π(s + κ)/2)·L(1 − s, χ̄)
ε(χ) = τ(χ) / (i^κ·√q),   τ(χ) = Σ_a χ(a)·e^(2πia/q)
```

where κ = 0 for even and 1 for odd characters. Trivial zeros at s = −κ − 2m come out
as exact zeros. An imprimitive χ is its primitive character times the Euler factors
(1 − χ*(p)p^(−s)) for the primes that divide q but not the conductor.

On the critical line, θ(t, χ) = (t/2)·log(q/π) + arg Γ((1/2 + κ + it)/2) − arg(ε)/2
makes Z(t, χ) = e^(iθ)·L(1/2 + it, χ) real, the analogue of Riemann-Siegel Z. For a
complex χ the zeros are not symmetric about t = 0; the zeros of χ below the axis are
the zeros of χ̄ above it:

![Zeros of all primitive L-functions mod 31](docs/l_zeros.png)

Sign changes of Z only find zeros on the line. `count_zeros_box` counts every zero
in −1/2 ≤ σ ≤ 2, |t| ≤ T with the argument principle: it follows arg L around the
rectangle, bisecting any step where the argument jumps by more than 0.5 rad, and
subtracts the trivial zeros inside the box. `verify_zeros` finds zeros on the line
and checks that both counts agree, which shows (up to floating-point error, not with
interval arithmetic) that every zero in the box is on the critical line and simple.
When the counts disagree it doubles the sampling density and tries again. That
happened for real: χ₅₉(9, ·) has two zeros at t = 14.636869 and 14.662140, 0.025
apart, with |Z| no bigger than 0.0011 between them, and the default grid stepped
over both. mpmath's `dirichlet` confirms the sign pattern.

Checks:

- L(s, χ) against mpmath's `dirichlet` at random points with −15 < σ < 15 and
  |t| < 80, for ten characters up to q = 40: relative error below 1e-12, worst
  5e-14.
- The first three zeros of χ₃(2, ·), χ₄(3, ·) and χ₅(2, ·) against LMFDB: within
  2e-15. The root number of χ₅(2, ·) matches LMFDB's root angle 0.0881040955873917.
- L(1 − k, χ) = −B_(k,χ)/k with generalized Bernoulli numbers computed from exact
  Bernoulli polynomials, k = 1 to 7 for seven characters.
- The class number formula L(1, χ_d) = 2πh/(w√|d|) for every fundamental d with
  −2000 < d < 0, with h counted from reduced binary quadratic forms: max error in h
  2e-14. For real fields L(1, χ_d)·√d/(2 log ε_d) with the fundamental unit from
  the continued fraction of (1 + √d)/2 or √(d/4) lands within 4e-15 of an integer
  for every fundamental d < 2000, and on the known class numbers 2, 2, 3, 3, 4 for
  d = 40, 60, 229, 316, 328.
- Orthogonality of the character table, the count of primitive characters
  Σ μ(q/d)φ(d), χ_q(m, ·)χ_q(n, ·) = χ_q(mn, ·), the symmetry χ_q(m, n) = χ_q(n, m),
  |τ(χ)| = √q, and ε = 1 for the quadratic characters tested (Gauss's sign of the
  quadratic Gauss sum).
- `verify_zeros` with T = 30 on every primitive character of conductor 3 to 100:
  1,815 L-functions, 80,838 zeros, all on the critical line (342 s).

### Large conductors

![Cost of one value of L(1/2 + it, χ) for both methods](docs/afe_cost.png)

Euler-Maclaurin sums about q·(10 + |t|/π) terms, which is fine for q in the
thousands and hopeless at q = 10⁵. `afe.py` uses a smoothed approximate functional
equation in the style of Rubinstein instead. For a primitive χ with parity κ, the
completed L-function is

```
Λ(s, χ) = (q/π)^z·Γ(z)·L(s, χ) = ε(χ)·Λ(1 − s, χ̄),   z = (s + κ)/2
```

and Λ is the Mellin transform of the theta series Σ χ(n)·n^κ·exp(−πn²x/q). Split
the integral at |x| = 1 along the ray x = r·e^(iφ) instead of the positive axis, and
apply the theta transformation to the inner half. Both pieces become incomplete
gamma functions:

```
Λ(s, χ) = Σ χ(n)·(q/π)^z·n^(−s)·Γ(z, w_n) + ε·Σ χ̄(n)·(q/π)^z'·n^(s−1)·Γ(z', w̄_n)
w_n = πn²·e^(iφ)/q,   z' = (1 − s + κ)/2
```

With φ = 0 this is the classical formula, and it breaks high up. On the critical
line Γ(z) is about e^(−π|t|/4), the two sums cancel down to that size, and dividing
by Γ(z) loses every digit by t ≈ 50. Rotating fixes the scale. Along the ray,
|u^(z−1)| carries a factor e^(−φ·Im z), so with φ = sign(t)·(π/2 − θ) no term is
more than about e^(|t|θ/2) times Γ(z). The code takes |t|θ/2 = 3, which gives up
three e-folds (1.3 digits) and sets θ = 6/|t|. The terms decay like
exp(−πn²·sin θ/q), and each sum stops once that drops below e^(−43):

```
N ≈ √(43q / (π·sin θ)) ≈ √(43q|t| / 6π)
```

That is O(√(q|t|)) terms where Euler-Maclaurin needs O(q|t|). Below |t| = 3.8 the
angle stays at π/2 (no rotation) and N ≈ 3.7√q.

`gammaincc` computes Γ(z, w)/Γ(z) for z with real part 1/4 or 3/4 and imaginary part
up to t/2, and w nearly imaginary. It uses the power series of γ(z, w) when
|w| < |z| + 1 and Legendre's continued fraction (modified Lentz, with converged
entries dropped from the working set) otherwise. The fraction takes from a few to
about a hundred steps over the range the sums need. Both forms are relative to
Γ(z), the scale of the answer, so 1 − γ/Γ never cancels anything that matters. One
more rewrite keeps phases small: n^(−s)·w_n^z = n^κ·(π/q)^z·e^(izφ), so the t·log n
phase of each term cancels exactly instead of being computed and subtracted.

`dirichlet_l` picks a method per point with a cost model fit to the timings above:
the AFE when 0 < σ < 1 and q·(10 + |t|/π) > 25000 + 100·N. At t = 10 the crossover
is near q = 2000. At q = 100003 one value costs 2.0 ms against 57 ms, and at
q = 10007, t = 1000 it is 11.5 ms against 98 ms. `method="em"` forces the old path.

Checks:

- `gammaincc` against mpmath's `gammainc` at 80 random points with |Im z| ≤ 300 and
  10⁻⁴ ≤ |w| ≤ 10³: within 1e-11 (relative to max(1, |Q|)). Also Q(1, w) = e^(−w)
  and the recurrence Q(z + 1, w) = Q(z, w) + w^z·e^(−w)/Γ(z + 1).
- L against mpmath's `dirichlet` at random points of the strip with |t| < 200, for
  nine characters up to q = 97: relative error below 1e-12.
- The value cannot depend on the rotation, and a wrong formula would. At q = 1009,
  allowing 1.5, 3 or 6 e-folds of loss (which changes θ and the number of terms)
  moves L by less than 2e-12 relative at t = 0, 2, 40 and 400.
- AFE against Euler-Maclaurin for χ_(−2003) at σ = 0.2, 0.5 and 0.9: within 1e-11.
  The zeros of χ₁₀₀₀₇(2, ·) below t = 6 found both ways agree to 1e-10, and at
  q = 100003, t = 10 the two methods differ by 7e-14.

### Low-lying zeros

![One-level density of quadratic characters and of all characters mod 2003](docs/low_lying.png)

Katz and Sarnak predict that the zeros of a family of L-functions close to s = 1/2
follow the eigenvalues of a random matrix group, and that the group depends on the
family's symmetry. All characters mod a prime p form a unitary family: nothing
special happens at the central point, and the density of zeros there is flat. The
real quadratic characters χ_d form a symplectic family, and their zeros are pushed
away from s = 1/2 with density

```
W_USp(x) = 1 − sin(2πx) / (2πx)
```

`low_zeros` unfolds each zero with the L-function's own smooth count,
x = (θ(γ, χ) − θ(0, χ))/π, so every family has mean density 1 regardless of its
conductor. That needs θ to increase from t = 0, which holds when
log(q/π) + ψ((1/2 + κ)/2) > 0, i.e. q ≥ 10 for odd and q ≥ 216 for even characters.
Below that the code refuses. Only zeros with γ > 0 count; for a complex χ the zeros
below the axis belong to χ̄, which is in the family too.

The chart uses the 2,433 fundamental discriminants with 2000 ≤ |d| ≤ 6000 (both
signs) and the 2,001 nontrivial characters mod 2003, 13,347 zeros in all. The run
took 8 minutes on one process before the approximate functional equation, and 30 s
on 16 processes now. Results:

- Unitary: flat. Every bin of width 0.1 on [0, 3] is within 0.14 of 1.
- Symplectic: 0.029 zeros per L-function per unit x in [0, 0.1], where W_USp
  averages 0.022, and 0.37 in [0.1, 0.2], where it averages 0.15. Only 4.0% of the quadratic L-functions have a
  zero below x = 0.2, against 19.7% mod 2003.
- The quadratic density also oscillates with period 1 and peaks near x = 1/2, 3/2
  and 5/2, well beyond the oscillation in W_USp. This is Gram's law for L-functions.
  A real primitive χ has ε = 1, so θ(0) = 0 and the n = 1 term of
  Z(t) = Re(e^(iθ)·Σ χ(n)n^(−1/2−it)) is cos θ(t), which vanishes exactly at
  half-integer x. For complex characters x is measured from θ(0) = −arg(ε)/2, a
  phase spread around the circle, so the bias averages out. The higher terms only
  win slowly as log q grows. Stopple's experiment with about a million quadratic
  characters at |d| > 10¹² found that lower-order terms still matter at that size.

![Quadratic one-level density at three conductor sizes](docs/low_lying_scales.png)

With the approximate functional equation the same experiment reaches conductors
near 10⁵. The 2,133 fundamental discriminants with 100000 ≤ |d| ≤ 103500 take 60 s
on 16 processes, and 2,000 characters mod 100003 take 42 s. The unitary family
stays flat (every bin within 0.17 of 1, 20.0% with a zero below x = 0.2). The
quadratic family moves toward W_USp, which averages 0.178 on (0, 0.3]:

| \|d\| | characters | zeros with x ≤ 3 | density on (0, 0.3] | excess over W_USp | first zero below 0.2 |
|---|---|---|---|---|---|
| 2000 to 6000 | 2,433 | 7,350 | 0.377 | 0.199 ± 0.023 | 4.0% |
| 20000 to 23500 | 2,126 | 6,394 | 0.354 | 0.177 ± 0.024 | 3.4% |
| 100000 to 103500 | 2,133 | 6,415 | 0.320 | 0.143 ± 0.022 | 2.8% |

The ± is the Poisson error of the zero count. The excess shrinks about like
1/log|d|: the weighted line through the origin, 1.68/log|d|, passes through all
three error bars. The Gram oscillation is still there at 10⁵, a little lower.

### Lower-order terms from the ratios conjecture

The dashed curves on the left of that chart are not fits. They come from the
L-functions ratios conjecture (Conrey, Farmer and Zirnbauer), which Conrey and
Snaith turned into a one-level density for quadratic characters that is exact up
to O(X^(−1/2+ε)), lower-order terms included. For one χ_d it says the zeros near the
central point have density, per unit x,

```
ρ_d(x) = 1 + Re B(t) / θ'(t),    θ(t) = πx
B(t) = ζ'/ζ(1 + 2it) + A'(it; it) − e^(−2iθ(t)) ζ(1 − 2it) A(−it; it)
A(−r; r) = Π_p (p² − p^(2r)) / (p² − 1)
A'(r; r) = Σ_p log p / ((p + 1)(p^(1+2r) − 1))
```

with θ = θ(t, χ_d), so the gamma factor and the parity κ come in through the same
θ the zeros were unfolded with. The products are over every prime, 2 included.
That relies on a property of the family which `test_family_matches_local_densities`
checks directly: among fundamental discriminants a fraction 1/(p + 1) is
divisible by p, and the rest split evenly between χ_d(p) = 1 and −1. Both poles
at t = 0 cancel, and expanding to first order gives ρ_d(0) = 0, which the code
reproduces to 1e-6 at x = 10⁻⁴.

`ratios.py` evaluates this from scratch. ζ'/ζ on the line Re s = 1 is the
derivative of the entire function (s − 1)ζ(s), taken by the trapezoidal rule on a
circle of radius 1/4 (32 points, within 1e-13 of mpmath right next to the pole).
The Euler products run over primes up to 10⁴, which moves the density by 6e-5
against 10⁵. The arithmetic depends only on t and the gamma part only on the
parity, so a whole family of 2,433 characters costs 0.5 s. The family
prediction is the average of ρ_d over the actual discriminants, each mapped
through its own θ.

Against the histograms (30 bins on [0, 3], Poisson errors):

| \|d\| | χ² against the ratios prediction | χ² against W_USp | excess on (0, 0.3], measured | predicted |
|---|---|---|---|---|
| 2000 to 6000 | 20.7 | 718.9 | 0.199 ± 0.023 | 0.194 |
| 20000 to 23500 | 15.6 | 410.8 | 0.177 ± 0.024 | 0.160 |
| 100000 to 103500 | 17.7 | 298.3 | 0.143 ± 0.022 | 0.136 |

So the bumps at half-integers, the shoulder near the origin and the slow approach
to W_USp are all lower-order terms of the ratios conjecture. Nothing in the data
at this size needs anything beyond them.

The leading correction has a closed form. As |d| → ∞ at fixed x,

```
ρ_d(x) = W_USp(x) + 2(γ − ζ'(2)/ζ(2))(1 − cos 2πx) / (log(|d|/π) + ψ(1/4 + κ/2)) + ...
```

where γ − ζ'(2)/ζ(2) = 1.1472 and the −ζ'(2)/ζ(2) is A'(0; 0). The 1 − cos 2πx
is the half-integer peak pattern from before. The denominator is where the slow
convergence comes from: ψ(1/4) = −4.23 makes it log|d| − 5.37 for even
characters, so at |d| = 4000 the "large" parameter is 2.9. The orange curve on
the right panel is the full prediction for d = ±D out to D = 10⁶⁰. It bends below
the fitted line; the excess times log D goes 1.51, 1.32, 1.23 at D = 10⁶, 10¹²,
10²⁴, on its way to 2 · 1.1472 · 0.4954 = 1.137, where 0.4954 is the average of
1 − cos 2πx over the window. `test_leading_correction_is_closed_form` checks the
expansion at |d| = 10²⁰⁰ for both parities.

### The lowest zero and Fredholm determinants

![First zeros of quadratic L-functions against USp and SO(even) laws](docs/first_zero.png)

The ratios conjecture as used above gives the one-level density, not where the
lowest zero sits. For that the random matrix side needs a gap probability: the
chance that no eigenangle lands in (0, s) is a Fredholm determinant
E(s) = det(I − K) on L²(0, s), with

```
USp, N → ∞:        K(x, y) = S(x − y) − S(x + y),    S(z) = sin(πz)/(πz)
SO(even), N → ∞:   K(x, y) = S(x − y) + S(x + y)
USp(2N):           K(x, y) = D(x − y) − D(x + y),    D(z) = sin((2N + 1)πz/2N) / (2N sin(πz/2N))
```

in units where the mean spacing is 1. `rmt.py` follows Bornemann: Gauss-Legendre
nodes on (0, s), and the determinant of the 40 × 40 Nyström matrix
δ_ij − √w_i K(x_i, x_j) √w_j. The convergence is exponential; 20 and 60 nodes
agree to 1e-14 at s = 2. The density of the lowest zero is −E'(s), and instead of
a finite difference it comes from the resolvent at the endpoint,
−E'(s) = E(s)·R(s, s), where R(·, s) solves the same Nyström system with K(·, s) as
the right-hand side. That matches a central difference to 1e-8.

Checks:

- USp(2) is SU(2), one angle with density (2/π)sin²θ, so E(s) = 1 − s + sin(2πs)/2π
  in closed form. SO(2) is uniform. Both to 1e-14.
- USp(4) and SO(4) against a direct 2-D integral of the Weyl density
  (cos θ₁ − cos θ₂)²·sin²θ₁ sin²θ₂ (or without the sines): 1e-13.
- The sine kernel on (−s, s) splits into its even and odd parts, so
  E_USp(s)·E_SO(s) must equal the GUE gap on an interval of length 2s, and that
  must match Mehta's series 1 − s + π²s⁴/36 − π⁴s⁶/675 + π⁶s⁸/17640.
- Near the origin the USp lowest-zero density is 2π²s²/3, SO(even) starts at 2,
  the densities integrate to 1, and USp(2N) approaches the limit like 1/N.

The data are nowhere near the N → ∞ law. The mean lowest zero is 0.594, 0.618 and
0.635 at the three conductor sizes, against 0.783 for USp. What does fit is
USp(2N) at small N: maximum likelihood over real N (the kernel formula is smooth
in N) gives N = 1.75, 2.05 and 2.30, with likelihood-ratio 95% intervals 1.70 to
1.85, 1.95 to 2.15 and 2.15 to 2.45. The fitted laws have means 0.598, 0.621
and 0.636, and on the log-scale gap probability they follow the data down to
P ≈ 10⁻². The effective size grows by about 0.17 per unit of log|d|, a third of
the slope in the usual dictionary 2N ≈ log|d|. This is the same lesson as the one-level density: at conductors
you can compute, the arithmetic still shapes the zeros near the central point,
and a small matrix mimics it better than the limiting law does.

### Elliptic curves

![Hardy Z of four elliptic curves with central zeros of order 0 to 3](docs/ec_z.png)

`elliptic.py` takes a Weierstrass equation y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆
and builds its L-function from point counts. For an odd prime p, completing the
square turns the curve into Y² = f(x) = 4x³ + b₂x² + 2b₄x + b₆, so
a_p = −Σ_x (f(x)/p), one Legendre symbol lookup per x. p = 2 is counted by hand.
The same sum is right at a bad prime as long as the model is minimal there, and
it gives a_p = 1 or −1 for split or non-split multiplicative reduction. The code
only accepts semistable curves: every prime dividing Δ must leave c₄ a unit, which
also makes the model minimal. Then the conductor is the radical of Δ and the root
number is w = −∏_{p|N} (−a_p). Additive reduction would need Tate's algorithm and
raises an error. The a_n come from multiplicativity and the Hecke recursion
a_{p^(k+1)} = a_p·a_{p^k} − p·a_{p^(k−1)}.

In the analytic normalization (center at s = 1/2), Λ(s) = A^s Γ(s + 1/2) L(s)
with A = √N/2π satisfies Λ(s) = wΛ(1 − s). Splitting the Mellin integral of the
modular form at a point y₀ and folding the inner half over with the Fricke
involution gives

```
Λ(s) = Σ b_n [ (A/n)^s Γ(s + 1/2, n·y₀/A) + w (A/n)^(1−s) Γ(3/2 − s, n/(y₀·A)) ],   b_n = a_n/√n
```

for any y₀ with Re y₀ > 0. The answer can't depend on y₀, which makes a sharp
test: with the right N and w, y₀ = 1, 1.3, 0.8 and e^(0.6i) agree to 1e-12, and
flipping the root number of a twist breaks the agreement at the 1e-3 level. On the
critical line with |y₀| = 1 the two halves are complex conjugates, so
Λ(1/2 + it) = S + wS̄, and the Hardy function is 2·Re(e^(iθ)S) for w = +1 or
2·Im(e^(iθ)S) for w = −1, with θ(t) = t·log A + arg Γ(1 + it). For |t| above 6/π
the code rotates y₀ toward i, the same trick as in the Dirichlet smoothed sums,
so cancellation never costs more than a factor e³.

At the center the incomplete gamma functions collapse to the classical series
L(1/2) = 2Σ (a_n/n)·e^(−n/A) when w = +1 and L'(1/2) = 2Σ (a_n/n)·E₁(n/A) when
w = −1. E₁ is my own: the power series below 1, Legendre's continued fraction
above.

Checks:

- a_n of 11a3 against the eta product q∏(1 − qⁿ)²(1 − q¹¹ⁿ)² up to n = 300, and
  the isogenous curve 11a1 gives the same series. 37a1 against its LMFDB
  q-expansion. Hasse's bound for every p < 3000, and a_{p²} = a_p² − p.
- L(11a, 1) = 0.2538418608559107 and L'(37a, 1) = 0.3059997738340523 to 1e-14.
  The central values of 389a1 (rank 2) and 5077a1 (rank 3) come out below 1e-13.
- The first three zeros of 11a and 37a against LMFDB to 1e-11, and Z(−t) = ±Z(t)
  with the sign of w.

### Quadratic twists of 11a

For a fundamental discriminant d prime to N, the twist E_d has
a_n(E_d) = χ_d(n)·a_n(E), conductor N·d² and root number w(E)·χ_d(−N), so the
twists split by the sign of their functional equation. Katz and Sarnak predict
even orthogonal statistics for the even ones and odd orthogonal statistics for
the odd ones, where the forced central zero is set aside. I took the curve 11a3,
y² + y = x³ − x², which is the one Dueñez, Huynh, Keating, Miller and Snaith
(DHKMS, [arXiv:1107.4426](https://arxiv.org/abs/1107.4426)) use, and every
fundamental d with 5 ≤ d ≤ 6000 prime to 11: 1,674 twists, 832 even and 842 odd.
119 of the even twists have L(E_d, 1/2) = 0 (rank 2), and two odd ones, d = 1393
and d = 5777, have a triple zero (rank 3). The separation is clean: the vanishing
values come out around 1e-15, the smallest non-zero |L'(E_d, 1/2)| among the odd
twists is 0.486, and the smallest non-zero L(E_d, 1/2) is 0.082.

Finding the low zeros uses a symmetry. Λ(1/2 + it)/√w is a real entire function
of t that is even or odd, so it equals P(t²) or t·P(t²). `low_zeros` interpolates
P at Chebyshev points in u = t² up to the height of the third unfolded zero,
doubling the degree until the last coefficients fall below 1e-13 of the largest,
and takes the real roots of the Chebyshev series. The central zero is divided out
exactly, there is no sign scan, and 25 evaluations of Λ are enough at conductor
4·10⁸. On five test twists the roots agree with an Illinois search on sign
changes of Z to 1e-12. Zeros are unfolded with x = θ(γ)/π, so the mean spacing is
1. The whole family took 343 s on 12 processes.

![Central values of the even twists of 11a3](docs/central_values.png)

The central values are discrete. Waldspurger's theorem, in the explicit form of
Kohnen and Zagier, says L(E_d, 1/2) = κ·c(d)²/√d where c(d) is an integer (a
coefficient of a weight 3/2 modular form), and for this family DHKMS give
κ = 6.346046521. Every one of the 832 even twists lands on one of the curves
κc²/√d: the largest distance of √(√d·L/κ) from an integer is 1.6e-7, and c(d)
goes up to 16. So L(E_d, 1/2) is either 0 or at least κ/√d, never in between.

![Lowest zeros of the twists against SO(2N) with and without the cutoff](docs/twist_zeros.png)

That gap is how DHKMS explain the repulsion S. J. Miller found in the lowest zero
of elliptic curve L-functions. In a random matrix model the value of the
characteristic polynomial at 1, Λ_A(1) = det(I − A), plays the part of L(1/2),
and a small central value goes with an eigenvalue close to 1. In SO(2N) nothing
stops det(I − A) from being tiny, and the lowest eigenvalue has density 2 at the
origin. The L-functions can't have tiny central values, so DHKMS drop the
matrices with Λ_A(1) < c·e^(−N/2). `rmt.so_angles` samples Haar SO(2N): QR of a
Gaussian matrix with the sign of diag(R) fixed, then one column flipped when the
determinant is −1. Its one-level density and lowest-eigenvalue law match the
Fredholm determinants above within sampling error. With N = 8, close to
log(√11·6000/2π) = 8.06, and the DHKMS constant c = 2.188, 80.7% of 300,000
matrices survive the cut (14.3% of the even twists vanish). Nothing is fitted to
my data.

| model for the 713 rank 0 even twists | χ² on 15 bins of (0, 1.5) | KS distance |
|---|---|---|
| SO(even), N → ∞ | 181.6 | |
| SO(16) | 164.6 | 0.181 |
| SO(2N), maximum likelihood N = 3.85 | 158.7 | |
| SO(16) with det(I − A) ≥ 0.0389 | 34.6 | 0.043 |

The KS distance of the excised model is under the 5% critical value of 0.051 for
713 points. Mean lowest zero: 0.394 for the twists, 0.408 excised, 0.342 for
plain SO(16). Shrinking N doesn't do it: the likelihood is flat (95% interval 2.65
to 6.85) and every SO(2N) still has density near 2 at the origin, which the data
don't.

### Odd twists and the central derivative

![Lowest zeros of the odd twists against excised SO(17), and χ² as a function of the cut](docs/odd_excision.png)

The odd twists have a forced zero at the center, and Katz-Sarnak say the next one
should follow the first nontrivial eigenvalue of SO(2N + 1), whose forced
eigenvalue at 1 plays the part of the central zero. The data disagree in the same
direction as on the even side: among the 840 rank 1 twists the lowest
non-central zero never falls below x = 0.216, while SO(17) puts 1.7% of its mass
below 0.2, and the histogram is narrower (standard deviation 0.251 against 0.325).

The even fix was a cut on Λ_A(1) at the level where L(E_d, 1/2) is discretized.
For SO(2N + 1) the characteristic polynomial vanishes at 1, so the natural
quantity is the derivative, and det(I − zA) = (1 − z)∏(1 − 2z cos θ_j + z²) gives
|Λ_A'(1)| = ∏(2 − 2 cos θ_j), the same product over the N free angles. On the
arithmetic side Gross and Zagier write L'(E_d, 1/2) as a constant times the height
of a Heegner point over √d, and heights of non-torsion points are bounded below,
so L'(E_d, 1/2) ≥ m/√d. That is the argument in Miller's paper on zero repulsion
([arXiv:1109.0224](https://arxiv.org/abs/1109.0224)). Here the smallest √d·L' is
m = 6.437, at d = 8.

To turn that floor into a matrix cut I followed DHKMS step by step. The CFKRS
recipe says the moments of the family are an arithmetic factor times the matrix
moments. `elliptic.arithmetic_factor` is the Euler product a_s(E) from their
section 5, which gives a_(−1/2) = 0.73269 for 11a with the primes below 2^14,
against their 0.732728078. `rmt.so_moment` gives E|Λ_A(1)|^s for SO(2N) and
E|Λ_A'(1)|^s for SO(2N + 1) in closed form from Selberg's integral (in
x = sin²(θ/2) the angles form a Jacobi ensemble with weight x^(∓1/2)(1 − x)^(−1/2)),
and `rmt.small_value_law` is the residue at the first pole:
P(Λ ≤ ρ) ≈ (R/α)ρ^α with α = 1/2 for SO(2N) and α = 3/2 for SO(2N + 1). That
exponent is why rank 3 is so rare among odd twists. The dictionary holds up on
this data, with N = log(√11·d/2π) rounded per twist:

| s | family mean | a_s(E) × matrix mean | ratio |
|---|---|---|---|
| even, L(E_d, 1/2), s = 1/2 | 0.9574 | 0.9845 | 0.972 |
| even, L(E_d, 1/2), s = 1 | 1.4812 | 1.4828 | 0.999 |
| odd, L'(E_d, 1/2), s = −1 | 0.1994 | 0.2375 | 0.840 |
| odd, L'(E_d, 1/2), s = −1/2 | 0.3992 | 0.3918 | 1.019 |
| odd, L'(E_d, 1/2), s = 1/2 | 3.0766 | 3.0108 | 1.022 |
| odd, L'(E_d, 1/2), s = 1 | 11.232 | 9.375 | 1.198 |

For d > 0 every odd twist has χ_d(11) = −1, and a_(−3/2)(E) = 0.2030 with that
sign. Matching densities near 0, as DHKMS do, scales the floor by a², and matching
probabilities scales it by a^(2/3). With their convention that e^(−N/2) stands for
d^(−1/2), the cut is |Λ_A'(1)| ≥ c·e^(−N/2) with c = 0.265 or c = 2.22.

Neither does anything. With N = 8 and 300,000 Haar matrices from SO(17)
(`rmt.so_angles(..., odd=True)`), c = 2.22 removes 0.13% of them. The arithmetic
floor sits deep in the tail where P(Λ' ≤ ρ) ~ ρ^(3/2) is tiny. On 18 bins of
(0, 1.8):

| model for the 840 rank 1 odd twists | kept | χ² | KS | mean |
|---|---|---|---|---|
| SO(17) | 100% | 148.4 | 0.157 | 0.781 |
| SO(17), c = 2.22 | 99.87% | 147.1 | 0.156 | 0.782 |
| SO(17), best c = 112 | 83.1% | 61.5 | 0.056 | 0.860 |
| data | | | | 0.840 |

A cut 50 times the arithmetic one removes 17% of the matrices and cuts χ² by more
than half, but the KS distance stays above the 5% critical value of 0.047, the
excised mean overshoots, and the histogram is still too wide. Smaller matrices
don't rescue it: the best cut gives χ² 74.6 for SO(5), 78.3 for SO(9) and 68.7 for
SO(13). So the odd twists repel their first zero much harder than the size of
L'(E_d, 1/2) explains, and part of the mismatch is in the bulk. DHKMS handle the
bulk on the even side with an effective N from the lower-order terms of the
one-level density (their r₁ = 2.86 for this curve). The odd family needs the same
computation before an excised model can be judged fairly.

### Chebyshev's bias

![The mod 4 prime race against 3,130 zeros of L(s, χ₋₄), and the first lead change at 26861](docs/race.png)

Chebyshev noticed in 1853 that primes ≡ 3 (mod 4) seem to outnumber primes ≡ 1.
Counting confirms it: π(x; 4, 3) ≥ π(x; 4, 1) for every x below 26861, and up to
10⁸ the 1 class leads at only 1,940 primes. Rubinstein and Sarnak explained why,
assuming GRH and the linear independence of the zero ordinates (LI). With
E(x) = φ(q)·(π(x; q, a) − π(x; q, b))·log x / √x, the explicit formula for each
L(s, χ) gives

E(x) = c(q, b) − c(q, a) − Σ_χ (χ̄(a) − χ̄(b)) Σ_γ x^(iγ) / (1/2 + iγ) + o(1),

where c(q, a) + 1 is the number of square roots of a mod q. The constant is the
whole bias: prime squares are all quadratic residues, so the residue class pays
for them and starts behind by c(q, b) − c(q, a). As log x runs, LI makes the
phases γ log x equidistribute on a torus, so E(x) behaves like
X = C + Σ r_γ cos θ_γ with independent uniform θ_γ and r_γ = 2|χ(a) − χ(b)|/|ρ|.
Its characteristic function is e^(iCt) ∏ J₀(r_γ t), and the share of log x where
a leads is

δ(q; a, b) = P(X > 0) = 1/2 + (1/π) ∫₀^∞ sin(Ct)/t · ∏ J₀(r_γ t) dt.

```python
from zeta_explorer.bias import race, prime_race

R = race(4, 3, 1)                      # zeros of L(s, χ₋₄) up to T = 1000
len(R.gamma), R.bias, R.variance       # 868, 2.0, 0.6223
R.density()                            # 0.9959279226422321
p, d = prime_race(4, 3, 1, 30_000)
p[d < 0]                               # [26861]
```

`race` collects the zeros below T of every character with χ(a) ≠ χ(b), one per
conjugate pair (a complex χ contributes its zeros on both sides of the real axis,
which are the positive zeros of χ and χ̄). Zeros above T are not dropped: their
total weight Σ 1/(1/4 + γ²) has a closed form, log(q/π) + ψ((1 + κ)/2) +
2 Re L'/L(1, χ), so the tail variance is exact and enters as a Gaussian factor in
the product. L'(1, χ) comes from a 48-point Cauchy integral around s = 1, and J₀
from the midpoint rule on (1/π)∫ cos(z sin θ) dθ, which matches mpmath to 1e-14
with about z/2 + 15 nodes. The integral is a trapezoid sum; the integrand is
analytic and even, so the only error is aliasing, the mass of X beyond 2π/h, and h
is chosen so that is 40 standard deviations away.

| race | zeros (T = 1000) | δ here | Fiorilli-Martin | Var X | log share, 10³ to 10⁸ |
|---|---:|---:|---:|---:|---:|
| 3 vs 2 mod 3 (2 ahead) | 823 | 0.9990632 | 0.999063 | 0.4529 | 1.00000 |
| 3 vs 1 mod 4 | 868 | 0.9959279 | 0.995928 | 0.6223 | 0.99894 |
| 2 vs 1 mod 5 | 2,712 | 0.9521746 | 0.952175 | 1.4391 | 0.97730 |
| 3 vs 1 mod 8 | 1,847 | 0.9995688 | 0.999569 | 1.5640 | 1.00000 |
| 5 vs 1 mod 8 | 1,956 | 0.9973946 | 0.997395 | 2.2060 | 1.00000 |
| 7 vs 1 mod 8 | 1,845 | 0.9989378 | 0.998939 | 1.8866 | 1.00000 |
| 5 vs 1 mod 12 | 1,864 | 0.9992059 | 0.999206 | 1.7736 | 1.00000 |
| 7 vs 1 mod 12 | 1,909 | 0.9986061 | 0.998606 | 1.9429 | 1.00000 |
| 11 vs 1 mod 12 | 1,691 | 0.9999766 | 0.999977 | 1.0752 | 1.00000 |
| 5 vs 1 mod 24 | 3,820 | 0.9999875 | 0.999988 | 3.9796 | 1.00000 |
| 11 vs 1 mod 24 | 3,823 | 0.9999826 | 0.999983 | 4.2257 | 1.00000 |
| 23 vs 1 mod 24 | 3,820 | 0.9998886 | 0.999889 | 5.5277 | 1.00000 |

All twelve agree with Table 8 of Fiorilli and Martin (2013) to the six digits they
print, and the tests check eight of them with T = 200, which already lands within
5e-7 of the printed values. Feuerverger and Martin (2000) give 0.952140 for the mod 5 race and say their
complex-character computation interpolated between grid points; the value here
matches the later 0.952175. Their table of Σ_{γ>0} 1/(1/4 + γ²) for χ₋₃, χ₋₄,
χ₋₈ and χ₁₂ is matched to all six digits. As an independent check on the Fourier
inversion, 400,000 draws of X with random phases land within 4σ of the computed
density, and the density integrates to 1 with mean C. The last column is what
the primes actually do up to 10⁸, which is far from the limit for the mod 5 race
(1 leads after 300,619 of the primes below 10⁸), while the races mod 3, 8, 12 and
24 show no reversal at all in that range.

![Limiting distributions of five races](docs/bias.png)

`Race.explicit(x)` rebuilds E(x) itself from the zeros, with Riemann's Li(x^ρ)
for each zero. The plain sum is off by about one prime at small x because the
explicit formula for Π(x, χ) has a constant, log L(0, χ) for odd χ and
log L'(0, χ) − γ − log log x for even χ (L(0, χ) = 0 there), and π needs the
Möbius sum over Π(x^(1/k), χ^k). Those k ≥ 2 terms only involve primes up to √x,
so they are counted directly. Two details matter for complex characters: the
branch of log L(0, χ) comes from following L(σ, χ) along the real axis from
σ = 4, and Li(x^ρ) has to be the integral from −∞, which is Ei(ρ log x) minus
iπ·sign(γ). With Ei alone the error flips sign as T crosses each zero, because
zeros above and below the axis no longer pair up. After both fixes the mean error
over 200 ≤ x ≤ 10⁵ is at most 0.013 in E for the mod 4, 5, 8, 12 and 24 races. With
the 3,130 zeros below T = 3000 the rebuilt curve goes negative for x between
about 26857 and 26876, around the actual first lead change at 26861 (the right panel
of the first chart).

### The de Bruijn-Newman heat flow

![Zeros of H_t near a close pair as t runs from 0.2 to −0.2](docs/heat_flow.png)

Following Pólya and de Bruijn, put

```
Φ(u) = Σ_{n≥1} (2π²n⁴e^(9u) − 3πn²e^(5u)) exp(−πn²e^(4u))
H_t(z) = ∫_0^∞ e^(tu²) Φ(u) cos(zu) du
```

Then H_0(z) = ξ(1/2 + iz/2)/8 and ∂H/∂t = −∂²H/∂z². The de Bruijn-Newman constant
Λ is the smallest t for which H_t has only real zeros, so RH is Λ ≤ 0. Rodgers and
Tao proved Λ ≥ 0 in 2018 and Polymath 15 brought the upper bound to 0.22.

The integral is easy to write and hard to evaluate: H_0(x) decays like e^(−πx/8)
while the integrand is of size 1, so at x = 100 the cosine integral cancels 17
digits. `heat.ht` moves the contour. Φ is analytic and super-exponentially small in
the strip |Im u| < π/8, so the integral over the real line can be taken along
Im u = θ = π/8 − δ, which pulls out a factor e^(−θz) and leaves cancellation of only
e^(δx). With δ = (π/8)·min(1, 4/x) that factor is at most e^(π/2). Φ is even and
real, so the half of the line with u < 0 is the conjugate of the other half, and
what is left is a trapezoid sum over u > 0. On a line inside an analytic strip the
trapezoid rule is exact up to an aliasing error of order e^(−2πδ/h), so a step of
δ/8 is plenty, and the Φ sum stops once πn²e^(4u)·sin 4δ passes 40. At x = 1130
that is 11,700 nodes and 50 terms.

Checks:

- H_0(z) = ξ(1/2 + iz/2)/8 to 2e-12 relative for real and complex z up to x = 400.
  At x = 1128, where H is 5.4e-191, the agreement is 2.6e-11.
- Against mpmath's quadrature of the defining integral itself at 30 digits, for
  t = −1, −0.3, 0 and 0.25: 1e-10 relative.
- The backward heat equation by finite differences, Φ(u) = Φ(−u) (which is
  Jacobi's theta identity), H_t(z̄) = conj H_t(z), and the zeros of H_0 at exactly
  2γ for the first five zeros of ζ.

The zeros move. Write H_t = c∏(1 − z²/x_j²) over its zeros ±x_j. Differentiating
H(x_k(t), t) = 0 and using the heat equation gives ẋ_k = H''/H' at the zero, and
the product turns that into

```
dx_k/dt = 2 Σ_{j≠k} 1/(x_k − x_j) + 2 Σ_j 1/(x_k + x_j)
```

so as t grows the zeros repel like charges, and running t backwards pulls close
pairs together. `heat.velocities` evaluates the sum from a list of zeros plus the
smooth density log(x/4π)/4π beyond the last one. For the seven zeros of H_0 with
1119 < x < 1139 it matches `heat.zero_velocity` (H''/H' from the integral, no
differencing) to 6.8e-7 using the first 4,000 zeros of ζ:

| x_k = 2γ_k | 1120.48 | 1125.12 | 1128.32 | 1129.01 | 1133.40 | 1135.46 | 1137.85 |
|---|---|---|---|---|---|---|---|
| dx_k/dt at t = 0 | 0.23730 | −1.31979 | −3.44762 | 2.00431 | −1.49289 | −0.68430 | −0.55859 |

The third and fourth columns are the closest pair below γ = 600, γ = 564.1609 and
564.5061. Backwards in t they collide. `heat.collision` follows the critical point
of H_t between them and solves for the t where its value crosses zero, which lands
at t_c = −0.0614041 and x_c = 1128.71107, where H and H' both vanish to 1e-13
relative. Near a double zero H_t ≈ (z − x_c)² − 2(t − t_c), which solves the heat
equation exactly, so the squared gap of the pair should be 8(t − t_c) on both
sides: real for t > t_c and −4(Im z)² after. At t_c ± 10⁻⁴ the ratios are 0.99995
and 1.00005, and the right panel shows how the other zeros bend the line further
out. The other zeros hold the pair apart a little: an isolated pair with the same
gap would collide sooner, at −Δx²/8 = −0.0596.

This is the mechanism behind the old lower bounds on Λ. For a Lehmer pair,
Csordas, Smith and Varga show the collision time cannot be below
λ = [(1 − 5Δ²g/4)^(4/5) − 1]/(8g), and `lehmer.newman_bound` gives λ = −0.0616036
for this pair (Δ²g = 0.238). The real collision is only 2.0·10⁻⁴ above it.

![H_t near x = 1129 and its zeros in the complex plane as t decreases](docs/heat_flow.gif)

`heat.zero_paths` tracks every zero in a window as t decreases: real zeros by
sign changes and Illinois steps, and when two real zeros disappear between steps
it starts Newton from their midpoint with an imaginary part of half their old gap,
then follows the complex zero by Newton. Over t ∈ [−0.2, 0.2] the other five zeros
in the window stay real and the pair reaches Im z = 0.545 at t = −0.2. Each of the
81 frames takes about 1.6 s.

### In the browser

`web/` is a static page with no build step and no dependencies. `web/js/zeta.js`
is a port of the Python numerics: Euler-Maclaurin with Backlund's stopping rule,
the reflection formula through log χ, Stirling's log Γ, θ, Gram points, and the
Riemann-Siegel formula with C0 to C4. The Bernoulli numbers are exact BigInt
fractions. The C_k polynomials come precomputed in `web/data/rs.json`, trimmed to
the terms that still matter at |z| ≤ 1/2 (39 to 45 coefficients instead of 90; the
dropped tail changes nothing above 1e-16). Z uses Euler-Maclaurin below t = 200
and Riemann-Siegel above, same as the Python side.

The phase portrait is drawn by up to eight module workers, each rendering a band
of rows, first in 8×8 blocks and then at full resolution. A 948×474 frame of the
default view takes about 1 s on one core (2.2 μs per pixel), so it arrives in
roughly 150 ms. While the new frame renders, the old one is stretched to the new
view so panning stays smooth. The Z explorer samples Z once per pixel, takes zeros
from a table of the first 10,000 and finds the rest by sign changes and Illinois
refinement, and colors a Gram point orange when (−1)ⁿZ(g_n) < 0. At t = 10⁹ each Z
costs about 0.37 ms (12,615 terms), so the sample count is capped by a budget of
4·10⁶ terms per frame. Above 10⁸ or so the phases t·log n lose about 1e-6 radians
in doubles; the windowed evaluator would fix that, but the curve is fine for
looking at. The explicit formula keeps a running sum per x, so moving the slider
up only adds the new zeros.

Dirichlet L-functions have their own port. `web/js/characters.js` builds Conrey
characters the same way as `characters.py` (a generator per odd prime that is
primitive mod p², the ±5^k split for powers of 2) and reads off order, parity,
conductor, the Gauss sum and the root number. `web/js/lfunction.js` evaluates
L(s, χ) as q^(−s) Σ χ(a) ζ(s, a/q) with the head summed directly up to N·q,
N = 10 + ⌊|t|/π⌋, and an Euler-Maclaurin tail per residue class. The tail's
w^(1−s)/(s − 1) term is written as −log w · expm1(u)/u with u = (1 − s) log w; since
Σ χ(a) = 0 the dropped 1/(s − 1) cancels, and the first version, which skipped that
step, lost 9e-8 to cancellation at s = 1 + 1e-9. The terms χ(n) n^(−σ) are cached per
character, so one Z(t, χ) on the critical line costs 12 μs for q = 13 at t = 20,
96 μs for q = 163 and 265 μs for q = 199 at t = 100. For q ≤ 200 and |t| ≤ 1000 that
is well inside the range where the Python dispatcher also picks Euler-Maclaurin
over the smoothed approximate functional equation.

The section plots Z(t, χ) for the chosen character and, below it, the zeros with
|t| < 20 of every primitive character mod q, one row each. Z(−t, χ) = Z(t, χ̄), so a
complex pair is computed once and the second row is the mirror image. The rows go to
a pool of up to four workers: mod 163 there are 161 primitive characters and 81
distinct rows, which take 6.3 s on one thread and arrive in about 2.3 s.

`python -m zeta_explorer.site` writes the three data files: the trimmed
Riemann-Siegel coefficients, the first 10,000 zeros, and the spacing histogram and
pair correlation of 40,000 zeros from n = 10¹² (computed with `zeros_fast`, 13 s
for all three). A test checks that the committed coefficients and the first 200
zeros still match what the code produces.

The JavaScript has its own tests under `web/test`, run with `node --test`. They
compare ζ at 18 points across the plane, Z at 13 heights, θ, Gram points and the
first 29 zeros with mpmath values stored in `reference.json` (written by
`reference.py`), check that the live zero finder agrees with the table, and check
the explicit formula against the Python version and against ψ(x) itself.

Those tests caught a real bug in the Python zero table. `zeros_between` stopped
the Illinois iteration once a step fell below 1e-12 relative, which is 1e-8 in
absolute terms at t = 10⁴. The stopping step is usually much larger than the
remaining error, but not always: 635 of the first 10,000 zeros were off by more
than 5e-10, the worst by 4.5e-9. `first_zeros` now refines to 1e-15 relative and
28 zeros sampled across the table agree with mpmath's `zetazero` to 3e-11. Just above t = 200 the
Riemann-Siegel truncation error (3e-9 in Z) still limits zeros to about 1.5e-9.

## Development

```bash
pip install -e ".[dev]"
pytest
ruff check .
cd web && npm test
python -m http.server -d web 8000
```

## License

MIT
