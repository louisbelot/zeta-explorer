# zeta-explorer devlog

## 2026-10-18: Kickoff, Riemann-Siegel numerics

- Scaffold: hatchling `src/` layout, ruff, GitHub Actions on Python 3.10 and 3.12.
  Dependencies: numpy, mpmath (reference values), matplotlib optional.
- `riemann_siegel.py`: theta from Stirling's series (plus its derivative), the
  C0 correction with its removable singularities handled, Z(t) from the
  Riemann-Siegel main sum + C0, and mpmath wrappers for reference values.
- `zeros.py`: sign-change scan + Illinois refinement, and Gram points via Newton.
- `plot.py`: Z(t) with zeros marked; rendered to `docs/z_function.png`.

Checks:

- theta matches mpmath's `siegeltheta` to 4e-13 from t = 10 to 1e5.
- All 29 zeros below t = 100 match `mpmath.zetazero` to 1e-9.
- Riemann-Siegel error / (2π/t)^(3/4) is flat at about 0.030 over [20, 5000],
  confirming the theoretical decay; the test uses 0.035 as the bound.
- Gram points g_0..g_2 match published values, and theta(g_n) = nπ up to n = 1e5.

Design notes:

- Web zeta visualizers are rare (mostly Wolfram or desktop apps), and implementing
  Riemann-Siegel from scratch is very doable, so that's the showpiece.
- First plot draft labeled the first three zeros, but the labels collided with
  the curve at every crossing; dropped them and listed the values in the README.

Next: higher Riemann-Siegel corrections C1 to C4, then mark
Gram points on the Z plot.

## 2026-10-18 (later): cleanup pass

Cut the module docstrings in `riemann_siegel.py` and `zeros.py` down to a line or
two, since the README already walks through the math, and dropped comments that
only restated the next line. The comments left in the Illinois loop and the Gram
point Newton step explain why those tricks work. Shortened two test names.

Also fixed a README mistake: Z at t = 10⁶ needs 398 terms, not 399
(floor(sqrt(10⁶/2π)) = floor(398.94)). No behavior changes, same 32 tests.

## 2026-10-20: Riemann-Siegel C1 to C4, own zeta, Gram's law, Turing's method

Built:

- `riemann_siegel.py`: theta and Z now take arrays. The main sum builds a
  chunked (t, n) phase matrix and masks n > N(t) per row. Z gets the corrections
  C0 to C4 (`z_rs(t, corrections=4)` is the default). Each C_k is a Taylor
  polynomial in p - 1/2, computed at import from a 110-digit series division of
  Psi(1/2 + z) = -cos(2pi z^2 - 5pi/8) / cos(2pi z). This replaced the old C0
  with its L'Hopital special case.
- `zeta.py`: exact Bernoulli numbers, complex log-gamma (shift past |z| = 15,
  then Stirling), Euler-Maclaurin with Backlund's error bound, and the
  functional equation for Re s < 0 with chi assembled in log space.
- `gram.py`: sign-safe Z at many points, Gram's law, Gram blocks, a Rosser's
  rule checker that refines sampling until it finds enough sign changes, and
  summary stats.
- `turing.py`: Brent's form of Turing's method with Lehman's constants.
  `verify_zeros(T)` proves N(g_n) = n + 1 for a good Gram point g_n >= T.
- Plots: Gram points on Z (bad ones flagged) and the Riemann-Siegel error per
  correction order. Both rendered and checked by eye.

Checks:

- C1 to C4 match mpmath's numerical derivatives of Psi to 1e-15.
- max |error| t^((2k+3)/4) over [200, 2000]: 0.119, 0.052, 0.0090, 0.029, 0.015
  for k = 0..4, about the same over [20, 200]. Past t ~ 3000, C3 and C4 hit a
  ~3e-11 floor from double-precision phases (θ - t log n is ~1e5 rad).
- zeta vs mpmath: relative error < 1e-12 on 300 random points in
  |Re s| <= 30, |Im s| <= 200. zeta(1 - 2k) matches -B_2k/2k to 1e-13. The
  Backlund bound holds and is within 50x of the true error in every test case.
- Gram's law first fails at n = 126 (matches Hutchinson). Up to n = 100,000,
  88.9% of Gram points are good, and the block lengths are 78,694 / 9,445 / 778 /
  19 / 1 for lengths 1 to 5. No Rosser violations.
- verify_zeros(100000): N(g_138068) = 138,069, with K = 2 blocks, in 9 s. That
  matches mpmath's nzeros. Smaller heights are cross-checked in the tests.

Design notes:

- Ran the corrections through a Taylor series instead of closed-form
  derivatives because Psi^(12) by hand is hopeless and Psi is entire. The
  series division is ill-conditioned (poles of 1/cos at |z| = 1/4), so it
  needs extended precision once, at import; evaluating the polynomials is
  plain doubles.
- I couldn't find Gabcke's constants online to cite, so the README reports the
  measured constants only.
- Lehman's integral bound is stated for t > 168 (Trudgian notes the 168 is
  probably a misprint and smaller would do), so verify_zeros never starts below
  168.
- A sign of Z is only trusted from z_rs when |Z| > 1e-6, roughly 100x the
  C4 error at t = 200. Smaller values get recomputed with mpmath.
- Slip: the `plots` commit went in with a failing test (my shell chain checked
  the exit code of `tail`, not pytest). Fixed in the next commit, `fix`.

Next: functional equation demos (they now just need a plot since
chi exists), then domain coloring of the strip with a vectorized zeta, then
spacing statistics vs GUE.

## 2026-10-22: zeta on arrays, phase portrait, continuation, zero spacings vs GUE

Built:

- `zeta.py`: `zeta`, `chi` and `loggamma` take arrays of any shape. Each point
  keeps its own N, the Dirichlet head is a chunked masked matrix, and the
  Euler-Maclaurin loop shrinks its index set as points meet their Backlund
  bound. Scalar calls go through the same code. Added `xi(s)`.
- `continuation.py`: partial sums of the Dirichlet series, the alternating eta
  series and the Euler-Maclaurin corrected sum.
- `zeros.py`: `illinois_many`, a vectorized Illinois method that retires each
  bracket once it converges (relative tolerance, since an absolute 1e-10 is
  below one ulp of t past 10⁶).
- `spacing.py`: `zeros_between(n_lo, n_hi)` brackets every zero between two good
  Gram points (refining Gram blocks until Rosser's count is met) and refines all
  brackets together; unfolding with θ/π; the GUE gap probability as a Fredholm
  determinant with Gauss-Legendre (Bornemann's method); the spacing density and
  CDF from a Chebyshev interpolant of E(s); Wigner surmise; Montgomery's pair
  correlation; an empirical pair correlation with no edge bias; KS distance.
- Plots: phase portrait of the strip, the spiral of ζ(1/2 + it), the three
  partial sums, and the spacing histogram with pair correlation. All four
  rendered and checked by eye; moved two legends and a label off the data.

Checks:

- Array zeta vs mpmath on 2,000 random points (|σ| ≤ 30, |t| ≤ 200): max
  relative error 2.9e-13. Array vs scalar calls agree to 1e-15 relative.
- ξ against mpmath to 1e-12, ξ(0) = ξ(1) = 1/2, ξ(2) = π/6, symmetry, real on
  the line.
- E(s) matches 1 − s + π²s⁴/36 − π⁴s⁶/675 for small s; 20 and 60 nodes agree to
  1e-12. p(s) agrees with 30-digit mpmath derivatives to about 3e-10, integrates
  to 1 with mean 1, variance 0.17999.
- 300 random 120×120 GUE matrices: KS distance to the sine-kernel CDF 0.0044,
  to Poisson 0.28.
- zeros_between(1000, 1060) gives zeros 1002 to 1061, matching `zetazero` to
  1e-8. Near n = 10⁷ it finds exactly 100,000 zeros between g_n and g_(n+100000)
  in 18 s.
- Spacings of 100,000 zeros at n = 10⁵, 10⁶, 10⁷: variance 0.1636, 0.1668,
  0.1693; KS distance to GUE 0.0174, 0.0138, 0.0110.
- Phase portrait grid (594,000 points): 1.05 s vectorized, 6.5 s for a loop over
  the old scalar zeta.

Design notes:

- The price of one code path: a single scalar zeta call went from about 11 µs to
  about 200 µs, mostly numpy overhead in the Euler-Maclaurin loop. Nothing hot
  calls zeta one point at a time, so I kept one implementation.
- E''(s) by Chebyshev differentiation loses about 6 digits near s = 0 (the
  derivative of a degree-80 series amplifies errors by deg² per order), which
  still leaves p(s) good to 1e-9. The CDF is clamped to 1 and the density to 0
  past s = 6, where E(6) ≈ 2e-20.
- The first draft used an absolute tolerance in the vectorized root finder and
  none of the brackets converged at t = 6·10⁵. Switched to relative.
- The zeros are more regular than GUE at these heights. The variance and KS
  numbers move toward GUE with height, matching the finite-height picture of
  Bogomolny and Keating (effective matrix size log(t/2π)/√(12Λ), about 3 at
  t = 5·10⁶). The dip in the pair correlation at 0.3 < u < 0.9 is about −0.013 at
  all three heights, so this range is too short to see it move.

Next: an animation sweeping σ across the strip, then the explicit formula
rebuilding π(x) from the zeros, then Lehmer pairs (the closest zeros found by
zeros_between are a natural place to start).

## 2026-10-23: explicit formula, Lehmer pairs, value distribution

Built:

- `spacing.first_zeros(n)`: the first n zeros, mpmath sign scan below the good
  Gram point g_79 and `zeros_between` above it, with a count check at the seam.
- `explicit.py`: complex Ei from scratch (series, E1 continued fraction with a
  depth that depends on |w|, asymptotic series past |w| = 40), li, Riemann's R
  from the Gram series, a sieve for π and ψ, Möbius, the tail integral by
  Gauss-Laguerre, and the explicit formulas for ψ(x), J(x) and π(x) from a list
  of zeros.
- `lehmer.py`: the Csordas-Smith-Varga sum g for every consecutive pair (exact
  over a window of neighbors and mirrored negative zeros, closed-form smooth
  density outside it), the Lehmer test Δ²g < 4/5, and the resulting lower bound
  on the de Bruijn-Newman constant.
- `values.py`: extrema of Z between zeros by vectorized golden section, running
  records, random samples of log|Z|, exact CUE cumulants from polygamma sums and
  the CUE density by Fourier inversion of the Keating-Snaith product.
- Plots: π(x) and ψ(x) − x rebuilt from 10, 100 and 1000 zeros, Lehmer's pair
  with a zoom, and the log|ζ| histogram against CUE and Selberg. All three
  rendered and checked by eye; the first drafts had a missing 1000-zero curve
  and colliding tick labels on the zoom, both fixed.

Checks:

- Ei vs mpmath on 2,851 random points: max relative error 1.5e-15. R(x) and li(x)
  to 1e-13, the tail integral to 1e-10 against mpmath quadrature.
- With 1000 zeros: π(x) within 0.033 on the half-integers up to 100.5 and rounds
  correctly on 1000.5 to 1099.5; ψ within 0.148. 100 zeros are 4 to 6 times worse.
  At a prime both formulas land on the midpoint of the jump.
- g against an independent formula from the Hadamard product of Ξ (Z', Z'', Z'''
  from mpmath at the zeros): relative error about 3e-6.
- 29 Lehmer pairs below t = 1000, matching Stopple (2017). Near t = 10⁶, 7,384 of
  114,660 pairs, against Stopple's 7,398 of 114,661. Lehmer's pair at 7005 gives
  λ = −7.112e-4, and the maximum of Z between its zeros sits at 7005.0818554, next
  to Stopple's ζ' zero at 7005.08185555. The van de Lune pair's published Δ gives
  −5.894e-9 against the published −5.895e-9.
- Z extrema have |Z'| < 1e-6 by mpmath; the record |Z| below 10⁵ is 28.216 at
  t = 77,403.722.
- CUE: cumulants vs mpmath polygamma to 1e-13; the density integrates to 1 and
  reproduces κ₂, κ₃; 20,000 Haar unitaries of size 12 fit the density within 0.015.
- log|Z| at T = 2π·e^N for N = 8 to 16 (200,000 samples each): variance 1.73 to
  2.09, CUE 1.83 to 2.18, Selberg 1.14 to 1.44.

Design notes:

- li(x^ρ) has to be Ei(ρ log x). The principal branch of li at the complex
  point x^ρ is wrong once γ log x > π, which is almost always.
- The continued fraction for E1 needs over 200 terms near the negative real
  axis at |z| = 10, and the series loses e^(0.3|w|) digits within π/4 of the
  positive real axis. Rather than paper over that wedge, `ei` raises there. The
  explicit formula only evaluates near the imaginary axis.
- The Lehmer bound is stated for zeros of H₀, which are 2γ. Plugging γ straight
  into the formula gives −1.78e-4 for Lehmer's pair instead of −7.11e-4. The
  1993 paper's van de Lune numbers (λ ≈ −Δ²/8 in x units) settled the
  convention. Stopple's paper writes the definition in γ, which is fine for Δ²g
  but not for λ.
- Negative zeros matter at the bottom: for the first pair they are 14% of g
  (0.0229 of 0.1644). They are summed exactly over the first 500 zeros.
- The zeta variance sits 0.07 to 0.10 below CUE at every height. That's
  presumably the arithmetic factor from small primes, which I haven't derived, so
  the README just reports it.
- An mpmath `zetazero(6709)` call ran for minutes, so the reference g uses
  `findroot` on `siegelz` started from our own zeros instead.

Next: an animation sweeping σ across the strip, then the 3D |ζ| landscape, then
a first static web page with precomputed zeros and the explicit-formula curves.

## 2026-10-23 (afternoon): fast multi-evaluation of Z, the 3D landscape, the σ sweep

Built:

- `multieval.py`: `nufft`, a type-1 nonuniform FFT by Gaussian gridding (2x
  oversampled grid, spreading width from the tolerance, deconvolution by
  e^(τk²)); `fft_size`, the next 5-smooth length; and `ZWindow`, a simplified
  Odlyzko-Schönhage evaluator. The window computes the Riemann-Siegel main sum on
  an evenly spaced grid with one NUFFT, interpolates between grid points with a
  Gaussian-windowed sinc, and adds the C0 to C4 remainder. `zeros_fast(n_lo, n_hi)`
  runs the Gram-block zero finder on top of it.
- `riemann_siegel.py`: split out `theta_tail` and `rs_remainder` so the window
  can reuse them.
- `spacing.py`: the Gram-block logic became `gram_block_zeros(g, n0, sample,
  refine)` so it works with any Z. It now handles Rosser violations by merging a
  short block with the previous block and then the following ones.
  `illinois_many` takes an absolute tolerance.
- Plots: `plot_landscape` (log|ζ| as a 3D surface), `sweep_sigma` (an animation
  of ζ(σ + it) as σ crosses the strip, saved as a GIF), and `plot_multieval`
  (accuracy and cost against the direct sum). All three rendered and checked by
  eye.

Checks:

- NUFFT against the direct sum: error below tol·Σ|c| for m = 64, 250, 1000 and
  tol = 1e-6, 1e-10, 1e-13.
- Window against mpmath `siegelz`: 2e-13, 8e-14, 3e-13, 2e-13, 1e-13 at
  t = 10⁴, 10⁶, 10⁸, 10¹⁰, 10¹². The direct sum at the same points: 1e-11, 3e-9,
  5e-7, 6e-5, 3.5e-3. Also exact on both sides of a jump in N at t = 2π·200².
- All 10⁴ zeros in Odlyzko's table for n = 10¹² + 1 to 10¹² + 10⁴: max difference
  2.8e-9, RMS 7e-10, in 2.5 s. The numbering matches the Gram-point count.
- Spacing statistics recomputed at n = 10⁷ with the window match the old row
  exactly (0.1693, 0.0110, −0.013), in 11 s instead of 18 s.
- Rosser violations: none in 100,000 Gram intervals from n = 10⁹, four from
  n = 10¹². A synthetic Z with a borrowed pair of zeros tests the merge, and one
  with a missing pair tests the error.

Benchmarks (this machine, 20,000 points per height, after the window is built):

| t | N | direct sum | window | window build |
|---|---|---|---|---|
| 10⁶ | 399 | 7.4 μs | 2.6 μs | 0.01 s |
| 10⁸ | 3,989 | 58 μs | 2.5 μs | 0.03 s |
| 10¹⁰ | 39,894 | 795 μs | 2.6 μs | 0.25 s |
| 10¹² | 206,393 | 7.7 ms | 2.6 μs | about 1.3 s |

Design notes:

- Precision was the bigger win. t·log n is ~10¹³ rad at t = 10¹², so a double
  phase is off by ~10⁻³. The window computes t₀ log n mod 2π and θ(t₀) once in
  mpmath and works in the offset u = t − t₀ from then on. θ(t₀ + u) − θ(t₀) is
  written with log1p so no term is large. Zeros and Gram points come back as
  offsets, since a double near 2.7·10¹¹ only resolves 6e-5.
- Interpolation error: with sampling at α = 2β, a sinc of bandwidth α and a
  Gaussian window with ε² = (α − β)/(Kδ), both the spectral leakage and the
  truncation are e^(−πK/4). K = 36 for 1e-12, and the bound
  e^(−9π) = 5e-13 sits just above the measured 1e-13 to 3e-13.
- Gaussian gridding needed w = 15 grid points on each side for 1e-13. The
  exponential-of-semicircle kernel FINUFFT uses would need fewer, but its
  Fourier transform has no closed form, and spreading is not the bottleneck
  here. The mpmath phases are: about 1 s of the 1.3 s build at 10¹².
- When N changes inside a window the grid uses the larger N and subtracts the
  last term directly where t is below 2πN². That beats splitting the window.
- zeros_fast's profile at 10¹² is now dominated by `np.polyval` in the
  degree-90 remainder polynomials, which get called on small arrays many times
  during root refinement.
- The first landscape drafts plotted |ζ| capped at 4 or 5. The pits at the zeros
  were hidden behind the front of the surface, and the cap made a flat yellow
  plateau. log₁₀|ζ| clipped to [−1.5, 1], with the zeros added to the t grid,
  fixed both. Viewing from higher up floated the zero markers off the surface,
  so the camera stays at 42°.
- The sweep uses σ = 0.5 + 0.4·sin³ so the animation pauses near the critical
  line; with a plain sine it rushed through σ = 1/2 in one frame.

Next: a static web page with precomputed data, then ζ in the browser with a
canvas domain coloring. Shrinking the remainder polynomials (degree 90 is far
more than doubles need) would speed up zeros_fast.

## 2026-10-26: the browser version

Built:

- `site.py`: writes `web/data/` (the C0 to C4 polynomials trimmed for |z| ≤ 1/2,
  the first 10,000 zeros, and spacing and pair-correlation histograms for 40,000
  zeros from n = 10¹²) with `python -m zeta_explorer.site`. 13 s.
- `web/js/zeta.js`: ζ(s) by Euler-Maclaurin plus reflection, log Γ, θ, Z by
  Euler-Maclaurin below t = 200 and Riemann-Siegel above, Gram points, Illinois,
  and a sign-change zero finder. `explicit.js`: sieve, ψ(x), and the explicit
  formula as a running sum over zeros.
- The page: a phase portrait rendered in bands by module workers (coarse pass,
  then full resolution), a Z(t) explorer with table and live zeros and Gram points
  colored by Gram's law, the explicit formula with a slider and a play button,
  and the spacing and pair-correlation charts at height 10¹². Light and dark
  themes, pan and pinch on touch.
- CI: a node job runs `npm test` in `web/`, and a Pages workflow deploys
  `index.html`, `style.css`, `js/` and `data/`.

Checks:

- 24 node tests. ζ against mpmath at 18 points from σ = −10 to 4 and |t| ≤ 100
  (1e-12 relative), Z at 13 heights up to 10⁶ (1e-12 below 200, 1e-8 above), θ and
  Gram points to 1e-12 relative, the first 29 zeros and zeros 1990 to 2000 to
  1e-10 and 1e-9, the table against the live finder near t = 9800, and ψ from 100
  zeros against the Python function to 1e-9.
- Python: the trimmed coefficients differ from the full degree-90 polynomials by
  less than 1e-16 on [−1/2, 1/2]; the committed data match the generator.
- Looked at the page in light and dark mode at desktop width and at 375 px.

Bug found on the way: the JS zero finder disagreed with the Python table by 2e-9
near t = 9800, and mpmath sided with the JS. `zeros_between` stops when the
Illinois step drops under rtol·t, which is 1e-8 at that height, and sometimes the
error at that point is still a few 1e-9. 635 of the first 10,000 table zeros were
off by more than 5e-10 (worst 4.5e-9). `first_zeros` now asks for 1e-15 relative,
and a sample of 28 zeros matches `zetazero` to 3e-11. The JS Illinois uses an
absolute tolerance plus a few ulps for the same reason. A regression test uses
three of the zeros that used to be wrong.

Design notes:

- Workers instead of WebGL. A float32 shader would get ζ at t = 50 to about four
  digits (t·log n is already 160 radians there), and doubles in workers are fast
  enough: 2.2 μs per pixel, about 150 ms for a 948×474 frame on 8 workers. Each
  worker gets one band at a time from a queue, and a new view replaces the queue,
  so dragging never piles up stale work.
- Scroll zoom only kicks in after the plot has been clicked (or with a pinch,
  which browsers send as ctrl+wheel). Otherwise scrolling down the page got stuck
  on the first canvas.
- The log n and n^(−1/2) tables in the Riemann-Siegel sum are cached and grown on
  demand. At t = 10⁹ Z is still 0.37 ms, so the explorer budgets 4·10⁶ terms per
  frame and samples less often up there.
- The JS keeps Riemann-Siegel from t = 200 like the Python side. Its error there is
  3e-9, which shows up as a 1.5e-9 shift in zero 81; moving the switch up to 1000
  would cut that to about 1e-10 for a few hundred extra terms per point.

Next: math notes for the page (how each piece is computed, with the formulas), a
GIF of the explicit-formula animation for the README, then a PyPI release.

## 2026-10-27: Dirichlet L-functions and low-lying zeros

Built:

- `characters.py`: Dirichlet characters with LMFDB's Conrey labels. Odd primes use
  the least primitive root mod p² (10 for p = 40487, where 5 fails mod p²), powers
  of 2 use ±5^a. Values are stored as integer numerators over the order, so real
  characters are exactly real and products are exact. Conductor, inducing primitive
  character, parity, Gauss sum, root number. Also the Kronecker symbol, fundamental
  discriminants, quadratic characters, class numbers of imaginary quadratic fields
  from reduced forms, and fundamental units of real quadratic fields from the
  continued fraction of (1 + √d)/2 or √(d/4).
- `lfunction.py`: L(s, χ) anywhere, by Euler-Maclaurin per residue class for σ ≥ 0
  and the asymmetric functional equation for σ < 0; imprimitive characters as the
  primitive one times Euler factors. θ(t, χ) and the Hardy Z function, zeros by sign
  changes plus Illinois, and `count_zeros_box`, an argument-principle count over
  −1/2 ≤ σ ≤ 2, |t| ≤ T. `verify_zeros` compares the two counts and refines the grid
  when they disagree.
- `lowlying.py`: zeros unfolded by θ(t, χ), one-level densities for quadratic
  characters and for all characters mod a prime, and W_USp.
- Plots: Z(t, χ) for three characters, the zeros of all primitive L-functions mod
  31 with conjugate pairs side by side, and the one-level density and first-zero
  histograms. All three rendered and checked by eye; the density chart needed its
  legend moved off the peaks.

Checks:

- Against mpmath's `dirichlet`: ten characters up to q = 40, random points with
  −15 < σ < 15 and |t| < 80, worst relative error 5e-14. At t = 2000 within 1e-10.
- Against LMFDB (queried live once, values hardcoded in the tests): first three
  zeros of 3.2, 4.3 and 5.2 within 2e-15, and the root angle of 5.2.
- Generalized Bernoulli numbers L(1 − k, χ) = −B_(k,χ)/k, exact trivial zeros, the
  functional equation with mpmath's log Γ, χ_q(m, ·)χ_q(n, ·) = χ_q(mn, ·) and the
  symmetry χ_q(m, n) = χ_q(n, m), character orthogonality, the primitive count
  Σ μ(q/d)φ(d), and Euler's criterion for the Kronecker symbol.
- Class number formula: imaginary fields for every fundamental −2000 < d < 0 (max
  error in h 2e-14), real fields land within 4e-15 of an integer for every
  fundamental d < 2000.
- `verify_zeros(χ, 30)` for all 1,815 primitive characters of conductor 3 to 100:
  80,838 zeros, every box count matched, 342 s.

Bug found on the way: the first run of that sweep reported four mismatches, two
conjugate pairs mod 59 and mod 67, with two more zeros in the box than on the line.
χ₅₉(9, ·) has zeros at 14.636869 and 14.662140, 0.025 apart, and |Z| reaches only
0.0011 between them; the grid at 1/8 of the mean spacing (0.11) stepped over both.
mpmath agrees on the signs. `verify_zeros` now doubles the density until the counts
match, and a test pins that pair.

Design notes:

- The integral term of each Hurwitz tail is (N + a/q)^(1−s)/(s − 1). The 1/(s − 1)
  parts cancel over a because Σχ(a) = 0, so the code writes the term as
  −log(w)·expm1(z)/z with z = (1 − s)·log w and never forms them. L(1, χ) needs no
  limit, and the finite difference of L at 1 ± 1e-9 matches mpmath's derivative to
  1e-6.
- I first assumed the inducing character of χ_q(n, ·) has label n mod f. That holds
  for conductors p but not for higher prime powers or powers of 2: with n = g^a and
  a = p^(e−f)·a', the inducing label is g^(a') mod p^f. χ₂₇(8, ·), (10), (17), (19)
  are counterexamples. The test now checks the defining property instead.
- Unfolding by θ only works when θ increases from 0. θ'(0) has the sign of
  log(q/π) + ψ((1/2 + κ)/2), and ψ(1/4) = −γ − π/2 − 3 log 2 ≈ −4.23, so even
  characters need q ≥ 216. `low_zeros` checks this up front.
- The cost is about q·(10 + |t|/π) terms per value, fine for q in the thousands at
  small t (0.07 s per character mod 2003, 0.13 s per quadratic character) but the wrong shape for
  large conductors.

Low-lying zeros (2,433 quadratic characters with 2000 ≤ |d| ≤ 6000, 2,001 characters
mod 2003, 8 minutes): the mod-2003 family is flat to within 0.14 in every bin, and
the quadratic family has 0.029 zeros per unit x in [0, 0.1] against 0.022 for W_USp.
4.0% of quadratic L-functions have a zero below x = 0.2, 19.7% mod 2003. The
quadratic density also oscillates with period 1, peaking near half-integers. That is
Gram's law for real characters: ε = 1 puts the zeros of the n = 1 term, cos θ, at
half-integer x, and complex characters average it out through their root numbers.

Next: an approximate functional equation so L costs O(√(q·t)), which would let the
family runs reach conductors near 10⁵, and the ratios-conjecture lower-order terms
drawn on the quadratic chart.

## 2026-10-29: Smoothed approximate functional equation

Built:

- `afe.py`: L(s, χ) for primitive χ inside the critical strip from the rotated
  theta integral, two sums of incomplete gamma terms with N ≈ √(43q|t|/6π) terms
  each. `gammaincc(z, w)` is the regularized upper incomplete gamma for complex z
  and w: power series below |w| = |z| + 1, Legendre's continued fraction (modified
  Lentz) above, with converged entries dropped from the working set.
- `dirichlet_l` now chooses per point between Euler-Maclaurin and the AFE with a
  cost model, and `method="em"` forces the old path. The Gauss sum is cached on the
  character, since θ(t, χ), Z and the AFE all need ε and it costs O(q).
- Plots: cost per value against q and t for both methods, and the quadratic
  one-level density at |d| ≈ 4·10³, 2·10⁴ and 10⁵ with the excess near the origin
  against 1/log|d|. Both rendered and checked by eye. The first version of the scale
  chart silently drew two of the three families (a two-color palette zipped with
  `strict=False`); the zip is strict now.

Checks:

- `gammaincc` against mpmath at 80 random points, |Im z| ≤ 300, |w| from 1e-4 to
  1e3: within 1e-11. Q(1, w) = e^(−w) and the recurrence in z.
- L against mpmath's `dirichlet` at random strip points, |t| < 200, nine characters
  up to q = 97: below 1e-12 relative. A wider run with 110 points had a worst error
  of 2e-13.
- Rotation independence: allowing 1.5, 3 or 6 e-folds of loss at q = 1009 changes L
  by less than 2e-12 at t up to 400. This is the check I trust most, because a
  sign or normalization mistake in the formula would make the answer depend on φ.
- AFE against Euler-Maclaurin at q = 2003 (1e-11), zeros of χ₁₀₀₀₇(2, ·) both ways
  (1e-10), and 7e-14 at q = 100003, t = 10.
- Rerunning the old 2000 ≤ |d| ≤ 6000 family through the new dispatch reproduces
  the previous numbers exactly: 0.029 and 0.37 zeros per unit x in the first two
  bins, 4.0% with a zero below x = 0.2.

Benchmarks (one value of L(1/2 + it, χ), 8 points averaged, this machine):

| | Euler-Maclaurin | AFE |
|---|---|---|
| q = 101, t = 10 | 0.11 ms | 1.08 ms |
| q = 1009, t = 10 | 0.61 ms | 1.21 ms |
| q = 10007, t = 10 | 5.45 ms | 1.45 ms |
| q = 100003, t = 10 | 57.3 ms | 1.99 ms |
| q = 10007, t = 1000 | 98.3 ms | 11.5 ms |

The AFE has a floor of about 1 ms from Python overhead per point (two kernels, two
log Γ), which is why small q stays on Euler-Maclaurin.

Design notes:

- The rotation angle comes from how much growth is allowed: |t|θ/2 = 3, so terms
  are at most e³ times the answer. Doubling the allowance would cut N by √2 and
  cost another 1.3 digits; 3 kept every check below 1e-12.
- n^(−s)·w_n^z = n^κ·(π/q)^z·e^(izφ). Writing the smooth part this way means no
  term ever computes the t·log n phase, which is the largest rounding error in the
  Euler-Maclaurin path.
- The two normalizations e^(−log Γ(z)) and the e^(−φ Im z) inside each sum are huge
  and tiny at large t (e^(±7800) at t = 10⁴), so they are combined in one exponent
  before calling exp.
- The AFE is only used for 0 < σ < 1. Near s = 0 and s = 1 log Γ of z or z' goes
  to a pole and the series part cancels; the strip edges were checked against mpmath
  (3e-15 near s = 0, 3e-13 near s = 1) but Euler-Maclaurin is better there anyway.
- LMFDB's API answered with a reCAPTCHA page today, so the large-conductor checks
  lean on mpmath, rotation independence and the cross-check between methods.

Low-lying zeros at larger conductors: 2,133 quadratic characters with
100000 ≤ |d| ≤ 103500 in 60 s and 2,000 characters mod 100003 in 42 s, both on 16
processes. The density on (0, 0.3] minus W_USp goes 0.199, 0.177, 0.143 (±0.023)
at |d| ≈ 4·10³, 2·10⁴, 10⁵, which a single 1.68/log|d| fits. The mod-p family is
flat to 0.17 in every bin.

Next: the ratios-conjecture lower-order terms for the quadratic family, drawn on
top of the scale chart, to see whether they account for the 1.68/log|d| excess.

## 2026-11-04: Ratios-conjecture lower-order terms and the lowest zero

Built:

- `ratios.py`: the Conrey-Snaith one-level density for quadratic characters,
  ρ_d(x) = 1 + Re B(t)/θ'(t) with B = ζ'/ζ(1 + 2it) + A'(it; it) − e^(−2iθ)ζ(1 − 2it)A(−it; it),
  evaluated per discriminant and averaged over the family. ζ'/ζ comes from the
  trapezoidal rule for the derivative of (s − 1)ζ(s) on a circle of radius 1/4,
  so it is accurate right next to the pole. A(−it; it) simplifies to
  Π (p² − p^(2it))/(p² − 1), summed as logs.
- `rmt.py`: Fredholm determinants with Bornemann's Gauss-Legendre Nyström method,
  limit kernels for the sine process, USp and SO(even), and the finite USp(2N) and
  SO(2N) kernels in unfolded units. The lowest-eigenvalue density comes from the
  resolvent at the endpoint, −E' = E·R(s, s), and `fit_symplectic_n` does a
  maximum-likelihood fit of N with a likelihood-ratio interval.
- `digamma` next to `loggamma` (same Stirling shift), and a public `primes`.
- Plots: the scale chart now carries the ratios prediction for the same characters
  and the predicted excess out to D = 10⁶⁰, and a new lowest-zero chart compares
  the histograms and the log-scale gap probability with USp, SO(even) and fitted
  USp(2N). Both rendered and checked by eye.

Checks:

- ζ'/ζ against mpmath at six points including 1 + 0.01i: within 1e-13 relative.
  The Laurent constant γ at u = 10⁻⁶.
- A'(0; 0) = −ζ'(2)/ζ(2) = 0.5699610 to 2e-5 with primes up to 10⁶ (the tail is
  about 1/P). A(−r; r) and A'(r; r) against the unsimplified product A_D(α; γ)
  from the paper and its numerical α-derivative.
- The family assumption behind the products: among fundamental discriminants
  with |d| ≤ 30000, the share divisible by p = 2, 3, 5, 7 is 1/(p + 1) to 0.01
  and χ_d(p) is balanced.
- I derived the first-order expansion, ρ_d = W_USp + 2(γ − ζ'(2)/ζ(2))(1 − cos 2πx)/(log(|d|/π) + ψ(1/4 + κ/2)),
  and the code matches it to 0.05 (times the log) at |d| = 10²⁰⁰ for both parities.
- Fredholm: USp(2) and SO(2) closed forms to 1e-14, USp(4) and SO(4) against a
  direct 2-D Weyl integral to 1e-13, the parity split E_USp·E_SO = E_GUE(2s) to
  1e-14, Mehta's series for the GUE gap through s⁸, 2π²s²/3 near the origin, and a
  synthetic USp(6) sample from which the fit recovers N = 3.

Results on the three families from the last entry (zeros recomputed today, same
counts as before: 7,350, 6,394, 6,415):

| \|d\| | χ² ratios | χ² W_USp | excess measured | excess predicted | N_eff (95%) |
|---|---|---|---|---|---|
| 2000 to 6000 | 20.7 | 718.9 | 0.199 ± 0.023 | 0.194 | 1.75 (1.70 to 1.85) |
| 20000 to 23500 | 15.6 | 410.8 | 0.177 ± 0.024 | 0.160 | 2.05 (1.95 to 2.15) |
| 100000 to 103500 | 17.7 | 298.3 | 0.143 ± 0.022 | 0.136 | 2.30 (2.15 to 2.45) |

χ² is over 30 bins with Poisson errors. The ratios prediction accounts for the
whole 1.68/log|d| excess, the half-integer bumps and the shoulder. The mean lowest
zero is 0.594, 0.618, 0.635, against 0.783 for the USp limit and 0.598, 0.621,
0.636 for the fitted USp(2N).

Design notes:

- My first t-grid search started at t = 1 and grew by 1.5×, which is already past
  x = 3 for |d| = 10¹⁰⁰, so almost every grid point was wasted and the large-|d|
  check looked like it diverged. θ is convex for t > 0 (Re ψ grows with |Im|), so
  Newton started from the bound πx/θ'(0) comes down monotonically onto the root.
- Everything arithmetic depends only on t and the gamma factor only on parity, so
  the per-character loop is a θ evaluation and an interpolation. Dropping the
  prime cutoff from 10⁵ to 10⁴ moved the density by 6e-5 and the cost from 5 s to
  0.5 s for 2,433 characters.
- The resolvent formula for −E' avoids a finite difference, and caching the
  Gauss-Legendre nodes made the N fit four times faster (leggauss solves an
  eigenproblem on every call).
- Fitting real N uses the kernel formula outside the integers. It is smooth there,
  and N must exceed the largest data point, since USp(2N) angles fill only x ≤ N.

Next: the even orthogonal side, either quadratic twists of an elliptic curve or
real quadratic characters with the central zero statistics, so the SO(even) curve
on the lowest-zero chart has data of its own.

## 2026-11-06: Elliptic curves and the even/odd twist families

Built:

- `elliptic.py`: `Curve` from Weierstrass coefficients (b and c invariants,
  discriminant, conductor and root number for semistable curves), a_p by counting
  points with a Legendre symbol table, a_n by multiplicativity and the Hecke
  recursion, cached in power-of-two blocks. `EllipticL(E, d)` is the L-function of
  the quadratic twist, with χ_d(n)·a_n, conductor N·d² and root number
  w(E)·χ_d(−N).
- `completed(s, L, y0)`: Λ(s) in the strip from the split Mellin integral with
  incomplete gamma functions, split point y0 free. `hardy_z` uses the fact that the
  two halves are conjugate on the critical line, so it needs one incomplete gamma
  sum per point, and rotates y0 for |t| > 6/π. `central` gives L(1/2) or L'(1/2)
  from the classical e^(−n/A) and E₁(n/A) series, with `exp1` written from scratch.
- `low_zeros` and `zeros_below`: Chebyshev interpolation of P(u), where
  Λ(1/2 + it)/√w = P(t²) or t·P(t²), roots of the Chebyshev series, central order
  from the central value.
- `rmt.py`: `odd_orthogonal_n` (the SO(2N + 1) kernel), `so_angles` (Haar SO(2N)
  by QR with the sign fix and a column flip), `char_poly_at_one` and
  `excised_first`, the DHKMS excised orthogonal ensemble. `fit_symplectic_n`
  became `fit_n` with the kernel family as a parameter.
- Plots: Z(t) for 11a3, 37a1, 389a1 and 5077a1 (central orders 0 to 3), the
  central values of the even twists on the curves κc²/√d, and the lowest-zero
  histograms of the even and odd twists against SO(even), fitted SO(2N), excised
  SO(16), SO(odd) and SO(17).

Checks:

- a_n of 11a3 equal the eta product q∏(1 − qⁿ)²(1 − q¹¹ⁿ)² to n = 300; 11a1
  (isogenous) gives the same series; 37a1 matches LMFDB; Hasse bound to 3000.
- L(11a, 1) and L'(37a, 1) to 1e-14 against Cremona's values. 389a1 and 5077a1
  give central values below 1e-13 with w = +1 and −1, as their ranks 2 and 3 need.
- Λ independent of the split point to 1e-12 for two curves and two twists, and a
  wrong root number breaks that by more than 1e-3.
- First three zeros of 11a and 37a against LMFDB to 1e-11. Chebyshev zeros of five
  twists against an Illinois sign scan of Z to 1e-12.
- Waldspurger: for the 52 even twists with d ≤ 400, √(√d·L(E_d, 1/2)/κ) is an
  integer to 1e-7, and the five zeros are d = 232, 265, 273, 344, 364.
- Haar SO(8): one-level density against (2N − 1)/2π + sin((2N − 1)θ)/(2π sin θ)
  within 4.5σ per bin, lowest eigenvalue of SO(6) against the Fredholm law by KS.
  SO(3) closed form 1 − s + sin(πs)/π to 1e-14.

Results for 11a3 twists with 5 ≤ d ≤ 6000 (1,674 curves, 832 even, 842 odd):
119 even twists vanish at the center, two odd ones (d = 1393, 5777) have a triple
zero. All 832 even central values sit on κc²/√d with κ = 6.346046521, c up to 16,
worst distance from an integer 1.6e-7. For the 713 rank 0 even twists, the lowest
zero has χ² 181.6 against SO(even), 164.6 against SO(16), 158.7 against the ML
SO(2N) (N = 3.85), and 34.6 against SO(16) with the DHKMS cutoff
det(I − A) ≥ 2.188·e^(−N/2), with no fitted parameter. KS distance 0.043 against
0.181 without the cut. The odd twists' first non-central zero has mean 0.840,
above the SO(odd) limit (0.783), and an SO(2N + 1) fit hits the boundary N > max x.

Design notes:

- Λ(1/2 + it) is entire of order one, so Chebyshev interpolation in u = t²
  converges geometrically and degree 24 was enough for every twist I tried, up to
  conductor 4·10⁸. Z itself carries 1/|Γ(1 + it)|, which has branch points at
  t = ±i and would slow the interpolation, so the interpolant is for Λ.
- Working in t² factors out the parity: the central zero of an odd twist is
  divided out exactly instead of being a root the solver has to find next to the
  ones I want.
- Point counting is O(p) per prime, which is the slow part: 3.6 s for all a_n up
  to 65,536 and 15 s to 131,072. It runs once per curve and the twists reuse it.
- The family run is 1,674 independent twists, so it went through a process pool:
  343 s on 12 processes, about 2.5 s of process time per twist.
- The quadratic Dirichlet experiment wanted USp(2N) with small N. Here plain
  SO(2N) at any N fails near the origin, and the fix is a cut on the central value
  instead of a smaller matrix.

Next: a model for the odd twists. Either an excised SO(2N + 1) with a cut on the
derivative of the characteristic polynomial at 1, or the Heegner point heights
behind L'(E_d, 1/2). After that, Dirichlet L-functions in the browser.

## 2026-11-07: Chebyshev's bias

Built `bias.py`: prime races mod q from the zeros of Dirichlet L-functions, after
Rubinstein and Sarnak.

- `race(q, a, b, T)` gathers the zeros below T of one character per conjugate
  pair with χ(a) ≠ χ(b) (both sides of the axis for complex χ), the amplitudes
  2|χ(a) − χ(b)|/|ρ|, the bias c(q, b) − c(q, a) from counting square roots, and
  the exact variance of the zeros above T.
- `zero_sum(χ)` = Σ_γ 1/(1/4 + γ²) = log(q/π) + ψ((1 + κ)/2) + 2 Re L'/L(1, χ),
  from the Hadamard product of Λ. L'(1, χ) by a Cauchy integral on |s − 1| = 0.4.
- `bessel_j0` written from the midpoint rule on the integral representation.
- `Race.density()` is Gil-Pelaez on the characteristic function
  e^(iCt) ∏ J₀(r t) · e^(−V t²/4); `pdf`, `sample` and `variance` alongside.
- `Race.explicit(x)` rebuilds φ(q)(π(x; q, a) − π(x; q, b)) log x/√x from the
  zeros with Riemann's Li(x^ρ), the constant log L(0, χ) or
  log L'(0, χ) − γ − log log x, the Euler factors at p | q for imprimitive
  characters, and the k ≥ 2 Möbius terms counted directly from primes up to √x.
- `prime_race` and `observed_density` for the actual counts.
- Plots: the mod 4 race against 3,130 zeros with a zoom on the first lead change,
  and the limiting densities of five races.

Checks:

- J₀ against mpmath to 1e-14 for 0 ≤ z ≤ 80.
- L'/L(1, χ₋₄) against γ + log 2π + 2 log(Γ(3/4)/Γ(1/4)) to 1e-13.
- Σ_{γ>0} 1/(1/4 + γ²) for χ₋₃, χ₋₄, χ₋₈, χ₁₂ against Feuerverger-Martin Table
  3.2 to all six digits, and the part above T = 300 against the zero density
  (log(qT/2π) + 1)/(2πT) to 2%.
- δ for twelve races mod 3, 4, 5, 8, 12, 24 against Fiorilli-Martin Table 8 to
  six digits at T = 1000, eight of them in the tests at T = 200.
- 400,000 Monte Carlo draws of X against the Fourier density (within 4σ), the
  pdf integrating to 1 with mean C, δ(b, a) = 1 − δ(a, b), and the four N vs S
  races mod 5 equal to 1e-10.
- First lead of 1 mod 4 at 26861 (Leech), no lead for 1 mod 3 below 10⁶.
- Explicit formula against the real race for 200 ≤ x ≤ 10⁵ at T = 200: mean
  error under 0.03 for mod 4, 5 and 12.

Results: δ(4; 3, 1) = 0.9959279, δ(3; 2, 1) = 0.9990632, δ(5; 2, 1) = 0.9521746
and the most lopsided two-way race in the table, δ(24; 5, 1) = 0.9999875. One
race mod 4 at T = 1000 takes about 0.8 s, mod 24 (3,820 zeros) about 4 s. Up to
10⁸, 1 mod 4 leads after 1,940 primes and the share of log x in [10³, 10⁸] where
3 leads is 0.99894.

Design notes:

- Feuerverger and Martin (2000) print 0.952140 for the mod 5 race; this code gives
  0.952175, which is Fiorilli and Martin's 2013 value. The older paper says it
  interpolated F_T between grid points for complex characters.
- The tail is exact in variance and Gaussian in shape. The next cumulant of the
  tail is of order Σ r⁴ ~ 1/T³, which is why T = 200 already gives six digits.
- The trapezoid rule on the even analytic integrand only errs by aliasing, so the
  step is set from the spread of X (period 40σ + |C|) and the sum stops once the
  characteristic function falls below 1e-17.
- The first version of `explicit` used Ei(ρ log x) and was off by ±0.96 in E for
  the mod 5 race, with the sign flipping with T. Ei carries iπ·sign(Im z), which
  cancels between γ and −γ for real characters but not for complex ones, whose
  zeros above and below the axis differ. Riemann's Li(x^ρ) is the integral from
  −∞ and drops that constant.

Next: Dirichlet L-functions in the browser, then races with three or more
residue classes (where δ depends on the order, not just the set).

## 2026-11-07: Dirichlet L-functions in the browser

New section on the web page: pick a modulus q ≤ 200 and a primitive character by
its Conrey label, see Z(t, χ) with its zeros, and below it the zeros with |t| < 20
of every primitive character mod q, one row per character. Clicking a row plots it.

- `web/js/characters.js`: Conrey characters as exponent tables over a common
  denominator, built prime power by prime power (the least generator that is also
  primitive mod p², and the ±5^k split mod 2^e). Order, parity, conductor, conjugate
  label n⁻¹ mod q, Gauss sum, root number ε = τ(χ)/(i^a √q).
- `web/js/lfunction.js`: L(s, χ) = q^(−s) Σ χ(a) ζ(s, a/q), head summed to N·q with
  N = 10 + ⌊|t|/π⌋ and an Euler-Maclaurin tail per residue. θ(t, χ), Z(t, χ), the
  zero finder (8 samples per mean gap, Illinois refinement) and the smooth count
  (θ(T) − θ(−T))/π.
- `web/js/lzeros-worker.js` and `dirichlet.js`: the rows go to up to four workers,
  one job per conjugate pair since Z(−t, χ) = Z(t, χ̄).

Checks:

- Character values against `characters.py` for 14 moduli up to 99, order, parity
  and conductor for ten labels against the LMFDB, multiplicativity over all pairs
  mod 63, 64, 99 and 120, the count of primitive characters for eleven moduli
  (including 0 for 30), |τ| = √q, ε = 1 for every real primitive character.
- L(s, χ) at 12 points with mpmath's `dirichlet` to 1e-12 relative, L(1, χ₋₄) = π/4
  to 1e-14, Z against the Python `l_z` to 1e-11.
- Im(e^(iθ) L(1/2 + it, χ)) = 0 to 1e-11 for every primitive character mod 5, 7, 16,
  21 and 37 at five heights, and Z(−t, χ) = Z(t, χ̄).
- First three zeros of χ₃(2, ·), χ₄(3, ·), χ₅(2, ·) against the LMFDB to 1e-10, and
  the exact zero counts from the Python argument-principle check (31, 24, 24, 25
  and 15 zeros for 5.2, 7.3, 13.2, 40.3 and 99.2).
- In the browser: q = 30 says there are no primitive characters, 13.5 shows 18
  zeros below 20 with the lowest at 2.195553, the same as the Python code.

Timings (Node 24, one thread): one Z(t, χ) costs 12 μs for q = 13 at t = 20, 96 μs
for q = 163, 265 μs for q = 199 at t = 100. All 81 distinct rows mod 163 take 6.3 s
on one thread and about 2.3 s in the page with four workers.

Design notes:

- The tail term w^(1−s)/(s − 1) is computed as −log w · expm1(u)/u with
  u = (1 − s) log w, the same trick as the Python version. The dropped 1/(s − 1)
  cancels because Σ χ(a) = 0. Without it L(1 + 1e-9, χ₋₄) was off by 9e-8.
- Caching χ(n) n^(−σ) per character took the critical-line evaluation from four
  transcendental calls per term to two, but the Euler-Maclaurin tail (10 to 20
  terms per residue) costs about as much as the head, so the real speedup for the
  rows came from the conjugate symmetry and the worker pool.
- `conreyGenerator(40487)` first returned 5: the check 5^(p−1) mod p² overflowed
  doubles at p² ≈ 1.6e9. 40487 is the classic prime where 5 is a primitive root
  mod p but not mod p², so the right answer is 10; powmod now runs in BigInt.

Next: odd quadratic twists against excised SO(2N+1), then elliptic curves in the
browser.

## 2026-11-07: Odd twists against excised SO(2N + 1), and the heat flow

Built:

- `rmt.so_angles(..., odd=True)`: Haar SO(2N + 1) by the same QR recipe, dropping
  the forced eigenvalue at 1. `excised_first` takes the same flag; for SO(2N + 1)
  the product ∏(2 − 2 cos θ_j) over the free angles is |Λ_A'(1)|.
- `rmt.so_moment` and `rmt.small_value_law`: E|Λ|^s in closed form for both
  parities from Selberg's integral, and the residue at the first pole, which gives
  P(Λ ≤ ρ) ≈ (R/α)ρ^α with α = 1/2 (SO(2N)) or 3/2 (SO(2N + 1)). `binned_chi2` for
  comparing samples on fixed bins.
- `elliptic.arithmetic_factor`: the CFKRS factor a_s(E) for prime conductor, as in
  section 5 of DHKMS.
- `heat.py`: Φ(u), H_t(z) and its first two z-derivatives on a shifted contour,
  real zeros, complex Newton, `collision` (the double zero where a pair leaves the
  line), `zero_paths` (every zero in a window as t decreases, complex ones
  included), `zero_velocity` (H''/H') and `velocities` (the particle sum over a
  list of zeros plus a smooth tail).
- Plots: the odd twists against excised SO(17) with χ² against the cut, the zeros
  of H_t near the pair γ = 564.161, 564.506 with the squared gap against
  8(t − t_c), and a GIF of H_t and its zeros in the complex plane.

Checks:

- a_(−1/2)(11a) = 0.73269 with the primes below 2^14, DHKMS print 0.732728078;
  a_0 = 1; the conductor factor between the two signs is exactly (12/10)^(3/2).
- SO(3) moments in closed form to 1e-13, E Λ_A(1) = 2 for every SO(2N), Monte
  Carlo moments for both parities within 4σ, the residue against 1/(4π) for SO(3)
  and against DHKMS's Barnes-G asymptotic at N = 100, the small-value count of
  60,000 SO(9) matrices within 4σ. One-level density and first-eigenvalue law of
  the SO(2N + 1) sampler against the kernel and the Fredholm determinant.
- CFKRS dictionary on the family (N rounded per twist): the first moment of the
  even central values is 1.4812 against 1.4828 predicted, and the ±1/2 moments of
  L'(E_d, 1/2) on the odd side come out 2% from a_s(E) times the SO(2N + 1) moment.
  The s = 1 and s = −1 moments are off by 20% and 16%.
- H_0 = ξ(1/2 + iz/2)/8 to 2e-12 for real and complex z (2.6e-11 at x = 1128,
  where H = 5.4e-191), the defining integral by mpmath quadrature at 30 digits to
  1e-10 for t from −1 to 0.25, the backward heat equation, Φ(u) = Φ(−u),
  conjugate symmetry, zeros of H_0 at 2γ.
- Collision of the pair: H and H' vanish at (t_c, x_c) to 1e-13 relative, squared
  gap over 8(t − t_c) equal to 1 within 5e-4 at distance 1e-3 and 5e-5 at 1e-4,
  the same for −4(Im z)² after the collision.
- Velocities from 4,000 zeros plus the smooth tail against H''/H' at the seven
  zeros in 1119 < x < 1139: 6.8e-7. H''/H' against finite differences of the zero
  positions in t: 1e-7.

Results: the 3,345 twists with 5 ≤ |d| ≤ 6000 took 591 s on 14 processes, and
the analysis uses the 1,674 with d > 0, as before. Smallest √d·L'(E_d, 1/2) is
6.437 at d = 8, a_(−3/2) = 0.2030, and the arithmetic cut is c = 0.265 (densities)
or 2.22 (probabilities) in |Λ_A'(1)| ≥ c·e^(−N/2). At N = 8 that removes 0.13% of
SO(17) and χ² goes from 148.4 to 147.1. The minimum of χ² over c is 61.5 at
c = 112, keeping 83% of the matrices, with KS 0.056, still above the critical
0.047. SO(5), SO(9) and SO(13) with their best cuts do worse. The heat flow
collision is at t_c = −0.0614041, x_c = 1128.71107; `lehmer.newman_bound` for the
same pair gives −0.0616036.

Design notes:

- The direct cosine integral loses e^(πx/8) to cancellation. Moving the contour to
  Im u = π/8 − δ costs only e^(δx), and the trapezoid rule on the shifted line has
  aliasing error e^(−2πδ/h), so h = δ/8 is safe. δ ∝ 1/x keeps the cancellation
  bounded, the price is node count ∝ x. At x = 1130 a point costs about 1 ms in a
  batch.
- Φ is real and even, so the integral over u < 0 is the conjugate of the one over
  u > 0. That halves the nodes and makes H exactly real on the real axis.
- The first zero scan found nothing near x = 1130: the sign test multiplied two
  values near 1e-191 and underflowed to zero. It compares signs now.
- Measuring a zero's velocity by differencing positions in t works, but the heat
  equation gives it exactly: ẋ = −H_t/H_x = H_xx/H_x. Same integral, one more
  power of iw.
- Matching the cut to the observed fraction of rank 3 twists (2 of 842) gives the
  same picture as the probability recipe: almost nothing is excised. Whatever
  repels the first odd zero, it is not the size of the central derivative.

Next: lower-order terms of the one-level density for the odd family, so that an
effective N can be computed instead of fitted. Then the Polymath 15 approximation
of H_t to reach Lehmer's pair at x ≈ 14010, and elliptic curves in the browser.

## 2026-11-08: Moments, the arithmetic factor, and Li's coefficients

Built:

- `moments.py`: the CFKRS moment polynomials (P_1 = x + 2γ in closed form, P_2 and
  P_3 as printed in the paper), their leading coefficient a_k·G(k + 1)²/G(2k + 1)
  from the Euler product, exact integrals of P_k(log(t/2π)) over blocks and against
  e^(−t/T), and `z_power_integrals`, which integrates Z^(2k) on long ranges through
  the Odlyzko-Schönhage windows with Gregory end corrections.
- `values.arithmetic_cumulants` and `hybrid_density`: the arithmetic factor of the
  moment conjecture read as a correction to the first four cumulants of log|ζ|, and
  the CUE density with that correction. The value distribution chart got a second
  panel with histogram minus CUE against the predicted correction.
- `li.py`: Li's coefficients λ_n two ways, as a sum over zeros with a smooth tail
  and a boundary term, and as Taylor coefficients of log ξ(1/(1 − z)) from a
  discrete Fourier transform on |z| = 1/2 at 40 digits.
- Charts: moments over windows of 10⁵ up to T = 10¹⁰ against the full polynomials
  and the leading term, λ_n up to 20,000 with the residual after Voros's trend.

Checks:

- Leading coefficients: 1/(2π²) for k = 2 and CFKRS's 5.708527034652e−6 for k = 3,
  both to 3e−11. a_2 = 6/π² to 1e−10.
- Polynomial integrals reproduce the conjecture columns of CFKRS's tables 5.1.2 and
  5.1.3 to every printed digit, and the smoothed one matches mpmath quadrature to
  1e−12.
- Gregory weights: (3/8, 7/6, 23/24) for three points, and exact on x^d, d < 8,
  for eight.
- Z^(2k) integrals against Gauss-Legendre on mpmath's siegelz: 1e−11 or better on
  3-unit blocks at 60000, at the record values near 50626 and 77403 (Z⁶ above
  10⁸), and across the seam at t = 1000 where the code switches from ζ to the
  windows. Window length and step don't move a long integral beyond 5e−12.
- Arithmetic κ₂ against ½Σ_p (Li₂(1/p) + log(1 − 1/p)) with mpmath's polylog, κ₃
  and κ₄ against cumulants of −log|1 − e^(iθ)/√p| by trapezoid in θ, all to
  1e−10. The hybrid density integrates to 1 within 1e−8 with the predicted variance
  and skew.
- λ_1 = 1 + γ/2 − ½ log 4π to 1e−15. The two routes to λ_n agree to 1.3e−12 for
  n ≤ 40 with a million zeros, and to 1e−6 with 1,000 zeros. Dropping half of the
  million zeros moves λ_n by 1.6e−7 at most for n up to 10⁵.

Results:

- The arithmetic factor contributes κ₂ = −0.0881, κ₃ = +0.2337, κ₄ = −0.0989. On the
  same 200,000-point samples as before (N = 8 to 16), CUE plus this correction gets
  the variance within 0.015 of the data at every height (the sampling error is
  0.008) and moves the predicted skew from −2.4 to −2.14..−2.19, measured
  −2.10..−2.20. That closes the gap the last value-distribution entry left open,
  with nothing fitted.
- Smoothed moments at T = 10⁴: 79499.9312639, 5088336.51209, 708967539.4 for k = 1,
  2, 3. CFKRS print 79499.9312635 for k = 1. Their k = 2 row has the two numbers
  the other way round: my conjecture integral is 5088332.55512, which they list as
  the measured value. Their sixth-moment blocks are 2e−7 to 1.1e−6 below mine,
  growing with k while k = 1 agrees to 12 digits.
- Over windows [T, T + 10⁵] up to 10¹⁰ the second moment stays within 8e−4 of P_1,
  the fourth within 0.5% up to 10⁹, and the sixth within 0.7% up to 3·10⁶ and then
  within 11%. The leading term alone is off by a factor 105 for k = 3 at 10³ and
  still 11 at 10¹⁰.
- The first 1,000,001 zeros (to g_(10⁶) ≈ 600270) in 148 s. λ_n is positive and
  increasing for every n ≤ 20,000 and at 120 log-spaced n up to 10⁵. The residual
  after (n/2)(log n − 1 + γ − log 2π) has standard deviation 0.43, 2.8, 4.4 and 7.5
  on [1, 100], [100, 1000], [1000, 5000], [5000, 20000].

Design notes:

- The trapezoid rule on Z^(2k) looked like it converged at second order, which
  made no sense for a band-limited integrand. It was all endpoint error: f' at the
  ends is in the hundreds. Gregory corrections derived from Euler-Maclaurin (solve
  Σ c_j j^d = B_(d+1)/(d+1) for odd d, 0 for even d, in exact rationals) fix it
  without touching the interior weights.
- The full arithmetic factor can't be inverted into a density on its own: its
  variance is negative, so the characteristic function grows like e^(0.044u²).
  Keeping its first four cumulants gives a transform that decays because κ₄ < 0.
  The first attempt multiplied the full quadrature-based factor into the CUE
  transform and the density integrated to 1.022.
- λ_n via mpmath.taylor fails (no ξ in mpmath, and finite differences at z = 0 hit
  the pole of ζ). The DFT on a circle avoids the point and gives all coefficients
  from one set of 160 samples. The phase of log ξ is unwrapped along the circle;
  the winding number is 0 because nothing vanishes inside.
- The tail of λ_n is the integral in φ = arctan(1/2t) with φ = φ_T y⁴, which turns
  the log singularity of θ' at φ = 0 into y³ log y for Gauss-Legendre.

Next: the Fyodorov-Hiary-Keating maximum of |ζ| in short intervals using the same
window machinery, then Jensen polynomials of ξ, and interval arithmetic for
certified zero counts.
