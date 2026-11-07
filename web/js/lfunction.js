import { B_OVER_FACT, M_MAX, TWO_PI, illinois, loggamma } from "./zeta.js"

const TAIL_TOL = 1e-18

function terms(chi, limit, sr) {
  let tab = chi._terms
  if (!tab || tab.limit < limit) {
    const cap = Math.max(limit, 2 * (tab?.limit ?? 0))
    const logn = []
    const cr = []
    const ci = []
    for (let n = 1; n < cap; n++) {
      if (chi.num[n % chi.q] < 0) continue
      const [re, im] = chi.value(n)
      logn.push(Math.log(n))
      cr.push(re)
      ci.push(im)
    }
    tab = chi._terms = {
      limit: cap,
      phi: chi.num.filter((k) => k >= 0).length,
      logn: Float64Array.from(logn),
      cr: Float64Array.from(cr),
      ci: Float64Array.from(ci),
      sr: NaN,
    }
  }
  if (tab.sr !== sr) {
    tab.sr = sr
    tab.ar = tab.cr.map((c, k) => c * Math.exp(-sr * tab.logn[k]))
    tab.ai = tab.ci.map((c, k) => c * Math.exp(-sr * tab.logn[k]))
  }
  return tab
}

export function expm1Ratio(x, y) {
  if (x === 0 && y === 0) return [1, 0]
  const s = Math.sin(y / 2)
  const er = Math.expm1(x) * Math.cos(y) - 2 * s * s
  const ei = Math.exp(x) * Math.sin(y)
  const n = x * x + y * y
  return [(er * x + ei * y) / n, (ei * x - er * y) / n]
}

export function emLength(t) {
  return 10 + Math.floor(Math.abs(t) / Math.PI)
}

export function lValue(sr, si, chi) {
  if (chi.isPrincipal || !chi.isPrimitive) {
    throw new RangeError("lValue needs a primitive non-principal character")
  }
  const q = chi.q
  const N = emLength(si)
  const cut = N * q
  const { logn, ar, ai, phi } = terms(chi, cut, sr)

  let hr = 0
  let hi = 0
  for (let k = 0; k < N * phi; k++) {
    const c = Math.cos(si * logn[k])
    const s = Math.sin(si * logn[k])
    hr += ar[k] * c + ai[k] * s
    hi += ai[k] * c - ar[k] * s
  }

  let tr = 0
  let ti = 0
  for (let a = 1; a < q; a++) {
    if (chi.num[a] < 0) continue
    const w = N + a / q
    const lw = Math.log(w)
    const m = Math.exp(-sr * lw)
    const wr = m * Math.cos(si * lw)
    const wi = -m * Math.sin(si * lw)
    const [d0, d1] = expm1Ratio((1 - sr) * lw, -si * lw)
    let xr = -lw * d0 + wr / 2
    let xi = -lw * d1 + wi / 2
    let p0 = sr
    let p1 = si
    let k0 = wr / w
    let k1 = wi / w
    for (let j = 1; j <= M_MAX; j++) {
      const er = B_OVER_FACT[j] * (p0 * k0 - p1 * k1)
      const ei = B_OVER_FACT[j] * (p0 * k1 + p1 * k0)
      xr += er
      xi += ei
      if (er * er + ei * ei < TAIL_TOL * TAIL_TOL) break
      const u = sr + 2 * j - 1
      const v = sr + 2 * j
      const f0 = u * v - si * si
      const f1 = si * (u + v)
      ;[p0, p1] = [p0 * f0 - p1 * f1, p0 * f1 + p1 * f0]
      k0 /= w * w
      k1 /= w * w
    }
    const [vr, vi] = chi.value(a)
    tr += vr * xr - vi * xi
    ti += vr * xi + vi * xr
  }
  const lq = Math.log(q)
  const m = Math.exp(-sr * lq)
  const gr = m * Math.cos(si * lq)
  const gi = -m * Math.sin(si * lq)
  return [hr + gr * tr - gi * ti, hi + gr * ti + gi * tr]
}

function epsArg(chi) {
  if (chi._epsArg === undefined) {
    const [re, im] = chi.rootNumber()
    chi._epsArg = Math.atan2(im, re)
  }
  return chi._epsArg
}

export function lTheta(t, chi) {
  const g = loggamma((0.5 + chi.parity) / 2, t / 2)[1]
  return (t / 2) * Math.log(chi.q / Math.PI) + g - epsArg(chi) / 2
}

export function lZ(t, chi) {
  const th = lTheta(t, chi)
  const [re, im] = lValue(0.5, t, chi)
  return Math.cos(th) * re - Math.sin(th) * im
}

export function lMeanSpacing(t, q) {
  return TWO_PI / Math.max(1, Math.log((q * Math.max(Math.abs(t), 1)) / TWO_PI))
}

export function lZeros(chi, a, b, perGap = 8) {
  const step = lMeanSpacing(Math.max(Math.abs(a), Math.abs(b)), chi.q) / perGap
  const n = Math.max(2, Math.ceil((b - a) / step))
  const Z = (t) => lZ(t, chi)
  const out = []
  let t0 = a
  let z0 = Z(a)
  for (let i = 1; i <= n; i++) {
    const t1 = a + ((b - a) * i) / n
    const z1 = Z(t1)
    if (z0 === 0) out.push(t0)
    else if (z0 * z1 < 0) out.push(illinois(Z, t0, t1, z0, z1))
    t0 = t1
    z0 = z1
  }
  return out
}

export function smoothCount(chi, T) {
  return (lTheta(T, chi) - lTheta(-T, chi)) / Math.PI
}
