import { factor, gcd } from "./characters.js"
import { illinois, loggamma } from "./zeta.js"

const CUT = 36
const LOSS = 3
const EPS = 2 ** -52
const TINY = 1e-300
const MAX_TERMS = 10000
const EULER_GAMMA = 0.5772156649015329

const cmul = (a, b, c, d) => [a * c - b * d, a * d + b * c]
const cdiv = (a, b, c, d) => {
  const m = c * c + d * d
  return [(a * c + b * d) / m, (b * c - a * d) / m]
}
const cexp = (a, b) => {
  const m = Math.exp(a)
  return [m * Math.cos(b), m * Math.sin(b)]
}
const mod = (a, m) => ((a % m) + m) % m

export function curve(a1, a2, a3, a4, a6) {
  const [A1, A2, A3, A4, A6] = [a1, a2, a3, a4, a6].map(BigInt)
  const b2 = A1 * A1 + 4n * A2
  const b4 = 2n * A4 + A1 * A3
  const b6 = A3 * A3 + 4n * A6
  const b8 = A1 * A1 * A6 + 4n * A2 * A6 - A1 * A3 * A4 + A2 * A3 * A3 - A4 * A4
  const c4 = b2 * b2 - 24n * b4
  const disc = -b2 * b2 * b8 - 8n * b4 ** 3n - 27n * b6 * b6 + 9n * b2 * b4 * b6
  if (disc === 0n) throw new RangeError("singular curve")
  const size = disc < 0n ? -disc : disc
  if (size > 10n ** 15n) throw new RangeError("discriminant too large to factor here")
  const bad = factor(Number(size)).map(([p]) => p)
  const additive = bad.filter((p) => c4 % BigInt(p) === 0n)
  if (additive.length) {
    throw new RangeError(`additive reduction at ${additive.join(", ")}: only semistable curves`)
  }
  const E = {
    a: [a1, a2, a3, a4, a6],
    b: [b2, b4, b6],
    disc: Number(disc),
    bad,
    N: bad.reduce((x, p) => x * p, 1),
    coef: new Float64Array([0, 1]),
  }
  E.w = -bad.reduce((x, p) => x * -ap(E, p), 1)
  return E
}

export function ap(E, p) {
  if (p === 2) {
    const [a1, a2, a3, a4, a6] = E.a
    let pts = 0
    for (let x = 0; x < 2; x++) {
      for (let y = 0; y < 2; y++) {
        pts += mod(y * y + a1 * x * y + a3 * y - x ** 3 - a2 * x * x - a4 * x - a6, 2) === 0
      }
    }
    return 2 - pts
  }
  const P = BigInt(p)
  const [b2, b4, b6] = E.b.map((b) => Number(((b % P) + P) % P))
  const square = new Int8Array(p).fill(-1)
  for (let x = 1; x < p; x++) square[(x * x) % p] = 1
  square[0] = 0
  let f = b6
  let d1 = (4 + b2 + 2 * b4) % p
  let d2 = (24 + 2 * b2) % p
  const d3 = 24 % p
  let s = 0
  for (let x = 0; x < p; x++) {
    s -= square[f]
    f += d1
    if (f >= p) f -= p
    d1 += d2
    if (d1 >= p) d1 -= p
    d2 += d3
    if (d2 >= p) d2 -= p
  }
  return s
}

export function coefficients(E, M) {
  if (E.coef.length > M) return E.coef
  const size = Math.max(M + 1, 2 * E.coef.length)
  const spf = new Int32Array(size)
  for (let p = 2; p < size; p++) {
    if (spf[p]) continue
    for (let q = p; q < size; q += p) if (!spf[q]) spf[q] = p
  }
  const a = new Float64Array(size)
  const rest = new Int32Array(size)
  const bad = new Set(E.bad)
  a[1] = 1
  rest[1] = 1
  for (let n = 2; n < size; n++) {
    const p = spf[n]
    const m = n / p
    rest[n] = m % p === 0 ? rest[m] : m
    if (rest[n] !== 1) {
      a[n] = a[n / rest[n]] * a[rest[n]]
    } else if (m === 1) {
      a[n] = n < E.coef.length ? E.coef[n] : ap(E, p)
    } else {
      a[n] = bad.has(p) ? a[p] * a[m] : a[p] * a[m] - p * a[m / p]
    }
  }
  E.coef = a
  return a
}

export function kronecker(d, n) {
  if (n === 0) return Math.abs(d) === 1 ? 1 : 0
  let k = 1
  while (n % 2 === 0) {
    n /= 2
    if (d % 2 === 0) return 0
    if (mod(d, 8) === 3 || mod(d, 8) === 5) k = -k
  }
  let a = mod(d, n)
  while (a !== 0) {
    while (a % 2 === 0) {
      a /= 2
      if (n % 8 === 3 || n % 8 === 5) k = -k
    }
    ;[a, n] = [n, a]
    if (a % 4 === 3 && n % 4 === 3) k = -k
    a %= n
  }
  return n === 1 ? k : 0
}

export function isFundamental(d) {
  if (!Number.isInteger(d) || d === 0 || d === 1) return false
  const squarefree = (m) => factor(Math.abs(m)).every(([, e]) => e === 1)
  if (mod(d, 4) === 1) return squarefree(d)
  if (mod(d, 4) !== 0) return false
  const m = d / 4
  return (mod(m, 4) === 2 || mod(m, 4) === 3) && squarefree(m)
}

export function twist(E, d = 1) {
  if (d !== 1 && !isFundamental(d)) throw new RangeError(`${d} is not a fundamental discriminant`)
  if (gcd(d, E.N) !== 1) throw new RangeError("the twist needs gcd(d, N) = 1")
  const q = Math.abs(d)
  const chi = Int8Array.from({ length: q }, (_, n) => (d === 1 ? 1 : kronecker(d, n)))
  const N = E.N * d * d
  const w = d === 1 ? E.w : E.w * chi[mod(-E.N, q)]
  return { E, d, N, w, chi, A: Math.sqrt(N) / (2 * Math.PI) }
}

export function an(L, M) {
  const a = coefficients(L.E, M)
  const q = L.chi.length
  return Float64Array.from({ length: M + 1 }, (_, n) => a[n] * L.chi[n % q])
}

function lentz(br, bi, zr, zi) {
  let [cr, ci] = [1 / TINY, 0]
  let [dr, di] = cdiv(1, 0, br, bi)
  let [hr, hi] = [dr, di]
  for (let k = 1; k < MAX_TERMS; k++) {
    const ar = -k * (k - zr)
    const ai = k * zi
    br += 2
    ;[dr, di] = cmul(ar, ai, dr, di)
    ;[dr, di] = cdiv(1, 0, dr + br, di + bi)
    const [qr, qi] = cdiv(ar, ai, cr, ci)
    ;[cr, ci] = [br + qr, bi + qi]
    const [er, ei] = cmul(cr, ci, dr, di)
    ;[hr, hi] = cmul(hr, hi, er, ei)
    if (Math.hypot(er - 1, ei) <= 2 * EPS) return [hr, hi]
  }
  throw new RangeError("Legendre's continued fraction did not converge")
}

export function e1(x) {
  if (!(x > 0)) throw new RangeError("e1 needs x > 0")
  if (x >= 1) return Math.exp(-x) * lentz(x + 1, 0, 0, 0)[0]
  let term = -x
  let total = -x
  for (let k = 2; k < 40; k++) {
    term = (-term * x) / k
    total += term / k
  }
  return -EULER_GAMMA - Math.log(x) - total
}

export function gammaQ(zr, zi, wr, wi, lg = loggamma(zr, zi)) {
  const lwr = 0.5 * Math.log(wr * wr + wi * wi)
  const lwi = Math.atan2(wi, wr)
  const pr = zr * lwr - zi * lwi - wr - lg[0]
  const pi = zr * lwi + zi * lwr - wi - lg[1]
  if (Math.hypot(wr, wi) < Math.hypot(zr, zi) + 1) {
    let [tr, ti] = cdiv(1, 0, zr, zi)
    let [sr, si] = [tr, ti]
    for (let k = 1; ; k++) {
      if (k === MAX_TERMS) throw new RangeError("incomplete gamma series did not converge")
      ;[tr, ti] = cmul(tr, ti, ...cdiv(wr, wi, zr + k, zi))
      sr += tr
      si += ti
      if (Math.hypot(tr, ti) <= EPS * Math.hypot(sr, si)) break
    }
    const [xr, xi] = cmul(...cexp(pr, pi), sr, si)
    return [1 - xr, -xi]
  }
  return cmul(...cexp(pr, pi), ...lentz(wr + 1 - zr, wi - zi, zr, zi))
}

export function central(L) {
  const M = Math.ceil(CUT * L.A) + 1
  const a = an(L, M)
  let s = 0
  for (let n = 1; n <= M; n++) {
    if (a[n]) s += (a[n] / n) * (L.w === 1 ? Math.exp(-n / L.A) : e1(n / L.A))
  }
  return 2 * s
}

export function centralOrder(L, tol = 1e-6) {
  const value = central(L)
  return { value, order: (L.w === -1 ? 1 : 0) + (Math.abs(value) < tol ? 2 : 0) }
}

export function ecTheta(t, L) {
  return t * Math.log(L.A) + loggamma(1, t)[1]
}

export function hardyZ(t, L) {
  const ang = (Math.abs(t) * Math.PI) / 2 <= LOSS ? Math.PI / 2 : LOSS / Math.abs(t)
  const phi = Math.sign(t) * (Math.PI / 2 - ang)
  const [yr, yi] = [Math.cos(phi) / L.A, Math.sin(phi) / L.A]
  const M = Math.ceil(((CUT + LOSS) * L.A) / Math.sin(ang)) + 1
  const a = an(L, M)
  const lg = loggamma(1, t)
  let sr = 0
  let si = 0
  for (let n = 1; n <= M; n++) {
    if (!a[n]) continue
    const [qr, qi] = gammaQ(1, t, n * yr, n * yi, lg)
    const [er, ei] = cmul(Math.cos(t * Math.log(n)), -Math.sin(t * Math.log(n)), qr, qi)
    sr += (a[n] / n) * er
    si += (a[n] / n) * ei
  }
  const th = ecTheta(t, L)
  const [vr, vi] = cmul(Math.cos(th), Math.sin(th), sr, si)
  return 2 * (L.w === 1 ? vr : vi)
}

export function meanSpacing(t, L) {
  return Math.PI / Math.max(1, Math.log(L.A * Math.max(t, 1)))
}

export function zerosBetween(L, a, b, perSpacing = 10) {
  const out = []
  let t = Math.max(a, 1e-9)
  let z = hardyZ(t, L)
  while (t < b) {
    const next = Math.min(b, t + meanSpacing(t, L) / perSpacing)
    const zn = hardyZ(next, L)
    if (z * zn < 0) out.push(illinois((s) => hardyZ(s, L), t, next, z, zn))
    else if (zn === 0) out.push(next)
    t = next
    z = zn
  }
  return out
}

function monomial(c, s, first = false) {
  if (c === 0) return ""
  const body = Math.abs(c) === 1 && s ? s : `${Math.abs(c)}${s}`
  if (first) return c < 0 ? `−${body}` : body
  return `${c < 0 ? " − " : " + "}${body}`
}

export function equation([a1, a2, a3, a4, a6]) {
  const left = `y²${monomial(a1, "xy")}${monomial(a3, "y")}`
  const right = `x³${monomial(a2, "x²")}${monomial(a4, "x")}${monomial(a6, "")}`
  return `${left} = ${right}`
}

export function zRange(L, terms = 9000) {
  return Math.min(40, Math.max(10, terms / (13 * L.A)))
}

export function twistRange(E, terms = 30000) {
  return Math.min(1000, Math.floor((2 * Math.PI * terms) / (36 * Math.sqrt(E.N))))
}
