export const TWO_PI = 2 * Math.PI
export const RS_FROM = 200

const LOG_2 = Math.log(2)
const LOG_PI = Math.log(Math.PI)
const HALF_LOG_TWO_PI = 0.5 * Math.log(TWO_PI)
const M_MAX = 30

function gcd(a, b) {
  if (a < 0n) a = -a
  if (b < 0n) b = -b
  while (b) [a, b] = [b, a % b]
  return a
}

function reduce([p, q]) {
  const g = gcd(p, q) || 1n
  return q < 0n ? [-p / g, -q / g] : [p / g, q / g]
}

export function bernoulli(n) {
  const B = [[1n, 1n]]
  for (let m = 1; m <= n; m++) {
    let p = 0n
    let q = 1n
    let c = 1n
    for (let k = 0; k < m; k++) {
      const [bp, bq] = B[k]
      ;[p, q] = reduce([p * bq + c * bp * q, q * bq])
      c = (c * BigInt(m + 1 - k)) / BigInt(k + 1)
    }
    B.push(reduce([-p, q * BigInt(m + 1)]))
  }
  return B
}

const B = bernoulli(2 * M_MAX + 2)
const ratio = ([p, q], d = 1n) => Number(p) / Number(q * d)

let fact = 1n
const B_OVER_FACT = [1]
for (let k = 1; k <= M_MAX + 1; k++) {
  fact *= BigInt(2 * k - 1) * BigInt(2 * k)
  B_OVER_FACT.push(ratio(B[2 * k], fact))
}
const STIRLING = Array.from({ length: 10 }, (_, i) => {
  const k = i + 1
  return ratio(B[2 * k], BigInt(2 * k * (2 * k - 1)))
})

function cmul(a, b, c, d) {
  return [a * c - b * d, a * d + b * c]
}

function cdiv(a, b, c, d) {
  const n = c * c + d * d
  return [(a * c + b * d) / n, (b * c - a * d) / n]
}

function clog(a, b) {
  return [0.5 * Math.log(a * a + b * b), Math.atan2(b, a)]
}

function cexp(a, b) {
  const m = Math.exp(a)
  return [m * Math.cos(b), m * Math.sin(b)]
}

export function loggamma(re, im) {
  if (!(re > 0)) throw new RangeError("loggamma needs Re z > 0")
  let sr = 0
  let si = 0
  while (re * re + im * im < 225) {
    const [lr, li] = clog(re, im)
    sr += lr
    si += li
    re += 1
  }
  const [lr, li] = clog(re, im)
  let [ir, ii] = cdiv(1, 0, re, im)
  const [i2r, i2i] = cmul(ir, ii, ir, ii)
  let ser = 0
  let sei = 0
  for (const c of STIRLING) {
    ser += c * ir
    sei += c * ii
    ;[ir, ii] = cmul(ir, ii, i2r, i2i)
  }
  const [ar, ai] = cmul(re - 0.5, im, lr, li)
  return [ar - re + HALF_LOG_TWO_PI + ser - sr, ai - im + sei - si]
}

export function eulerMaclaurin(sr, si, N, M = M_MAX, tol = 1e-16) {
  let vr = 0
  let vi = 0
  for (let n = 1; n < N; n++) {
    const L = Math.log(n)
    const m = Math.exp(-sr * L)
    vr += m * Math.cos(si * L)
    vi -= m * Math.sin(si * L)
  }
  const LN = Math.log(N)
  const mN = Math.exp(-sr * LN)
  const nr = mN * Math.cos(si * LN)
  const ni = -mN * Math.sin(si * LN)
  const [hr, hi] = cdiv(N * nr, N * ni, sr - 1, si)
  vr += hr + nr / 2
  vi += hi + ni / 2
  let [pr, pi] = [sr, si]
  let [wr, wi] = [nr / N, ni / N]
  let bound = 0
  for (let k = 1; k <= M + 1; k++) {
    const [tr, ti] = cmul(pr, pi, wr, wi).map((x) => x * B_OVER_FACT[k])
    const kr = sr + 2 * k - 1
    bound = (Math.hypot(kr, si) / kr) * Math.hypot(tr, ti)
    if (k > M || bound <= tol * Math.max(1, Math.hypot(vr, vi))) break
    vr += tr
    vi += ti
    ;[pr, pi] = cmul(pr, pi, ...cmul(kr, si, kr + 1, si))
    wr /= N * N
    wi /= N * N
  }
  return { re: vr, im: vi, bound }
}

function zetaRight(sr, si) {
  const { re, im } = eulerMaclaurin(sr, si, 10 + Math.floor(Math.abs(si) / Math.PI))
  return [re, im]
}

function logSin(a, b) {
  if (b >= 0) {
    const [er, ei] = cexp(-2 * b, 2 * a)
    const [lr, li] = clog(...cdiv(er - 1, ei, 0, 2))
    return [lr + b, li - a]
  }
  const [er, ei] = cexp(2 * b, -2 * a)
  const [lr, li] = clog(...cdiv(1 - er, -ei, 0, 2))
  return [lr - b, li + a]
}

export function logChi(sr, si) {
  if (!(sr < 1)) throw new RangeError("chi is implemented for Re s < 1")
  const [sinr, sini] = logSin((Math.PI * sr) / 2, (Math.PI * si) / 2)
  const [gr, gi] = loggamma(1 - sr, -si)
  return [sr * LOG_2 + (sr - 1) * LOG_PI + sinr + gr, si * LOG_2 + si * LOG_PI + sini + gi]
}

export function zeta(sr, si = 0) {
  if (sr === 1 && si === 0) return [Infinity, 0]
  if (sr >= 0) return zetaRight(sr, si)
  if (si === 0 && sr % 2 === 0) return [0, 0]
  const [cr, ci] = logChi(sr, si)
  const [zr, zi] = zetaRight(1 - sr, -si)
  const [wr, wi] = cmul(...cexp(cr, ci), zr, zi)
  return si === 0 ? [wr, 0] : [wr, wi]
}

export function theta(t) {
  return (
    (t / 2) * Math.log(t / TWO_PI) - t / 2 - Math.PI / 8 +
    1 / (48 * t) + 7 / (5760 * t ** 3) + 31 / (80640 * t ** 5) + 127 / (430080 * t ** 7)
  )
}

export function thetaPrime(t) {
  return (
    0.5 * Math.log(t / TWO_PI) -
    1 / (48 * t ** 2) - 7 / (1920 * t ** 4) - 155 / (80640 * t ** 6) - 889 / (430080 * t ** 8)
  )
}

export function thetaExact(t) {
  return loggamma(0.25, t / 2)[1] - (t / 2) * LOG_PI
}

export function zExact(t) {
  t = Math.abs(t)
  const th = thetaExact(t)
  const [zr, zi] = zeta(0.5, t)
  return Math.cos(th) * zr - Math.sin(th) * zi
}

function horner(c, x) {
  let y = 0
  for (let i = c.length - 1; i >= 0; i--) y = y * x + c[i]
  return y
}

export function rsRemainder(t, rs, corrections = 4) {
  const a = Math.sqrt(t / TWO_PI)
  const N = Math.floor(a)
  const z = a - N - 0.5
  let r = 0
  let ak = 1
  for (let k = 0; k <= corrections; k++) {
    r += horner(rs[k], z) / ak
    ak *= a
  }
  return ((N % 2 === 1 ? 1 : -1) * r) / Math.sqrt(a)
}

let LOG_N = new Float64Array(0)
let INV_SQRT_N = new Float64Array(0)

function grow(N) {
  if (N <= LOG_N.length) return
  const m = Math.max(N, 2 * LOG_N.length)
  LOG_N = Float64Array.from({ length: m }, (_, i) => Math.log(i + 1))
  INV_SQRT_N = Float64Array.from({ length: m }, (_, i) => 1 / Math.sqrt(i + 1))
}

export function zRs(t, rs, corrections = 4) {
  t = Math.abs(t)
  const N = Math.floor(Math.sqrt(t / TWO_PI))
  grow(N)
  const th = theta(t)
  let sum = 0
  for (let k = 0; k < N; k++) sum += Math.cos(th - t * LOG_N[k]) * INV_SQRT_N[k]
  return 2 * sum + rsRemainder(t, rs, corrections)
}

export function zFunction(rs) {
  return (t) => (Math.abs(t) < RS_FROM ? zExact(t) : zRs(t, rs))
}

export function gramPoint(n, tol = 1e-13) {
  if (n < 0) throw new RangeError("Gram points need n >= 0")
  const target = n * Math.PI
  let t = 20
  while (theta(t) < target) t *= 2
  for (let i = 0; i < 100; i++) {
    const dt = (theta(t) - target) / thetaPrime(t)
    t -= dt
    if (Math.abs(dt) < tol * t) return t
  }
  throw new Error(`Newton didn't converge for g_${n}`)
}

export function gramPoints(a, b) {
  const lo = Math.max(0, Math.ceil(theta(Math.max(a, 18)) / Math.PI))
  const out = []
  for (let n = lo; ; n++) {
    const g = gramPoint(n)
    if (g > b) return out
    if (g >= a) out.push([n, g])
  }
}

export function illinois(f, a, b, fa, fb, tol = 1e-12) {
  if (fa * fb > 0) throw new RangeError("root is not bracketed")
  for (let i = 0; i < 200; i++) {
    const c = b - (fb * (b - a)) / (fb - fa)
    const fc = f(c)
    if (fc === 0) return c
    if (fc * fb < 0) {
      a = b
      fa = fb
    } else {
      fa /= 2
    }
    b = c
    fb = fc
    if (Math.abs(b - a) < tol + 4 * Number.EPSILON * Math.abs(b)) return b
  }
  return b
}

export function meanSpacing(t) {
  return TWO_PI / Math.log(Math.max(t, 20) / TWO_PI)
}

export function findZeros(Z, a, b, perSpacing = 12) {
  const step = meanSpacing(b) / perSpacing
  const n = Math.max(2, Math.ceil((b - a) / step))
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
