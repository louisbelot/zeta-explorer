const LOG_2PI = Math.log(2 * Math.PI)

export function sieve(n) {
  const composite = new Uint8Array(n + 1)
  const primes = []
  for (let p = 2; p <= n; p++) {
    if (composite[p]) continue
    primes.push(p)
    for (let q = p * p; q <= n; q += p) composite[q] = 1
  }
  return primes
}

export function primePowers(n) {
  const out = []
  for (const p of sieve(n)) {
    const L = Math.log(p)
    for (let q = p; q <= n; q *= p) out.push([q, L])
  }
  return out.sort((a, b) => a[0] - b[0])
}

export function chebyshevPsi(xs) {
  const jumps = primePowers(Math.floor(Math.max(...xs)))
  const order = xs.map((_, i) => i).sort((i, j) => xs[i] - xs[j])
  const out = new Float64Array(xs.length)
  let k = 0
  let psi = 0
  for (const i of order) {
    while (k < jumps.length && jumps[k][0] <= xs[i]) psi += jumps[k++][1]
    out[i] = psi
  }
  return out
}

export function smoothPart(x) {
  return x - LOG_2PI - 0.5 * Math.log1p(-1 / (x * x))
}

export function addZeros(acc, xs, gamma, from, to) {
  for (let i = 0; i < xs.length; i++) {
    const L = Math.log(xs[i])
    const r = Math.sqrt(xs[i])
    let s = 0
    for (let k = from; k < to; k++) {
      const g = gamma[k]
      const a = g * L
      s += (0.5 * Math.cos(a) + g * Math.sin(a)) / (0.25 + g * g)
    }
    acc[i] -= 2 * r * s
  }
  return acc
}

export function psiFromZeros(xs, gamma, K = gamma.length) {
  const acc = Float64Array.from(xs, smoothPart)
  return addZeros(acc, xs, gamma, 0, K)
}
