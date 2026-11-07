export function gcd(a, b) {
  while (b) [a, b] = [b, a % b]
  return Math.abs(a)
}

const lcm = (a, b) => (a / gcd(a, b)) * b

const powmod = (b, e, m) => {
  let r = 1n
  ;[b, e, m] = [BigInt(b) % BigInt(m), BigInt(e), BigInt(m)]
  for (; e > 0n; e >>= 1n) {
    if (e & 1n) r = (r * b) % m
    b = (b * b) % m
  }
  return Number(r)
}

export function factor(n) {
  if (!(Number.isInteger(n) && n >= 1)) throw new RangeError("n must be an integer >= 1")
  const out = []
  for (let p = 2; p * p <= n; p += p === 2 ? 1 : 2) {
    if (n % p) continue
    let e = 0
    while (n % p === 0) {
      n /= p
      e++
    }
    out.push([p, e])
  }
  if (n > 1) out.push([n, 1])
  return out
}

export function conreyGenerator(p) {
  const primes = factor(p - 1).map(([r]) => r)
  for (let g = 2; g < p; g++) {
    if (primes.some((r) => powmod(g, (p - 1) / r, p) === 1)) continue
    if (powmod(g, p - 1, p * p) !== 1) return g
  }
  throw new RangeError(`${p} is not an odd prime`)
}

function primePowerPart(p, e, n) {
  const pe = p ** e
  if (p === 2) {
    if (e === 1) return { num: new Int32Array(2), den: 1 }
    const den = 2 ** (e - 1)
    const eps = new Int8Array(pe)
    const idx = new Int32Array(pe)
    for (let k = 0, x = 1; k < Math.max(1, pe / 4); k++, x = (x * 5) % pe) {
      eps[x] = 1
      idx[x] = k
      eps[pe - x] = -1
      idx[pe - x] = k
    }
    const num = new Int32Array(pe)
    for (let m = 1; m < pe; m += 2) {
      const flip = eps[m] === -1 && eps[n] === -1 ? den / 2 : 0
      num[m] = (flip + 2 * idx[m] * idx[n]) % den
    }
    return { num, den }
  }
  const phi = pe - pe / p
  const idx = new Int32Array(pe)
  const g = conreyGenerator(p)
  for (let k = 0, x = 1; k < phi; k++, x = (x * g) % pe) idx[x] = k
  return { num: idx.map((k) => (k * idx[n]) % phi), den: phi }
}

export class Character {
  constructor(q, num, den, n = null) {
    let g = den
    for (const k of num) if (k >= 0) g = gcd(g, k)
    this.q = q
    this.n = n
    this.den = den / g
    this.num = num.map((k) => (k >= 0 ? k / g : -1))
  }

  get order() {
    return this.den
  }

  get isReal() {
    return this.den <= 2
  }

  get isPrincipal() {
    return this.den === 1
  }

  get parity() {
    return this.q < 3 || this.num[this.q - 1] === 0 ? 0 : 1
  }

  value(m) {
    const k = this.num[((m % this.q) + this.q) % this.q]
    if (k < 0) return [0, 0]
    return rootOfUnity(k, this.den)
  }

  get conductor() {
    this._conductor ??= this.#conductor()
    return this._conductor
  }

  #conductor() {
    for (let f = 1; f <= this.q; f++) {
      if (this.q % f) continue
      let trivial = true
      for (let m = 1 % f; m < this.q && trivial; m += f) trivial = this.num[m] <= 0
      if (trivial) return f
    }
  }

  get isPrimitive() {
    return this.conductor === this.q
  }

  conj() {
    const num = this.num.map((k) => (k >= 0 ? (this.den - k) % this.den : -1))
    return new Character(this.q, num, this.den, this.n === null ? null : inverse(this.n, this.q))
  }

  gaussSum() {
    let re = 0
    let im = 0
    for (let a = 1; a < this.q; a++) {
      const k = this.num[a]
      if (k < 0) continue
      const [cr, ci] = rootOfUnity(k, this.den)
      const [er, ei] = rootOfUnity(a, this.q)
      re += cr * er - ci * ei
      im += cr * ei + ci * er
    }
    return [re, im]
  }

  rootNumber() {
    if (!this.isPrimitive) throw new RangeError("the root number needs a primitive character")
    const [re, im] = this.gaussSum()
    const r = Math.sqrt(this.q)
    return this.parity ? [im / r, -re / r] : [re / r, im / r]
  }
}

function inverse(n, q) {
  for (let m = 1; m <= q; m++) if ((m * n) % q === 1 % q) return m
}

export function rootOfUnity(k, n) {
  k %= n
  if ((4 * k) % n === 0) return [[1, 0], [0, 1], [-1, 0], [0, -1]][(4 * k) / n]
  const a = (2 * Math.PI * k) / n
  return [Math.cos(a), Math.sin(a)]
}

export function conrey(q, n) {
  if (!(Number.isInteger(q) && q >= 1 && n >= 1 && n <= Math.max(q, 1) && gcd(n, q) === 1)) {
    throw new RangeError("need q >= 1 and 1 <= n <= q with gcd(n, q) = 1")
  }
  const ps = factor(q)
  const parts = ps.map(([p, e]) => primePowerPart(p, e, n % p ** e))
  const den = parts.reduce((d, part) => lcm(d, part.den), 1)
  const num = new Int32Array(q)
  for (let m = 0; m < q; m++) {
    if (gcd(m, q) !== 1) {
      num[m] = -1
      continue
    }
    let k = 0
    ps.forEach(([p, e], i) => {
      k += parts[i].num[m % p ** e] * (den / parts[i].den)
    })
    num[m] = k % den
  }
  return new Character(q, num, den, n % q || q)
}

export function primitiveCharacters(q) {
  const out = []
  for (let n = 1; n < Math.max(q, 2); n++) {
    if (gcd(n, q) !== 1) continue
    const chi = conrey(q, n)
    if (chi.isPrimitive) out.push(chi)
  }
  return out
}

export function conjugatePairs(chars) {
  const seen = new Set()
  const out = []
  for (const chi of chars) {
    if (seen.has(chi.n)) continue
    const bar = chi.isReal ? null : chi.conj().n
    seen.add(chi.n).add(bar)
    out.push([chi.n, bar])
  }
  return out
}
