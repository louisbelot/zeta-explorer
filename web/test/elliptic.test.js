import assert from "node:assert/strict"
import { readFileSync } from "node:fs"
import { test } from "node:test"

import {
  an,
  central,
  centralOrder,
  coefficients,
  curve,
  e1,
  equation,
  hardyZ,
  isFundamental,
  kronecker,
  twist,
  twistRange,
  zerosBetween,
  zRange,
} from "../js/elliptic.js"

const ref = JSON.parse(readFileSync(new URL("./reference.json", import.meta.url)))
const PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]

const close = (got, want, tol, what) =>
  assert.ok(Math.abs(got - want) <= tol, `${what}: got ${got}, want ${want}`)

const curves = ref.curves.map(([label, a]) => [label, curve(...a)])

test("conductor, root number and a_p match the python version", () => {
  ref.curves.forEach(([label, , N, w, aps], i) => {
    const E = curves[i][1]
    assert.equal(E.N, N, label)
    assert.equal(E.w, w, label)
    const a = coefficients(E, 60)
    assert.deepEqual(PRIMES.map((p) => a[p]), aps, label)
  })
})

test("a_n is multiplicative with the Hecke recursion at good primes", () => {
  for (const [label, E] of curves) {
    const a = coefficients(E, 2000)
    for (const [m, n] of [[4, 9], [8, 25], [7, 11], [27, 16], [13, 125]]) {
      assert.equal(a[m * n], a[m] * a[n], `${label}: a(${m}·${n})`)
    }
    for (const p of [3, 7, 13, 29]) {
      if (E.bad.includes(p)) continue
      assert.equal(a[p * p], a[p] * a[p] - p, `${label}: a(${p}²)`)
      assert.ok(Math.abs(a[p]) <= 2 * Math.sqrt(p), `${label}: Hasse at ${p}`)
    }
  }
})

test("central values match the python version and the LMFDB", () => {
  ref.curves.forEach(([label, , , , , value, rank], i) => {
    const { value: v, order } = centralOrder(twist(curves[i][1]))
    close(v, value, 1e-12, `central ${label}`)
    assert.equal(order, rank, label)
  })
  close(central(twist(curves[0][1])), 0.253841860855911, 1e-14, "L(11a, 1)")
  close(central(twist(curves[1][1])), 0.305999773834052, 1e-14, "L'(37a, 1)")
})

test("Hardy Z matches the python version", () => {
  ref.curves.forEach(([label, , , , , , , , zs], i) => {
    const L = twist(curves[i][1])
    ref.hardy_at.forEach((t, j) => close(hardyZ(t, L), zs[j], 1e-11, `Z(${t}, ${label})`))
  })
})

test("zeros below 20 match the python version", () => {
  ref.curves.forEach(([label, , , , , , , zeros], i) => {
    const got = zerosBetween(twist(curves[i][1]), 0, 20)
    assert.equal(got.length, zeros.length, label)
    got.forEach((t, j) => close(t, zeros[j], 1e-9, `${label} zero ${j + 1}`))
  })
})

test("twists match the python version", () => {
  for (const [i, d, N, w, value] of ref.twists) {
    const L = twist(curves[i][1], d)
    assert.equal(L.N, N, `d = ${d}`)
    assert.equal(L.w, w, `d = ${d}`)
    close(central(L), value, 1e-11 * Math.max(1, Math.abs(value)), `central at d = ${d}`)
  }
})

test("a twist's coefficients are a_n times the Kronecker symbol", () => {
  const E = curves[0][1]
  const L = twist(E, -7)
  const a = coefficients(E, 100)
  an(L, 100).forEach((v, n) => assert.equal(v, a[n] * kronecker(-7, n)))
})

test("Kronecker symbols", () => {
  assert.deepEqual([1, 2, 3, 4, 5].map((n) => kronecker(-4, n)), [1, 0, -1, 0, 1])
  assert.deepEqual([1, 2, 3, 4, 5, 6, 7].map((n) => kronecker(5, n)), [1, -1, -1, 1, 0, 1, -1])
  assert.equal(kronecker(8, 3), -1)
  assert.equal(kronecker(-3, 0), 0)
  for (const d of [-23, -20, 12, 13, 21, -47]) {
    for (let m = 1; m < 40; m++) {
      assert.equal(kronecker(d, m + Math.abs(d)), kronecker(d, m), `period of (${d}/·)`)
      for (const n of [3, 4, 10]) {
        assert.ok(kronecker(d, m * n) === kronecker(d, m) * kronecker(d, n), `(${d}/${m}·${n})`)
      }
    }
  }
})

test("fundamental discriminants up to 20", () => {
  const found = []
  for (let d = -20; d <= 20; d++) if (isFundamental(d)) found.push(d)
  assert.deepEqual(found, [-20, -19, -15, -11, -8, -7, -4, -3, 5, 8, 12, 13, 17])
})

test("E1 against tabulated values", () => {
  close(e1(0.5), 0.5597735947761608, 1e-15, "E1(0.5)")
  close(e1(2), 0.04890051070806112, 1e-16, "E1(2)")
  close(e1(10), 4.156968929685324e-6, 1e-20, "E1(10)")
})

test("curves that are not semistable or not smooth are refused", () => {
  assert.throws(() => curve(0, 0, 0, 0, 0), /singular/)
  assert.throws(() => curve(0, 0, 0, 0, 1), /additive/)
  assert.throws(() => twist(curves[0][1], -11), /gcd/)
  assert.throws(() => twist(curves[0][1], 9), /fundamental/)
})

test("equations print like the LMFDB", () => {
  assert.equal(equation([0, -1, 1, -10, -20]), "y² + y = x³ − x² − 10x − 20")
  assert.equal(equation([1, 0, 0, -1, 0]), "y² + xy = x³ − x")
  assert.equal(equation([0, 0, 1, -7, 6]), "y² + y = x³ − 7x + 6")
})

test("plot ranges shrink with the conductor", () => {
  const [, E11] = curves[0]
  const [, E5077] = curves[3]
  assert.equal(zRange(twist(E11)), 40)
  assert.ok(zRange(twist(E11, -1003)) === 10)
  assert.equal(twistRange(E11), 1000)
  assert.equal(twistRange(E5077), 73)
})
