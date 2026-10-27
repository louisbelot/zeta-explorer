import assert from "node:assert/strict"
import { readFileSync } from "node:fs"
import { test } from "node:test"

import {
  bernoulli,
  eulerMaclaurin,
  findZeros,
  gramPoint,
  gramPoints,
  loggamma,
  theta,
  thetaExact,
  zeta,
  zExact,
  zFunction,
  zRs,
} from "../js/zeta.js"

const load = (path) => JSON.parse(readFileSync(new URL(path, import.meta.url)))
const ref = load("./reference.json")
const rs = load("../data/rs.json")
const zeros = load("../data/zeros.json").gamma
const Z = zFunction(rs)

function close(got, want, tol, what) {
  assert.ok(Math.abs(got - want) <= tol, `${what}: got ${got}, want ${want}`)
}

test("bernoulli numbers", () => {
  const B = bernoulli(12)
  assert.deepEqual(B[1], [-1n, 2n])
  assert.deepEqual(B[3], [0n, 1n])
  assert.deepEqual(B[12], [-691n, 2730n])
})

test("loggamma matches known values", () => {
  close(loggamma(1, 0)[0], 0, 1e-15, "log G(1)")
  close(loggamma(0.5, 0)[0], 0.5 * Math.log(Math.PI), 1e-14, "log G(1/2)")
  close(loggamma(10, 0)[0], Math.log(362880), 1e-13, "log G(10)")
  const [re, im] = loggamma(1, 1)
  close(re, -0.6509231993018563, 1e-14, "Re log G(1+i)")
  close(im, -0.3016403204675331, 1e-14, "Im log G(1+i)")
})

test("zeta matches mpmath over the plane", () => {
  for (const [s, t, re, im] of ref.zeta) {
    const [zr, zi] = zeta(s, t)
    const tol = 1e-12 * Math.max(1, Math.hypot(re, im))
    close(zr, re, tol, `Re zeta(${s} + ${t}i)`)
    close(zi, im, tol, `Im zeta(${s} + ${t}i)`)
  }
})

test("zeta is real on the real axis and has a pole at 1", () => {
  assert.equal(zeta(-3)[1], 0)
  assert.equal(zeta(-6)[0], 0)
  assert.equal(zeta(1)[0], Infinity)
  close(zeta(1 + 2 ** -27)[0], 2 ** 27 + 0.5772156649, 1e-6, "zeta near the pole")
})

test("Euler-Maclaurin bound covers the error", () => {
  const exact = ref.zeta.find(([s, t]) => s === 0.5 && t === 100)
  for (const M of [2, 4, 8]) {
    const { re, im, bound } = eulerMaclaurin(0.5, 100, 40, M, 0)
    assert.ok(Math.hypot(re - exact[2], im - exact[3]) <= bound + 1e-14, `M = ${M}`)
  }
})

test("theta matches mpmath", () => {
  for (const [t, want] of ref.theta) {
    const got = t < 200 ? thetaExact(t) : theta(t)
    close(got, want, 2e-15 * Math.max(1, Math.abs(want)) + 1e-14, `theta(${t})`)
  }
  close(theta(50), thetaExact(50), 1e-12, "series against loggamma at t = 50")
})

test("Z matches mpmath siegelz", () => {
  for (const [t, want] of ref.Z) {
    const tol = t < 200 ? 1e-12 : 1e-8
    close(Z(t), want, tol, `Z(${t})`)
  }
})

test("Riemann-Siegel and Euler-Maclaurin agree at the switch", () => {
  for (const t of [200, 260.5, 420.25]) close(zRs(t, rs), zExact(t), 1e-8, `Z(${t})`)
})

test("Z is even", () => {
  for (const t of [3, 40, 333.3]) assert.equal(Z(-t), Z(t))
})

test("Gram points match mpmath", () => {
  for (const [n, want] of ref.gram) close(gramPoint(n), want, 1e-12 * want, `g_${n}`)
  const g = gramPoints(20, 57)
  assert.deepEqual(g.map(([n]) => n), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
  for (const [n, t] of g) close(theta(t), n * Math.PI, 1e-11, `theta(g_${n})`)
})

test("zero finder recovers the first zeros", () => {
  const found = findZeros(Z, 10, 101.3)
  assert.equal(found.length, 29)
  found.forEach((t, i) => close(t, ref.zeros["1-29"][i], 1e-10, `zero ${i + 1}`))
})

test("zero finder works above the switch to Riemann-Siegel", () => {
  const want = ref.zeros["1990-2000"]
  const found = findZeros(Z, want[0] - 0.3, want.at(-1) + 0.3)
  assert.equal(found.length, want.length)
  found.forEach((t, i) => close(t, want[i], 1e-9, `zero ${1990 + i}`))
})

test("precomputed zeros agree with the live finder", () => {
  const found = findZeros(Z, 9800, zeros.at(-1) + 0.01)
  const table = zeros.filter((t) => t >= 9800)
  assert.equal(found.length, table.length)
  found.forEach((t, i) => close(t, table[i], 1e-9, `table zero near ${t}`))
})
