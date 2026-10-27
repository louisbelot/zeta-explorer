import assert from "node:assert/strict"
import { readFileSync } from "node:fs"
import { test } from "node:test"

import { addZeros, chebyshevPsi, psiFromZeros, sieve, smoothPart } from "../js/explicit.js"

const load = (path) => JSON.parse(readFileSync(new URL(path, import.meta.url)))
const ref = load("./reference.json")
const gamma = load("../data/zeros.json").gamma

test("sieve", () => {
  assert.deepEqual(sieve(30), [2, 3, 5, 7, 11, 13, 17, 19, 23, 29])
  assert.equal(sieve(10000).length, 1229)
})

test("chebyshev psi at prime powers", () => {
  const [a, b, c] = chebyshevPsi([8, 9.5, 1])
  assert.ok(Math.abs(a - (3 * Math.log(2) + Math.log(3) + Math.log(5) + Math.log(7))) < 1e-12)
  assert.ok(Math.abs(b - (a + Math.log(3))) < 1e-12)
  assert.equal(c, 0)
})

test("explicit formula matches the Python version", () => {
  const xs = ref.psi100.map(([x]) => x)
  const ours = psiFromZeros(xs, gamma, 100)
  ref.psi100.forEach(([x, want], i) => assert.ok(Math.abs(ours[i] - want) < 1e-9, `x = ${x}`))
})

test("adding zeros in steps is the same as all at once", () => {
  const xs = [3.3, 17.5, 64.2]
  const acc = Float64Array.from(xs, smoothPart)
  for (const [a, b] of [[0, 7], [7, 50], [50, 300]]) addZeros(acc, xs, gamma, a, b)
  const once = psiFromZeros(xs, gamma, 300)
  xs.forEach((_, i) => assert.ok(Math.abs(acc[i] - once[i]) < 1e-11))
})

test("10000 zeros pin psi between prime powers", () => {
  const xs = [10.5, 30.5, 50.5, 70.5, 99.5]
  const exact = chebyshevPsi(xs)
  const approx = psiFromZeros(xs, gamma)
  xs.forEach((x, i) => assert.ok(Math.abs(approx[i] - exact[i]) < 0.02, `x = ${x}`))
})
