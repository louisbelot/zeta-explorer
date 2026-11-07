import assert from "node:assert/strict"
import { readFileSync } from "node:fs"
import { test } from "node:test"

import {
  conjugatePairs,
  conrey,
  conreyGenerator,
  factor,
  primitiveCharacters,
} from "../js/characters.js"

const ref = JSON.parse(readFileSync(new URL("./reference.json", import.meta.url)))

const close = (got, want, tol, what) =>
  assert.ok(Math.abs(got - want) <= tol, `${what}: got ${got}, want ${want}`)

test("factor", () => {
  assert.deepEqual(factor(1), [])
  assert.deepEqual(factor(360), [[2, 3], [3, 2], [5, 1]])
  assert.deepEqual(factor(9973), [[9973, 1]])
})

test("conrey generators", () => {
  assert.equal(conreyGenerator(3), 2)
  assert.equal(conreyGenerator(7), 3)
  assert.equal(conreyGenerator(29), 2)
  assert.equal(conreyGenerator(40487), 10)
})

test("order, parity and conductor match lmfdb", () => {
  const rows = [
    [7, 6, 2, 1, 7], [7, 2, 3, 0, 7], [7, 3, 6, 1, 7], [8, 3, 2, 1, 8], [8, 5, 2, 0, 8],
    [5, 2, 4, 1, 5], [9, 4, 3, 0, 9], [12, 5, 2, 1, 3], [16, 9, 2, 0, 8], [32, 17, 2, 0, 8],
  ]
  for (const [q, n, order, parity, conductor] of rows) {
    const chi = conrey(q, n)
    const got = [chi.order, chi.parity, chi.conductor]
    assert.deepEqual(got, [order, parity, conductor], `${q}.${n}`)
  }
})

test("values match the python characters", () => {
  for (const [q, n, values] of ref.conrey) {
    const chi = conrey(q, n)
    values.forEach(([re, im], m) => {
      const [vr, vi] = chi.value(m)
      close(vr, re, 1e-14, `Re chi_${q}(${n}, ${m})`)
      close(vi, im, 1e-14, `Im chi_${q}(${n}, ${m})`)
    })
  }
})

test("characters are multiplicative", () => {
  for (const [q, n] of [[63, 5], [64, 3], [99, 2], [120, 7]]) {
    const chi = conrey(q, n)
    for (let a = 1; a < q; a++) {
      for (let b = 1; b < q; b++) {
        const [ar, ai] = chi.value(a)
        const [br, bi] = chi.value(b)
        const [cr, ci] = chi.value(a * b)
        close(ar * br - ai * bi, cr, 1e-12, `${q}.${n} at ${a}*${b}`)
        close(ar * bi + ai * br, ci, 1e-12, `${q}.${n} at ${a}*${b}`)
      }
    }
  }
})

test("number of primitive characters", () => {
  const counts = { 3: 1, 4: 1, 5: 3, 8: 2, 9: 4, 12: 1, 15: 3, 16: 4, 30: 0, 36: 4, 100: 16 }
  for (const [q, k] of Object.entries(counts)) {
    assert.equal(primitiveCharacters(Number(q)).length, k, `q = ${q}`)
  }
})

test("gauss sums have modulus sqrt q", () => {
  for (const q of [3, 5, 8, 11, 16, 21, 25, 27, 40, 99]) {
    for (const chi of primitiveCharacters(q)) {
      close(Math.hypot(...chi.gaussSum()), Math.sqrt(q), 1e-11, `tau(${q}.${chi.n})`)
      close(Math.hypot(...chi.rootNumber()), 1, 1e-12, `eps(${q}.${chi.n})`)
    }
  }
})

test("root number of 5.2 matches lmfdb", () => {
  const [re, im] = conrey(5, 2).rootNumber()
  const a = 2 * Math.PI * 0.0881040955873917
  close(re, Math.cos(a), 1e-13, "Re eps")
  close(im, Math.sin(a), 1e-13, "Im eps")
})

test("real primitive characters have root number 1", () => {
  for (const q of [3, 4, 5, 8, 13, 24, 40, 101]) {
    for (const chi of primitiveCharacters(q).filter((c) => c.isReal)) {
      const [re, im] = chi.rootNumber()
      close(re, 1, 1e-12, `${q}.${chi.n}`)
      close(im, 0, 1e-12, `${q}.${chi.n}`)
    }
  }
})

test("conjugate is the inverse label", () => {
  const chi = conrey(5, 2)
  const bar = chi.conj()
  assert.equal(bar.n, 3)
  assert.deepEqual([...bar.num], [...conrey(5, 3).num])
})

test("bad labels throw", () => {
  assert.throws(() => conrey(12, 4), RangeError)
  assert.throws(() => conrey(7, 8), RangeError)
  assert.throws(() => conrey(12, 5).rootNumber(), RangeError)
})

test("conjugate pairs cover every character once", () => {
  const chars = primitiveCharacters(13)
  const pairs = conjugatePairs(chars)
  assert.equal(pairs.length, 6)
  assert.deepEqual(pairs.filter(([, m]) => m === null), [[12, null]])
  assert.deepEqual(pairs.flat().filter((n) => n !== null).sort((a, b) => a - b),
    chars.map((c) => c.n))
  assert.deepEqual(pairs[0], [2, 7])
})
