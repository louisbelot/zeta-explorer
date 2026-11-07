import assert from "node:assert/strict"
import { readFileSync } from "node:fs"
import { test } from "node:test"

import { conrey, primitiveCharacters } from "../js/characters.js"
import { lTheta, lValue, lZ, lZeros, smoothCount } from "../js/lfunction.js"

const ref = JSON.parse(readFileSync(new URL("./reference.json", import.meta.url)))

const close = (got, want, tol, what) =>
  assert.ok(Math.abs(got - want) <= tol, `${what}: got ${got}, want ${want}`)

test("L matches mpmath", () => {
  for (const [q, n, s, t, re, im] of ref.dirichlet) {
    const [lr, li] = lValue(s, t, conrey(q, n))
    const tol = 1e-12 * Math.max(1, Math.hypot(re, im))
    close(lr, re, tol, `Re L(${s} + ${t}i, ${q}.${n})`)
    close(li, im, tol, `Im L(${s} + ${t}i, ${q}.${n})`)
  }
})

test("L at s = 1 for chi_4 is pi / 4", () => {
  const [re, im] = lValue(1, 0, conrey(4, 3))
  close(re, Math.PI / 4, 1e-14, "Re L(1)")
  close(im, 0, 1e-15, "Im L(1)")
})

test("Z matches the python version", () => {
  for (const [q, n, ts, zs] of ref.lz) {
    const chi = conrey(q, n)
    ts.forEach((t, i) => close(lZ(t, chi), zs[i], 1e-11, `Z(${t}, ${q}.${n})`))
  }
})

test("e^(i theta) L is real on the critical line", () => {
  for (const q of [5, 7, 16, 21, 37]) {
    for (const chi of primitiveCharacters(q)) {
      for (const t of [-41.3, -3.2, 0.7, 12.6, 88.8]) {
        const th = lTheta(t, chi)
        const [re, im] = lValue(0.5, t, chi)
        close(Math.sin(th) * re + Math.cos(th) * im, 0, 1e-11, `Im Z(${t}, ${q}.${chi.n})`)
      }
    }
  }
})

test("Z(-t, chi) = Z(t, conj chi)", () => {
  for (const [q, n] of [[5, 2], [7, 3], [13, 2], [36, 11]]) {
    const chi = conrey(q, n)
    const bar = chi.conj()
    for (const t of [2.5, 17.1, 60.4]) close(lZ(-t, chi), lZ(t, bar), 1e-11, `${q}.${n} at ${t}`)
  }
})

test("first zeros match lmfdb", () => {
  const lmfdb = {
    "3.2": [8.03973715568146668, 11.24920620777293525, 15.70461917672162557],
    "4.3": [6.02094890469759665, 10.24377030416655455, 12.98809801231242251],
    "5.2": [6.18357819545085391, 8.45722917442323072, 12.67494641701135578],
  }
  for (const [label, want] of Object.entries(lmfdb)) {
    const [q, n] = label.split(".").map(Number)
    const got = lZeros(conrey(q, n), 0.5, 16)
    want.forEach((g, i) => close(got[i], g, 1e-10, `zero ${i + 1} of ${label}`))
  }
})

test("zero counts match the argument principle", () => {
  for (const [q, n, T, count, zeros] of ref.lzeros) {
    const got = lZeros(conrey(q, n), -T, T)
    assert.equal(got.length, count, `${q}.${n} up to ${T}`)
    got.forEach((g, i) => close(g, zeros[i], 1e-9, `${q}.${n} zero ${i}`))
  }
})

test("zero counts follow theta", () => {
  for (const [q, n, T] of [[7, 3, 60], [24, 5, 40], [61, 2, 30], [97, 5, 25]]) {
    const chi = conrey(q, n)
    const got = lZeros(chi, -T, T).length
    close(got, smoothCount(chi, T), 2, `${q}.${n} up to ${T}`)
  }
})

test("principal and imprimitive characters are rejected", () => {
  assert.throws(() => lValue(0.5, 1, conrey(5, 1)), RangeError)
  assert.throws(() => lValue(0.5, 1, conrey(12, 5)), RangeError)
})
