import assert from "node:assert/strict"
import { test } from "node:test"

import { phaseRgb } from "../js/color.js"
import { renderBand } from "../js/domain-worker.js"

const rgb = (re, im) => {
  const out = new Uint8ClampedArray(4)
  phaseRgb(re, im, out, 0)
  return [...out]
}

test("hue follows the argument", () => {
  const [r, g, b] = rgb(1, 0)
  assert.ok(r > g && r > b && g === b)
  const [r2, g2, b2] = rgb(-1, 0)
  assert.ok(g2 > r2 && b2 > r2)
})

test("shading steps at each doubling of the modulus", () => {
  assert.deepEqual(rgb(1, 0), rgb(2, 0))
  assert.deepEqual(rgb(1.5, 0), rgb(3, 0))
  assert.ok(rgb(1.9, 0)[0] > rgb(1.1, 0)[0])
})

test("infinite modulus is drawn at full value", () => {
  assert.equal(rgb(Infinity, 0)[0], 255)
})

test("band renderer fills every pixel and blocks", () => {
  const view = { tMin: 0, tMax: 30, sMin: -2, sMax: 3 }
  const fine = renderBand({ width: 12, height: 8, row0: 0, rows: 8, view, step: 1 })
  const coarse = renderBand({ width: 12, height: 8, row0: 0, rows: 8, view, step: 4 })
  for (const px of [fine, coarse]) for (let i = 3; i < px.length; i += 4) assert.equal(px[i], 255)
  assert.deepEqual(coarse.slice(0, 4), coarse.slice(12, 16))
})
