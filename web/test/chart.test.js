import assert from "node:assert/strict"
import { test } from "node:test"

import { decimals, niceTicks, tickFormat } from "../js/chart.js"

test("ticks land on round steps inside the range", () => {
  assert.deepEqual(niceTicks(0, 10, 5), [0, 2, 4, 6, 8, 10])
  assert.deepEqual(niceTicks(-0.3, 0.8, 5), [-0.25, 0, 0.25, 0.5, 0.75])
  const t = niceTicks(1e9 + 0.1, 1e9 + 3.9, 4)
  assert.deepEqual(t.map((v) => v - 1e9), [1, 2, 3])
})

test("labels carry just enough digits", () => {
  assert.equal(decimals(2), 0)
  assert.equal(decimals(0.25), 2)
  assert.equal(decimals(2.5), 1)
  assert.deepEqual(niceTicks(-0.3, 0.8, 5).map(tickFormat(niceTicks(-0.3, 0.8, 5))),
    ["-0.25", "0.00", "0.25", "0.50", "0.75"])
})
