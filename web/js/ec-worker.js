import { gcd } from "./characters.js"
import { centralOrder, curve, hardyZ, isFundamental, meanSpacing, twist } from "./elliptic.js"
import { illinois } from "./zeta.js"

const BATCH = 25

function zJob({ id, a, d, T }) {
  const L = twist(curve(...a), d)
  const n = Math.max(400, Math.ceil((8 * T) / meanSpacing(T, L)))
  const ts = Float64Array.from({ length: n + 1 }, (_, i) => (T * i) / n)
  const zs = ts.map((t) => hardyZ(t, L))
  const f = (t) => hardyZ(t, L)
  const zeros = []
  for (let i = 1; i < n; i++) {
    if (zs[i] * zs[i + 1] < 0) zeros.push(illinois(f, ts[i], ts[i + 1], zs[i], zs[i + 1]))
  }
  const { value, order } = centralOrder(L)
  self.postMessage({ kind: "z", id, ts, zs, zeros, value, order }, [ts.buffer, zs.buffer])
}

function twistJob({ id, a, D }) {
  const E = curve(...a)
  let batch = []
  let odd = 0
  const flush = (done) => {
    self.postMessage({ kind: "twists", id, points: batch, odd, done })
    batch = []
  }
  for (let m = 3; m <= D; m++) {
    for (const d of [m, -m]) {
      if (!isFundamental(d) || gcd(d, E.N) !== 1) continue
      const L = twist(E, d)
      if (L.w === -1) odd++
      else batch.push([d, centralOrder(L).value])
      if (batch.length === BATCH) flush(false)
    }
  }
  flush(true)
}

self.onmessage = ({ data }) => (data.kind === "z" ? zJob : twistJob)(data)
