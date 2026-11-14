import { Chart, niceTicks } from "./chart.js"
import { curve, equation, twist, twistRange, zRange } from "./elliptic.js"

const SUB = "₀₁₂₃₄₅₆₇₈₉"
const sub = (d) => [...String(d)].map((c) => (c === "-" ? "₋" : SUB[c])).join("")

function spawn(onmessage) {
  const w = new Worker(new URL("./ec-worker.js", import.meta.url), { type: "module" })
  w.onmessage = ({ data }) => onmessage(data)
  return w
}

function centralText(order, value, d) {
  const L = d === 1 ? "L(E, 1)" : `L(E ⊗ χ${sub(d)}, 1)`
  if (order === 0) return `${L} = ${value.toFixed(10)}`
  if (order === 1) return `${L.replace("L", "L′")} = ${value.toFixed(10)}`
  return `${L} vanishes to order ${order}`
}

export function curvePlots(el) {
  const zChart = new Chart(el.plot, { l: 44, r: 12, t: 14, b: 30 })
  const tChart = new Chart(el.twists, { l: 52, r: 12, t: 14, b: 30 })
  let E = null
  let d = 1
  let zJob = 0
  let tJob = 0
  let z = null
  let points = []
  let odd = 0
  let D = 0
  let twistsDone = false

  const zWorker = spawn((data) => {
    if (data.id !== zJob) return
    z = data
    el.info.textContent = header() + `   central order ${z.order}   ${centralText(z.order, z.value, d)}`
    drawZ()
  })
  const tWorker = spawn((data) => {
    if (data.id !== tJob) return
    points.push(...data.points)
    odd = data.odd
    twistsDone = data.done
    drawTwists()
  })

  const header = () => {
    const L = twist(E, d)
    const sign = L.w === 1 ? "+1" : "−1"
    const base = d === 1 ? equation(E.a) : `${equation(E.a)}, twisted by χ${sub(d)}`
    return `${base}   N = ${L.N}   ε = ${sign}`
  }

  function setCurve(a) {
    try {
      E = curve(...a)
    } catch (err) {
      el.info.textContent = err.message
      return
    }
    el.a.value = a.join(", ")
    d = 1
    el.d.value = "1"
    tJob++
    points = []
    odd = 0
    twistsDone = false
    D = twistRange(E)
    tWorker.postMessage({ kind: "twists", id: tJob, a, D })
    setTwist(1)
  }

  function setTwist(next) {
    if (!E) return
    let L
    try {
      L = twist(E, next)
    } catch (err) {
      el.readout.textContent = err.message
      return
    }
    d = next
    el.d.value = String(d)
    z = null
    zJob++
    el.info.textContent = header()
    el.readout.textContent = "computing Z(t)…"
    zWorker.postMessage({ kind: "z", id: zJob, a: E.a, d, T: zRange(L) })
    drawZ()
    drawTwists()
  }

  function drawZ() {
    zChart.resize()
    zChart.clear()
    if (!z) return
    const { ts, zs, zeros } = z
    const T = ts.at(-1)
    let top = 0
    for (const v of zs) top = Math.max(top, Math.abs(v))
    top = 1.12 * Math.max(top, 1e-12)
    zChart.setView(0, T, -top, top)
    const xTicks = niceTicks(0, T, Math.max(3, Math.floor(zChart.plotWidth / 90)))
    const yTicks = niceTicks(-top, top, 5)
    zChart.axes({ xTicks, yTicks, xLabel: "t", yLabel: "Z(t)" })
    const { colors } = zChart
    zChart.line([0, T], [0, 0], { color: colors.muted, width: 1 })
    zChart.line(ts, zs, { color: colors.series, width: 1.6 })
    zChart.dots(zeros, zeros.map(() => 0), { color: colors.ink, r: 3 })
    if (z.order > 0) zChart.dots([0], [0], { color: colors.accent, r: 3 + 1.5 * z.order })
    el.readout.textContent =
      `${zeros.length} zeros with 0 < t < ${T.toFixed(1)}` +
      (zeros.length ? `, the first at ${zeros[0].toFixed(8)}` : "")
  }

  function drawTwists() {
    tChart.resize()
    tChart.clear()
    if (!E) return
    let top = 0
    for (const [, v] of points) top = Math.max(top, v)
    top = 1.1 * Math.max(top, 1e-9)
    tChart.setView(-D, D, -0.04 * top, top)
    const xTicks = niceTicks(-D, D, Math.max(3, Math.floor(tChart.plotWidth / 90)))
    const yTicks = niceTicks(0, top, 4)
    tChart.axes({ xTicks, yTicks, xLabel: "d", yLabel: "L(1/2, E ⊗ χ_d)" })
    const { colors } = tChart
    const zero = points.filter(([, v]) => Math.abs(v) < 1e-6)
    const rest = points.filter(([, v]) => Math.abs(v) >= 1e-6)
    tChart.dots(rest.map(([x]) => x), rest.map(([, v]) => v), { color: colors.series, r: 2.2 })
    tChart.dots(zero.map(([x]) => x), zero.map(() => 0), { color: colors.accent, r: 3.2 })
    const here = points.find(([x]) => x === d)
    if (here) tChart.dots([here[0]], [here[1]], { color: colors.ink, r: 4 })
    if (!el.twistsReadout.dataset.hover) {
      el.twistsReadout.textContent =
        `${points.length} even twists with |d| ≤ ${D}${twistsDone ? "" : " so far"}, ` +
        `${zero.length} with a vanishing central value, and ${odd} odd ones (ε = −1)`
    }
  }

  function nearest(e) {
    const rect = el.twists.getBoundingClientRect()
    const x = tChart.invX(e.clientX - rect.left)
    let best = null
    for (const p of points) if (!best || Math.abs(p[0] - x) < Math.abs(best[0] - x)) best = p
    return best && Math.abs(tChart.X(best[0]) - (e.clientX - rect.left)) < 10 ? best : null
  }

  el.twists.addEventListener("pointermove", (e) => {
    const p = nearest(e)
    if (!p) return
    el.twistsReadout.dataset.hover = "1"
    el.twistsReadout.textContent = `d = ${p[0]}   L(1/2, E ⊗ χ${sub(p[0])}) = ${p[1].toFixed(8)}`
  })
  el.twists.addEventListener("pointerleave", () => {
    delete el.twistsReadout.dataset.hover
    drawTwists()
  })
  el.twists.addEventListener("click", (e) => {
    const p = nearest(e)
    if (p) setTwist(p[0])
  })
  el.a.addEventListener("change", () => {
    const a = el.a.value.split(/[\s,[\]]+/).filter(Boolean).map(Number)
    if (a.length === 5 && a.every(Number.isInteger)) setCurve(a)
    else el.info.textContent = "enter five integers a1, a2, a3, a4, a6"
  })
  el.d.addEventListener("change", () => setTwist(Math.round(Number(el.d.value)) || 1))
  for (const button of el.section.querySelectorAll("button[data-a]")) {
    button.addEventListener("click", () => setCurve(button.dataset.a.split(",").map(Number)))
  }

  setCurve(el.a.value.split(",").map(Number))
  return () => {
    drawZ()
    drawTwists()
  }
}
