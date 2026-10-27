import { Chart, decimals, niceTicks } from "./chart.js"
import { frameLoop, panZoom } from "./gestures.js"
import { findZeros, gramPoints, meanSpacing } from "./zeta.js"

const T_MAX = 1e9 + 1e3
const MAX_SPACINGS = 160
const MIN_SPAN = 0.05
const TERM_BUDGET = 4e6

export function zExplorer(canvas, Z, gamma, readout, controls) {
  const chart = new Chart(canvas, { l: 44, r: 12, t: 14, b: 30 })
  const tableEnd = gamma.at(-1)
  let view = [0, 60]
  let marks = []

  function clamp([a, b]) {
    const w = Math.max(MIN_SPAN, Math.min(b - a, MAX_SPACINGS * meanSpacing(b)))
    a = Math.min(Math.max(0, a), T_MAX - w)
    return [a, a + w]
  }

  function zerosIn(a, b) {
    const lo = lowerBound(gamma, a)
    const hi = lowerBound(gamma, Math.min(b, tableEnd + 1e-9))
    const out = []
    for (let i = lo; i < hi; i++) out.push({ t: gamma[i], n: i + 1 })
    const from = Math.max(a, tableEnd + 1e-6)
    if (b > from) for (const t of findZeros(Z, from, b)) out.push({ t })
    return out
  }

  function draw() {
    chart.resize()
    const [a, b] = view
    const terms = Math.sqrt(b / (2 * Math.PI))
    const n = Math.round(Math.min(chart.plotWidth, Math.max(300, TERM_BUDGET / terms)))
    const ts = Array.from({ length: n + 1 }, (_, i) => a + ((b - a) * i) / n)
    const zs = ts.map(Z)
    const top = Math.max(1, 1.12 * Math.max(...zs.map(Math.abs)))
    chart.setView(a, b, -top, top)

    const base = a >= 1e5 ? Math.floor(a) : 0
    const xTicks = niceTicks(a, b, Math.max(3, Math.floor(chart.plotWidth / 90)))
    const d = xTicks.length > 1 ? decimals(xTicks[1] - xTicks[0]) : 0
    chart.clear()
    chart.axes({
      xTicks,
      xFormat: base ? (v) => `+${(v - base).toFixed(d)}` : (v) => v.toFixed(d),
      xLabel: base ? `t − ${base.toLocaleString("en-US")}` : "t",
      yLabel: "Z(t)",
    })
    const { colors, ctx } = chart
    chart.line([a, b], [0, 0], { color: colors.muted, width: 1 })
    chart.line(ts, zs, { color: colors.series, width: 1.6 })

    const y0 = chart.Y(0)
    chart.clip(() => {
      for (const [k, g] of gramPoints(a, b)) {
        const good = (k % 2 === 0 ? 1 : -1) * Z(g) > 0
        ctx.strokeStyle = good ? colors.muted : colors.accent
        ctx.lineWidth = good ? 1 : 2
        ctx.beginPath()
        ctx.moveTo(chart.X(g), y0 - 7)
        ctx.lineTo(chart.X(g), y0 + 7)
        ctx.stroke()
      }
    })
    marks = zerosIn(a, b)
    chart.dots(marks.map((m) => m.t), marks.map(() => 0), { color: colors.ink, r: 3 })
  }

  const schedule = frameLoop(draw)

  function show(next) {
    view = clamp(next)
    schedule()
  }

  panZoom(canvas, {
    pan(dx) {
      const shift = (dx * (view[1] - view[0]) * chart.w) / chart.plotWidth
      show([view[0] - shift, view[1] - shift])
    },
    zoom(f, fx) {
      const t = chart.invX(fx * chart.w)
      show([t + (view[0] - t) * f, t + (view[1] - t) * f])
    },
    hover(fx) {
      if (fx === null) {
        readout.textContent = " "
        return
      }
      const t = chart.invX(fx * chart.w)
      const dist = (m) => Math.abs(m.t - t)
      const near = marks.reduce((m, z) => (!m || dist(z) < dist(m) ? z : m), null)
      let text = `t = ${t.toFixed(4)}   Z = ${Z(t).toFixed(5)}`
      if (near && Math.abs(chart.X(near.t) - fx * chart.w) < 12) {
        const label = near.n ? `zero #${near.n.toLocaleString("en-US")}` : "zero"
        text += `   ${label} at ${near.t.toFixed(6)}`
      }
      readout.textContent = text
    },
  })

  for (const button of controls.querySelectorAll("button[data-t]")) {
    button.addEventListener("click", () => {
      const t = Number(button.dataset.t)
      show([t, t + Number(button.dataset.w)])
    })
  }
  const go = controls.querySelector("#z-go")
  go.addEventListener("change", () => {
    const t = Number(go.value)
    if (!Number.isFinite(t)) return
    const w = 12 * meanSpacing(t)
    show([t - w / 2, t + w / 2])
  })

  return schedule
}

function lowerBound(xs, v) {
  let lo = 0
  let hi = xs.length
  while (lo < hi) {
    const mid = (lo + hi) >> 1
    if (xs[mid] < v) lo = mid + 1
    else hi = mid
  }
  return lo
}
