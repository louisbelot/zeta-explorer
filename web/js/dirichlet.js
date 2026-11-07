import { Chart, decimals, niceTicks } from "./chart.js"
import { conjugatePairs, primitiveCharacters } from "./characters.js"
import { frameLoop, panZoom } from "./gestures.js"
import { emLength, lMeanSpacing, lZ, lZeros } from "./lfunction.js"

const Q_MAX = 200
const ROWS_T = 20
const T_MAX = 1000
const MAX_SPACINGS = 120
const MIN_SPAN = 0.05
const TERM_BUDGET = 3e6
const SUB = "₀₁₂₃₄₅₆₇₈₉"

const label = (chi) => `χ${[...String(chi.q)].map((c) => SUB[c]).join("")}(${chi.n}, ·)`

function describe(chi) {
  const kind = chi.order === 2 ? "real" : `order ${chi.order}`
  return `${kind}, ${chi.parity ? "odd" : "even"}`
}

function rootNumberText(chi) {
  const [re, im] = chi.rootNumber()
  if (Math.abs(im) < 1e-12) return re > 0 ? "ε = 1" : "ε = −1"
  const turn = (Math.atan2(im, re) / (2 * Math.PI) + 1) % 1
  return `ε = exp(2πi · ${turn.toFixed(6)})`
}

export function dirichletPlots(el) {
  const chart = new Chart(el.plot, { l: 44, r: 12, t: 14, b: 30 })
  const rowChart = new Chart(el.rows, { l: 44, r: 12, t: 10, b: 30 })
  const workers = Array.from({ length: Math.min(4, navigator.hardwareConcurrency || 2) }, () => {
    const w = new Worker(new URL("./lzeros-worker.js", import.meta.url), { type: "module" })
    w.onmessage = ({ data }) => {
      w.busy = false
      if (data.id === job) {
        rowZeros.set(data.n, data.zeros)
        const m = partner.get(data.n)
        if (m !== null) rowZeros.set(m, data.zeros.map((t) => -t).reverse())
        scheduleRows()
      }
      pump()
    }
    return w
  })

  let chars = []
  let chi = null
  let view = [-30, 30]
  let marks = []
  let job = 0
  let queue = []
  let rowZeros = new Map()
  let partner = new Map()

  function pump() {
    for (const w of workers) {
      if (w.busy || queue.length === 0) continue
      w.busy = true
      w.postMessage(queue.shift())
    }
  }

  function setModulus(q, n) {
    chars = primitiveCharacters(q)
    el.select.replaceChildren(
      ...chars.map((c) => new Option(`${c.n}: ${describe(c)}`, String(c.n))),
    )
    el.select.disabled = chars.length === 0
    const height = Math.min(480, Math.max(120, 40 + 18 * chars.length))
    el.rows.parentElement.style.height = `${height}px`

    const pairs = conjugatePairs(chars)
    partner = new Map(pairs.flatMap(([a, b]) => (b === null ? [[a, null]] : [[a, b], [b, a]])))
    rowZeros = new Map()
    job++
    queue = pairs.map(([a]) => ({ id: job, q, n: a, T: ROWS_T }))
    pump()
    pick(chars.find((c) => c.n === n) ?? chars[0] ?? null)
  }

  function pick(next) {
    chi = next
    if (chi) el.select.value = String(chi.n)
    view = [-30, 30]
    el.info.textContent = chi
      ? `${label(chi)}   ${describe(chi)}   ${rootNumberText(chi)}   ` +
        (chi.isReal ? "self-dual" : `dual ${label(chi.conj())}`)
      : `no primitive characters mod ${el.q.value}: a modulus ≡ 2 mod 4 never has one`
    schedule()
    scheduleRows()
  }

  function clamp([a, b]) {
    const spacing = lMeanSpacing(Math.max(Math.abs(a), Math.abs(b)), chi.q)
    const w = Math.max(MIN_SPAN, Math.min(b - a, MAX_SPACINGS * spacing))
    a = Math.min(Math.max(-T_MAX, a), T_MAX - w)
    return [a, a + w]
  }

  function draw() {
    chart.resize()
    chart.clear()
    if (!chi) return
    const [a, b] = view
    const cost = chi.q * (emLength(Math.max(Math.abs(a), Math.abs(b))) + 15)
    const n = Math.round(Math.min(chart.plotWidth, Math.max(200, TERM_BUDGET / cost)))
    const ts = Array.from({ length: n + 1 }, (_, i) => a + ((b - a) * i) / n)
    const zs = ts.map((t) => lZ(t, chi))
    const top = Math.max(1, 1.12 * Math.max(...zs.map(Math.abs)))
    chart.setView(a, b, -top, top)
    const xTicks = niceTicks(a, b, Math.max(3, Math.floor(chart.plotWidth / 90)))
    const d = xTicks.length > 1 ? decimals(xTicks[1] - xTicks[0]) : 0
    chart.axes({ xTicks, xFormat: (v) => v.toFixed(d), xLabel: "t", yLabel: `Z(t, ${label(chi)})` })
    const { colors } = chart
    chart.line([a, b], [0, 0], { color: colors.muted, width: 1 })
    chart.line(ts, zs, { color: colors.series, width: 1.6 })
    marks = lZeros(chi, a, b)
    chart.dots(marks, marks.map(() => 0), { color: colors.ink, r: 3 })
  }

  function rowGeometry() {
    const { pad, h } = rowChart
    const step = (h - pad.t - pad.b) / Math.max(1, chars.length)
    return { step, y: (i) => pad.t + step * (i + 0.5) }
  }

  function drawRows() {
    rowChart.resize()
    rowChart.clear()
    if (!chars.length) return
    rowChart.setView(-ROWS_T, ROWS_T, 0, 1)
    rowChart.axes({ yTicks: [], xLabel: "t" })
    const { ctx, colors } = rowChart
    const { step, y } = rowGeometry()
    const r = Math.max(1.2, Math.min(3, step / 3))
    ctx.save()
    ctx.font = `10px ${colors.mono}`
    ctx.textAlign = "right"
    ctx.textBaseline = "middle"
    chars.forEach((c, i) => {
      const selected = c === chi
      if (selected) {
        ctx.fillStyle = colors.soft
        ctx.globalAlpha = 0.35
        ctx.fillRect(rowChart.pad.l, y(i) - step / 2, rowChart.plotWidth, step)
        ctx.globalAlpha = 1
      }
      if (step >= 9) {
        ctx.fillStyle = selected ? colors.ink : colors.muted
        ctx.fillText(String(c.n), rowChart.pad.l - 6, y(i))
      }
      const zeros = rowZeros.get(c.n)
      if (!zeros) return
      ctx.fillStyle = selected ? colors.series : colors.muted
      for (const t of zeros) {
        ctx.beginPath()
        ctx.arc(rowChart.X(t), y(i), selected ? r + 0.6 : r, 0, 2 * Math.PI)
        ctx.fill()
      }
    })
    ctx.restore()
    const done = chars.filter((c) => rowZeros.has(c.n)).length
    if (done < chars.length) {
      el.rowsReadout.textContent = `computing zeros: ${done} of ${chars.length} characters`
    } else if (!el.rowsReadout.dataset.hover) {
      el.rowsReadout.textContent = " "
    }
  }

  function rowAt(e) {
    const rect = el.rows.getBoundingClientRect()
    const { step } = rowGeometry()
    const i = Math.floor((e.clientY - rect.top - rowChart.pad.t) / step)
    return chars[i] ?? null
  }

  const schedule = frameLoop(draw)
  const scheduleRows = frameLoop(drawRows)

  function show(next) {
    view = clamp(next)
    schedule()
  }

  panZoom(el.plot, {
    pan(dx) {
      const shift = (dx * (view[1] - view[0]) * chart.w) / chart.plotWidth
      show([view[0] - shift, view[1] - shift])
    },
    zoom(f, fx) {
      const t = chart.invX(fx * chart.w)
      show([t + (view[0] - t) * f, t + (view[1] - t) * f])
    },
    hover(fx) {
      if (fx === null || !chi) {
        el.readout.textContent = " "
        return
      }
      const t = chart.invX(fx * chart.w)
      let text = `t = ${t.toFixed(4)}   Z = ${lZ(t, chi).toFixed(5)}`
      const dist = (z) => Math.abs(z - t)
      const near = marks.reduce((m, z) => (m === null || dist(z) < dist(m) ? z : m), null)
      if (near !== null && Math.abs(chart.X(near) - fx * chart.w) < 12) {
        text += `   zero at ${near.toFixed(8)}`
      }
      el.readout.textContent = text
    },
  })

  el.rows.addEventListener("click", (e) => {
    const c = rowAt(e)
    if (c && c !== chi) pick(c)
  })
  el.rows.addEventListener("pointermove", (e) => {
    const c = rowAt(e)
    const zeros = c && rowZeros.get(c.n)
    if (!zeros) return
    const low = zeros.reduce((m, t) => (Math.abs(t) < Math.abs(m) ? t : m), Infinity)
    el.rowsReadout.dataset.hover = "1"
    el.rowsReadout.textContent =
      `${label(c)}   ${describe(c)}   ${zeros.length} zeros with |t| < ${ROWS_T}` +
      (zeros.length ? `, lowest at ${low.toFixed(6)}` : "")
  })
  el.rows.addEventListener("pointerleave", () => {
    delete el.rowsReadout.dataset.hover
    scheduleRows()
  })

  const readQ = () => Math.min(Q_MAX, Math.max(3, Math.round(Number(el.q.value) || 3)))
  el.q.addEventListener("change", () => {
    el.q.value = readQ()
    setModulus(readQ())
  })
  el.select.addEventListener("change", () => {
    pick(chars.find((c) => c.n === Number(el.select.value)))
  })
  for (const button of el.section.querySelectorAll("button[data-q]")) {
    button.addEventListener("click", () => {
      el.q.value = button.dataset.q
      setModulus(Number(button.dataset.q), Number(button.dataset.n))
    })
  }

  setModulus(readQ(), Number(el.select.dataset.n))
  return () => {
    schedule()
    scheduleRows()
  }
}
