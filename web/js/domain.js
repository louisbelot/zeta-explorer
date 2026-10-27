import { niceTicks, palette, tickFormat } from "./chart.js"
import { frameLoop, panZoom } from "./gestures.js"
import { zeta } from "./zeta.js"

const HOME = { tMin: -4, tMax: 56, sMin: -6, sMax: 7 }
const T_LIMIT = 1500
const S_LIMIT = 40
const PASSES = [8, 1]

export function domainPlot(img, overlay, readout, reset) {
  const ctx = img.getContext("2d")
  const snap = document.createElement("canvas")
  const workers = Array.from({ length: Math.min(8, navigator.hardwareConcurrency || 4) }, () => {
    const w = new Worker(new URL("./domain-worker.js", import.meta.url), { type: "module" })
    w.onmessage = ({ data }) => {
      w.busy = false
      if (data.id === job.id) {
        ctx.putImageData(new ImageData(data.px, job.width, data.rows), 0, data.row0)
        if (--job.left === 0) keep()
      }
      pump()
    }
    return w
  })

  let view = { ...HOME }
  let job = { id: 0 }
  let queue = []
  let snapView = null

  function pump() {
    for (const w of workers) {
      if (w.busy || queue.length === 0) continue
      w.busy = true
      w.postMessage(queue.shift())
    }
  }

  function keep() {
    snap.width = img.width
    snap.height = img.height
    snap.getContext("2d").drawImage(img, 0, 0)
    snapView = { ...job.view }
  }

  function preview() {
    if (!snapView) return
    const { width, height } = img
    const sx = (width / (view.tMax - view.tMin)) * (snapView.tMax - snapView.tMin)
    const sy = (height / (view.sMax - view.sMin)) * (snapView.sMax - snapView.sMin)
    const x = ((snapView.tMin - view.tMin) / (view.tMax - view.tMin)) * width
    const y = ((view.sMax - snapView.sMax) / (view.sMax - view.sMin)) * height
    ctx.clearRect(0, 0, width, height)
    ctx.drawImage(snap, x, y, sx, sy)
  }

  function render() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2)
    const { width: cw, height: ch } = img.getBoundingClientRect()
    const width = Math.round(cw * dpr)
    const height = Math.round(ch * dpr)
    if (width === 0 || height === 0) return
    if (img.width !== width || img.height !== height) {
      img.width = width
      img.height = height
      snapView = null
    }
    preview()
    const id = job.id + 1
    const rows = Math.ceil(height / (workers.length * 2) / 8) * 8
    queue = []
    for (const step of PASSES) {
      for (let row0 = 0; row0 < height; row0 += rows) {
        queue.push({ id, width, height, row0, rows: Math.min(rows, height - row0), view, step })
      }
    }
    job = { id, width, view: { ...view }, left: queue.length }
    pump()
    drawOverlay()
  }

  function drawOverlay() {
    const dpr = window.devicePixelRatio || 1
    const { width: w, height: h } = overlay.getBoundingClientRect()
    overlay.width = Math.round(w * dpr)
    overlay.height = Math.round(h * dpr)
    const o = overlay.getContext("2d")
    o.setTransform(dpr, 0, 0, dpr, 0, 0)
    const colors = palette()
    const X = (t) => ((t - view.tMin) / (view.tMax - view.tMin)) * w
    const Y = (s) => ((view.sMax - s) / (view.sMax - view.sMin)) * h
    o.strokeStyle = "rgba(255, 255, 255, 0.75)"
    for (const [s, dash] of [[0, []], [1, []], [0.5, [5, 5]]]) {
      o.setLineDash(dash)
      o.lineWidth = s === 0.5 ? 1.4 : 0.8
      o.beginPath()
      o.moveTo(0, Y(s))
      o.lineTo(w, Y(s))
      o.stroke()
    }
    o.setLineDash([])
    o.font = `11px ${colors.mono}`
    const label = (text, x, y, align, base) => {
      o.textAlign = align
      o.textBaseline = base
      o.lineWidth = 3
      o.strokeStyle = "rgba(0, 0, 0, 0.55)"
      o.strokeText(text, x, y)
      o.fillStyle = "#fff"
      o.fillText(text, x, y)
    }
    const tt = niceTicks(view.tMin, view.tMax, Math.max(3, Math.floor(w / 90)))
    const tf = tickFormat(tt)
    for (const t of tt) {
      if (X(t) > 72 && X(t) < w - 24) label(tf(t), X(t), h - 6, "center", "bottom")
    }
    const st = niceTicks(view.sMin, view.sMax, 5)
    const sf = tickFormat(st)
    for (const s of st) {
      if (Y(s) > 24 && Y(s) < h - 24) label(`σ = ${sf(s)}`, 6, Y(s), "left", "middle")
    }
    label("t = Im s", w - 6, 6, "right", "top")
  }

  const schedule = frameLoop(render)

  function clampView() {
    const tw = Math.min(view.tMax - view.tMin, T_LIMIT)
    const sw = Math.min(view.sMax - view.sMin, S_LIMIT)
    view.tMin = Math.min(Math.max(view.tMin, -T_LIMIT), T_LIMIT - tw)
    view.tMax = view.tMin + tw
    view.sMin = Math.min(Math.max(view.sMin, -S_LIMIT), S_LIMIT - sw)
    view.sMax = view.sMin + sw
  }

  panZoom(overlay, {
    pan(dx, dy) {
      const tw = view.tMax - view.tMin
      const sw = view.sMax - view.sMin
      view.tMin -= dx * tw
      view.tMax -= dx * tw
      view.sMin += dy * sw
      view.sMax += dy * sw
      clampView()
      schedule()
    },
    zoom(f, fx, fy) {
      const t = view.tMin + fx * (view.tMax - view.tMin)
      const s = view.sMax - fy * (view.sMax - view.sMin)
      if (f < 1 && view.tMax - view.tMin < 1e-3) return
      view.tMin = t + (view.tMin - t) * f
      view.tMax = t + (view.tMax - t) * f
      view.sMin = s + (view.sMin - s) * f
      view.sMax = s + (view.sMax - s) * f
      clampView()
      schedule()
    },
    hover(fx, fy) {
      if (fx === null) {
        readout.textContent = " "
        return
      }
      const t = view.tMin + fx * (view.tMax - view.tMin)
      const s = view.sMax - fy * (view.sMax - view.sMin)
      const [re, im] = zeta(s, t)
      const sign = (v) => (v < 0 ? "−" : "+")
      readout.textContent =
        `s = ${s.toFixed(3)} ${sign(t)} ${Math.abs(t).toFixed(3)}i   ` +
        `ζ(s) = ${fmt(re)} ${sign(im)} ${fmt(Math.abs(im))}i`
    },
  })

  reset.addEventListener("click", () => {
    view = { ...HOME }
    schedule()
  })

  return schedule
}

function fmt(v) {
  return Math.abs(v) >= 1e5 || (v !== 0 && Math.abs(v) < 1e-3) ? v.toExponential(3) : v.toFixed(4)
}
