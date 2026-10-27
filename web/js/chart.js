export function niceTicks(lo, hi, target = 6) {
  const span = hi - lo
  if (!(span > 0)) return [lo]
  const raw = span / target
  const mag = 10 ** Math.floor(Math.log10(raw))
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw)
  const first = Math.ceil(lo / step)
  const out = []
  for (let i = first; i * step <= hi + step * 1e-9; i++) out.push(i * step + 0)
  return out
}

export function decimals(step) {
  const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9))
  const m = step * 10 ** d
  return Math.abs(m - Math.round(m)) > 1e-6 ? d + 1 : d
}

export function tickFormat(ticks) {
  const d = ticks.length > 1 ? decimals(ticks[1] - ticks[0]) : 0
  return (v) => v.toFixed(d)
}

export function palette() {
  const css = getComputedStyle(document.documentElement)
  const get = (name) => css.getPropertyValue(name).trim()
  return {
    ink: get("--ink"),
    muted: get("--ink-2"),
    grid: get("--grid"),
    surface: get("--surface"),
    series: get("--series"),
    accent: get("--accent"),
    soft: get("--series-soft"),
    mono: get("--mono"),
  }
}

export class Chart {
  constructor(canvas, pad = { l: 48, r: 12, t: 12, b: 32 }) {
    this.canvas = canvas
    this.ctx = canvas.getContext("2d")
    this.pad = pad
    this.view = [0, 1, 0, 1]
  }

  resize() {
    const dpr = window.devicePixelRatio || 1
    const { width, height } = this.canvas.getBoundingClientRect()
    this.w = width
    this.h = height
    this.canvas.width = Math.round(width * dpr)
    this.canvas.height = Math.round(height * dpr)
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
    this.colors = palette()
  }

  setView(x0, x1, y0, y1) {
    this.view = [x0, x1, y0, y1]
  }

  get plotWidth() {
    return this.w - this.pad.l - this.pad.r
  }

  X(x) {
    const [x0, x1] = this.view
    return this.pad.l + ((x - x0) / (x1 - x0)) * this.plotWidth
  }

  Y(y) {
    const [, , y0, y1] = this.view
    return this.pad.t + ((y1 - y) / (y1 - y0)) * (this.h - this.pad.t - this.pad.b)
  }

  invX(px) {
    const [x0, x1] = this.view
    return x0 + ((px - this.pad.l) / this.plotWidth) * (x1 - x0)
  }

  clear() {
    this.ctx.clearRect(0, 0, this.w, this.h)
  }

  axes({ xTicks, yTicks, xLabel = "", yLabel = "", xFormat } = {}) {
    const { ctx, colors } = this
    const [x0, x1, y0, y1] = this.view
    xTicks ??= niceTicks(x0, x1, Math.max(3, Math.floor(this.plotWidth / 90)))
    yTicks ??= niceTicks(y0, y1, 5)
    xFormat ??= tickFormat(xTicks)
    const yFormat = tickFormat(yTicks)
    ctx.save()
    ctx.lineWidth = 1
    ctx.strokeStyle = colors.grid
    ctx.fillStyle = colors.muted
    ctx.font = `11px ${colors.mono}`
    ctx.textAlign = "center"
    ctx.textBaseline = "top"
    const labelLeft = this.w - this.pad.r - ctx.measureText(xLabel).width - 10
    for (const v of xTicks) {
      const x = Math.round(this.X(v)) + 0.5
      ctx.beginPath()
      ctx.moveTo(x, this.pad.t)
      ctx.lineTo(x, this.h - this.pad.b)
      ctx.stroke()
      const text = xFormat(v)
      if (x + ctx.measureText(text).width / 2 > labelLeft) continue
      ctx.fillText(text, x, this.h - this.pad.b + 6)
    }
    ctx.textAlign = "right"
    ctx.textBaseline = "middle"
    for (const v of yTicks) {
      const y = Math.round(this.Y(v)) + 0.5
      ctx.beginPath()
      ctx.moveTo(this.pad.l, y)
      ctx.lineTo(this.w - this.pad.r, y)
      ctx.stroke()
      ctx.fillText(yFormat(v), this.pad.l - 6, y)
    }
    ctx.fillStyle = colors.ink
    ctx.textAlign = "right"
    ctx.textBaseline = "top"
    ctx.fillText(xLabel, this.w - this.pad.r, this.h - this.pad.b + 6)
    ctx.textAlign = "left"
    ctx.textBaseline = "top"
    ctx.fillText(yLabel, this.pad.l + 4, this.pad.t + 2)
    ctx.restore()
  }

  clip(fn) {
    const { ctx } = this
    ctx.save()
    ctx.beginPath()
    ctx.rect(this.pad.l, this.pad.t, this.plotWidth, this.h - this.pad.t - this.pad.b)
    ctx.clip()
    fn(ctx)
    ctx.restore()
  }

  line(xs, ys, { color, width = 1.5, dash = [] } = {}) {
    this.clip((ctx) => {
      ctx.strokeStyle = color
      ctx.lineWidth = width
      ctx.setLineDash(dash)
      ctx.lineJoin = "round"
      ctx.beginPath()
      for (let i = 0; i < xs.length; i++) {
        const x = this.X(xs[i])
        const y = this.Y(ys[i])
        if (i === 0) ctx.moveTo(x, y)
        else ctx.lineTo(x, y)
      }
      ctx.stroke()
    })
  }

  steps(xs, ys, style) {
    const sx = []
    const sy = []
    for (let i = 0; i < xs.length; i++) {
      if (i > 0) {
        sx.push(xs[i])
        sy.push(ys[i - 1])
      }
      sx.push(xs[i])
      sy.push(ys[i])
    }
    this.line(sx, sy, style)
  }

  bars(edges, heights, color) {
    this.clip((ctx) => {
      ctx.fillStyle = color
      for (let i = 0; i < heights.length; i++) {
        const x0 = this.X(edges[i])
        const x1 = this.X(edges[i + 1])
        const y = this.Y(heights[i])
        ctx.fillRect(x0 + 0.5, y, Math.max(0, x1 - x0 - 1), this.Y(0) - y)
      }
    })
  }

  dots(xs, ys, { color, r = 2.5 } = {}) {
    this.clip((ctx) => {
      ctx.fillStyle = color
      for (let i = 0; i < xs.length; i++) {
        ctx.beginPath()
        ctx.arc(this.X(xs[i]), this.Y(ys[i]), r, 0, 2 * Math.PI)
        ctx.fill()
      }
    })
  }
}
