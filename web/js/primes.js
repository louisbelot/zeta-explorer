import { Chart } from "./chart.js"
import { addZeros, chebyshevPsi, smoothPart } from "./explicit.js"

const SAMPLES = 1600
const PLAY_MS = 12000

export function primesPlot(canvas, gamma, { slider, count, play, xmax }) {
  const chart = new Chart(canvas, { l: 44, r: 12, t: 14, b: 30 })
  const top = Number(slider.max)
  const total = gamma.length.toLocaleString("en-US")
  let xs = []
  let exact = []
  let acc = null
  let K = 0
  let timer = null

  const zerosFor = (v) => (v <= 0 ? 0 : Math.round(gamma.length ** (v / top)))

  function setRange() {
    const hi = Number(xmax.value)
    xs = Array.from({ length: SAMPLES }, (_, i) => 1.5 + ((hi - 1.5) * i) / (SAMPLES - 1))
    exact = chebyshevPsi(xs)
    acc = Float64Array.from(xs, smoothPart)
    K = 0
  }

  function draw() {
    const want = zerosFor(Number(slider.value))
    if (want < K) {
      acc = Float64Array.from(xs, smoothPart)
      K = 0
    }
    addZeros(acc, xs, gamma, K, want)
    K = want
    count.textContent = `${K.toLocaleString("en-US")} of ${total} zeros`

    chart.resize()
    const hi = xs.at(-1)
    chart.setView(0, hi, -0.06 * hi, 1.1 * hi)
    chart.clear()
    chart.axes({ xLabel: "x", yLabel: "ψ(x)" })
    chart.steps(xs, exact, { color: chart.colors.ink, width: 1.3 })
    chart.line(xs, acc, { color: chart.colors.series, width: 1.4 })
  }

  function stop() {
    cancelAnimationFrame(timer)
    timer = null
    play.textContent = "play"
  }

  play.addEventListener("click", () => {
    if (timer !== null) return stop()
    if (Number(slider.value) >= top) slider.value = 0
    const start = performance.now() - (Number(slider.value) / top) * PLAY_MS
    play.textContent = "pause"
    const tick = (now) => {
      slider.value = Math.min(top, Math.round(((now - start) / PLAY_MS) * top))
      draw()
      if (Number(slider.value) >= top) stop()
      else timer = requestAnimationFrame(tick)
    }
    timer = requestAnimationFrame(tick)
  })
  slider.addEventListener("input", () => {
    stop()
    draw()
  })
  xmax.addEventListener("change", () => {
    setRange()
    draw()
  })

  setRange()
  return draw
}
