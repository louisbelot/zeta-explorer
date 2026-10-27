import { Chart } from "./chart.js"

export function statsPlots(spacingCanvas, pairCanvas, intro, data) {
  const [m, e] = data.t0.toExponential(2).split("e+")
  const t0 = `${m}·10${[...e].map((c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[c]).join("")}`
  intro.textContent +=
    ` The ${data.count.toLocaleString("en-US")} zeros below start at n = 10¹² (t ≈ ${t0}).` +
    ` Their spacings sit ${data.ks.toFixed(4)} from the GUE law in Kolmogorov-Smirnov distance.`

  const spacing = new Chart(spacingCanvas, { l: 40, r: 10, t: 12, b: 30 })
  const pair = new Chart(pairCanvas, { l: 40, r: 10, t: 12, b: 30 })

  return () => {
    const sp = data.spacing
    spacing.resize()
    spacing.setView(0, sp.s.at(-1), 0, 1.15)
    spacing.clear()
    spacing.axes({ xLabel: "s", yLabel: "density" })
    spacing.bars(sp.edges, sp.density, spacing.colors.soft)
    spacing.line(sp.s, sp.poisson, { color: spacing.colors.muted, width: 1.4, dash: [2, 3] })
    spacing.line(sp.s, sp.wigner, { color: spacing.colors.muted, width: 1.4, dash: [6, 4] })
    spacing.line(sp.s, sp.gue, { color: spacing.colors.series, width: 2 })

    const pc = data.pair
    pair.resize()
    pair.setView(0, pc.uu.at(-1), 0, 1.3)
    pair.clear()
    pair.axes({ xLabel: "u", yLabel: "R₂(u)" })
    pair.line(pc.uu, pc.montgomery, { color: pair.colors.series, width: 2 })
    pair.dots(pc.u, pc.r, { color: pair.colors.ink, r: 2.4 })
  }
}
