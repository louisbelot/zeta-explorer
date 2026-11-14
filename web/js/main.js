import { curvePlots } from "./curves.js"
import { dirichletPlots } from "./dirichlet.js"
import { domainPlot } from "./domain.js"
import { zExplorer } from "./explorer.js"
import { primesPlot } from "./primes.js"
import { statsPlots } from "./stats.js"
import { zFunction } from "./zeta.js"

const $ = (id) => document.getElementById(id)
const load = (name) => fetch(new URL(`../data/${name}`, import.meta.url)).then((r) => r.json())

const [rs, zeros, spacing] = await Promise.all(["rs.json", "zeros.json", "spacing.json"].map(load))

const redraw = [
  domainPlot($("domain-img"), $("domain-overlay"), $("domain-readout"), $("domain-reset")),
  zExplorer($("z-plot"), zFunction(rs), zeros.gamma, $("z-readout"), $("zeros")),
  primesPlot($("psi-plot"), zeros.gamma, {
    slider: $("psi-k"),
    count: $("psi-count"),
    play: $("psi-play"),
    xmax: $("psi-x"),
  }),
  dirichletPlots({
    section: $("dirichlet"),
    q: $("l-q"),
    select: $("l-chi"),
    info: $("l-info"),
    plot: $("l-plot"),
    readout: $("l-readout"),
    rows: $("l-rows"),
    rowsReadout: $("l-rows-readout"),
  }),
  curvePlots({
    section: $("elliptic"),
    a: $("ec-a"),
    d: $("ec-d"),
    info: $("ec-info"),
    plot: $("ec-plot"),
    readout: $("ec-readout"),
    twists: $("ec-twists"),
    twistsReadout: $("ec-twists-readout"),
  }),
  statsPlots($("spacing-plot"), $("pair-plot"), $("stats-intro"), spacing),
]

const all = () => redraw.forEach((f) => f())
new ResizeObserver(all).observe(document.querySelector("main"))
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", all)
