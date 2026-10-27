import { domainPlot } from "./domain.js"
import { zExplorer } from "./explorer.js"
import { zFunction } from "./zeta.js"

const $ = (id) => document.getElementById(id)
const load = (name) => fetch(new URL(`../data/${name}`, import.meta.url)).then((r) => r.json())

const [rs, zeros] = await Promise.all(["rs.json", "zeros.json"].map(load))

const redraw = [
  domainPlot($("domain-img"), $("domain-overlay"), $("domain-readout"), $("domain-reset")),
  zExplorer($("z-plot"), zFunction(rs), zeros.gamma, $("z-readout"), $("zeros")),
]

const all = () => redraw.forEach((f) => f())
new ResizeObserver(all).observe(document.querySelector("main"))
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", all)
