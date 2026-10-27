import { domainPlot } from "./domain.js"

const $ = (id) => document.getElementById(id)

const redraw = [
  domainPlot($("domain-img"), $("domain-overlay"), $("domain-readout"), $("domain-reset")),
]

const all = () => redraw.forEach((f) => f())
new ResizeObserver(all).observe(document.querySelector("main"))
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", all)
