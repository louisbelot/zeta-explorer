import { phaseRgb } from "./color.js"
import { zeta } from "./zeta.js"

export function renderBand({ width, height, row0, rows, view, step }) {
  const { tMin, tMax, sMin, sMax } = view
  const px = new Uint8ClampedArray(width * rows * 4)
  const cell = new Uint8ClampedArray(4)
  for (let j = 0; j < rows; j += step) {
    const sigma = sMax - ((row0 + j + step / 2) / height) * (sMax - sMin)
    for (let i = 0; i < width; i += step) {
      const t = tMin + ((i + step / 2) / width) * (tMax - tMin)
      const [re, im] = zeta(sigma, t)
      phaseRgb(re, im, cell, 0)
      for (let y = j; y < Math.min(j + step, rows); y++) {
        for (let x = i; x < Math.min(i + step, width); x++) px.set(cell, 4 * (y * width + x))
      }
    }
  }
  return px
}

if (typeof WorkerGlobalScope !== "undefined" && self instanceof WorkerGlobalScope) {
  self.onmessage = ({ data }) => {
    const px = renderBand(data)
    self.postMessage({ id: data.id, row0: data.row0, rows: data.rows, px }, [px.buffer])
  }
}
