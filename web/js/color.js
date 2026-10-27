const SAT = 0.78

export function phaseRgb(re, im, out, o) {
  let h = Math.atan2(im, re) / (2 * Math.PI)
  h -= Math.floor(h)
  const b = Math.log2(Math.hypot(re, im))
  const v = Number.isFinite(b) ? 0.62 + 0.38 * (b - Math.floor(b)) : 1
  const k = Math.floor(h * 6)
  const f = h * 6 - k
  const p = v * (1 - SAT)
  const q = v * (1 - f * SAT)
  const t = v * (1 - (1 - f) * SAT)
  const [r, g, bl] = [
    [v, t, p],
    [q, v, p],
    [p, v, t],
    [p, q, v],
    [t, p, v],
    [v, p, q],
  ][k % 6]
  out[o] = 255 * r
  out[o + 1] = 255 * g
  out[o + 2] = 255 * bl
  out[o + 3] = 255
}
