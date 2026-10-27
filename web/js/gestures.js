export function panZoom(el, { pan, zoom, hover }) {
  const pointers = new Map()
  let pinch = 0
  let active = false

  const local = (e) => {
    const r = el.getBoundingClientRect()
    return [(e.clientX - r.left) / r.width, (e.clientY - r.top) / r.height]
  }
  const spread = () => {
    const [p, q] = [...pointers.values()]
    return Math.hypot(p.clientX - q.clientX, p.clientY - q.clientY)
  }

  el.addEventListener("pointerdown", (e) => {
    active = true
    el.setPointerCapture(e.pointerId)
    pointers.set(e.pointerId, e)
    if (pointers.size === 2) pinch = spread()
  })

  el.addEventListener("pointermove", (e) => {
    const last = pointers.get(e.pointerId)
    if (!last) {
      hover?.(...local(e))
      return
    }
    pointers.set(e.pointerId, e)
    const r = el.getBoundingClientRect()
    if (pointers.size === 1) {
      pan((e.clientX - last.clientX) / r.width, (e.clientY - last.clientY) / r.height)
    } else if (pointers.size === 2) {
      const d = spread()
      const [p, q] = [...pointers.values()]
      const mid = { clientX: (p.clientX + q.clientX) / 2, clientY: (p.clientY + q.clientY) / 2 }
      if (pinch > 0 && d > 0) zoom(pinch / d, ...local(mid))
      pinch = d
    }
  })

  const release = (e) => {
    pointers.delete(e.pointerId)
    pinch = 0
  }
  el.addEventListener("pointerup", release)
  el.addEventListener("pointercancel", release)
  el.addEventListener("pointerleave", () => {
    active = false
    hover?.(null, null)
  })

  el.addEventListener(
    "wheel",
    (e) => {
      if (!active && !e.ctrlKey) return
      e.preventDefault()
      zoom(Math.exp(e.deltaY * (e.deltaMode === 1 ? 0.05 : 0.0015)), ...local(e))
    },
    { passive: false },
  )
}

export function frameLoop(draw) {
  let queued = false
  return () => {
    if (queued) return
    queued = true
    requestAnimationFrame(() => {
      queued = false
      draw()
    })
  }
}
