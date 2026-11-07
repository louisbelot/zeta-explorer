import { conrey } from "./characters.js"
import { lZeros } from "./lfunction.js"

self.onmessage = ({ data: { id, q, n, T } }) => {
  const zeros = Float64Array.from(lZeros(conrey(q, n), -T, T))
  self.postMessage({ id, n, zeros }, [zeros.buffer])
}
