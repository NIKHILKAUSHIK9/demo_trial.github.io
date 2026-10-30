// StarNest <-> Streamlit bridge (runs inside the Streamlit component iframe, BEFORE app.js).
// Streamlit Community Cloud exposes one port, so the page cannot call a FastAPI server. Instead every
// fetch('<API_BASE>/...') is serialised, sent to the Streamlit script with setComponentValue, executed there
// in-process on the FastAPI app, and the reply comes back through the component args.
(function () {
  const API = window.STARNEST_API_BASE;                       // fake origin, never really fetched
  const queue = []; let busy = false; let cur = null;
  const uid = Math.random().toString(36).slice(2, 10); let n = 0;
  const send = (type, extra) => window.parent.postMessage(Object.assign({ isStreamlitMessage: true, type, apiVersion: 1 }, extra), '*');
  const b64 = (buf) => { let s = '', a = new Uint8Array(buf); for (let i = 0; i < a.length; i += 32768) s += String.fromCharCode.apply(null, a.subarray(i, i + 32768)); return btoa(s); };
  const unb64 = (t) => { const s = atob(t), a = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) a[i] = s.charCodeAt(i); return a; };

  function pump() {
    if (busy || !queue.length) return;
    busy = true; cur = queue.shift();
    if (cur.aborted) { busy = false; cur = null; return pump(); }
    sendPart(0);
  }
  // the browser -> Streamlit websocket limits one message (~25 MB): big bodies (CSV uploads) go in parts
  const PART = 2000000;
  function sendPart(i) {
    const r = cur.req, parts = Math.max(1, Math.ceil(r.body.length / PART));
    send('streamlit:setComponentValue', { dataType: 'json', value: { id: r.id, part: i, parts, method: r.method, path: r.path, headers: r.headers, body: r.body.slice(i * PART, (i + 1) * PART) } });
  }
  window.addEventListener('message', (ev) => {
    const d = ev.data || {};
    if (d.type !== 'streamlit:render') return;
    const r = d.args && d.args.resp;
    if (cur && r && r.id === cur.req.id && r.ack !== undefined) { sendPart(r.ack + 1); return; }
    if (cur && r && r.id === cur.req.id) {
      const c = cur; cur = null; busy = false;
      c.resolve(new Response(unb64(r.body), { status: r.status, headers: r.headers }));
      pump();
    }
  });

  const realFetch = window.fetch.bind(window);
  window.fetch = async (u, o) => {
    const url = typeof u === 'string' ? u : (u && u.url) || '';
    if (!url.startsWith(API)) return realFetch(u, o);
    const req = new Request(url, o || {});
    const body = ['GET', 'HEAD'].includes(req.method) ? '' : b64(await req.arrayBuffer());
    const headers = {}; req.headers.forEach((v, k) => { headers[k] = v; });
    return new Promise((resolve, reject) => {
      const item = { resolve, aborted: false, req: { id: uid + '-' + (++n), method: req.method, path: url.slice(API.length), headers, body } };
      const sig = o && o.signal;
      if (sig) sig.addEventListener('abort', () => { item.aborted = true; reject(new DOMException('Aborted', 'AbortError')); if (cur === item) { cur = null; busy = false; pump(); } });
      queue.push(item); pump();
    });
  };

  send('streamlit:componentReady', {});
  const resize = () => send('streamlit:setFrameHeight', { height: window.STARNEST_FRAME_HEIGHT || 1000 });
  resize(); window.addEventListener('resize', resize);
})();
