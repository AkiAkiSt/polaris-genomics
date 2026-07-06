// POLARIS — shared utilities. d3 is a global (loaded via CDN before this module).
export const V = 'r2-3';                       // bump when data/JSON change (cache-bust)
export const C = {
  navy:'#0b1f3a', ink:'#12263f', slate:'#5b6b7b', mute:'#8595a4',
  teal:'#1f9e89', tealD:'#137a68', amber:'#e08a00', amberL:'#f0a202',
  crimson:'#b3123a', violet:'#6a4c93', grass:'#3a8d5b', line:'#e4ecf4', line2:'#d7e2ee',
  paper2:'#f4f7fb'
};

export async function load(name){
  const r = await fetch(`data/${name}.json?v=${V}`);
  if(!r.ok) throw new Error(`fetch ${name}: ${r.status}`);
  return r.json();
}
export async function loadAll(names){
  const out = {};
  await Promise.all(names.map(async n => { out[n] = await load(n); }));
  return out;
}

export const fmt = (x, d=2) => (x==null || isNaN(x)) ? '—' : (+x).toFixed(d);
export const pct = (x, d=0) => (x==null || isNaN(x)) ? '—' : (100*x).toFixed(d)+'%';
export const sci = x => (x==null) ? '—' : (x < 1e-3 ? (+x).toExponential(0).replace('e','×10^') : (+x).toFixed(3));

// tooltip
const tipEl = () => document.getElementById('tip');
export function tip(html, ev){
  const t = tipEl(); if(!t) return;
  t.innerHTML = html; t.classList.add('show');
  const pad = 14, w = t.offsetWidth, h = t.offsetHeight;
  let x = ev.clientX + pad, y = ev.clientY + pad;
  if(x + w > innerWidth - 8) x = ev.clientX - w - pad;
  if(y + h > innerHeight - 8) y = ev.clientY - h - pad;
  t.style.left = x+'px'; t.style.top = y+'px';
}
export function untip(){ const t = tipEl(); if(t) t.classList.remove('show'); }

// scroll-reveal
export function observeReveal(){
  const io = new IntersectionObserver((es)=>es.forEach(e=>{
    if(e.isIntersecting){ e.target.classList.add('in'); io.unobserve(e.target); }
  }), {threshold:0.12, rootMargin:'0px 0px -8% 0px'});
  document.querySelectorAll('.reveal').forEach(el=>io.observe(el));
}

// nav shadow + active section
export function navScroll(){
  const nav = document.getElementById('nav');
  const links = [...document.querySelectorAll('nav .links a')];
  const secs = links.map(a=>document.querySelector(a.getAttribute('href'))).filter(Boolean);
  const onScroll = () => {
    nav.classList.toggle('scrolled', scrollY > 8);
    let cur = null;
    for(const s of secs){ if(s.getBoundingClientRect().top < 120) cur = s.id; }
    links.forEach(a=>a.classList.toggle('active', a.getAttribute('href')==='#'+cur));
  };
  addEventListener('scroll', onScroll, {passive:true}); onScroll();
}

// tiny Gaussian KDE for smooth distribution curves
export function kde(samples, {bw=null, n=64, lo=null, hi=null}={}){
  const s = samples.filter(v=>v!=null && !isNaN(v));
  if(s.length < 2) return {pts:[], max:0};
  const mn = lo==null ? Math.min(...s) : lo, mx = hi==null ? Math.max(...s) : hi;
  const mean = s.reduce((a,b)=>a+b,0)/s.length;
  const sd = Math.sqrt(s.reduce((a,b)=>a+(b-mean)**2,0)/s.length) || 1;
  const h = bw || 1.06*sd*Math.pow(s.length,-0.2) || 0.1;
  const pts = d3.range(n).map(i=>{
    const x = mn + (mx-mn)*i/(n-1);
    let y = 0; for(const v of s) y += Math.exp(-0.5*((x-v)/h)**2);
    return {x, y: y/(s.length*h*Math.sqrt(2*Math.PI))};
  });
  return {pts, max:d3.max(pts,p=>p.y), n:s.length};
}

// responsive svg helper: returns {svg,g,w,h,inner}
export function svgFrame(sel, {h=200, m={t:16,r:16,b:28,l:34}}={}){
  const host = d3.select(sel); host.selectAll('*').remove();
  const w = host.node().clientWidth || 460;
  const svg = host.append('svg').attr('width',w).attr('height',h).attr('viewBox',`0 0 ${w} ${h}`);
  const g = svg.append('g').attr('transform',`translate(${m.l},${m.t})`);
  return {svg, g, w, h, iw:w-m.l-m.r, ih:h-m.t-m.b, m};
}
export const debounce = (fn,ms=180)=>{let t;return(...a)=>{clearTimeout(t);t=setTimeout(()=>fn(...a),ms)}};
