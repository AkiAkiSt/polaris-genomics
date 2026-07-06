// POLARIS — D3 charts (findings distributions + validation). d3 is global.
import {C, fmt, pct, tip, untip, kde, svgFrame} from './util.js';

// Two overlapping density curves: causal (teal) vs passenger (grey)
export function densityCompare(sel, d, {xlabel, reading, flip=false}={}){
  const F = svgFrame(sel, {h:210, m:{t:14,r:14,b:42,l:14}});
  if(!d || (!d.causal?.length && !d.passenger?.length)){ F.g.append('text').text('no data').attr('fill',C.mute); return; }
  const all = [...d.causal, ...d.passenger];
  const lo = d3.min(all), hi = d3.max(all);
  const kC = kde(d.causal,{lo,hi}), kP = kde(d.passenger,{lo,hi});
  const x = d3.scaleLinear().domain([lo,hi]).range([0,F.iw]);
  const y = d3.scaleLinear().domain([0, Math.max(kC.max,kP.max)*1.12]).range([F.ih,0]);
  const area = d3.area().x(p=>x(p.x)).y0(F.ih).y1(p=>y(p.y)).curve(d3.curveBasis);
  const line = d3.line().x(p=>x(p.x)).y(p=>y(p.y)).curve(d3.curveBasis);
  // passenger behind
  F.g.append('path').datum(kP.pts).attr('d',area).attr('fill',C.slate).attr('opacity',.16);
  F.g.append('path').datum(kP.pts).attr('d',line).attr('fill','none').attr('stroke',C.slate).attr('stroke-width',1.4).attr('opacity',.6);
  F.g.append('path').datum(kC.pts).attr('d',area).attr('fill',C.teal).attr('opacity',.20);
  F.g.append('path').datum(kC.pts).attr('d',line).attr('fill','none').attr('stroke',C.teal).attr('stroke-width',2.4);
  // medians
  const med = a=>{const s=[...a].sort((p,q)=>p-q);return s[Math.floor(s.length/2)];};
  [[med(d.causal),C.teal],[med(d.passenger),C.slate]].forEach(([m,col])=>{
    F.g.append('line').attr('x1',x(m)).attr('x2',x(m)).attr('y1',F.ih).attr('y2',F.ih-8).attr('stroke',col).attr('stroke-width',2);
  });
  F.g.append('g').attr('class','axis').attr('transform',`translate(0,${F.ih})`).call(d3.axisBottom(x).ticks(5));
  F.g.append('text').attr('x',F.iw/2).attr('y',F.ih+34).attr('text-anchor','middle').attr('fill',C.slate).attr('font-size',11).text(xlabel);
  // legend
  const lg = F.g.append('g').attr('font-size',10.5).attr('transform',`translate(${F.iw-4},2)`).attr('text-anchor','end');
  lg.append('text').attr('fill',C.tealD).attr('font-weight',600).text('● fine-mapped causal');
  lg.append('text').attr('y',14).attr('fill',C.slate).text('● LD passenger');
  if(reading){ F.svg.append('text').attr('x',14).attr('y',18).attr('fill',flip?C.crimson:C.tealD).attr('font-size',11.5).attr('font-weight',600).text(reading); }
}

// grouped bars: POLARIS vs nearest, per disease
export function geneBars(sel, meta){
  const F = svgFrame(sel, {h:230, m:{t:16,r:12,b:40,l:40}});
  const dz = [['T2D',meta.T2D],['CAD',meta.CAD]];
  const groups = dz.map(([k,d])=>({k, POLARIS:d.gene_top1, nearest:d.gene_nearest}));
  const x0 = d3.scaleBand().domain(groups.map(g=>g.k)).range([0,F.iw]).padding(.34);
  const x1 = d3.scaleBand().domain(['POLARIS','nearest']).range([0,x0.bandwidth()]).padding(.18);
  const y = d3.scaleLinear().domain([0,1]).range([F.ih,0]);
  F.g.append('g').attr('class','axis').call(d3.axisLeft(y).ticks(5).tickFormat(d3.format('.0%')));
  F.g.append('g').attr('class','axis').attr('transform',`translate(0,${F.ih})`).call(d3.axisBottom(x0));
  const col = {POLARIS:C.teal, nearest:C.slate};
  groups.forEach(g=>{
    ['POLARIS','nearest'].forEach(key=>{
      const val=g[key];
      F.g.append('rect').attr('x',x0(g.k)+x1(key)).attr('y',y(val)).attr('width',x1.bandwidth())
        .attr('height',F.ih-y(val)).attr('fill',col[key]).attr('rx',3)
        .on('mousemove',e=>tip(`<b>${g.k}</b> · ${key==='POLARIS'?'POLARIS (L2G linkage)':'nearest gene'}<br>top-1 = <b>${pct(val)}</b>`,e)).on('mouseleave',untip);
      F.g.append('text').attr('x',x0(g.k)+x1(key)+x1.bandwidth()/2).attr('y',y(val)-5).attr('text-anchor','middle')
        .attr('font-size',11).attr('font-weight',600).attr('fill',col[key]).text(pct(val));
    });
  });
  const lg=F.svg.append('g').attr('font-size',11).attr('transform',`translate(${F.w-12},14)`).attr('text-anchor','end');
  lg.append('text').attr('fill',C.tealD).attr('font-weight',600).text('■ POLARIS');
  lg.append('text').attr('y',14).attr('fill',C.slate).text('■ nearest gene');
}

// calibration reliability curve
export function calibCurve(sel, cal){
  const F = svgFrame(sel, {h:230, m:{t:16,r:16,b:38,l:40}});
  const x=d3.scaleLinear().domain([0,1]).range([0,F.iw]), y=d3.scaleLinear().domain([0,1]).range([F.ih,0]);
  F.g.append('line').attr('x1',0).attr('y1',F.ih).attr('x2',F.iw).attr('y2',0).attr('stroke',C.line2).attr('stroke-dasharray','4 4');
  F.g.append('g').attr('class','axis').call(d3.axisLeft(y).ticks(5).tickFormat(d3.format('.0%')));
  F.g.append('g').attr('class','axis').attr('transform',`translate(0,${F.ih})`).call(d3.axisBottom(x).ticks(5).tickFormat(d3.format('.0%')));
  const pts=cal.points||[];
  const line=d3.line().x(p=>x(p.pred)).y(p=>y(p.obs)).curve(d3.curveMonotoneX);
  F.g.append('path').datum(pts).attr('d',line).attr('fill','none').attr('stroke',C.teal).attr('stroke-width',2.6);
  F.g.selectAll('circle').data(pts).join('circle').attr('cx',p=>x(p.pred)).attr('cy',p=>y(p.obs))
    .attr('r',p=>4+Math.sqrt(p.n)).attr('fill',C.teal).attr('opacity',.85)
    .on('mousemove',(e,p)=>tip(`predicted <b>${pct(p.pred)}</b> → observed <b>${pct(p.obs)}</b><br>${p.n} variants`,e)).on('mouseleave',untip);
  F.g.append('text').attr('x',F.iw).attr('y',F.ih-6).attr('text-anchor','end').attr('font-size',10.5).attr('fill',C.mute).text('perfect calibration = diagonal');
  F.g.append('text').attr('x',F.iw/2).attr('y',F.ih+32).attr('text-anchor','middle').attr('font-size',11).attr('fill',C.slate).text('predicted P(causal)');
}

// reliability weights with bootstrap CIs
export function relBars(sel, rel){
  const ch=[...rel.channels].sort((a,b)=>a.weight-b.weight);
  const F = svgFrame(sel, {h:230, m:{t:12,r:60,b:26,l:140}});
  const y=d3.scaleBand().domain(ch.map(c=>c.label)).range([F.ih,0]).padding(.3);
  const x=d3.scaleLinear().domain([0,d3.max(ch,c=>c.ci_hi)*1.05||1]).range([0,F.iw]);
  F.g.append('g').attr('class','axis').attr('transform',`translate(0,${F.ih})`).call(d3.axisBottom(x).ticks(4));
  F.g.append('g').attr('class','axis').call(d3.axisLeft(y)).selectAll('text').attr('font-size',10.5);
  ch.forEach(c=>{
    const yc=y(c.label)+y.bandwidth()/2, top=c.weight===ch[ch.length-1].weight;
    F.g.append('rect').attr('x',0).attr('y',y(c.label)).attr('height',y.bandwidth()).attr('width',x(c.weight))
      .attr('fill',top?C.teal:C.slate).attr('rx',3)
      .on('mousemove',e=>tip(`<b>${c.label}</b><br>weight ${fmt(c.weight)} · d′ ${fmt(c.dprime)}<br>95% CI [${fmt(c.ci_lo)}, ${fmt(c.ci_hi)}]`,e)).on('mouseleave',untip);
    F.g.append('line').attr('x1',x(c.ci_lo)).attr('x2',x(c.ci_hi)).attr('y1',yc).attr('y2',yc).attr('stroke',C.ink).attr('stroke-width',1.4);
    ['ci_lo','ci_hi'].forEach(k=>F.g.append('line').attr('x1',x(c[k])).attr('x2',x(c[k])).attr('y1',yc-4).attr('y2',yc+4).attr('stroke',C.ink).attr('stroke-width',1.4));
    F.g.append('text').attr('x',x(c.ci_hi)+6).attr('y',yc+3.5).attr('font-size',11).attr('font-weight',600).attr('fill',top?C.tealD:C.slate).text(fmt(c.weight));
  });
}

// Marchenko–Pastur spectrum
export function mpSpectrum(sel, mp){
  const F = svgFrame(sel, {h:230, m:{t:16,r:16,b:38,l:34}});
  const ev=mp.evals, x=d3.scaleLinear().domain([0,d3.max(ev)*1.08]).range([0,F.iw]);
  // MP density
  const g=mp.gamma, lm=(1-Math.sqrt(g))**2, lp=mp.edge;
  const dens=d3.range(80).map(i=>{const v=0.02+(d3.max(ev)*1.05)*i/79;
    const y=(v>lm&&v<lp)?Math.sqrt((lp-v)*(v-lm))/(2*Math.PI*g*v):0; return {v,y};});
  const y=d3.scaleLinear().domain([0,d3.max(dens,d=>d.y)*1.1||1]).range([F.ih,0]);
  F.g.append('path').datum(dens).attr('d',d3.area().x(d=>x(d.v)).y0(F.ih).y1(d=>y(d.y)).curve(d3.curveBasis)).attr('fill',C.slate).attr('opacity',.14);
  F.g.append('path').datum(dens).attr('d',d3.line().x(d=>x(d.v)).y(d=>y(d.y)).curve(d3.curveBasis)).attr('fill','none').attr('stroke',C.slate).attr('stroke-width',1.4);
  F.g.append('line').attr('x1',x(mp.edge)).attr('x2',x(mp.edge)).attr('y1',0).attr('y2',F.ih).attr('stroke',C.crimson).attr('stroke-dasharray','4 3').attr('stroke-width',1.4);
  F.g.append('text').attr('x',x(mp.edge)).attr('y',-4).attr('text-anchor','middle').attr('font-size',10).attr('fill',C.crimson).text('noise edge');
  ev.forEach(e=>{const sig=e>mp.edge;
    F.g.append('line').attr('x1',x(e)).attr('x2',x(e)).attr('y1',F.ih).attr('y2',sig?F.ih*0.15:F.ih*0.6)
      .attr('stroke',sig?C.amber:C.slate).attr('stroke-width',sig?3:1.6).attr('opacity',sig?1:.6)
      .on('mousemove',ev2=>tip(`eigenvalue <b>${fmt(e)}</b>${sig?' — <b>shared signal</b>':' — noise-level'}`,ev2)).on('mouseleave',untip);
  });
  F.g.append('g').attr('class','axis').attr('transform',`translate(0,${F.ih})`).call(d3.axisBottom(x).ticks(5));
  F.g.append('text').attr('x',F.iw/2).attr('y',F.ih+32).attr('text-anchor','middle').attr('font-size',11).attr('fill',C.slate).text('eigenvalue of channel-correlation matrix');
  F.svg.append('text').attr('x',F.w-14).attr('y',18).attr('text-anchor','end').attr('font-size',11.5).attr('font-weight',600).attr('fill',C.amber).text(`${mp.n_signal}/${mp.n_channels} channels share signal`);
}
