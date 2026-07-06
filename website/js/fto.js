// POLARIS — FTO rs1421085 interactive vignette (locus map · sequence logo · animated ODE)
import {C, fmt} from './util.js';

const BASE = ['A','C','G','T'];
const BASECOL = {A:'#4bd48f', C:'#5aa9e6', G:'#f0b34a', T:'#ff6f86'};
let F;

export function renderFto(f){
  if(!f || !f.motif){ return; }
  F = f;
  const app = document.getElementById('fto-app');
  const g = f.genes;
  const near = g.find(x=>x.role==='nearest') || g.find(x=>x.nearest);
  const l2g = g.find(x=>x.gene==='FTO');
  const truth = g.find(x=>x.gene==='IRX3');
  app.innerHTML = `
    <div class="fto-locus">
      <h3>The locus — half a megabase of ambiguity</h3>
      <div class="sub">Where does rs1421085 act? The three methods point to three different genes.</div>
      <div id="ftoMap"></div>
      <div class="answers3">
        <div class="ans near"><div class="who">nearest gene</div><div class="g">${near?.gene||'RPGRIP1L'}</div><div class="d">${kb(near?.dist)} away · <span style="color:#ff9db0">wrong</span></div></div>
        <div class="ans l2g"><div class="who">L2G / proximal tools</div><div class="g">FTO</div><div class="d">${kb(l2g?.dist)} away · <span style="color:#ff9db0">wrong</span></div></div>
        <div class="ans truth"><div class="who">experimental truth</div><div class="g">IRX3</div><div class="d">${kb(truth?.dist)} away · <span style="color:#ffd98a">distal</span></div></div>
      </div>
    </div>
    <div class="fto-grid">
      <div class="fto-panel">
        <h3>The mechanism — a broken ARID5B grip</h3>
        <div class="sub">POLARIS scores every motif transparently. The risk allele destroys an ARID5B site.</div>
        <div id="ftoLogo"></div>
        <div id="ftoSeq"></div>
      </div>
      <div class="fto-panel">
        <h3>The consequence — a gene turned up <button class="replay" id="ftoReplay">▶ replay</button></h3>
        <div class="sub">Feed the binding change into a gene-regulatory ODE. The risk allele derepresses IRX3.</div>
        <div class="occ" id="ftoOcc"></div>
        <div id="ftoOde"></div>
        <div class="concl">Risk allele → ARID5B falls off → enhancer derepressed → <b>IRX3/IRX5 ↑</b> → less adipocyte browning.
          Direction is robust across model parameters, matching Claussnitzer <i>et&nbsp;al.</i> 2015.</div>
      </div>
    </div>`;
  locusMap('#ftoMap', f);
  motifLogo('#ftoLogo', f);
  seqSite('#ftoSeq', f);
  occBars('#ftoOcc', f);
  odePlot('#ftoOde', f);
  document.getElementById('ftoReplay').addEventListener('click', ()=>{ occBars('#ftoOcc', f, true); odePlot('#ftoOde', f, true); });
}
const kb = d => d==null ? '' : (d>=1000 ? (d/1000).toFixed(0)+' kb' : d+' bp');

function frame(sel, h, m={t:14,r:18,b:26,l:18}){
  const host = d3.select(sel); host.selectAll('*').remove();
  const w = host.node().clientWidth || 460;
  const svg = host.append('svg').attr('width',w).attr('height',h);
  return {svg, g:svg.append('g').attr('transform',`translate(${m.l},${m.t})`), iw:w-m.l-m.r, ih:h-m.t-m.b};
}

function locusMap(sel, f){
  const A = frame(sel, 128, {t:30,r:26,b:34,l:26});
  const genes = f.genes, vpos = f.variant.pos;
  const xs = d3.extent([...genes.map(g=>g.tss), vpos]);
  const pad = (xs[1]-xs[0])*0.06;
  const x = d3.scaleLinear().domain([xs[0]-pad, xs[1]+pad]).range([0, A.iw]);
  const y = A.ih*0.55;
  A.g.append('line').attr('x1',0).attr('x2',A.iw).attr('y1',y).attr('y2',y).attr('stroke','rgba(255,255,255,.25)').attr('stroke-width',1.5);
  // genes
  genes.forEach(gn=>{
    const gx=x(gn.tss), key=gn.gene;
    const col = key==='IRX3'?'#f0b34a':key==='FTO'?'#3fd0b6':(gn.nearest?'#a9c0d6':'rgba(180,200,220,.5)');
    const big = ['IRX3','FTO','RPGRIP1L'].includes(key);
    A.g.append('circle').attr('cx',gx).attr('cy',y).attr('r',big?5:3).attr('fill',col);
    if(big){ A.g.append('text').attr('x',gx).attr('y',y-12).attr('text-anchor','middle').attr('font-size',12).attr('font-weight',700).attr('fill',col).text(key); }
  });
  // variant
  const vx=x(vpos);
  A.g.append('path').attr('d',`M${vx-6},${y+14} L${vx+6},${y+14} L${vx},${y+4} Z`).attr('fill','#ffd98a');
  A.g.append('text').attr('x',vx).attr('y',y+30).attr('text-anchor','middle').attr('font-size',11).attr('fill','#ffd98a').text('rs1421085');
  // distance arc variant->IRX3
  const irx=genes.find(g=>g.gene==='IRX3');
  if(irx){ const ix=x(irx.tss);
    A.g.append('path').attr('d',`M${vx},${y-4} Q${(vx+ix)/2},${y-30} ${ix},${y-4}`).attr('fill','none').attr('stroke','#f0b34a').attr('stroke-width',1.3).attr('stroke-dasharray','4 3').attr('opacity',.8);
    A.g.append('text').attr('x',(vx+ix)/2).attr('y',y-32).attr('text-anchor','middle').attr('font-size',10.5).attr('fill','#f0b34a').text('acts 520 kb away → IRX3');
  }
}

// D3 information-content sequence logo (getBBox-scaled glyphs)
function motifLogo(sel, f){
  const m = f.motif, w = m.w;
  const A = frame(sel, 150, {t:12,r:10,b:24,l:26});
  const cw = A.iw / w, maxIC = 2;
  const y = d3.scaleLinear().domain([0,maxIC]).range([A.ih,0]);
  // y axis (bits)
  A.g.append('line').attr('x1',0).attr('x2',0).attr('y1',0).attr('y2',A.ih).attr('stroke','rgba(255,255,255,.2)');
  [0,1,2].forEach(t=>{A.g.append('text').attr('x',-6).attr('y',y(t)+3).attr('text-anchor','end').attr('font-size',9).attr('fill','#9fbad3').text(t);});
  A.g.append('text').attr('transform',`translate(-20,${A.ih/2}) rotate(-90)`).attr('text-anchor','middle').attr('font-size',10).attr('fill','#9fbad3').text('bits');
  for(let i=0;i<w;i++){
    const col=[0,1,2,3].map(b=>({b:BASE[b], p:m.prob[b][i]})).sort((a,b)=>a.p-b.p);
    let yb=A.ih;
    col.forEach(({b,p})=>{
      const hgt=p*m.ic[i]/maxIC*A.ih; if(hgt<0.8) return;
      const t=A.g.append('text').text(b).attr('font-family','Inter, sans-serif').attr('font-weight',800)
        .attr('font-size',100).attr('fill',BASECOL[b]).attr('x',0).attr('y',0);
      const bb=t.node().getBBox();
      const sx=(cw*0.92)/bb.width, sy=hgt/bb.height;
      t.attr('transform',`translate(${i*cw+cw*0.04 - bb.x*sx}, ${yb - hgt - bb.y*sy}) scale(${sx},${sy})`);
      yb-=hgt;
    });
    if(i===m.var_offset){ A.g.append('rect').attr('x',i*cw).attr('y',0).attr('width',cw).attr('height',A.ih)
      .attr('fill','none').attr('stroke','#ffd98a').attr('stroke-width',1.5).attr('rx',3).attr('opacity',.9); }
  }
}

function seqSite(sel, f){
  const m=f.motif; const host=document.querySelector(sel);
  const cells=(seq,varIdx,alt)=>[...seq].map((ch,i)=>
    `<div class="b base${ch} ${i===varIdx?'var':''}">${ch}</div>`).join('');
  host.innerHTML = `
    <div style="display:flex;align-items:center;gap:.8rem;flex-wrap:wrap">
      <div><div style="font-size:.72rem;color:#9fbad3">reference (T)</div><div class="seqrow">${cells(m.site_ref,m.var_offset)}</div></div>
      <div style="font-size:1.4rem;color:#9fbad3">→</div>
      <div><div style="font-size:.72rem;color:#ff9db0">risk allele (C)</div><div class="seqrow">${cells(m.site_alt,m.var_offset)}</div></div>
    </div>
    <div class="llrline">match score <span class="big" style="color:#4bd48f">${fmt(m.llr_ref,1)}</span>
      <span style="color:#9fbad3">→</span> <span class="big" style="color:#ff6f86">${fmt(m.llr_alt,1)}</span> bits
      <span style="margin-left:.4rem;color:#ffd98a;font-weight:600">motif abolished (Δ ${fmt(m.dLLR,1)})</span></div>`;
}

function occBars(sel, f, animate){
  const o=f.ode; const host=document.querySelector(sel);
  host.innerHTML = `
    <div class="o"><div class="lab">ARID5B occupancy · reference</div><div class="track"><i style="width:0;background:#3fd0b6" data-w="${(o.occ_ref*100).toFixed(0)}%"></i></div><div style="font-size:.72rem;color:#7fd8c8">bound (${fmt(o.occ_ref,2)})</div></div>
    <div class="o"><div class="lab">ARID5B occupancy · risk</div><div class="track"><i style="width:0;background:#ff6f86" data-w="${(o.occ_alt*100).toFixed(0)}%"></i></div><div style="font-size:.72rem;color:#ff9db0">lost (${fmt(o.occ_alt,2)})</div></div>`;
  requestAnimationFrame(()=>host.querySelectorAll('i').forEach(i=>setTimeout(()=>i.style.width=i.dataset.w,40)));
}

function odePlot(sel, f, animate){
  const o=f.ode, A=frame(sel, 200, {t:14,r:16,b:30,l:36});
  const x=d3.scaleLinear().domain([0,d3.max(o.t)]).range([0,A.iw]);
  const ymax=d3.max([...o.x_ref,...o.x_alt])*1.1;
  const y=d3.scaleLinear().domain([0,ymax]).range([A.ih,0]);
  A.g.append('g').attr('transform',`translate(0,${A.ih})`).call(d3.axisBottom(x).ticks(5)).selectAll('text,line,path').attr('stroke','rgba(255,255,255,.25)').attr('fill','#9fbad3');
  A.g.append('g').call(d3.axisLeft(y).ticks(4)).selectAll('text,line,path').attr('stroke','rgba(255,255,255,.25)').attr('fill','#9fbad3');
  A.g.append('text').attr('x',A.iw/2).attr('y',A.ih+26).attr('text-anchor','middle').attr('font-size',10.5).attr('fill','#9fbad3').text('time (a.u.)');
  A.g.append('text').attr('transform',`translate(-26,${A.ih/2}) rotate(-90)`).attr('text-anchor','middle').attr('font-size',10.5).attr('fill','#9fbad3').text('IRX3 expression');
  const line=d3.line().x((d,i)=>x(o.t[i])).y(d=>y(d)).curve(d3.curveMonotoneX);
  const draw=(data,col,lab,yl)=>{
    const p=A.g.append('path').datum(data).attr('d',line).attr('fill','none').attr('stroke',col).attr('stroke-width',2.6);
    const L=p.node().getTotalLength();
    p.attr('stroke-dasharray',`${L} ${L}`).attr('stroke-dashoffset',L).transition().duration(1500).ease(d3.easeCubicOut).attr('stroke-dashoffset',0);
    A.g.append('text').attr('x',A.iw-2).attr('y',y(data[data.length-1])+yl).attr('text-anchor','end').attr('font-size',11).attr('font-weight',600).attr('fill',col).text(lab);
  };
  draw(o.x_ref, '#3fd0b6', 'reference', 12);
  draw(o.x_alt, '#ff6f86', 'risk allele', -6);
}
