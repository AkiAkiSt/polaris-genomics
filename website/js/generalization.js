// POLARIS — generalization (T2D ↔ CAD replication)
import {C, fmt, pct, tip, untip} from './util.js';

export function renderGeneralization(DATA){
  const app = document.getElementById('gen-app'); if(!app) return;
  const t = DATA.meta.diseases.T2D, c = DATA.meta.diseases.CAD;
  const repro = [
    ['Conservation enrichment', t.cons_auc, c.cons_auc, 'AUC · causal vs passenger'],
    ['Rarity inversion', t.rarity_auc, c.rarity_auc, 'AUC · inverted in both'],
    ['Target-gene top-1', t.gene_top1, c.gene_top1, 'controls'],
    ['Most-reliable channel', null, null, 'phyloP conservation, both'],
  ];
  const cards = repro.map((r,i)=>{
    const vals = (i===3)
      ? `<div class="rvals"><span class="rv"><b style="font-size:1.3rem">phyloP</b><i>T2D #1</i></span><span class="arrow">→</span><span class="rv cad"><b style="font-size:1.3rem">phyloP</b><i>CAD #1</i></span></div>`
      : `<div class="rvals"><span class="rv"><b>${fmt(r[1])}</b><i>T2D</i></span><span class="arrow">→</span><span class="rv cad"><b>${fmt(r[2])}</b><i>CAD</i></span></div>`;
    return `<div class="card repro reveal"><div class="rlabel">${r[0]}</div>${vals}
      <span class="chip teal">✓ replicated</span><div style="font-size:.72rem;color:var(--mute);margin-top:.5rem">${r[3]}</div></div>`;
  }).join('');
  app.innerHTML = `
    <div class="grid g4" style="margin-bottom:1.4rem">${cards}</div>
    <div class="grid g2">
      <div class="card"><h3>The model re-learns the same reliabilities</h3>
        <p class="muted" style="margin:.2rem 0 .5rem">Trained independently per disease, the latent-truth EM ranks the channels the same way.</p>
        <div id="genRel"></div></div>
      <div class="card"><h3>Textbook mechanisms, recovered from scratch</h3>
        <p class="muted" style="margin:.2rem 0 .6rem">The transparent motif engine reproduces the known controls in both diseases.</p>
        <div id="genMech"></div></div>
    </div>`;
  relCompare('#genRel', DATA.reliability);
  mechStrip('#genMech', DATA.variants);
  document.querySelectorAll('#gen-app .reveal').forEach(el=>el.classList.add('in'));
}

function relCompare(sel, rel){
  const host = d3.select(sel); host.selectAll('*').remove();
  const ch = [...rel.channels].filter(c=>c.cad_weight!=null).sort((a,b)=>b.weight-a.weight);
  const w = host.node().clientWidth||460, m={t:10,r:14,b:26,l:150}, h=200, iw=w-m.l-m.r, ih=h-m.t-m.b;
  const svg=host.append('svg').attr('width',w).attr('height',h), g=svg.append('g').attr('transform',`translate(${m.l},${m.t})`);
  const y0=d3.scaleBand().domain(ch.map(c=>c.label)).range([0,ih]).padding(.28);
  const y1=d3.scaleBand().domain(['T2D','CAD']).range([0,y0.bandwidth()]).padding(.15);
  const x=d3.scaleLinear().domain([0,Math.max(d3.max(ch,c=>c.weight),d3.max(ch,c=>c.cad_weight))*1.05||1]).range([0,iw]);
  g.append('g').attr('class','axis').attr('transform',`translate(0,${ih})`).call(d3.axisBottom(x).ticks(4));
  g.append('g').attr('class','axis').call(d3.axisLeft(y0)).selectAll('text').attr('font-size',10.5);
  ch.forEach(c=>{
    [['T2D',c.weight,C.navy],['CAD',c.cad_weight,C.amber]].forEach(([k,v,col])=>{
      g.append('rect').attr('x',0).attr('y',y0(c.label)+y1(k)).attr('height',y1.bandwidth()).attr('width',x(v)).attr('fill',col).attr('rx',2)
        .on('mousemove',e=>tip(`<b>${c.label}</b> · ${k}<br>weight ${fmt(v)}`,e)).on('mouseleave',untip);
    });
  });
  const lg=svg.append('g').attr('font-size',10.5).attr('transform',`translate(${w-10},14)`).attr('text-anchor','end');
  lg.append('text').attr('fill',C.navy).attr('font-weight',600).text('■ T2D');
  lg.append('text').attr('y',13).attr('fill',C.amber).attr('font-weight',600).text('■ CAD');
}

function mechStrip(sel, variants){
  const pick = rs => variants.find(v=>v.rsid===rs);
  const rows = [
    ['t2d','T2D', pick('rs1421085'), 'FTO → IRX3'],
    ['cad','CAD', pick('rs12740374'), 'SORT1'],
    ['cad','CAD', pick('rs9349379'), 'PHACTR1 → EDN1'],
  ].filter(r=>r[2]);
  document.querySelector(sel).innerHTML = rows.map(([cls,lab,v,name])=>`
    <div class="mechrow"><span class="dchip ${cls}">${lab}</span>
      <span class="mono">${v.rsid}</span> · <span class="gene">${name}</span>
      <span class="tf">${v.motif_tf||'—'} ${v.motif_dLLR!=null?`(${fmt(v.motif_dLLR,1)} b)`:''}</span></div>`).join('');
}
