// POLARIS — interactive explorer (master-detail, evidence triangulation, honesty layer)
import {C, fmt, pct, tip, untip} from './util.js';

let DATA, state = {disease:'T2D', locus:null, vid:null, q:''};
const byLocus = {};            // disease -> locus -> [variants]
const clamp01 = x => Math.max(0, Math.min(1, x));
const nz = (x,d=0) => (x==null||isNaN(x)) ? d : x;

export function renderExplorer(data){
  DATA = data;
  for(const v of data.variants){ (byLocus[v.disease] ??= {})[v.locus] ??= []; byLocus[v.disease][v.locus].push(v); }
  for(const dis in byLocus) for(const L in byLocus[dis]) byLocus[dis][L].sort((a,b)=>nz(b.POLARIS)-nz(a.POLARIS));
  const app = document.getElementById('explorer-app');
  app.innerHTML = `
    <div class="exp-top">
      <div class="toggle" id="expDis">
        <button data-d="T2D" class="on">Type 2 diabetes</button>
        <button data-d="CAD">Coronary artery disease</button>
      </div>
      <input class="exp-search" id="expSearch" placeholder="search gene or rsID…">
      <span class="count" id="expCount"></span>
    </div>
    <div class="exp">
      <div class="exp-list" id="expList"></div>
      <div class="exp-detail" id="expDetail"></div>
    </div>`;
  app.querySelector('#expDis').addEventListener('click', e=>{
    const b = e.target.closest('button'); if(!b) return;
    app.querySelectorAll('#expDis button').forEach(x=>x.classList.toggle('on', x===b));
    state.disease = b.dataset.d; state.locus = null; buildList(); autoSelect();
  });
  app.querySelector('#expSearch').addEventListener('input', e=>{ state.q = e.target.value.toLowerCase(); buildList(); });
  buildList(); autoSelect();
}

function lociFor(dis){
  return DATA.loci.filter(l=>l.disease===dis).sort((a,b)=>nz(b.top_POLARIS)-nz(a.top_POLARIS));
}
const byCount = dis => DATA.variants.filter(v=>v.disease===dis).length;

function buildList(){
  const list = document.getElementById('expList');
  let loci = lociFor(state.disease);
  if(state.q){
    loci = loci.filter(l => (l.locus||'').toLowerCase().includes(state.q)
      || (l.polaris_gene||'').toLowerCase().includes(state.q)
      || (l.known_effector||'').toLowerCase().includes(state.q)
      || (l.top_rsid||'').toLowerCase().includes(state.q));
  }
  const maxS = Math.max(...lociFor(state.disease).map(l=>nz(l.top_POLARIS)), 0.01);
  document.getElementById('expCount').textContent = `${loci.length} loci · ${byCount(state.disease)} variants`;
  list.innerHTML = loci.map(l=>{
    const hit = l.gene_hit;
    const badge = hit===true ? `<span class="badge hit">✓ ${l.polaris_gene}</span>`
      : hit===false ? `<span class="badge miss">✗ ${l.polaris_gene}</span>`
      : `<span class="badge unk">${l.polaris_gene||'—'}</span>`;
    const flag = l.distal_flag ? `<span style="color:var(--amber)">⚑ distal</span>` : '';
    return `<div class="exp-row ${l.locus===state.locus?'sel':''}" data-loc="${l.locus}">
      <div class="rmain"><div class="loc">${l.locus}</div>
        <div class="sub">${badge}${flag}<span class="mono">${l.top_rsid}</span> · ${l.n_variants} variant${l.n_variants>1?'s':''}</div></div>
      <div class="scorebar"><div class="sbar"><i style="width:${100*clamp01(nz(l.top_POLARIS)/maxS)}%"></i></div>
        <div class="sval">${fmt(l.top_POLARIS)}</div></div>
    </div>`;
  }).join('') || `<div style="padding:1.4rem;color:var(--slate)">No loci match “${state.q}”.</div>`;
  list.querySelectorAll('.exp-row').forEach(r=>r.addEventListener('click',()=>selectLocus(r.dataset.loc)));
}

function autoSelect(){
  const loci = lociFor(state.disease);
  if(!loci.length) return;
  const marquee = {T2D:'TCF7L2', CAD:'SORT1'}[state.disease];
  const pick = loci.find(l=>l.locus===marquee) || loci[0];
  selectLocus(pick.locus);
}

function selectLocus(locus){
  state.locus = locus;
  document.querySelectorAll('#expList .exp-row').forEach(r=>r.classList.toggle('sel', r.dataset.loc===locus));
  const vs = (byLocus[state.disease][locus]||[]);
  state.vid = vs[0]?.vid;
  renderDetail();
}

function factorStrength(v){
  return [nz(v.pip)>0.5, nz(v.P_functional)>0.5, nz(v.l2g_score)>0.5].filter(Boolean).length;
}
function confidence(v){
  if(v.distal_flag) return ['flagged · distal','warn'];
  if(nz(v.pip)>0.9 && nz(v.l2g_score)>0.5) return ['high','ok'];
  if(nz(v.pip)>0.5) return ['moderate','info'];
  return ['low · uncertain','warn'];
}

function renderDetail(){
  const vs = byLocus[state.disease][state.locus] || [];
  const v = vs.find(x=>x.vid===state.vid) || vs[0];
  const d = document.getElementById('expDetail');
  if(!v){ d.innerHTML = '<p class="muted">Select a locus.</p>'; return; }
  const [conf, confc] = confidence(v);
  const known = v.known_effector, gene = v.polaris_gene;
  const hit = (known && gene) ? known===gene : null;
  const gcell = (who, g, cls) => `<div class="gcell ${cls}"><div class="who">${who}</div><div class="g">${g||'—'}</div></div>`;
  const g3 = `<div class="genes3">
    ${gcell('nearest gene', v.nearest_gene, v.nearest_gene && known && v.nearest_gene!==known ? 'disagree':'')}
    ${gcell('POLARIS (via L2G)', gene, hit===true?'agree':hit===false?'disagree':'')}
    ${gcell('known effector', known||'not established', known?'truth':'')}
  </div>`;
  const bars = [
    ['Fine-mapping (PIP)', nz(v.pip), C.navy, `posterior this is THE causal variant = ${fmt(v.pip)}`, null],
    ['Molecular function', nz(v.P_functional), C.teal, `relative molecular-evidence score = ${fmt(v.P_functional)}`, null],
    ['Conservation (phyloP)', clamp01((nz(v.phylop_max)+1)/7), C.teal, `phyloP = ${fmt(v.phylop_max)} (mammalian constraint)`, fmt(v.phylop_max)],
    ['TF-motif disruption', clamp01(Math.abs(nz(v.motif_dLLR))/15), C.teal,
      v.motif_tf ? `${v.motif_dLLR<0?'abolishes':'alters'} a ${v.motif_tf} motif, ΔLLR ${fmt(v.motif_dLLR,1)} bits` : 'no significant motif effect',
      v.motif_tf ? `${fmt(v.motif_dLLR,1)}b` : 'n/a'],
    ['Regulatory element', v.in_ccre?0.85:0.06, C.teal, v.ccre?`in ENCODE cCRE: ${v.ccre}`:'not in an annotated cCRE', v.ccre?'yes':'no'],
    ['Gene linkage (L2G)', nz(v.l2g_score), C.violet, `Open Targets locus-to-gene score = ${fmt(v.l2g_score)}`, null],
  ];
  const barHtml = bars.map(([name,val,col,tipTxt,label])=>{
    const na = (name.includes('motif') && !v.motif_tf);
    const shown = label!=null ? label : pct(val,0);
    return `<div class="evbar ${na?'na':''}" data-tip="${(tipTxt||'').replace(/"/g,'&quot;')}">
      <div class="name">${name}</div>
      <div class="track"><i style="width:${100*clamp01(val)}%;background:${col}"></i></div>
      <div class="num">${shown}</div></div>`;
  }).join('');
  const flags = [];
  if(hit===true) flags.push(['ok','✓ matches the known effector gene']);
  if(hit===false) flags.push(['err','✗ disagrees with the known effector ('+known+')']);
  if(v.distal_flag) flags.push(['warn','⚑ distal-regulation signal — the proximal gene call is untrustworthy here']);
  if(v.coding) flags.push(['info','coding lead variant (kept as contrast)']);
  if(nz(v.pip)<0.5) flags.push(['warn','low fine-mapping certainty — many variants share the signal']);
  flags.push(['info','hypothesis — a molecular effect is not proof of disease causation']);
  const conv = factorStrength(v);

  d.innerHTML = `
    <div class="dt-head">
      <div><div class="id">${v.rsid||v.vid}</div>
        <div class="pos">${v.locus} locus · chr${v.chrom}:${(v.pos&&v.pos.toLocaleString)?v.pos.toLocaleString():v.pos} · ${v.ref}→${v.alt}</div></div>
      <div class="dt-gauge"><div class="v" style="color:${confc==='ok'?C.grass:confc==='warn'?C.amber:C.navy}">${fmt(v.POLARIS)}</div>
        <div class="l">POLARIS score</div>
        <div class="chip ${confc==='ok'?'teal':confc==='warn'?'amber':'slate'}" style="margin-top:.3rem">${conf}</div></div>
    </div>
    ${vs.length>1 ? `<div class="csset" id="csSet">${vs.slice(0,14).map(x=>
      `<span class="cv ${x.vid===v.vid?'on':''}" data-vid="${x.vid}" title="PIP ${fmt(x.pip)}">${x.rsid||x.vid.split('_').slice(0,2).join(':')}</span>`).join('')}
      ${vs.length>14?`<span class="cv" style="border:0;cursor:default">+${vs.length-14} more</span>`:''}</div>
      <div style="font-size:.72rem;color:var(--mute);margin-top:.2rem">credible set — the variants that plausibly carry the signal (click to compare)</div>`:''}

    <div class="dt-sec"><div class="h">Which gene? — all three answers, disagreements shown</div>${g3}</div>

    <div class="dt-sec"><div class="h">Evidence — ${conv}/3 factors converge
      <span style="float:right;font-weight:400;text-transform:none;letter-spacing:0">
        <span style="color:${C.navy}">■</span> statistical · <span style="color:${C.teal}">■</span> molecular · <span style="color:${C.violet}">■</span> linkage</span></div>
      ${barHtml}</div>

    <div class="dt-sec"><div class="h">Transparent mechanism</div>
      <div class="mechbox">${mechHtml(v)}</div></div>

    <div class="flags">${flags.map(f=>`<span class="flag ${f[0]}">${f[1]}</span>`).join('')}</div>
    <div class="dt-note">POLARIS ranks hypotheses to prioritize experiments; it does not diagnose. A variant reaches the top of its locus by converging across independent evidence — not from any single score.</div>`;

  d.querySelectorAll('#csSet .cv[data-vid]').forEach(el=>el.addEventListener('click',()=>{ state.vid = el.dataset.vid; renderDetail(); }));
  d.querySelectorAll('.evbar').forEach(el=>{
    el.addEventListener('mousemove',e=>tip(el.dataset.tip,e)); el.addEventListener('mouseleave',untip);
    const i = el.querySelector('i'); const w=i.style.width; i.style.width='0'; requestAnimationFrame(()=>setTimeout(()=>i.style.width=w,30));
  });
}

function mechHtml(v){
  if(!v.motif_tf) return `No transparent TF-motif mechanism detected. Statistical and linkage evidence still nominate <span class="mono">${v.polaris_gene||'—'}</span>; the regulatory mechanism is not captured by current motif annotations.`;
  const verb = v.motif_dLLR < -2 ? 'abolishes' : 'weakens';
  const cc = v.ccre ? ` within a <span class="mono">${v.ccre}</span> element` : '';
  const p = v.motif_pref!=null ? (+v.motif_pref).toExponential(0) : '—';
  return `The <span class="mono">${v.alt}</span> allele ${verb} a <span class="mono">${v.motif_tf}</span> transcription-factor binding site (ΔLLR <b>${fmt(v.motif_dLLR,1)} bits</b>, p<sub>ref</sub> ${p})${cc}, plausibly changing regulation of <span class="mono">${v.polaris_gene||'its target'}</span>.`;
}
