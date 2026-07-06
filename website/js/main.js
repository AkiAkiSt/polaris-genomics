// POLARIS — orchestrator
import {C, fmt, pct, sci, loadAll, observeReveal, navScroll, debounce} from './util.js';
import * as charts from './charts.js';
import {renderExplorer} from './explorer.js';
import {renderFto} from './fto.js';
import {renderGeneralization} from './generalization.js';

let DATA = null;

async function boot(){
  navScroll(); observeReveal(); heroSky();
  try{
    DATA = await loadAll(['meta','variants','loci','reliability','calibration','mp','findings','fto','coloc']);
  }catch(e){ console.error(e); document.getElementById('heroStats').innerHTML =
      `<div style="color:#f0a">Data failed to load: ${e.message}</div>`; return; }
  renderHero(DATA.meta);
  renderStaticCharts();
  renderExplorer(DATA);
  renderFto(DATA.fto);
  renderGeneralization(DATA);
  addEventListener('resize', debounce(()=>{ renderStaticCharts(); renderGeneralization(DATA); }, 200));
}

function renderHero(meta){
  const t = meta.diseases.T2D;
  const stats = [
    [pct(t.gene_top1,0), 'target-gene top-1'],
    [meta.diseases.T2D.n_variants + meta.diseases.CAD.n_variants, 'variants triangulated'],
    [fmt(t.ece,3), 'calibration error (ECE)'],
    ['2', 'diseases · findings replicate'],
  ];
  document.getElementById('heroStats').innerHTML = stats.map(s=>
    `<div class="hstat"><div class="v">${s[0]}</div><div class="l">${s[1]}</div></div>`).join('');
}

function renderStaticCharts(){
  const {findings, meta, reliability, calibration, mp} = DATA;
  charts.densityCompare('#consChart', findings.T2D.phylop,
    {xlabel:'conservation (phyloP, higher = more conserved)', reading:`causal are more conserved · AUC ${fmt(findings.T2D.cons_auc)}`});
  charts.densityCompare('#rarityChart', findings.T2D.rarity,
    {xlabel:'rarity (−log₁₀ allele freq, higher = rarer)', flip:true, reading:`causal are more COMMON · AUC ${fmt(findings.T2D.rarity_auc)} (inverted)`});
  charts.geneBars('#geneChart', meta.diseases);
  charts.calibCurve('#calibChart', calibration);
  charts.relBars('#relChart', reliability);
  charts.mpSpectrum('#mpChart', mp);
  document.getElementById('eceChip').textContent = `ECE ${fmt(calibration.ece,3)}`;
  document.getElementById('permChip').textContent = `perm p = ${fmt(reliability.perm_p,3)}`;
}

/* hero starfield with periodic "triangulation" onto Polaris */
function heroSky(){
  const cv = document.getElementById('sky'); if(!cv) return;
  const ctx = cv.getContext('2d'); let W,H,stars,polaris,dpr;
  function size(){ dpr=Math.min(devicePixelRatio||1,2); const r=cv.getBoundingClientRect();
    W=cv.width=r.width*dpr; H=cv.height=r.height*dpr; init(); }
  function init(){
    const n=Math.min(140, Math.round(W*H/26000));
    stars=d3.range(n).map(()=>({x:Math.random()*W,y:Math.random()*H,r:(Math.random()*1.4+.3)*dpr,
      tw:Math.random()*Math.PI*2, sp:.4+Math.random()*.8, vx:(Math.random()-.5)*.04*dpr}));
    polaris={x:W*(0.6+Math.random()*0.25), y:H*(0.28+Math.random()*0.3), r:2.6*dpr};
  }
  let t=0, tri=null, triT=0;
  function frame(){
    t+=0.016; ctx.clearRect(0,0,W,H);
    // triangulation event every ~6s
    if(!tri && Math.random()<0.006){ const pick=()=>stars[Math.floor(Math.random()*stars.length)];
      tri={a:pick(),b:pick(),c:pick(),life:0}; triT=0; }
    if(tri){ tri.life+=0.016; const k=Math.sin(Math.min(tri.life/1.6,1)*Math.PI); // ease in/out
      ctx.strokeStyle=`rgba(127,216,200,${0.28*k})`; ctx.lineWidth=1*dpr;
      [tri.a,tri.b,tri.c].forEach(s=>{ctx.beginPath();ctx.moveTo(s.x,s.y);ctx.lineTo(polaris.x,polaris.y);ctx.stroke();});
      if(tri.life>1.9) tri=null; }
    stars.forEach(s=>{ s.x+=s.vx; if(s.x<0)s.x=W; if(s.x>W)s.x=0; s.tw+=0.02*s.sp;
      const a=0.35+0.4*Math.sin(s.tw);
      ctx.beginPath(); ctx.arc(s.x,s.y,s.r,0,7); ctx.fillStyle=`rgba(200,224,244,${a})`; ctx.fill(); });
    // polaris glow
    const gl=ctx.createRadialGradient(polaris.x,polaris.y,0,polaris.x,polaris.y,26*dpr);
    gl.addColorStop(0,'rgba(127,216,200,.9)'); gl.addColorStop(1,'rgba(127,216,200,0)');
    ctx.fillStyle=gl; ctx.beginPath(); ctx.arc(polaris.x,polaris.y,26*dpr,0,7); ctx.fill();
    ctx.fillStyle='#eafaf6'; ctx.beginPath(); ctx.arc(polaris.x,polaris.y,polaris.r,0,7); ctx.fill();
    requestAnimationFrame(frame);
  }
  size(); addEventListener('resize', debounce(size,200)); requestAnimationFrame(frame);
}

boot();
