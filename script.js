'use strict';

// ─── PARAMETERS (MATLAB convention: Vrest=0) ─────────────
const P={gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:6.3,blockNa:false,blockK:false};
// Default: continuous Iext=10, 120 ms window
const STIM={type:'continuous',I0:10,tStart:0,dur:100,freq:40};
let DT=0.01;
const STORE_EVERY=2, MAX_HIST=15000;
const SPEED_LEVELS=[10,25,50,100,200,350,500,800];
let speedIdx=2, winMs=120;
let viewStart=0, viewEnd=winMs; // CC plot view state

// ─── SIM STATE ────────────────────────────────────────────
let t=0, y=[0,0,0,0], stepCnt=0;
let running=false, animId=null;
let viewLocked=false;
let currentV=0;

// History arrays (rolling: capped at MAX_HIST)
const hist={t:[],V:[],m:[],h:[],n:[],INa:[],IK:[],IL:[],Iext:[]};

// AP detection
let apSpikes=0, apLastSpike=-1e9, apPeakV=-Infinity, apSpikeTs=[];
const AP_THRESH=-20; // mV — absolute scale, crosses on AP upstroke

// Phase portrait
const TRAIL_LEN=6000;
const phaseV=new Float32Array(TRAIL_LEN), phaseN=new Float32Array(TRAIL_LEN);
let phasePtr=0, phaseCount=0;

// Auto-pause: triggered 20 ms after voltage becomes negligible post-stimulus
const STAB_THRESH=1.0;  // mV  — negligible change threshold
const STAB_DUR=50;      // ms  — must be stable this long
let stabTimer=0, postStimStart=-1;

// ─── RATE FUNCTIONS — absolute voltage scale (Vrest = −65 mV) ──
// Internally shifts V to HH convention: V' = V + 65
function am(V){const x=(V+65)-25;return Math.abs(x)<1e-7?1.0:-0.1*x/(Math.exp(-x/10)-1);}
function bm(V){return 4*Math.exp(-(V+65)/18);}
function ah(V){return 0.07*Math.exp(-(V+65)/20);}
function bh(V){return 1/(Math.exp(-((V+65)-30)/10)+1);}
function an(V){const x=(V+65)-10;return Math.abs(x)<1e-7?0.1:-0.01*x/(Math.exp(-x/10)-1);}
function bn(V){return 0.125*Math.exp(-(V+65)/80);}
function phi(){return Math.pow(3,(P.T-6.3)/10);}
function ss(V){return[am(V)/(am(V)+bm(V)),ah(V)/(ah(V)+bh(V)),an(V)/(an(V)+bn(V))];}

// ─── STIMULUS ─────────────────────────────────────────────
function iext(t){
  const{type,I0,tStart,dur,freq}=STIM, tEnd=tStart+dur;
  if(t<tStart) return 0;
  if(type==='continuous') return I0;                          // no end — always on
  if(type==='step')  return t<=tEnd ? I0 : 0;
  if(type==='ramp')  return t>tEnd  ? 0  : I0*(t-tStart)/dur;
  if(type==='pulse'){
    if(t>tEnd) return 0;
    const p=1000/freq, pw=Math.min(1,p*0.4);
    return (t-tStart)%p < pw ? I0 : 0;
  }
  if(type==='sine'){
    if(t>tEnd) return 0;
    return I0*Math.sin(2*Math.PI*freq*(t-tStart)/1000);
  }
  return 0;
}

// ─── ODE ─────────────────────────────────────────────────
function dydt(y,t){
  const[V,m,h,n]=y, ph=phi();
  const INa=P.blockNa?0:P.gNa*m*m*m*h*(V-P.ENa);
  const IK =P.blockK ?0:P.gK *n*n*n*n*(V-P.EK);
  const IL =P.gL*(V-P.EL), Ie=iext(t);
  return[(Ie-INa-IK-IL)/P.Cm, ph*(am(V)*(1-m)-bm(V)*m), ph*(ah(V)*(1-h)-bh(V)*h), ph*(an(V)*(1-n)-bn(V)*n)];
}

// ─── RK4 ─────────────────────────────────────────────────
function rk4(){
  const dt=DT;
  const k1=dydt(y,t);
  const y2=y.map((v,i)=>v+0.5*dt*k1[i]);
  const k2=dydt(y2,t+0.5*dt);
  const y3=y.map((v,i)=>v+0.5*dt*k2[i]);
  const k3=dydt(y3,t+0.5*dt);
  const y4=y.map((v,i)=>v+dt*k3[i]);
  const k4=dydt(y4,t+dt);
  y=y.map((v,i)=>v+dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]));
  t+=dt; stepCnt++;
}

// ─── STORE & DETECT ──────────────────────────────────────
function storeAndDetect(){
  const[V,m,h,n]=y, Ie=iext(t);
  const INa=P.blockNa?0:P.gNa*m*m*m*h*(V-P.ENa);
  const IK =P.blockK ?0:P.gK *n*n*n*n*(V-P.EK);
  const IL =P.gL*(V-P.EL);
  hist.t.push(t);hist.V.push(V);hist.m.push(m);hist.h.push(h);
  hist.n.push(n);hist.INa.push(INa);hist.IK.push(IK);hist.IL.push(IL);hist.Iext.push(Ie);
  // Periodic memory trim: every 500 stored points, discard old data
  // For continuous mode keep 2 full display cycles; otherwise keep MAX_HIST
  if(hist.t.length % 500 === 0 && hist.t.length > MAX_HIST){
    // Keep last 2 display windows (ECG needs current + previous cycle)
    const keepMs = STIM.type==='continuous' ? winMs*2+20 : winMs*3;
    const keepPts = Math.ceil(keepMs / (DT*STORE_EVERY));
    const drop = Math.max(0, hist.t.length - keepPts);
    if(drop > 0){
      for(const k of Object.keys(hist))hist[k].splice(0,drop);
    }
  }
  phaseV[phasePtr]=V;phaseN[phasePtr]=n;phasePtr=(phasePtr+1)%TRAIL_LEN;
  if(phaseCount<TRAIL_LEN)phaseCount++;

  // Spike detection
  const prev=hist.V.length>=2?hist.V[hist.V.length-2]:V;
  if(prev<AP_THRESH&&V>=AP_THRESH&&t-apLastSpike>3){
    apSpikes++;if(apSpikes>1)apSpikeTs.push(t-apLastSpike);
    apLastSpike=t;apPeakV=V;
  }
  if(V>apPeakV&&apSpikes>0)apPeakV=V;

  // Stability / auto-pause: only fires when stimulus is NOT continuous
  // AND after the stimulus window has ended AND V has settled
  const isContinuous = STIM.type==='continuous';
  if(!isContinuous){
    const stimDone = t > STIM.tStart + STIM.dur + 2;
    if(stimDone){
      if(postStimStart<0) postStimStart=t;
      if(Math.abs(V+65)<STAB_THRESH) stabTimer+=DT; else stabTimer=0;
    }
  }
}

function shouldAutoPause(){return P.autoPause && STIM.type!=='continuous' && stabTimer>=STAB_DUR;}

// ─── INIT ─────────────────────────────────────────────────
function initState(){
  y=[-65,0,0,0]; t=0; stepCnt=0;
  for(const k of Object.keys(hist))hist[k]=[];
  apSpikes=0;apLastSpike=-1e9;apPeakV=-Infinity;apSpikeTs=[];
  phasePtr=0;phaseCount=0;stabTimer=0;postStimStart=-1;
}

// ─── SIMULATION LOOP ──────────────────────────────────────
function animate(){
  if(!running)return;
  const steps=SPEED_LEVELS[speedIdx];
  for(let i=0;i<steps;i++){
    rk4();
    if(stepCnt%STORE_EVERY===0)storeAndDetect();
    if(shouldAutoPause()){
      // Lock view to show full trace from t=0
      viewLocked=true;
      stopSim();
      const sd=document.getElementById('sim-dot');if(sd)sd.classList.add('paused');
      const st=document.getElementById('sim-status-text');if(st)st.textContent='Stable — auto-paused';
      renderAll();updateStats();return;
    }
  }
  renderAll();updateStats();
  animId=requestAnimationFrame(animate);
}

function startSim(){
  if(running)return;
  stabTimer=0;postStimStart=-1;  // clear stale pause state
  viewLocked=false;
  running=true;
  const sd=document.getElementById('sim-dot');if(sd)sd.classList.remove('paused');
  const st=document.getElementById('sim-status-text');if(st)st.textContent='Running';
  animId=requestAnimationFrame(animate);
}

function stopSim(){
  running=false;
  if(animId){cancelAnimationFrame(animId);animId=null;}
  const sd=document.getElementById('sim-dot');if(sd)sd.classList.add('paused');
  const st=document.getElementById('sim-status-text');
  if(st&&st.textContent==='Running')st.textContent='Paused';
}

function resetSim(){
  stopSim();initState();
  viewLocked=false;
  renderAll();updateStats();
  const st=document.getElementById('sim-status-text');if(st)st.textContent='Ready';
}

// ─── CANVAS SETUP ─────────────────────────────────────────
const PAD={top:18,right:16,bottom:40,left:58};

function setupCanvas(id){
  const c=document.getElementById(id);
  if(!c)return{ctx:null,W:0,H:0};
  if(!c.dataset.lh)c.dataset.lh=c.getAttribute('height')||'160';
  const minH=parseInt(c.dataset.lh);
  const W=Math.max(100,c.parentElement.clientWidth||600);
  const H=Math.max(minH,c.clientHeight||minH);
  const dpr=window.devicePixelRatio||1;
  c.width=Math.round(W*dpr);c.height=Math.round(H*dpr);
  c.style.width=W+'px';c.style.height=H+'px';
  const ctx=c.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);
  return{ctx,W,H};
}

// ─── DRAW AXES ────────────────────────────────────────────
function drawAxes(ctx,W,H,yMin,yMax,xMin,xMax,yLabel){
  if(!ctx)return;
  const p=PAD,pw=W-p.left-p.right,ph=H-p.top-p.bottom;
  ctx.clearRect(0,0,W,H);
  ctx.fillStyle='#f9fafb';ctx.fillRect(0,0,W,H);
  ctx.fillStyle='#fff';ctx.fillRect(p.left,p.top,pw,ph);
  ctx.strokeStyle='#f3f4f6';ctx.lineWidth=1;
  for(let i=0;i<=5;i++){ctx.beginPath();ctx.moveTo(p.left,p.top+ph*i/5);ctx.lineTo(p.left+pw,p.top+ph*i/5);ctx.stroke();}
  for(let i=0;i<=6;i++){ctx.beginPath();ctx.moveTo(p.left+pw*i/6,p.top);ctx.lineTo(p.left+pw*i/6,p.top+ph);ctx.stroke();}
  ctx.strokeStyle='#d1d5db';ctx.lineWidth=1;ctx.strokeRect(p.left,p.top,pw,ph);
  if(yMin<0&&yMax>0){
    const y0=p.top+ph*(1-(0-yMin)/(yMax-yMin));
    ctx.strokeStyle='rgba(0,0,0,.1)';ctx.lineWidth=1;ctx.setLineDash([3,3]);
    ctx.beginPath();ctx.moveTo(p.left,y0);ctx.lineTo(p.left+pw,y0);ctx.stroke();ctx.setLineDash([]);
  }
  const font="bold 12px 'DM Mono',monospace";
  ctx.fillStyle='#4a5568';ctx.font=font;ctx.textAlign='right';
  for(let i=0;i<=5;i++){
    const val=yMin+(yMax-yMin)*(1-i/5);
    ctx.fillText(val.toFixed(Math.abs(val)>=10||val===0?0:1),p.left-4,p.top+ph*i/5+4);
  }
  ctx.textAlign='center';
  for(let i=0;i<=6;i++){
    // X axis shows time within the current window (0..winMs)
    ctx.fillText(Math.round(xMin+(xMax-xMin)*i/6),p.left+pw*i/6,H-8);
  }
  ctx.save();ctx.translate(10,p.top+ph/2);ctx.rotate(-Math.PI/2);
  ctx.textAlign='center';ctx.fillStyle='#718096';ctx.font="bold 12px 'DM Sans',sans-serif";
  ctx.fillText(yLabel,0,0);ctx.restore();
}

// ─── ECG-STYLE PLOT (heart rate monitor — wraps at winMs) ──────────────
// Maps time data to X position within a fixed [0, winMs] window.
// The current "write head" (t % winMs) sweeps right; old data shows faded.
function plotECG(ctx,ts,vs,W,H,vMin,vMax,color,lw){
  if(!ctx||!ts||ts.length<2)return;
  const p=PAD,pw=W-p.left-p.right,ph=H-p.top-p.bottom;
  // Current cycle base
  const tCur=ts[ts.length-1];
  const cycleBase=Math.floor(tCur/winMs)*winMs; // start of current page
  const headX=p.left+pw*((tCur%winMs)/winMs);  // current write-head X

  // ── clip to plot area ──
  ctx.save();ctx.beginPath();ctx.rect(p.left,p.top,pw,ph);ctx.clip();

  // Draw the "eraser" strip just ahead of the write head
  const eraseW=pw*0.04; // ~4% of width
  ctx.fillStyle='#fff';
  ctx.fillRect(headX,p.top,eraseW,ph);

  // Helper: draw a segment of time data with a given alpha
  function drawSegment(fromT,toT,alpha){
    ctx.strokeStyle=color;ctx.lineWidth=lw||2;ctx.globalAlpha=alpha;
    ctx.beginPath();let first=true;
    for(let i=0;i<ts.length;i++){
      if(ts[i]<fromT-DT||ts[i]>toT+DT) continue;
      // x within window: ((ts[i] - cycleBase) % winMs) / winMs * pw
      // but for previous cycle data: offset by -winMs
      let tOffset=(ts[i]-cycleBase);
      if(tOffset<0)tOffset+=winMs; // wrap previous cycle
      const x=p.left+pw*(tOffset/winMs);
      const yc=p.top+ph*(1-(vs[i]-vMin)/(vMax-vMin));
      const ycc=Math.max(p.top-2,Math.min(p.top+ph+2,yc));
      if(first||Math.abs(x-(p.left+pw*((ts[i-1]-cycleBase<0?ts[i-1]-cycleBase+winMs:ts[i-1]-cycleBase)/winMs)))>pw*0.5){
        ctx.moveTo(x,ycc);first=false;
      }else{
        ctx.lineTo(x,ycc);
      }
    }
    ctx.stroke();ctx.globalAlpha=1;
  }

  // Previous cycle (faded) — data from [cycleBase-winMs, cycleBase]
  if(tCur>=winMs) drawSegment(cycleBase-winMs, cycleBase, 0.28);
  // Current cycle (full brightness) — data from [cycleBase, tCur]
  drawSegment(cycleBase, tCur, 1.0);

  // Write-head cursor line
  ctx.strokeStyle='rgba(11,122,110,0.6)';ctx.lineWidth=1.5;ctx.setLineDash([3,3]);
  ctx.beginPath();ctx.moveTo(headX,p.top);ctx.lineTo(headX,p.top+ph);ctx.stroke();
  ctx.setLineDash([]);

  ctx.restore();
}

// ─── RENDER PLOTS ─────────────────────────────────────────
function renderAll(){
  if(!hist.t.length)return;
  currentV=hist.V[hist.V.length-1];
  renderV();renderGates();renderCurrents();renderIext_plot();
  renderSS();renderTau();renderPhase();
}

function renderV(){
  const{ctx,W,H}=setupCanvas('cv-V');if(!ctx)return;
  drawAxes(ctx,W,H,-85,60,0,winMs,'mV');
  plotECG(ctx,hist.t,hist.V,W,H,-85,60,'#0b7a6e',2);
}

function renderGates(){
  const{ctx,W,H}=setupCanvas('cv-gates');if(!ctx)return;
  drawAxes(ctx,W,H,0,1,0,winMs,'');
  plotECG(ctx,hist.t,hist.m,W,H,0,1,'#dc2626',1.5);
  plotECG(ctx,hist.t,hist.h,W,H,0,1,'#f59e0b',1.5);
  plotECG(ctx,hist.t,hist.n,W,H,0,1,'#7c3aed',1.5);
}

function renderCurrents(){
  const{ctx,W,H}=setupCanvas('cv-currents');if(!ctx)return;
  const tCur=hist.t[hist.t.length-1];
  const cycleBase=Math.floor(tCur/winMs)*winMs;
  let mxI=5;
  for(let i=0;i<hist.t.length;i++){
    if(hist.t[i]<cycleBase-winMs-0.5)continue;
    const mx=Math.max(Math.abs(hist.INa[i]),Math.abs(hist.IK[i]),Math.abs(hist.IL[i]));
    if(mx>mxI)mxI=mx;
  }
  mxI=Math.max(5,mxI*1.12);
  drawAxes(ctx,W,H,-mxI,mxI,0,winMs,'µA/cm²');
  plotECG(ctx,hist.t,hist.INa,W,H,-mxI,mxI,'#dc2626',1.5);
  plotECG(ctx,hist.t,hist.IK, W,H,-mxI,mxI,'#1d4ed8',1.5);
  plotECG(ctx,hist.t,hist.IL, W,H,-mxI,mxI,'#15803d',1);
}

function renderIext_plot(){
  const{ctx,W,H}=setupCanvas('cv-Iext');if(!ctx)return;
  const maxI=Math.max(Math.abs(STIM.I0),1);
  const yMin=STIM.I0<0?-maxI*1.2:-maxI*0.15, yMax=maxI*1.2;
  drawAxes(ctx,W,H,yMin,yMax,0,winMs,'µA/cm²');
  if(hist.t.length>1){
    plotECG(ctx,hist.t,hist.Iext,W,H,yMin,yMax,'#0b7a6e',2);
  } else {
    // Preview
    const p=PAD,pw=W-p.left-p.right,ph=H-p.top-p.bottom;
    const step=0.5, pts=[],is=[];
    for(let tt=0;tt<=winMs;tt+=step){pts.push(tt);is.push(iext(tt));}
    // Draw as normal trace
    ctx.save();ctx.strokeStyle='rgba(11,122,110,.5)';ctx.lineWidth=1.5;ctx.setLineDash([4,3]);
    ctx.beginPath();let first=true;
    for(let i=0;i<pts.length;i++){
      const x=p.left+pw*(pts[i]/winMs);
      const yc=p.top+ph*(1-(is[i]-yMin)/(yMax-yMin));
      first?ctx.moveTo(x,yc):ctx.lineTo(x,yc);first=false;
    }
    ctx.stroke();ctx.setLineDash([]);ctx.restore();
  }
}

function renderSS(){
  const{ctx,W,H}=setupCanvas('cv-ss');if(!ctx)return;
  const Vs=[];for(let V=-20;V<=130;V+=0.5)Vs.push(V);
  const p=PAD,pw=W-p.left-p.right,ph=H-p.top-p.bottom;
  drawAxes(ctx,W,H,0,1,-85,65,'x∞');
  const cols=['#dc2626','#f59e0b','#7c3aed'],lbls=['m∞','h∞','n∞'];
  for(let c=0;c<3;c++){
    ctx.beginPath();ctx.strokeStyle=cols[c];ctx.lineWidth=2;
    Vs.forEach((V,i)=>{const s=ss(V);const x=p.left+pw*(V+20)/150;const y=p.top+ph*(1-s[c]);i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);});
    ctx.stroke();
  }
  const xV=p.left+pw*(currentV+20)/150;
  ctx.strokeStyle='rgba(11,122,110,.5)';ctx.lineWidth=1;ctx.setLineDash([3,3]);
  ctx.beginPath();ctx.moveTo(xV,p.top);ctx.lineTo(xV,p.top+ph);ctx.stroke();ctx.setLineDash([]);
  [[30,0],[-10,1],[10,2]].forEach(([V,c])=>{const s=ss(V);ctx.fillStyle=cols[c];ctx.font="bold 11px 'DM Mono',monospace";ctx.textAlign='left';ctx.fillText(lbls[c],p.left+pw*(V+20)/150+3,p.top+ph*(1-s[c])-3);});
}

function renderTau(){
  const{ctx,W,H}=setupCanvas('cv-tau');if(!ctx)return;
  const ph2=phi();
  const Vs=[];for(let V=-20;V<=130;V+=0.5)Vs.push(V);
  const taus=Vs.map(V=>[1/((am(V)+bm(V))*ph2),1/((ah(V)+bh(V))*ph2),1/((an(V)+bn(V))*ph2)]);
  const mxT=Math.max(...taus.flat())*1.1;
  const p=PAD,pw=W-p.left-p.right,pph=H-p.top-p.bottom;
  drawAxes(ctx,W,H,0,mxT,-85,65,'τ (ms)');
  const cols=['#dc2626','#f59e0b','#7c3aed'],lbls=['τm','τh','τn'];
  for(let c=0;c<3;c++){
    ctx.beginPath();ctx.strokeStyle=cols[c];ctx.lineWidth=2;
    Vs.forEach((V,i)=>{const x=p.left+pw*(V+20)/150;const y=p.top+pph*(1-taus[i][c]/mxT);i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);});
    ctx.stroke();
  }
  [[35,0],[-10,1],[15,2]].forEach(([V,c])=>{const idx=Math.round((V+20)/0.5);if(idx<0||idx>=taus.length)return;const τ=taus[idx][c];ctx.fillStyle=cols[c];ctx.font="bold 11px 'DM Mono',monospace";ctx.textAlign='left';ctx.fillText(lbls[c],p.left+pw*(V+20)/150+3,p.top+pph*(1-τ/mxT)-3);});
  const xV=p.left+pw*(currentV+20)/150;ctx.strokeStyle='rgba(11,122,110,.5)';ctx.lineWidth=1;ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(xV,p.top);ctx.lineTo(xV,p.top+pph);ctx.stroke();ctx.setLineDash([]);
}

function renderPhase(){
  const{ctx,W,H}=setupCanvas('cv-phase');if(!ctx)return;
  const p={top:20,right:20,bottom:36,left:54},pw=W-p.left-p.right,ph=H-p.top-p.bottom;
  ctx.clearRect(0,0,W,H);ctx.fillStyle='#f9fafb';ctx.fillRect(0,0,W,H);ctx.fillStyle='#fff';ctx.fillRect(p.left,p.top,pw,ph);
  ctx.strokeStyle='#f3f4f6';ctx.lineWidth=1;
  for(let i=0;i<=5;i++){ctx.beginPath();ctx.moveTo(p.left,p.top+ph*i/5);ctx.lineTo(p.left+pw,p.top+ph*i/5);ctx.stroke();ctx.beginPath();ctx.moveTo(p.left+pw*i/5,p.top);ctx.lineTo(p.left+pw*i/5,p.top+ph);ctx.stroke();}
  ctx.strokeStyle='#d1d5db';ctx.lineWidth=1;ctx.strokeRect(p.left,p.top,pw,ph);
  const vMin=-85,vMax=65,nMin=0,nMax=1;
  const xp=(n)=>p.left+pw*(n-nMin)/(nMax-nMin), yp=(v)=>p.top+ph*(1-(v-vMin)/(vMax-vMin));

  // ── Draw trajectory FIRST so annotations layer on top ────────────
  const cnt=Math.min(phaseCount,TRAIL_LEN);
  if(cnt>2){
    for(let i=1;i<cnt;i++){
      const a=(phasePtr-cnt+i-1+TRAIL_LEN)%TRAIL_LEN,b=(phasePtr-cnt+i+TRAIL_LEN)%TRAIL_LEN;
      ctx.strokeStyle=`rgba(11,122,110,${Math.pow(i/cnt,0.6)*0.9})`;ctx.lineWidth=1.5;
      ctx.beginPath();ctx.moveTo(xp(phaseN[a]),yp(phaseV[a]));ctx.lineTo(xp(phaseN[b]),yp(phaseV[b]));ctx.stroke();
    }
    const ci=(phasePtr-1+TRAIL_LEN)%TRAIL_LEN;
    // Current position dot — pulsing larger dot at simulation tip
    ctx.fillStyle='rgba(11,122,110,0.9)';
    ctx.beginPath();ctx.arc(xp(phaseN[ci]),yp(phaseV[ci]),6,0,Math.PI*2);ctx.fill();
    ctx.strokeStyle='#fff';ctx.lineWidth=1.5;ctx.stroke();
  } else {
    ctx.fillStyle='#9ca3af';ctx.font="12px 'DM Sans',sans-serif";ctx.textAlign='center';
    ctx.fillText('Run the simulation to see the phase portrait',p.left+pw/2,p.top+ph/2-8);
    ctx.fillText('or click "Phase Demo" preset',p.left+pw/2,p.top+ph/2+12);
  }

  // ── Compute annotation positions from live simulation data ────────
  // Rest = fixed point of HH system (analytical, never changes)
  const restV=-65, restN=0.317; // n∞(-65 mV) = 0.317
  let threshV=-55, threshN=0.34;
  let peakV=20, peakN=0.42;

  if(hist.V && hist.V.length>10){

    // Peak: from spike detection
    if(apPeakV>-Infinity && apPeakV>-40){
      peakV=apPeakV;
      // Find n at the time of peak V
      if(apLastSpike>0 && hist.t.length>0){
        let bestIdx=0, bestDist=Infinity;
        for(let i=0;i<hist.t.length;i++){
          const d=Math.abs(hist.t[i]-apLastSpike);
          if(d<bestDist){ bestDist=d; bestIdx=i; }
        }
        // Find actual peak index near apLastSpike (within ±2ms)
        let pkIdx=bestIdx;
        for(let i=Math.max(0,bestIdx-200);i<Math.min(hist.t.length,bestIdx+200);i++){
          if(hist.V[i]>hist.V[pkIdx]) pkIdx=i;
        }
        peakN=hist.n[pkIdx];
        peakV=hist.V[pkIdx];
      }
    }

    // Threshold: find point with maximum dV/dt on most recent upstroke
    if(hist.V.length>5){
      let maxDV=0, threshIdx=-1;
      for(let i=1;i<hist.V.length;i++){
        const dv=hist.V[i]-hist.V[i-1];
        if(dv>maxDV && hist.V[i]>-65 && hist.V[i]<0){
          maxDV=dv; threshIdx=i;
        }
      }
      if(threshIdx>=0){ threshV=hist.V[threshIdx]; threshN=hist.n[threshIdx]; }
    }
  }

  // ── Draw annotations ON TOP of trajectory ─────────────────────────
  const annotations=[
    {V:restV,  n:restN,   label:'Rest',      align:'left',  color:'#22c55e'},
    {V:threshV,n:threshN, label:'Threshold', align:'left',  color:'#f59e0b'},
    {V:peakV,  n:peakN,   label:'AP Peak',   align:'center', color:'#ef4444'},
  ];
  ctx.font="bold 11px 'DM Sans',sans-serif";
  for(const{V,n,label,align,color} of annotations){
    const ax=xp(n), ay=yp(V);
    if(ax<p.left-5||ax>p.left+pw+5||ay<p.top-5||ay>p.top+ph+5) continue;
    // Marker: square for Rest (fixed point), circle for Threshold & Peak
    ctx.fillStyle=color;
    if(label==='Rest'){
      ctx.fillRect(ax-4,ay-4,8,8);
      ctx.strokeStyle='#fff';ctx.lineWidth=1.5;ctx.strokeRect(ax-4,ay-4,8,8);
    } else {
      ctx.beginPath();ctx.arc(ax,ay,4,0,Math.PI*2);ctx.fill();
      ctx.strokeStyle='#fff';ctx.lineWidth=1.5;ctx.stroke();
    }
    // Label
    ctx.fillStyle=color;
    ctx.textAlign=align||'center';
    const ox=align==='left'?8:align==='right'?-8:0;
    ctx.fillText(label,ax+ox,ay-7);
  }

  // ── Axis labels ───────────────────────────────────────────────────
  ctx.fillStyle='#6b7280';ctx.font="bold 12px 'DM Mono',monospace";ctx.textAlign='right';
  for(let i=0;i<=5;i++)ctx.fillText((vMin+(vMax-vMin)*(1-i/5)).toFixed(0),p.left-4,p.top+ph*i/5+4);
  ctx.textAlign='center';ctx.font="bold 12px 'DM Mono',monospace";
  for(let i=0;i<=5;i++)ctx.fillText((nMin+(nMax-nMin)*i/5).toFixed(2),p.left+pw*i/5,H-8);
  ctx.save();ctx.translate(12,p.top+ph/2);ctx.rotate(-Math.PI/2);
  ctx.textAlign='center';ctx.fillStyle='#4a5568';ctx.font="bold 12px 'DM Sans',sans-serif";
  ctx.fillText('V (mV)',0,0);ctx.restore();
  ctx.fillStyle='#4a5568';ctx.font="bold 12px 'DM Sans',sans-serif";
  ctx.textAlign='center';ctx.fillText('n (K⁺ activation)',p.left+pw/2,H-2);
}

// ─── STATS ────────────────────────────────────────────────
function updateStats(){
  const el=(id,v)=>{const e=document.getElementById(id);if(e)e.textContent=v;};
  el('stat-spikes',apSpikes);
  el('stat-peak',apSpikes>0?apPeakV.toFixed(1):'—');
  const meanISI=apSpikeTs.length>0?apSpikeTs.reduce((a,b)=>a+b,0)/apSpikeTs.length:0;
  el('stat-freq',meanISI>0?(1000/meanISI).toFixed(1):'—');
  el('stat-t',t.toFixed(0));
  el('stat-phi',phi().toFixed(2)+'×');
}

// ─── UI HELPERS ───────────────────────────────────────────
function setSpeed(v){
  speedIdx=parseInt(v)-1;
  const labels=['⅛×','¼×','½×','1×','2×','4×','7×','14×'];
  const el=document.getElementById('sl-speed-val');if(el)el.textContent=labels[speedIdx]||v+'×';
}

function updateWindow(){
  // Called when winMs slider changes; refit ECG view to new window size
  if(!viewLocked){ viewStart=Math.max(0,t-winMs); viewEnd=Math.max(winMs,t); }
}

function slUpdate(param,val,valId){
  const dp=param==='EL'?1:param==='gL'?2:param==='Cm'?1:param==='T'?1:0;
  const numVal=parseFloat(val);P[param]=numVal;
  const spanEl=document.getElementById(valId);if(spanEl)spanEl.textContent=numVal.toFixed(dp);
  const numEl=document.getElementById(valId.replace('-val','-num'));if(numEl)numEl.value=numVal.toFixed(dp);
  if(param==='T'){const e=document.getElementById('stat-phi');if(e)e.textContent=phi().toFixed(2)+'×';renderSS();renderTau();}
  if(param==='win'){winMs=numVal;if(!running){renderAll();}}
}


function slUpdateDT(val){
  const v=parseFloat(val);
  if(isNaN(v)||v<=0)return;
  DT=Math.min(0.115,Math.max(0.005,v));
  const dp=3;
  const spanEl=document.getElementById('sl-dt-val');if(spanEl)spanEl.textContent=DT.toFixed(dp);
  const numEl=document.getElementById('sl-dt-num');if(numEl)numEl.value=DT.toFixed(dp);
  // Update subtitle to reflect current dt
  const sub=document.querySelector('.sim-subtitle');
  if(sub) sub.textContent=`4th-order Runge–Kutta · dt = ${DT.toFixed(3)} ms · Squid giant axon · V_rest = −65 mV`;
  // Warn if dt is large — likely to be numerically unstable
  const st=document.getElementById('sim-status-text');
  if(st && DT>0.1) st.textContent='⚠ Large dt — may be unstable';
}

function updateStim(){
  STIM.type=document.getElementById('stim-type').value;
  STIM.I0=parseFloat(document.getElementById('sl-I0').value);
  STIM.tStart=parseFloat(document.getElementById('sl-tStart').value);
  STIM.dur=parseFloat(document.getElementById('sl-dur').value);
  STIM.freq=parseFloat(document.getElementById('sl-freq').value);
  const showFreq=['pulse','sine'].includes(STIM.type);
  const frEl=document.getElementById('freq-row');
  if(frEl){ if(!showFreq) frEl.classList.add('disabled'); else frEl.classList.remove('disabled'); }
  const isCont=STIM.type==='continuous';
  const durRow=document.getElementById('dur-row');
  if(durRow){ if(isCont) durRow.classList.add('disabled'); else durRow.classList.remove('disabled'); }
  const t0=STIM.tStart,t1=STIM.tStart+STIM.dur;
  const d0=document.getElementById('stim-t0'),d1=document.getElementById('stim-t1');
  if(d0)d0.textContent=t0.toFixed(0);if(d1)d1.textContent=t1.toFixed(0);
  if(!running)renderIext_plot();
}

function clearPhase(){phasePtr=0;phaseCount=0;renderPhase();}

// ─── PRESETS ──────────────────────────────────────────────
const PRESETS={
  standard:  {gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:6.3,I0:10,tStart:0,dur:100,type:'continuous',freq:40,win:120,speed:3,autoStart:true},
  step100:   {gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:6.3,I0:10,tStart:0,dur:100,type:'step',    freq:40,win:140,speed:3},
  body:      {gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:37, I0:10,tStart:0,dur:20, type:'continuous',freq:40,win:25, speed:2},
  burst:     {gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:6.3,I0:20,tStart:0,dur:200,type:'step',    freq:40,win:220,speed:4},
  subthresh: {gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:6.3,I0:5, tStart:0,dur:100,type:'step',    freq:40,win:140,speed:3},
  phasedemo: {gNa:120,gK:36,gL:0.3,ENa:50,EK:-77,EL:-54.4,Cm:1,T:6.3,I0:12,tStart:0,dur:150,type:'continuous',freq:40,win:160,speed:6,autoStart:true,keepPhase:true},
};

function applyPreset(name){
  const pr=PRESETS[name];if(!pr)return;
  P.gNa=pr.gNa;P.gK=pr.gK;P.gL=pr.gL;P.ENa=pr.ENa;P.EK=pr.EK;P.EL=pr.EL;P.Cm=pr.Cm;P.T=pr.T;
  STIM.I0=pr.I0;STIM.tStart=pr.tStart;STIM.dur=pr.dur;STIM.type=pr.type;STIM.freq=pr.freq;
  winMs=pr.win;speedIdx=(pr.speed||3)-1;
  const sl=(id,val,dp)=>{
    const e=document.getElementById(id);if(e)e.value=val;
    const d=document.getElementById(id+'-val');if(d)d.textContent=parseFloat(val).toFixed(dp??0);
    const n=document.getElementById(id+'-num');if(n)n.value=parseFloat(val).toFixed(dp??0);
  };
  sl('sl-gNa',pr.gNa);sl('sl-gK',pr.gK);sl('sl-gL',pr.gL,2);
  sl('sl-ENa',pr.ENa);sl('sl-EK',pr.EK);sl('sl-EL',pr.EL,1);
  sl('sl-Cm',pr.Cm,1);sl('sl-T',pr.T,1);
  // Reset dt to default when applying a preset
  const dtPreset=pr.dt||0.01;
  DT=dtPreset;
  const dtSl=document.getElementById('sl-dt');if(dtSl)dtSl.value=dtPreset;
  const dtNum=document.getElementById('sl-dt-num');if(dtNum)dtNum.value=dtPreset.toFixed(3);
  sl('sl-I0',pr.I0);sl('sl-tStart',pr.tStart);sl('sl-dur',pr.dur);sl('sl-freq',pr.freq);sl('sl-win',pr.win);
  const spEl=document.getElementById('sl-speed');if(spEl)spEl.value=pr.speed||3;
  const labels=['⅛×','¼×','½×','1×','2×','4×','7×','14×'];
  const spV=document.getElementById('sl-speed-val');if(spV)spV.textContent=labels[speedIdx]||'—';
  const stEl=document.getElementById('stim-type');if(stEl)stEl.value=pr.type;
  // Show/hide dur row based on type
  const isCont2=pr.type==='continuous';
  const durRow=document.getElementById('dur-row');
  if(durRow){ if(isCont2) durRow.classList.add('disabled'); else durRow.classList.remove('disabled'); }
  const showFreq=['pulse','sine'].includes(pr.type);
  const frEl=document.getElementById('freq-row');
  if(frEl){ if(!showFreq) frEl.classList.add('disabled'); else frEl.classList.remove('disabled'); }
  if(pr.keepPhase){stopSim();initState();}else{resetSim();}
  renderSS();renderTau();
  if(pr.autoStart)setTimeout(startSim,100);
}

function applyDefaults(){applyPreset('standard');}

// ─── INITIAL SIMULATION ───────────────────────────────────
function runInitialSimulation(){
  if(running)return;
  DT=0.01; // reset to stable default for initial simulation
  initState();
  // Run 12000 steps = 120 ms full window (complete stimulus trace)
  for(let i=0;i<12000;i++){rk4();if(stepCnt%STORE_EVERY===0)storeAndDetect();}
  renderAll();updateStats();
  const st=document.getElementById('sim-status-text');if(st)st.textContent='Ready — click ▶ Run';
}


// ─── LOAD ─────────────────────────────────────────────────
window.addEventListener('load',()=>{
  initState();
  const labels=['⅛×','¼×','½×','1×','2×','4×','7×','14×'];
  const spV=document.getElementById('sl-speed-val');if(spV)spV.textContent=labels[speedIdx];
  const dtSlEl=document.getElementById('sl-dt');if(dtSlEl)dtSlEl.value='0.01';
  const durRowInit=document.getElementById('dur-row');
  if(durRowInit) durRowInit.classList.add('disabled');
  const freqRowInit=document.getElementById('freq-row');
  if(freqRowInit) freqRowInit.classList.add('disabled');
  const dtNumEl=document.getElementById('sl-dt-num');if(dtNumEl)dtNumEl.value='0.010';
  renderSS();renderTau();renderPhase();renderIext_plot();updateStats();
  runInitialSimulation();
  vcInitInteraction(); // attach zoom/pan to VC canvases on page load
});

window.addEventListener('resize',()=>{
  if(hist.t.length>0)renderAll();else{renderSS();renderTau();renderPhase();}
});

// ════════════════════════════════════════════════════════════
// ════════════════════════════════════════════════════════════════


// ─── plotTrace: generic straight-line plot for voltage clamp ─────────────
// (plotECG is used for current-clamp ECG wrap; this is a simple x-y trace)
function plotTrace(ctx, xs, ys, W, H, yMin, yMax, xMin, xMax, color, lw) {
  if (!ctx || !xs || xs.length < 2) return;
  const p = PAD, pw = W-p.left-p.right, ph = H-p.top-p.bottom;
  ctx.save();
  ctx.beginPath(); ctx.rect(p.left, p.top, pw, ph); ctx.clip();
  ctx.strokeStyle = color; ctx.lineWidth = lw || 1.5;
  ctx.beginPath();
  let first = true;
  for (let i = 0; i < xs.length; i++) {
    if (xs[i] < xMin - 0.01 || xs[i] > xMax + 0.01) continue;
    const x = p.left + pw * ((xs[i] - xMin) / (xMax - xMin));
    const y = p.top  + ph * (1 - (ys[i] - yMin) / (yMax - yMin));
    if (first) { ctx.moveTo(x, y); first = false; }
    else         ctx.lineTo(x, y);
  }
  ctx.stroke();
  ctx.restore();
}

// ════════════════════════════════════════════════════════════════
// VOLTAGE CLAMP ENGINE — clean rewrite, verified against Python sim
// ════════════════════════════════════════════════════════════════

// ── State ────────────────────────────────────────────────────────
const vcP = {
  gNa:120, gK:36, gL:0.3, ENa:50, EK:-77, EL:-54.4,
  blockNa:false, blockK:false,
  Vhold:-65, Vcmd:20, tStart:10, dur:40, T:6.3,
};
const vcDT = 0.01;
let vcT = 0;
let vcY = [-65, 0.053, 0.596, 0.318]; // [V, m, h, n] at rest (V=-65 mV)
let vcRunning = false, vcAnimId = null;
const vcHist = {t:[],V:[],m:[],h:[],n:[],INa:[],IK:[],IL:[]};
const vcMAX = 12000; // 120ms @ dt=0.01, every step stored

// I–V accumulated data
const ivData = [];
// Which curves to show on I-V plot
const vcIVShow = {INa: true, IK: true, Itot: true};
const vcIShow  = {INa: true, IK: true, IL: true};   // I(t) plot toggles
const vcGShow  = {m: true, h: true, n: true};        // Gates plot toggles
// Persistent run history — each entry is {t,V,INa,IK,IL,m,h,n,vcmd,vhold}
const vcAllRuns = [];

// ── Voltage waveform ─────────────────────────────────────────────
function vcVoltage(t) {
  const tEnd = vcP.tStart + vcP.dur;
  if (t < vcP.tStart) return vcP.Vhold;
  if (t <= tEnd)       return vcP.Vcmd;
  return vcP.Vhold;
}

// ── Gate ODE (V is clamped — only m, h, n evolve) ────────────────
function phiVC() { return Math.pow(3, (vcP.T - 6.3) / 10); }

function vcDydt(y, t) {
  const V = vcVoltage(t);
  const [, m, h, n] = y;
  const ph = phiVC();
  return [
    0,
    ph * (am(V)*(1-m) - bm(V)*m),
    ph * (ah(V)*(1-h) - bh(V)*h),
    ph * (an(V)*(1-n) - bn(V)*n),
  ];
}

// ── RK4 ──────────────────────────────────────────────────────────
function vcRK4() {
  const dt = vcDT;
  const k1 = vcDydt(vcY, vcT);
  const y2 = vcY.map((v,i) => v + 0.5*dt*k1[i]);
  const k2 = vcDydt(y2, vcT + 0.5*dt);
  const y3 = vcY.map((v,i) => v + 0.5*dt*k2[i]);
  const k3 = vcDydt(y3, vcT + 0.5*dt);
  const y4 = vcY.map((v,i) => v + dt*k3[i]);
  const k4 = vcDydt(y4, vcT + dt);
  vcY = vcY.map((v,i) => v + (dt/6)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
  vcY[0] = vcVoltage(vcT + dt); // enforce clamp
  vcT += dt;
}

// ── Store ─────────────────────────────────────────────────────────
function vcStore() {
  const [V, m, h, n] = vcY;
  const INa = vcP.blockNa ? 0 : vcP.gNa * m*m*m * h * (V - vcP.ENa);
  const IK  = vcP.blockK  ? 0 : vcP.gK  * n*n*n*n   * (V - vcP.EK);
  const IL  = vcP.gL * (V - vcP.EL);
  vcHist.t.push(vcT);  vcHist.V.push(V);
  vcHist.m.push(m);    vcHist.h.push(h);   vcHist.n.push(n);
  vcHist.INa.push(INa);vcHist.IK.push(IK); vcHist.IL.push(IL);
  if (vcHist.t.length > vcMAX) {
    const drop = vcHist.t.length - vcMAX;
    for (const k of Object.keys(vcHist)) vcHist[k].splice(0, drop);
  }
}

// ── Duration of one clamp run ────────────────────────────────────
function vcTotalDur() { return vcP.tStart + vcP.dur + 10; }

// ── Animation ─────────────────────────────────────────────────────
const VC_STEPS_PER_FRAME = 50;

function vcAnimate() {
  if (!vcRunning) return;
  for (let i = 0; i < VC_STEPS_PER_FRAME; i++) {
    vcRK4();
    vcStore();
    if (vcT >= vcTotalDur()) {
      vcStop();
      vcSaveRun();
      vcRecordIV();
      vcDrawAll();
      return;
    }
  }
  vcDrawAll();
  vcAnimId = requestAnimationFrame(vcAnimate);
}

function vcRun() {
  if (vcRunning) return;
  // Sync params from sliders before starting
  vcSyncSliders();
  // Re-initialize gates at current Vhold
  vcInitGates();
  vcRunning = true;
  const st = document.getElementById('sim-status-text');
  if (st) st.textContent = 'VC Running…';
  vcAnimId = requestAnimationFrame(vcAnimate);
}

function vcStop() {
  vcRunning = false;
  if (vcAnimId) { cancelAnimationFrame(vcAnimId); vcAnimId = null; }
  const st = document.getElementById('sim-status-text');
  if (st && st.textContent === 'VC Running…') st.textContent = 'VC Paused';
}

function vcReset() {
  vcStop();
  vcSyncSliders();
  vcInitGates();
  for (const k of Object.keys(vcHist)) vcHist[k] = [];
  vcAllRuns.length = 0; // clear all trace history
  vcT = 0;
  vcDrawAll();
  const st = document.getElementById('sim-status-text');
  if (st) st.textContent = 'VC Ready — click ▶ Run Clamp';
}

// ── Sync all VC slider values into vcP ──────────────────────────
function vcSyncSliders() {
  const map = {
    'vc-vhold': 'Vhold', 'vc-vcmd': 'Vcmd',
    'vc-tpre': 'tStart', 'vc-dur': 'dur',
    'vc-gna': 'gNa',    'vc-gk': 'gK',  'vc-gl': 'gL',
    'vc-ena': 'ENa',    'vc-ek': 'EK',  'vc-el': 'EL',  'vc-T': 'T',
  };
  for (const [id, field] of Object.entries(map)) {
    const el = document.getElementById(id);
    if (el) {
      const v = parseFloat(el.value);
      if (!isNaN(v)) vcP[field] = v;
    }
  }
  const naEl = document.getElementById('vc-block-na');
  const kEl  = document.getElementById('vc-block-k');
  if (naEl) vcP.blockNa = naEl.checked;
  if (kEl)  vcP.blockK  = kEl.checked;
}

// ── Initialize gates at steady state for current Vhold ──────────
function vcInitGates() {
  const V0 = vcP.Vhold;
  const m0 = am(V0)/(am(V0)+bm(V0));
  const h0 = ah(V0)/(ah(V0)+bh(V0));
  const n0 = an(V0)/(an(V0)+bn(V0));
  vcY = [V0, m0, h0, n0];
  vcT = 0;
}

// ── UI helpers ───────────────────────────────────────────────────
function vcSlider(param, val) {
  const v = parseFloat(val);
  if (isNaN(v)) return;
  const map = {
    vhold:'Vhold', vcmd:'Vcmd', tpre:'tStart', dur:'dur',
    gna:'gNa', gk:'gK', gl:'gL', ena:'ENa', ek:'EK', el:'EL', T:'T',
  };
  if (map[param] !== undefined) vcP[map[param]] = v;
  const numEl = document.getElementById(`vc-${param}-num`);
  if (numEl) numEl.value = v;
}

function syncVCNum(sliderId, val) {
  const el = document.getElementById(sliderId);
  if (el) { el.value = val; el.dispatchEvent(new Event('input')); }
}

// ── I–V curve ────────────────────────────────────────────────────

// Save completed run to persistent history
function vcSaveRun() {
  if (!vcHist.t.length) return;
  vcAllRuns.push({
    t:    [...vcHist.t],
    V:    [...vcHist.V],
    INa:  [...vcHist.INa],
    IK:   [...vcHist.IK],
    IL:   [...vcHist.IL],
    m:    [...vcHist.m],
    h:    [...vcHist.h],
    n:    [...vcHist.n],
    vcmd: vcP.Vcmd,
    vhold: vcP.Vhold,
  });
  // Keep last 20 runs max
  if (vcAllRuns.length > 20) vcAllRuns.shift();
}

function vcRecordIV() {
  if (!vcHist.t.length) return;
  let peakINa = 0, peakIK = 0, peakItot = 0;
  for (let i = 0; i < vcHist.t.length; i++) {
    const t = vcHist.t[i];
    if (t < vcP.tStart || t > vcP.tStart + vcP.dur) continue;
    const itot = vcHist.INa[i] + vcHist.IK[i] + vcHist.IL[i];
    if (vcHist.INa[i] < peakINa) peakINa = vcHist.INa[i];
    if (vcHist.IK[i]  > peakIK)  peakIK  = vcHist.IK[i];
    if (Math.abs(itot) > Math.abs(peakItot)) peakItot = itot;
  }
  const Vc = vcP.Vcmd;
  const existing = ivData.findIndex(d => d.vcmd === Vc);
  const entry = {vcmd:Vc, peakINa, peakIK, peakItot};
  if (existing >= 0) ivData[existing] = entry; else ivData.push(entry);
  ivData.sort((a,b) => a.vcmd - b.vcmd);
  vcUpdateIVTable();
}

function vcUpdateIVTable() {
  const tbody = document.getElementById('iv-tbody');
  if (!tbody) return;
  if (!ivData.length) {
    tbody.innerHTML = '<tr><td colspan="4" style="text-align:center;color:var(--dim);padding:8px">Run clamp at different V<sub>cmd</sub> to populate</td></tr>';
    return;
  }
  tbody.innerHTML = ivData.map(d =>
    `<tr><td>${d.vcmd}</td><td style="color:#dc2626">${d.peakINa.toFixed(2)}</td><td style="color:#1d4ed8">${d.peakIK.toFixed(2)}</td><td>${d.peakItot.toFixed(2)}</td></tr>`
  ).join('');
}

function vcClearIV() {
  ivData.length = 0;
  vcUpdateIVTable();
  vcDrawIV();
}

// ── Drawing ───────────────────────────────────────────────────────


// ── Voltage Sweep ────────────────────────────────────────────────
// Runs multiple clamps sequentially at evenly-spaced Vcmd values.
let vcSweepQueue = [];
let vcSweepIdx = 0;

function vcStartSweep() {
  if (vcRunning) return;
  vcSyncSliders();
  // Read sweep range from UI
  const vFrom  = parseFloat(document.getElementById('vc-sweep-from')?.value  ?? -20);
  const vTo    = parseFloat(document.getElementById('vc-sweep-to')?.value    ?? 100);
  const vStep  = parseFloat(document.getElementById('vc-sweep-step')?.value  ?? 10);
  if (isNaN(vFrom) || isNaN(vTo) || isNaN(vStep) || vStep <= 0) return;
  vcSweepQueue = [];
  const sign = vTo >= vFrom ? 1 : -1;
  for (let v = vFrom; sign*(v - vTo) <= 1e-6; v += sign * Math.abs(vStep)) {
    vcSweepQueue.push(Math.round(v * 100) / 100);
  }
  if (!vcSweepQueue.length) return;
  vcSweepIdx = 0;
  const prog = document.getElementById('vc-sweep-prog');
  if (prog) { prog.textContent = `0 / ${vcSweepQueue.length}`; prog.style.display=''; }
  vcRunNextSweep();
}

function vcRunNextSweep() {
  if (vcSweepIdx >= vcSweepQueue.length) {
    // All done — draw final result and update UI
    vcDrawAll();
    const prog = document.getElementById('vc-sweep-prog');
    if (prog) prog.textContent = 'Done ✓';
    const st = document.getElementById('sim-status-text');
    if (st) st.textContent = `Sweep complete — ${vcSweepQueue.length} steps`;
    return;
  }
  const vcmd = vcSweepQueue[vcSweepIdx];
  vcP.Vcmd = vcmd;
  const slEl = document.getElementById('vc-vcmd');   if (slEl)  slEl.value  = vcmd;
  const numEl = document.getElementById('vc-vcmd-num'); if (numEl) numEl.value = vcmd;

  // Run this voltage step SYNCHRONOUSLY — no animation, just crunch numbers
  vcInitGates();
  for (const k of Object.keys(vcHist)) vcHist[k] = [];
  vcT = 0;
  const nSteps = Math.ceil(vcTotalDur() / vcDT);
  for (let i = 0; i < nSteps; i++) { vcRK4(); vcStore(); }
  vcSaveRun();
  vcRecordIV();

  vcSweepIdx++;
  const prog = document.getElementById('vc-sweep-prog');
  if (prog) prog.textContent = `${vcSweepIdx} / ${vcSweepQueue.length}`;

  // Yield to browser every 4 steps so UI stays responsive and shows progress
  if (vcSweepIdx < vcSweepQueue.length) {
    if (vcSweepIdx % 4 === 0) {
      vcDrawAll(); // repaint to show accumulated traces so far
      setTimeout(vcRunNextSweep, 0); // yield, then continue
    } else {
      vcRunNextSweep(); // next step immediately (no yield)
    }
  } else {
    vcRunNextSweep(); // triggers "done" branch above
  }
}

// vcSweepAnimate no longer needed — kept as no-op for safety
function vcSweepAnimate() {}

function vcStopSweep() {
  vcSweepQueue = []; vcSweepIdx = 0;
  vcStop();
  const prog = document.getElementById('vc-sweep-prog');
  if (prog) { prog.textContent = 'Stopped'; }
}


// ─── VC Canvas Pan / Zoom state ────────────────────────────────
// Each VC canvas has an independent view: {xMin,xMax,yMin,yMax} or null=auto
const vcView = {'vc-cv-V':null,'vc-cv-I':null,'vc-cv-gates':null,'vc-cv-IV':null};
const vcDrag = {}; // tracks active drag per canvas id

function vcResetView(id) {
  vcView[id] = null;
  vcDrawAll();
}

function vcResetAllViews() {
  for (const k of Object.keys(vcView)) vcView[k] = null;
  vcDrawAll();
}

// Attach wheel + drag pan to a canvas
function vcAttachInteraction(id) {
  const el = document.getElementById(id);
  if (!el || el._vcInteract) return;
  el._vcInteract = true;

  // Helper: get cursor position in CSS pixels relative to canvas top-left
  function cursorPos(e) {
    const r = el.getBoundingClientRect();
    return {x: e.clientX - r.left, y: e.clientY - r.top, W: r.width, H: r.height};
  }

  // Helper: is cursor inside the plot area?
  function inPlot(x, y, W, H) {
    const p = PAD;
    return x >= p.left && x <= W - p.right && y >= p.top && y <= H - p.bottom;
  }

  // Wheel = zoom centered on cursor (uses CSS pixel coords throughout)
  el.addEventListener('wheel', function(e) {
    e.preventDefault();
    const {x, y, W, H} = cursorPos(e);
    if (!W || !inPlot(x, y, W, H)) return;
    const p = PAD;
    const pw = W - p.left - p.right, ph = H - p.top - p.bottom;
    const view = vcGetView(id);
    const {xMin, xMax, yMin, yMax} = view;
    // Fraction within plot area (0-1)
    const fx = (x - p.left) / pw;
    const fy = 1 - (y - p.top) / ph;
    // Data coordinate under cursor
    const dataX = xMin + fx * (xMax - xMin);
    const dataY = yMin + fy * (yMax - yMin);
    // Zoom factor
    const f = e.deltaY > 0 ? 1.3 : 1 / 1.3;
    vcView[id] = {
      xMin: dataX - fx * (xMax - xMin) * f,
      xMax: dataX + (1 - fx) * (xMax - xMin) * f,
      yMin: dataY - fy * (yMax - yMin) * f,
      yMax: dataY + (1 - fy) * (yMax - yMin) * f,
    };
    vcDrawAll();
  }, {passive: false});

  // Drag = pan
  let drag = null;
  el.addEventListener('mousedown', function(e) {
    if (e.button !== 0) return;
    const {x, y, W, H} = cursorPos(e);
    if (!W || !inPlot(x, y, W, H)) return;
    drag = {x0: e.clientX, y0: e.clientY, view: {...vcGetView(id)}, W, H};
    el.style.cursor = 'grabbing';
    e.preventDefault();
  });
  el.addEventListener('mousemove', function(e) {
    if (!drag) return;
    const p = PAD;
    const pw = drag.W - p.left - p.right;
    const ph = drag.H - p.top - p.bottom;
    const {xMin, xMax, yMin, yMax} = drag.view;
    const dx = -(e.clientX - drag.x0) / pw * (xMax - xMin);
    const dy =  (e.clientY - drag.y0) / ph * (yMax - yMin);
    vcView[id] = {xMin: xMin + dx, xMax: xMax + dx, yMin: yMin + dy, yMax: yMax + dy};
    vcDrawAll();
  });
  function endDrag() { drag = null; el.style.cursor = 'default'; }
  el.addEventListener('mouseup', endDrag);
  el.addEventListener('mouseleave', endDrag);
}

// Get current view or compute auto view for a given canvas
function vcGetView(id) {
  if (vcView[id]) return vcView[id];
  // Compute auto view
  const tMax = vcTotalDur();
  if (id === 'vc-cv-V') {
    const lo = Math.min(vcP.Vhold, vcP.Vcmd)-10, hi = Math.max(vcP.Vhold, vcP.Vcmd)+30;
    return {xMin:0, xMax:tMax, yMin:Math.min(-10,lo), yMax:Math.max(130,hi)};
  }
  if (id === 'vc-cv-I') {
    let mxI=2;
    for (let i=0;i<vcHist.INa.length;i++) {
      const t=vcHist.t[i];
      if (t<vcP.tStart-1||t>vcP.tStart+vcP.dur+1) continue;
      const mx=Math.max(Math.abs(vcHist.INa[i]),Math.abs(vcHist.IK[i]),Math.abs(vcHist.IL[i]));
      if (mx>mxI) mxI=mx;
    }
    mxI=Math.max(2, mxI*1.2);
    return {xMin:0, xMax:tMax, yMin:-mxI, yMax:mxI};
  }
  if (id === 'vc-cv-gates') return {xMin:0, xMax:tMax, yMin:0, yMax:1};
  if (id === 'vc-cv-IV') {
    let vMin=-20,vMax=130,iMin=-5,iMax=5;
    if (ivData.length>0) {
      const allV=ivData.map(d=>d.vcmd);
      const allI=ivData.flatMap(d=>[d.peakINa,d.peakIK,d.peakItot]);
      const vSpan=Math.max(...allV)-Math.min(...allV);
      const pad=vSpan>0?vSpan*0.1:10;
      vMin=Math.floor((Math.min(...allV)-pad)/10)*10;
      vMax=Math.ceil((Math.max(...allV)+pad)/10)*10;
      iMin=Math.min(-1,Math.min(...allI)*1.2);
      iMax=Math.max(1,Math.max(...allI)*1.2);
    }
    return {xMin:vMin, xMax:vMax, yMin:iMin, yMax:iMax};
  }
  return {xMin:0, xMax:tMax, yMin:-10, yMax:10};
}

// Call once when switching to VC tab
function vcInitInteraction() {
  ['vc-cv-V','vc-cv-I','vc-cv-gates','vc-cv-IV'].forEach(vcAttachInteraction);
}

// Draw all previous runs faintly behind the current run
function vcDrawPrevRuns(ctx, W, H, yMin, yMax, xMin, xMax, field, color) {
  if (!vcAllRuns.length) return;
  const p = PAD, pw = W-p.left-p.right, ph = H-p.top-p.bottom;
  ctx.save();
  // Clip to plot area so traces never bleed into axis labels
  ctx.beginPath(); ctx.rect(p.left, p.top, pw, ph); ctx.clip();
  ctx.globalAlpha = 0.18;
  ctx.strokeStyle = color;
  ctx.lineWidth = 1;
  for (const run of vcAllRuns) {
    const xs = run.t, ys = run[field];
    ctx.beginPath();
    let first = true;
    for (let i = 0; i < xs.length; i++) {
      const x = p.left + pw * ((xs[i] - xMin) / (xMax - xMin));
      const y = p.top  + ph * (1 - (ys[i] - yMin) / (yMax - yMin));
      if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }
  ctx.globalAlpha = 1;
  ctx.restore();
}



function vcToggleI(key) {
  const naEl = document.getElementById('i-show-na');
  const kEl  = document.getElementById('i-show-k');
  const lEl  = document.getElementById('i-show-l');
  vcIShow.INa = naEl ? naEl.checked : true;
  vcIShow.IK  = kEl  ? kEl.checked  : true;
  vcIShow.IL  = lEl  ? lEl.checked  : true;
  vcView['vc-cv-I'] = null; // reset Y auto-fit to visible curves
  vcDrawI();
}

function vcToggleG(key) {
  const mEl = document.getElementById('g-show-m');
  const hEl = document.getElementById('g-show-h');
  const nEl = document.getElementById('g-show-n');
  vcGShow.m = mEl ? mEl.checked : true;
  vcGShow.h = hEl ? hEl.checked : true;
  vcGShow.n = nEl ? nEl.checked : true;
  vcDrawGates();
}

function vcToggleIV(key) {
  // Sync all three checkbox states directly from DOM → no flip desync
  const naEl  = document.getElementById('iv-show-na');
  const kEl   = document.getElementById('iv-show-k');
  const totEl = document.getElementById('iv-show-tot');
  vcIVShow.INa  = naEl  ? naEl.checked  : true;
  vcIVShow.IK   = kEl   ? kEl.checked   : true;
  vcIVShow.Itot = totEl ? totEl.checked : true;
  // Reset Y axis to re-fit to visible curves
  vcView['vc-cv-IV'] = null;
  vcDrawIV();
}

function vcDrawAll() {
  vcDrawV(); vcDrawI(); vcDrawGates(); vcDrawIV();
}

function vcDrawV() {
  const {ctx,W,H} = setupCanvas('vc-cv-V'); if (!ctx) return;
  const tMax = vcTotalDur();
  const {xMin,xMax,yMin,yMax} = vcGetView('vc-cv-V');
  drawAxes(ctx, W, H, yMin, yMax, xMin, xMax, 'mV');
  // Step onset marker — only if visible in current view
  const p = PAD, pw = W-p.left-p.right, ph = H-p.top-p.bottom;
  const xS = p.left + pw*((vcP.tStart - xMin)/(xMax - xMin));
  if (xS > p.left && xS < p.left+pw) {
    ctx.save();
    ctx.beginPath(); ctx.rect(p.left,p.top,pw,ph); ctx.clip();
    ctx.strokeStyle='rgba(220,38,38,.5)'; ctx.lineWidth=1; ctx.setLineDash([4,3]);
    ctx.beginPath(); ctx.moveTo(xS,p.top); ctx.lineTo(xS,p.top+ph); ctx.stroke();
    ctx.setLineDash([]); ctx.restore();
  }
  vcDrawPrevRuns(ctx, W, H, yMin, yMax, xMin, xMax, 'V', '#0b7a6e');
  if (vcHist.t.length > 1)
    plotTrace(ctx, vcHist.t, vcHist.V, W, H, yMin, yMax, xMin, xMax, '#0b7a6e', 2);
}

function vcDrawI() {
  const {ctx,W,H} = setupCanvas('vc-cv-I'); if (!ctx) return;
  const tMax = vcTotalDur();
  const {xMin,xMax,yMin,yMax} = vcGetView('vc-cv-I');
  drawAxes(ctx, W, H, yMin, yMax, xMin, xMax, 'µA/cm²');
  if (vcIShow.INa) vcDrawPrevRuns(ctx, W, H, yMin, yMax, xMin, xMax, 'INa', '#dc2626');
  if (vcIShow.IK)  vcDrawPrevRuns(ctx, W, H, yMin, yMax, xMin, xMax, 'IK',  '#1d4ed8');
  if (vcHist.t.length > 1) {
    if (vcIShow.INa) plotTrace(ctx, vcHist.t, vcHist.INa, W, H, yMin, yMax, xMin, xMax, '#dc2626', 2);
    if (vcIShow.IK)  plotTrace(ctx, vcHist.t, vcHist.IK,  W, H, yMin, yMax, xMin, xMax, '#1d4ed8', 1.5);
    if (vcIShow.IL)  plotTrace(ctx, vcHist.t, vcHist.IL,  W, H, yMin, yMax, xMin, xMax, '#15803d', 1);
  }
}

function vcDrawGates() {
  const {ctx,W,H} = setupCanvas('vc-cv-gates'); if (!ctx) return;
  const tMax = vcTotalDur();
  const {xMin,xMax,yMin,yMax} = vcGetView('vc-cv-gates');
  drawAxes(ctx, W, H, yMin, yMax, xMin, xMax, '');
  if (vcGShow.m) vcDrawPrevRuns(ctx, W, H, yMin, yMax, xMin, xMax, 'm', '#dc2626');
  if (vcGShow.n) vcDrawPrevRuns(ctx, W, H, yMin, yMax, xMin, xMax, 'n', '#7c3aed');
  if (vcHist.t.length > 1) {
    if (vcGShow.m) plotTrace(ctx, vcHist.t, vcHist.m, W, H, yMin, yMax, xMin, xMax, '#dc2626', 1.5);
    if (vcGShow.h) plotTrace(ctx, vcHist.t, vcHist.h, W, H, yMin, yMax, xMin, xMax, '#f59e0b', 1.5);
    if (vcGShow.n) plotTrace(ctx, vcHist.t, vcHist.n, W, H, yMin, yMax, xMin, xMax, '#7c3aed', 1.5);
  }
}

function vcDrawIV() {
  const {ctx,W,H} = setupCanvas('vc-cv-IV'); if (!ctx) return;
  const p = PAD, pw = W-p.left-p.right, ph = H-p.top-p.bottom;
  // Use vcGetView — respects pan/zoom state; auto-fits when null
  const {xMin:vMin, xMax:vMax, yMin:iMin, yMax:iMax} = vcGetView('vc-cv-IV');
  drawAxes(ctx, W, H, iMin, iMax, vMin, vMax, 'µA/cm²');
  // Zero lines (clipped to plot area)
  ctx.save();
  ctx.beginPath(); ctx.rect(p.left, p.top, pw, ph); ctx.clip();
  const y0 = p.top + ph*(1-(0-iMin)/(iMax-iMin));
  const x0 = p.left + pw*(0-vMin)/(vMax-vMin);
  ctx.strokeStyle='rgba(0,0,0,.15)'; ctx.lineWidth=1;
  if (y0 >= p.top && y0 <= p.top+ph) {
    ctx.beginPath(); ctx.moveTo(p.left,y0); ctx.lineTo(p.left+pw,y0); ctx.stroke();
  }
  if (x0 >= p.left && x0 <= p.left+pw) {
    ctx.beginPath(); ctx.moveTo(x0,p.top); ctx.lineTo(x0,p.top+ph); ctx.stroke();
  }
  if (!ivData.length) {
    ctx.fillStyle='#9ca3af'; ctx.font="11px 'DM Sans',sans-serif"; ctx.textAlign='center';
    ctx.fillText('Run clamp at multiple Vcmd values to build the I–V curve', p.left+pw/2, p.top+ph/2-8);
    ctx.fillText('Use channel blockers (TTX/TEA) to isolate INa or IK', p.left+pw/2, p.top+ph/2+10);
    ctx.restore();
    return;
  }
  const xp = v => p.left + pw*(v-vMin)/(vMax-vMin);
  const yp = i => p.top  + ph*(1-(i-iMin)/(iMax-iMin));
  const drawLine = (arr, color) => {
    ctx.strokeStyle=color; ctx.fillStyle=color; ctx.lineWidth=2;
    ctx.beginPath();
    arr.forEach(({x,y},idx) => idx===0 ? ctx.moveTo(x,y) : ctx.lineTo(x,y));
    ctx.stroke();
    arr.forEach(({x,y}) => { ctx.beginPath(); ctx.arc(x,y,4,0,Math.PI*2); ctx.fill(); });
  };
  if (vcIVShow.INa)  drawLine(ivData.map(d=>({x:xp(d.vcmd),y:yp(d.peakINa)})),  '#dc2626');
  if (vcIVShow.IK)   drawLine(ivData.map(d=>({x:xp(d.vcmd),y:yp(d.peakIK)})),   '#1d4ed8');
  if (vcIVShow.Itot) drawLine(ivData.map(d=>({x:xp(d.vcmd),y:yp(d.peakItot)})), '#374151');
  ctx.restore();
}

// ── Tab switching ────────────────────────────────────────────────
function switchSimTab(tab) {
  ['cc','vc'].forEach(t => {
    const tb = document.getElementById(`tab-${t}`);
    const pn = document.getElementById(`panel-${t}`);
    if (tb) tb.classList.toggle('active', t === tab);
    if (pn) pn.classList.toggle('active', t === tab);
  });
  if (tab === 'vc') {
    vcSyncSliders();
    vcInitGates();
    for (const k of Object.keys(vcHist)) vcHist[k] = [];
    vcAllRuns.length = 0;
    vcT = 0;
    // Force sweep defaults in case of cached HTML
    const sfEl=document.getElementById('vc-sweep-from');  if(sfEl&&!sfEl._init){sfEl.value='-70';sfEl._init=true;}
    const stEl=document.getElementById('vc-sweep-to');    if(stEl&&!stEl._init){stEl.value='50'; stEl._init=true;}
    const ssEl=document.getElementById('vc-sweep-step'); if(ssEl&&!ssEl._init){ssEl.value='5';  ssEl._init=true;}
    vcInitInteraction(); // attach zoom/pan listeners (safe to call multiple times)
    vcDrawAll();
    const st = document.getElementById('sim-status-text');
    if (st) st.textContent = 'VC Ready — click ▶ Run Clamp';
  }
}

// ─── SYNC: slider ↔ number input ───────────────────────
function syncSlider(slId, val){
  const numEl = document.getElementById(slId+'-num');
  // Don't overwrite if the number field is currently focused (user is typing)
  if(numEl && document.activeElement !== numEl){
    const dp = numEl.step && numEl.step.includes('.') ?
      numEl.step.split('.')[1].length : 0;
    numEl.value = parseFloat(val).toFixed(dp);
  }
}
function syncNum(slId, val){
  const slEl = document.getElementById(slId);
  if(!slEl) return;
  const numVal = parseFloat(val);
  if(isNaN(numVal)) return;
  const mn = parseFloat(slEl.min), mx = parseFloat(slEl.max);
  const clamped = Math.min(mx, Math.max(mn, numVal));
  slEl.value = clamped;
  // Trigger the slider's oninput
  slEl.dispatchEvent(new Event('input', {bubbles:true}));
}

// ─── DEFAULTS ────────────────────────────────────────────
const DEFAULTS = {
  gNa:120, gK:36, gL:0.3, ENa:50, EK:-77, EL:-54.387, Cm:1, T:6.3,
  I0:10, tStart:10, dur:40, type:'step', freq:40, win:160, speed:3
};
function applyDefaults(){
  applyPreset('standard');
}

// ─── ZOOM / PAN / HOME ───────────────────────────────────
// Drag to pan, scroll to zoom, Home resets to auto-follow
let _dragActive=false, _dragX0=0, _dragV0=0, _dragV1=0;

function attachPlotInteraction(){
  const plotIds=['cv-V','cv-gates','cv-currents'];
  plotIds.forEach(id=>{
    const c=document.getElementById(id);
    if(!c) return;
    c.style.cursor='grab';

    // Mouse drag — pan
    c.addEventListener('mousedown',e=>{
      _dragActive=true; _dragX0=e.clientX;
      _dragV0=viewStart; _dragV1=viewEnd;
      c.style.cursor='grabbing'; e.preventDefault();
    });
    c.addEventListener('mousemove',e=>{
      if(!_dragActive) return;
      const pxPerMs=(c.clientWidth||600)/(viewEnd-viewStart);
      const dMs=(e.clientX-_dragX0)/pxPerMs;
      const span=_dragV1-_dragV0;
      const newStart=_dragV0-dMs;
      viewStart=Math.max(0,newStart);
      viewEnd=viewStart+span;
      viewLocked=true;
      renderAll();
    });
    const stopDrag=()=>{ _dragActive=false; c.style.cursor='grab'; };
    c.addEventListener('mouseup',stopDrag);
    c.addEventListener('mouseleave',stopDrag);

    // Scroll — zoom time axis
    c.addEventListener('wheel',e=>{
      e.preventDefault();
      const span=viewEnd-viewStart;
      const factor=e.deltaY>0?1.2:0.833;
      const newSpan=Math.max(20,Math.min(2000,span*factor));
      const center=(viewStart+viewEnd)/2;
      viewStart=Math.max(0,center-newSpan/2);
      viewEnd=viewStart+newSpan;
      winMs=newSpan;
      viewLocked=true;
      // Update window slider
      const slW=document.getElementById('sl-win');
      if(slW){ slW.value=Math.round(newSpan); syncSlider('sl-win',Math.round(newSpan)); }
      renderAll();
    },{passive:false});

    // Touch drag
    let _tx0=0,_tv0=0,_tv1=0;
    c.addEventListener('touchstart',e=>{
      _tx0=e.touches[0].clientX; _tv0=viewStart; _tv1=viewEnd;
    },{passive:true});
    c.addEventListener('touchmove',e=>{
      const pxPerMs=(c.clientWidth||600)/(viewEnd-viewStart);
      const dMs=(e.touches[0].clientX-_tx0)/pxPerMs;
      const span=_tv1-_tv0;
      viewStart=Math.max(0,_tv0-dMs);
      viewEnd=viewStart+span;
      viewLocked=true;
      renderAll();
    },{passive:true});
  });

  // Double-click any plot → home
  plotIds.forEach(id=>{
    const c=document.getElementById(id);
    if(c) c.addEventListener('dblclick',()=>homeView());
  });
}

function homeView(){
  viewLocked=false;
  if(hist.t.length){
    const tEnd=hist.t[hist.t.length-1];
    viewStart=Math.max(0,tEnd-winMs);
    viewEnd=Math.max(winMs,tEnd);
  } else {
    viewStart=0; viewEnd=winMs;
  }
  renderAll();
}

window.addEventListener('DOMContentLoaded',()=>{
  attachPlotInteraction();
});

(function(){
  var FONTS=['','font-serif','font-mono'];
  var LS=window.localStorage;

  function lsGet(k){try{return LS?LS.getItem(k):null;}catch(e){return null;}}
  function lsSet(k,v){try{if(LS)LS.setItem(k,v);}catch(e){}}

  function applyPrefs(){
    var dark=lsGet('acc_dark')==='1';
    var sz=lsGet('acc_sz')||'md';
    var fi=parseInt(lsGet('acc_fi')||'0');
    var b=document.body;
    // Theme
    b.classList.toggle('dark-mode',dark);
    b.classList.toggle('light-mode',!dark);
    // Size
    b.classList.remove('size-sm','size-md','size-lg');
    b.classList.add('size-'+sz);
    // Font
    b.classList.remove('font-serif','font-mono');
    if(FONTS[fi]) b.classList.add(FONTS[fi]);
    // Sync controls
    var cb=document.getElementById('acc-dark');
    if(cb) cb.checked=dark;
    ['sm','md','lg'].forEach(function(s){
      var el=document.getElementById('acc-'+s);
      if(el) el.classList.toggle('active',s===sz);
    });
    [0,1,2].forEach(function(i){
      var el=document.getElementById('acc-f'+i);
      if(el) el.classList.toggle('active',i===fi);
    });
  }

  window.accToggle=function(){
    var p=document.getElementById('acc-panel');
    if(p) p.classList.toggle('open');
  };
  window.accDark=function(cb){lsSet('acc_dark',cb.checked?'1':'0');applyPrefs();};
  window.accSize=function(s){lsSet('acc_sz',s);applyPrefs();};
  window.accFont=function(i){lsSet('acc_fi',i);applyPrefs();};

  // Close on outside click
  document.addEventListener('click',function(e){
    var p=document.getElementById('acc-panel');
    var b2=document.getElementById('acc-btn');
    if(p&&b2&&p.classList.contains('open')&&!p.contains(e.target)&&e.target!==b2)
      p.classList.remove('open');
  });

  // Apply on load (before DOMContentLoaded to avoid flash)
  applyPrefs();
  document.addEventListener('DOMContentLoaded',applyPrefs);
})();