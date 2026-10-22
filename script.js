/* ================= STATE (3D) ================= */
const MAXS=10000,MAXN=MAXS*8;
let N=4000,CORE_B=2400;const G=0.8,SOFT=16;
const cfg={theta:.8,baseDT:.4,dt:.4,grid:false,approx:true,slow:true,brute:false,follow:-1,zoom:1,tz:1,cx:0,cy:0,cz:0,yaw:.6,pitch:.55,manual:false,userRot:false};
const mk=()=>new Float32Array(MAXS);
const px=mk(),py=mk(),pz=mk(),vx=mk(),vy=mk(),vz=mk(),ax=mk(),ay=mk(),az=mk(),mass=mk();
const nCx=new Float32Array(MAXN),nCy=new Float32Array(MAXN),nCz=new Float32Array(MAXN),nHalf=new Float32Array(MAXN),
      nChild=new Int32Array(MAXN*8),nPart=new Int32Array(MAXN),nMass=new Float32Array(MAXN),
      nComX=new Float32Array(MAXN),nComY=new Float32Array(MAXN),nComZ=new Float32Array(MAXN);
let nodeCount=0,checks=0,coreDist=520,minDist=9999,frames=0,impacted=false,shock=null,shake=0,simTime=0;
const hist=[];
const RM=matchMedia('(prefers-reduced-motion:reduce)').matches;
const approx=[],lists=[[],[],[],[],[],[]];

/* ================= OCTREE ================= */
function alloc(cx,cy,cz,h){
  if(nodeCount>=MAXN)return -1;
  const id=nodeCount++;
  nCx[id]=cx;nCy[id]=cy;nCz[id]=cz;nHalf[id]=h;nPart[id]=-1;nMass[id]=0;
  nChild.fill(-1,id*8,id*8+8);
  return id;
}
const oct=(n,x,y,z)=>(x>=nCx[n]?1:0)+(y>=nCy[n]?2:0)+(z>=nCz[n]?4:0);
function child(n,o){
  const h=nHalf[n]*.5;
  const id=alloc(nCx[n]+(o&1?h:-h),nCy[n]+(o&2?h:-h),nCz[n]+(o&4?h:-h),h);
  if(id!==-1)nChild[n*8+o]=id;
  return id;
}
function insert(n,p,d){
  if(n===-1)return;
  if(nPart[n]===-1){nPart[n]=p;return}
  if(d>18)return;
  if(nPart[n]>=0){
    const old=nPart[n];nPart[n]=-2;
    insert(child(n,oct(n,px[old],py[old],pz[old])),old,d+1);
  }
  const o=oct(n,px[p],py[p],pz[p]);
  let c=nChild[n*8+o];
  if(c===-1)c=child(n,o);
  insert(c,p,d+1);
}
function build(){
  let x0=1e9,x1=-1e9,y0=1e9,y1=-1e9,z0=1e9,z1=-1e9;
  for(let i=0;i<N;i++){
    if(px[i]<x0)x0=px[i];if(px[i]>x1)x1=px[i];
    if(py[i]<y0)y0=py[i];if(py[i]>y1)y1=py[i];
    if(pz[i]<z0)z0=pz[i];if(pz[i]>z1)z1=pz[i];
  }
  nodeCount=0;
  const half=Math.max(x1-x0,y1-y0,z1-z0)*.5*1.01+1;
  const root=alloc((x0+x1)/2,(y0+y1)/2,(z0+z1)/2,half);
  for(let i=0;i<N;i++)insert(root,i,0);
  calcMass(root);
}
function calcMass(n){
  let m=0,mx=0,my=0,mz=0;
  if(nPart[n]>=0){const p=nPart[n];m=mass[p];mx=px[p]*m;my=py[p]*m;mz=pz[p]*m}
  else if(nPart[n]===-2){
    for(let o=0;o<8;o++){
      const c=nChild[n*8+o];
      if(c!==-1){calcMass(c);const cm=nMass[c];m+=cm;mx+=nComX[c]*cm;my+=nComY[c]*cm;mz+=nComZ[c]*cm}
    }
  }
  nMass[n]=m;
  if(m>0){nComX[n]=mx/m;nComY[n]=my/m;nComZ[n]=mz/m}
}

/* ================= PHYSICS ================= */
const acc={x:0,y:0,z:0};
function walk(n,p,xi,yi,zi,rec){
  const m=nMass[n];
  if(m===0||nPart[n]===p)return;
  const dx=nComX[n]-xi,dy=nComY[n]-yi,dz=nComZ[n]-zi;
  const d2=dx*dx+dy*dy+dz*dz+SOFT,d=Math.sqrt(d2);
  const leaf=nPart[n]>=0;
  if(leaf||(nHalf[n]*2/d<cfg.theta)){
    const s=G*m/(d2*d);
    acc.x+=dx*s;acc.y+=dy*s;acc.z+=dz*s;checks++;
    if(rec&&!leaf)approx.push(n);
  }else{
    for(let o=0;o<8;o++){const c=nChild[n*8+o];if(c!==-1)walk(c,p,xi,yi,zi,rec)}
  }
}
function forces(){
  checks=0;approx.length=0;
  if(cfg.brute){
    for(let i=0;i<N;i++){
      let fx=0,fy=0,fz=0;const xi=px[i],yi=py[i],zi=pz[i];
      for(let j=0;j<N;j++){
        if(i===j)continue;
        const dx=px[j]-xi,dy=py[j]-yi,dz=pz[j]-zi,d2=dx*dx+dy*dy+dz*dz+SOFT,inv=1/Math.sqrt(d2);
        const s=G*mass[j]*inv*inv*inv;fx+=dx*s;fy+=dy*s;fz+=dz*s;checks++;
      }
      ax[i]=fx;ay[i]=fy;az[i]=fz;
    }
    return;
  }
  for(let i=0;i<N;i++){
    acc.x=0;acc.y=0;acc.z=0;
    walk(0,i,px[i],py[i],pz[i],i===cfg.follow&&cfg.approx);
    ax[i]=acc.x;ay[i]=acc.y;az[i]=acc.z;
  }
}
let rx=0,ry=0,rz=0;
function rot(x,y,z,ca,sa,cb,sb){
  const y1=y*ca-z*sa,z1=y*sa+z*ca;
  rx=x*cb+z1*sb;ry=y1;rz=-x*sb+z1*cb;
}
function init(){
  CORE_B=Math.floor(N*.6);
  for(const l of lists)l.length=0;
  const T=[[.9,.3],[-.6,-.5]].map(([a,b])=>[Math.cos(a),Math.sin(a),Math.cos(b),Math.sin(b)]);
  for(let i=0;i<N;i++){
    const A=i<CORE_B,t=T[A?0:1];
    const gx=A?-260:260,gy=A?-90:90,gz=A?-60:60,bx=A?1.15:-1.15,by=A?.35:-.35,bz=A?.18:-.18,cm=A?6500:4200;
    if(i===0||i===CORE_B){
      px[i]=gx;py[i]=gy;pz[i]=gz;vx[i]=bx;vy[i]=by;vz[i]=bz;mass[i]=cm;
    }else{
      const span=A?210:150,f=Math.pow(Math.random(),.85),r=18+f*span;
      // two spiral arms plus a smooth disc population
      const arm=Math.random()<.5?0:Math.PI;
      const a=Math.random()<.72?arm-r*.026+(Math.random()+Math.random()-1)*.7:Math.random()*6.2832;
      const bulge=r<70?(70-r)*.45:0;
      const h=(Math.random()+Math.random()-1)*(5+.09*r+bulge);
      rot(Math.cos(a)*r,Math.sin(a)*r,h,t[0],t[1],t[2],t[3]);
      px[i]=gx+rx;py[i]=gy+ry;pz[i]=gz+rz;
      const sp=Math.sqrt(G*cm/Math.sqrt(r*r+SOFT));
      rot(-Math.sin(a)*sp,Math.cos(a)*sp,0,t[0],t[1],t[2],t[3]);
      vx[i]=bx+rx;vy[i]=by+ry;vz[i]=bz+rz;
      mass[i]=1+Math.random()*2;
      lists[(f<.33?0:f<.66?1:2)+(A?0:3)].push(i);
    }
  }
  cfg.follow=-1;cfg.tz=1;cfg.manual=false;cfg.dt=cfg.baseDT;minDist=9999;
  cfg.cx=0;cfg.cy=0;cfg.cz=0;cfg.yaw=.6;cfg.pitch=.55;cfg.userRot=false;
  frames=0;capKey='';impacted=false;shock=null;shake=0;simTime=0;hist.length=0;
  build();forces();
}
function step(){
  coreDist=Math.hypot(px[0]-px[CORE_B],py[0]-py[CORE_B],pz[0]-pz[CORE_B]);
  if(coreDist<minDist)minDist=coreDist;
  let t=cfg.baseDT;
  if(story.on&&story.beat>=0&&BEATS[story.beat].dt)t=BEATS[story.beat].dt;
  if(cfg.slow&&coreDist<280)t=.08;
  if(cfg.follow!==-1)t=Math.min(t,.14);
  cfg.dt+=(t-cfg.dt)*.06;
  const dt=cfg.dt,h=dt*.5;simTime+=dt;
  for(let i=0;i<N;i++){
    vx[i]+=ax[i]*h;vy[i]+=ay[i]*h;vz[i]+=az[i]*h;
    px[i]+=vx[i]*dt;py[i]+=vy[i]*dt;pz[i]+=vz[i]*dt;
  }
  build();forces();
  for(let i=0;i<N;i++){vx[i]+=ax[i]*h;vy[i]+=ay[i]*h;vz[i]+=az[i]*h}
}

/* ================= BACKDROP + SPRITES ================= */
const bg=document.getElementById('bg'),bctx=bg.getContext('2d');
const cv=document.getElementById('sim'),ctx=cv.getContext('2d');
function drawBg(){
  bg.width=innerWidth;bg.height=innerHeight;
  const w=bg.width,h=bg.height;
  bctx.fillStyle='#060712';bctx.fillRect(0,0,w,h);
  [[.15,.25,.55,'rgba(80,50,220,.20)'],[.85,.75,.6,'rgba(0,150,220,.14)'],[.65,.1,.35,'rgba(210,60,130,.10)']].forEach(([x,y,r,c])=>{
    const g=bctx.createRadialGradient(w*x,h*y,0,w*x,h*y,Math.max(w,h)*r);
    g.addColorStop(0,c);g.addColorStop(1,'rgba(0,0,0,0)');bctx.fillStyle=g;bctx.fillRect(0,0,w,h);
  });
  for(let i=0;i<450;i++){
    bctx.fillStyle='rgba(255,255,255,'+(Math.random()*.6+.1)+')';
    bctx.fillRect(Math.random()*w,Math.random()*h,Math.random()<.9?1:2,1.5);
  }
}
function resize(){cv.width=innerWidth;cv.height=innerHeight;drawBg()}
addEventListener('resize',resize);resize();
function glowSprite(c1,c2){
  const s=document.createElement('canvas');s.width=s.height=128;
  const g=s.getContext('2d'),r=g.createRadialGradient(64,64,0,64,64,64);
  r.addColorStop(0,'#fff');r.addColorStop(.12,c1);r.addColorStop(.4,c2);r.addColorStop(1,'rgba(0,0,0,0)');
  g.fillStyle=r;g.fillRect(0,0,128,128);return s;
}
const glowA=glowSprite('rgba(160,235,255,.95)','rgba(40,150,255,.35)'),glowB=glowSprite('rgba(255,225,150,.95)','rgba(255,90,90,.35)');
const COLORS=['#cfeeff','#4aa8ff','#a07bff','#fff1c2','#ffa24a','#ff5f9a'];
const DOTS=COLORS.map(c=>{
  const s=document.createElement('canvas');s.width=s.height=32;
  const g=s.getContext('2d'),r=g.createRadialGradient(16,16,0,16,16,16);
  r.addColorStop(0,'#fff');r.addColorStop(.28,c);r.addColorStop(1,'rgba(0,0,0,0)');
  g.fillStyle=r;g.fillRect(0,0,32,32);return s;
});

/* ================= 3D CAMERA / INPUT ================= */
const scaleOf=()=>Math.min(cv.width/900,cv.height/620)*cfg.zoom;
const PD=1100;
let cYa=1,sYa=0,cPi=1,sPi=0,PS=1,PW=0,PH=0,SX=0,SY=0,SK=0;
function setCam(s,W,H,ox,oy){
  cYa=Math.cos(cfg.yaw);sYa=Math.sin(cfg.yaw);cPi=Math.cos(cfg.pitch);sPi=Math.sin(cfg.pitch);
  PS=s;PW=W/2+ox;PH=H/2+oy;
}
function proj(x,y,z){
  const dx=x-cfg.cx,dy=y-cfg.cy,dz=z-cfg.cz;
  const x1=dx*cYa+dz*sYa,z1=-dx*sYa+dz*cYa;
  const y2=dy*cPi-z1*sPi,z2=dy*sPi+z1*cPi;
  const den=PD+z2;SK=den>40?PD/den:0;
  SX=PW+x1*SK*PS;SY=PH+y2*SK*PS;
}
let drag=null,moved=false;
cv.addEventListener('pointerdown',e=>{drag={sx:e.clientX,sy:e.clientY,x:e.clientX,y:e.clientY};moved=false});
addEventListener('pointermove',e=>{
  if(!drag)return;
  if(Math.abs(e.clientX-drag.sx)+Math.abs(e.clientY-drag.sy)>5)moved=true;
  if(moved&&!story.on){
    cfg.yaw+=(e.clientX-drag.x)*.006;
    cfg.pitch=Math.max(-1.45,Math.min(1.45,cfg.pitch+(e.clientY-drag.y)*.006));
    cfg.userRot=true;
  }
  drag.x=e.clientX;drag.y=e.clientY;
});
addEventListener('pointerup',()=>{drag=null});
cv.addEventListener('click',e=>{
  if(moved||story.on)return;
  $('set').classList.remove('open');$('learn').classList.remove('open');
  setCam(scaleOf(),cv.width,cv.height,0,0);
  let best=-1,bd=24*24;
  for(let i=0;i<N;i++){
    proj(px[i],py[i],pz[i]);if(!SK)continue;
    const d=(SX-e.clientX)*(SX-e.clientX)+(SY-e.clientY)*(SY-e.clientY);
    if(d<bd){bd=d;best=i}
  }
  cfg.follow=best;cfg.tz=best!==-1?2.2:1;cfg.manual=best!==-1;
});
addEventListener('wheel',e=>{if(story.on)return;cfg.manual=true;cfg.tz=Math.max(.4,Math.min(5,cfg.tz-e.deltaY*.001))},{passive:true});

/* ================= AUDIO: original cinematic score (no files) =================
   A low pulsing arpeggio over a deep pad and sub-bass, in the spirit of
   modern film scores. Procedurally generated, so nothing is downloaded. */
let ac=null,master,bus,arpFlt,padFlt,bassBus,on=false,timer=null,nextT=0,stepI=0,tn=0;
const CH=[[110,164.81,220,261.63,329.63],[87.31,130.81,174.61,220,261.63],[130.81,196,261.63,329.63,392],[98,146.83,196,246.94,293.66]];
const ARP=[0,1,2,3,4,3,2,1],STEP=.34;
function impulse(sec){
  const n=Math.floor(ac.sampleRate*sec),b=ac.createBuffer(2,n,ac.sampleRate);
  for(let c=0;c<2;c++){const d=b.getChannelData(c);for(let i=0;i<n;i++)d[i]=(Math.random()*2-1)*Math.pow(1-i/n,2.8)}
  return b;
}
function tone(f,t,dur,type,vol,dest,att){
  const o=ac.createOscillator(),g=ac.createGain();
  o.type=type;o.frequency.value=f;
  g.gain.setValueAtTime(0.0001,t);g.gain.linearRampToValueAtTime(vol,t+att);g.gain.exponentialRampToValueAtTime(.0001,t+dur);
  o.connect(g);g.connect(dest);o.start(t);o.stop(t+dur+.1);
}
function pad(f,t,dur){
  const o=ac.createOscillator(),g=ac.createGain();
  o.type='triangle';o.frequency.value=f;o.detune.value=(Math.random()-.5)*8;
  g.gain.setValueAtTime(0.0001,t);g.gain.linearRampToValueAtTime(.05,t+2);
  g.gain.setValueAtTime(.05,t+dur-2.2);g.gain.linearRampToValueAtTime(.0001,t+dur);
  o.connect(g);g.connect(padFlt);o.start(t);o.stop(t+dur+.1);
}
function sched(){
  while(nextT<ac.currentTime+.3){
    const ch=CH[(stepI>>4)%4],k=stepI%8;
    if(stepI%16===0){ch.slice(0,4).forEach(f=>pad(f,nextT,STEP*16+2));tone(ch[0]/2,nextT,3,'sine',.5,bassBus,.03)}
    if(k===0)tone(ch[0]/2,nextT,.9,'sine',.38,bassBus,.01);
    tone(ch[ARP[k]],nextT,.7,'triangle',.07+tn*.06,arpFlt,.015);
    stepI++;nextT+=STEP;
  }
}
function audio(state){
  const b=$('audio');
  if(state&&!on){
    if(!ac){
      ac=new(window.AudioContext||window.webkitAudioContext)();
      master=ac.createGain();master.gain.value=0;master.connect(ac.destination);
      bus=ac.createGain();bus.connect(master);
      const conv=ac.createConvolver();conv.buffer=impulse(3);
      const send=ac.createGain();send.gain.value=.6;bus.connect(send);send.connect(conv);conv.connect(master);
      arpFlt=ac.createBiquadFilter();arpFlt.type='lowpass';arpFlt.frequency.value=900;arpFlt.Q.value=.6;arpFlt.connect(bus);
      padFlt=ac.createBiquadFilter();padFlt.type='lowpass';padFlt.frequency.value=700;padFlt.connect(bus);
      bassBus=ac.createGain();bassBus.connect(master);
    }
    ac.resume();nextT=ac.currentTime+.1;timer=setInterval(sched,100);
    master.gain.setTargetAtTime(.3,ac.currentTime,1.5);
    on=true;b.classList.add('on');
  }else if(!state&&on){
    clearInterval(timer);master.gain.setTargetAtTime(0,ac.currentTime,.5);
    setTimeout(()=>{if(!on)ac.suspend()},2200);
    on=false;b.classList.remove('on');
  }
}
function tension(){
  if(!on)return;
  tn=Math.max(0,Math.min(1,1-coreDist/500));
  arpFlt.frequency.setTargetAtTime(800+tn*2200,ac.currentTime,.3);
  master.gain.setTargetAtTime(.28+tn*.1,ac.currentTime,.6);
}

/* ================= CAPTIONS ================= */
let capKey='';
function caption(){
  let k,t,p;
  if(cfg.follow!==-1){k='f'+cfg.follow;t='Following a star';p='It skips most stars by treating distant groups as one point (green boxes).'}
  else if(cfg.brute){k='b';t='Brute force';p='Every star checks every other star. Watch the frame rate drop.'}
  else if(frames<260&&coreDist>450){k='0';t='Two galaxies';p='Drag to rotate. Click any star to follow it.'}
  else if(coreDist>300&&minDist>300){k='1';t='Approach';p='Gravity pulls the two cores together.'}
  else if(coreDist>=120&&minDist>=120){k='2';t='First contact';p='The galaxies overlap and time slows down.'}
  else if(coreDist<120){k='3';t='Merger';p='The two cores spiral around each other.'}
  else{k='4';t='Tidal tails';p='Stars are flung out into long arcs.'}
  if(k===capKey)return;
  capKey=k;
  const el=$('capin');el.classList.add('h');
  setTimeout(()=>{el.innerHTML='<b>'+t+'</b>'+p;el.classList.remove('h')},400);
}

/* ================= UI ================= */
function $(id){return document.getElementById(id)}
const fmt=v=>(+v).toLocaleString();
function setStars(v,apply){
  N=+v;['stars','stars0'].forEach(id=>$(id).value=v);
  $('starv').textContent=$('starv0').textContent=fmt(v);
  if(apply)init();
}
['stars','stars0'].forEach(id=>{
  $(id).oninput=e=>{$('starv').textContent=$('starv0').textContent=fmt(e.target.value)};
  $(id).onchange=e=>setStars(e.target.value,true);
});
$('theta').oninput=e=>{cfg.theta=+e.target.value;$('thl').textContent=cfg.theta.toFixed(2)};
$('grid').onchange=e=>cfg.grid=e.target.checked;
$('slow').onchange=e=>cfg.slow=e.target.checked;
$('brute').onchange=e=>cfg.brute=e.target.checked;
$('audio').onclick=()=>audio(!on);
$('reset').onclick=init;
$('gear').onclick=()=>{$('learn').classList.remove('open');$('set').classList.toggle('open')};
$('help').onclick=()=>{$('set').classList.remove('open');updateLearn();$('learn').classList.toggle('open')};
$('whyAcc').onclick=e=>{e.preventDefault();$('set').classList.remove('open');updateLearn();$('learn').classList.add('open')};
function updateLearn(){
  $('l-n').textContent=fmt(N*(N-1));$('l-ops').textContent=fmt(checks);
  $('l-x').textContent=(N*(N-1)/Math.max(1,checks)).toFixed(0)+'x';
  $('l-acc').textContent=cfg.theta===0?'Right now θ = 0, so no box is ever merged. The result is exact.'
    :'Right now θ = '+cfg.theta.toFixed(2)+': a box 100 wide is merged only when it is more than '+Math.round(100/cfg.theta)+' away.';
}

/* ---------- intro story ---------- */
const story={on:false,beat:-1,t0:0};
const speedNow=()=>(N*(N-1)/Math.max(1,checks)).toFixed(0);
const BEATS=[
 {t:'These are two galaxies.',p:'Each dot is a star. The two bright points are their heavy cores.',min:5500,zoom:.85,dt:.05},
 {t:'Gravity pulls them together.',p:'Every star attracts every other star, so the galaxies fall toward each other.',min:4500,max:18000,until:()=>coreDist<300,zoom:1},
 {t:'Now they collide.',p:'Time slows down so you can watch the stars being thrown around.',min:6500,max:18000,until:()=>coreDist<140,zoom:1.5},
 {t:'Here is the hard part.',p:()=>fmt(N*(N-1))+' pairs of stars pull on each other every frame. Barnes-Hut groups faraway stars and does about '+speedNow()+'x less work.',min:7000,zoom:1.5},
 {t:'Now it is your turn.',p:'Drag to rotate, change the number of stars, adjust accuracy, or click a star to follow it.',min:5000,zoom:1.3}
];
function showBeat(){
  const b=BEATS[story.beat],el=$('story');
  el.classList.remove('show');
  setTimeout(()=>{
    if(!story.on)return;
    $('sh').textContent=b.t;$('sp').textContent=typeof b.p==='function'?b.p():b.p;el.classList.add('show');
  },650);
  cfg.tz=b.zoom;story.t0=performance.now();
  if(story.beat===BEATS.length-1)document.body.classList.add('run');
}
function nextBeat(){story.beat++;if(story.beat>=BEATS.length)endStory();else showBeat()}
function storyTick(now){
  const b=BEATS[story.beat],el=now-story.t0;
  if(el<b.min)return;
  if(!b.until||b.until()||el>=b.max)nextBeat();
}
function endStory(){
  story.on=false;$('story').classList.remove('show');
  document.body.classList.remove('story');document.body.classList.add('run');cfg.manual=false;
}
function playStory(){
  $('intro').classList.add('gone');document.body.classList.remove('run');document.body.classList.add('story');
  init();cfg.zoom=.85;story.on=true;story.beat=-1;
  if($('snd').checked)audio(true);
  nextBeat();
}
function skipToSim(){
  story.on=false;$('story').classList.remove('show');
  $('intro').classList.add('gone');document.body.classList.remove('story');document.body.classList.add('run');
  init();if($('snd').checked)audio(true);
}
function leaveToIntro(){
  story.on=false;$('story').classList.remove('show');
  document.body.classList.remove('story','run');$('intro').classList.remove('gone');
  $('set').classList.remove('open');$('learn').classList.remove('open');
}
$('begin').onclick=playStory;
$('skipIntro').onclick=skipToSim;$('skip').onclick=skipToSim;
$('brand').onclick=leaveToIntro;

/* ================= LOOP ================= */
init();
let last=performance.now(),fpsAvg=60,physMs=0;
$('lhead').onclick=()=>$('live').classList.toggle('min');
function drawSpark(){
  const c=$('spark'),g=c.getContext('2d'),w=c.width,h=c.height;
  g.clearRect(0,0,w,h);if(hist.length<2)return;
  g.beginPath();
  hist.forEach((d,i)=>{const x=i/119*w,y=h-2-Math.min(1,d/600)*(h-4);i?g.lineTo(x,y):g.moveTo(x,y)});
  g.strokeStyle='#5ce1ff';g.lineWidth=1.5;g.stroke();
}
function updateLive(){
  $('v-time').textContent=simTime.toFixed(0)+' t';
  const r=cfg.dt/cfg.baseDT;$('v-ts').textContent=r.toFixed(2)+'x'+(r<.6?' slow-mo':'');
  $('v-dist').textContent=coreDist.toFixed(0)+' u';
  const dx=px[0]-px[CORE_B],dy=py[0]-py[CORE_B],dz=pz[0]-pz[CORE_B];
  const cl=-((vx[0]-vx[CORE_B])*dx+(vy[0]-vy[CORE_B])*dy+(vz[0]-vz[CORE_B])*dz)/Math.max(1e-6,coreDist);
  $('v-close').textContent=(cl>=0?'closing ':'separating ')+Math.abs(cl).toFixed(2)+' u/t';
  $('v-ops').textContent=fmt(checks);$('v-naive').textContent=fmt(N*(N-1));
  $('v-save').textContent=(N*(N-1)/Math.max(1,checks)).toFixed(0)+'x';
  $('v-ms').textContent=physMs.toFixed(1)+' ms';$('v-fps').textContent=Math.round(fpsAvg)+' fps';
  hist.push(coreDist);if(hist.length>120)hist.shift();drawSpark();
  const f=cfg.follow,box=$('v-follow');
  if(f===-1){box.style.display='none';return}
  box.style.display='block';$('v-fid').textContent='#'+f;
  $('v-fspd').textContent=Math.hypot(vx[f],vy[f],vz[f]).toFixed(2)+' u/t';
  $('v-ffor').textContent=(mass[f]*Math.hypot(ax[f],ay[f],az[f])).toFixed(3)+' m·u/t²';
  $('v-fdist').textContent=Math.min(Math.hypot(px[f]-px[0],py[f]-py[0],pz[f]-pz[0]),Math.hypot(px[f]-px[CORE_B],py[f]-py[CORE_B],pz[f]-pz[CORE_B])).toFixed(0)+' u';
  $('v-fgrp').textContent=cfg.brute?'none':fmt(approx.length);
}
const EDGES=new Uint8Array([0,1,2,3,4,5,6,7,0,2,1,3,4,6,5,7,0,4,1,5,2,6,3,7]);
const CXs=new Float32Array(8),CYs=new Float32Array(8);
function cube(n){
  const h=nHalf[n];
  for(let c=0;c<8;c++){
    proj(nCx[n]+(c&1?h:-h),nCy[n]+(c&2?h:-h),nCz[n]+(c&4?h:-h));
    if(!SK)return;CXs[c]=SX;CYs[c]=SY;
  }
  for(let e=0;e<24;e+=2){ctx.moveTo(CXs[EDGES[e]],CYs[EDGES[e]]);ctx.lineTo(CXs[EDGES[e+1]],CYs[EDGES[e+1]])}
}
function loop(now){
  fpsAvg+=(1000/Math.max(1,now-last)-fpsAvg)*.1;last=now;
  const tp=performance.now();step();physMs=performance.now()-tp;
  frames++;tension();if(story.on)storyTick(performance.now());

  if(cfg.follow!==-1){
    const f=cfg.follow;cfg.cx+=(px[f]-cfg.cx)*.12;cfg.cy+=(py[f]-cfg.cy)*.12;cfg.cz+=(pz[f]-cfg.cz)*.12;
  }else{
    cfg.cx+=((px[0]+px[CORE_B])/2-cfg.cx)*.04;
    cfg.cy+=((py[0]+py[CORE_B])/2-cfg.cy)*.04;
    cfg.cz+=((pz[0]+pz[CORE_B])/2-cfg.cz)*.04;
    if(!cfg.manual&&!story.on)cfg.tz=coreDist<300?1.5:1;
  }
  cfg.zoom+=(cfg.tz-cfg.zoom)*.06;
  if(!(drag&&moved)&&!RM)cfg.yaw+=story.on?.0025:(cfg.userRot?0:.0009);

  const s=scaleOf(),W=cv.width,H=cv.height;
  if(!impacted&&coreDist<130){
    impacted=true;shake=RM?0:9;
    shock={x:(px[0]+px[CORE_B])/2,y:(py[0]+py[CORE_B])/2,z:(pz[0]+pz[CORE_B])/2,r:0};
    const f=$('flash');f.style.animation='none';void f.offsetWidth;f.style.animation='flash 1.4s ease-out';
  }
  const shx=(Math.random()-.5)*shake,shy=(Math.random()-.5)*shake;shake*=.93;
  setCam(s,W,H,shx,shy);

  ctx.globalCompositeOperation='destination-out';
  ctx.globalAlpha=1;ctx.fillStyle='rgba(0,0,0,.3)';ctx.fillRect(0,0,W,H);
  ctx.globalCompositeOperation='lighter';

  if(cfg.grid&&!cfg.brute){
    ctx.strokeStyle='rgba(120,160,255,.22)';ctx.lineWidth=1;ctx.beginPath();
    let drawn=0;
    for(let i=0;i<nodeCount&&drawn<700;i++)if(nPart[i]===-2&&nHalf[i]*2*s>10){cube(i);drawn++}
    ctx.stroke();
  }
  if(cfg.follow!==-1&&cfg.approx&&!cfg.brute){
    const f=cfg.follow;proj(px[f],py[f],pz[f]);const sx0=SX,sy0=SY;
    const cap=Math.min(approx.length,600);
    ctx.strokeStyle='rgba(109,255,168,.55)';ctx.lineWidth=1.2;ctx.beginPath();
    for(let k=0;k<cap;k++)cube(approx[k]);
    ctx.stroke();
    ctx.strokeStyle='rgba(109,255,168,.2)';ctx.beginPath();
    for(let k=0;k<cap;k++){const n=approx[k];proj(nComX[n],nComY[n],nComZ[n]);if(SK){ctx.moveTo(sx0,sy0);ctx.lineTo(SX,SY)}}
    ctx.stroke();
  }
  if(shock){
    shock.r+=5;const a=Math.max(0,1-shock.r/650);
    if(a<=0)shock=null;
    else{
      proj(shock.x,shock.y,shock.z);
      if(SK){
        const R=shock.r*SK*s,cx=SX,cy=SY;
        ctx.lineWidth=4+a*6;ctx.strokeStyle='rgba(150,190,255,'+(a*.55)+')';
        ctx.beginPath();ctx.arc(cx,cy,R,0,6.2832);ctx.stroke();
        if(shock.r>60){ctx.lineWidth=2;ctx.strokeStyle='rgba(255,120,190,'+(a*.4)+')';
          ctx.beginPath();ctx.arc(cx,cy,(shock.r-60)*SK*s,0,6.2832);ctx.stroke()}
      }
    }
  }

  // soft glowing stars, dimmer and smaller when far from the camera
  const base=4.4*Math.sqrt(s)*Math.min(1.8,Math.sqrt(4000/N));
  for(let b=0;b<6;b++){
    const L=lists[b],sp=DOTS[b];
    for(let k=0;k<L.length;k++){
      const i=L[k];proj(px[i],py[i],pz[i]);if(!SK)continue;
      const z=Math.min(28,Math.max(3,base*SK));
      ctx.globalAlpha=Math.max(.2,Math.min(1,.5+(SK-.75)*1.3));
      ctx.drawImage(sp,SX-z/2,SY-z/2,z,z);
    }
  }
  ctx.globalAlpha=1;
  const gs=170*Math.sqrt(s)*(1+.07*Math.sin(frames*.06));
  proj(px[0],py[0],pz[0]);if(SK){const g=gs*SK;ctx.drawImage(glowA,SX-g/2,SY-g/2,g,g)}
  proj(px[CORE_B],py[CORE_B],pz[CORE_B]);if(SK){const g=gs*SK;ctx.drawImage(glowB,SX-g/2,SY-g/2,g,g)}

  if(cfg.follow!==-1){
    const f=cfg.follow;proj(px[f],py[f],pz[f]);
    if(SK){
      ctx.globalCompositeOperation='source-over';
      ctx.strokeStyle='#ff5c8a';ctx.lineWidth=2;ctx.beginPath();
      ctx.arc(SX,SY,9+2*Math.sin(frames*.15),0,6.2832);ctx.stroke();
      ctx.strokeStyle='rgba(255,92,138,'+(.5-.5*((frames%50)/50))+')';ctx.beginPath();
      ctx.arc(SX,SY,9+(frames%50)*.6,0,6.2832);ctx.stroke();
    }
  }

  if(frames%8===0){
    updateLive();caption();
    if($('learn').classList.contains('open'))updateLearn();
  }
  requestAnimationFrame(loop);
}
requestAnimationFrame(loop);
