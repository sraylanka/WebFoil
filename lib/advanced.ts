import { validateGeometry } from './geometry';
import { linear, naca, solve, type Point } from './foil';

export type Element = { points: Point[] };
export type Surface = Point & { cp: number; velocity: number; signedVelocity: number; side: string; element: number };
type Panel = { a: Point; b: Point; x: number; y: number; l: number; tx: number; ty: number; nx: number; ny: number; element: number; side: string };
export type Flow = { cl: number; cd: number; cm: number; data: Surface[]; panels: Panel[]; strengths: number[]; alpha: number; ground?: number; flux: number[] };

// Exact constant-source line integral. Vortex influence is the source vector rotated +90°.
function influence(p: Panel, x: number, y: number): Point {
  const dx = x - p.a.x, dy = y - p.a.y, X = dx * p.tx + dy * p.ty, Y = -dx * p.ty + dy * p.tx;
  const r1 = Math.max(X * X + Y * Y, 1e-24), r2 = Math.max((X - p.l) ** 2 + Y * Y, 1e-24);
  const u = Math.log(r1 / r2) / (4 * Math.PI);
  let angle = Math.atan2(Y, X - p.l) - Math.atan2(Y, X);
  if (angle > Math.PI) angle -= 2 * Math.PI;
  if (angle < -Math.PI) angle += 2 * Math.PI;
  const v = angle / (2 * Math.PI);
  return { x: u * p.tx - v * p.ty, y: u * p.ty + v * p.tx };
}
function image(p: Panel): Panel {
  return { ...p, a: { x: p.a.x, y: -p.a.y }, b: { x: p.b.x, y: -p.b.y }, y: -p.y, ty: -p.ty, ny: -p.ny };
}
function vectors(p: Panel, x: number, y: number, ground: boolean, self = false) {
  let s = self ? { x: .5 * p.nx, y: .5 * p.ny } : influence(p, x, y);
  let v = self ? { x: .5 * p.tx, y: .5 * p.ty } : { x: -s.y, y: s.x };
  if (ground) {
    const q = influence(image(p), x, y);
    s = { x: s.x + q.x, y: s.y + q.y };
    // Mirroring preserves source strength and reverses circulation.
    v = { x: v.x + q.y, y: v.y - q.x };
  }
  return { s, v };
}

export function solveElements(elements: Element[], alpha: number, mach = 0, ground?: number): Flow {
  if (!elements.length || elements.length > 3) throw Error('Use one to three elements.');
  if (!Number.isFinite(alpha) || Math.abs(alpha) > 20 || !Number.isFinite(mach) || mach < 0 || mach > .5) throw Error('Use |α| ≤ 20° and 0 ≤ Mach ≤ 0.5.');
  if (ground !== undefined && (!Number.isFinite(ground) || ground <= 0)) throw Error('Ground clearance must be positive.');
  const ps: Panel[] = [], ranges: number[][] = [];
  for (const [element, { points }] of elements.entries()) {
    validateGeometry(points);
    const le = points.reduce((best, p, i) => p.x < points[best].x ? i : best, 0), indices: number[] = [];
    for (let i = 0; i < points.length - 1; i++) {
      const a = { ...points[i], y: points[i].y + (ground ?? 0) }, b = { ...points[i + 1], y: points[i + 1].y + (ground ?? 0) };
      const l = Math.hypot(b.x - a.x, b.y - a.y);
      if (l < 1e-8 || (ground !== undefined && Math.min(a.y, b.y) <= 0)) throw Error('Duplicate points or geometry intersects the ground.');
      const tx = (b.x - a.x) / l, ty = (b.y - a.y) / l;
      indices.push(ps.length); ps.push({ a, b, l, tx, ty, nx: ty, ny: -tx, x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, element, side: i < le ? 'upper' : 'lower' });
    }
    ranges.push(indices);
  }
  const enclosed=(p:Point,contour:Point[])=>{let found=false;for(let i=0;i<contour.length;i++){const a=contour[i],b=contour[(i+1)%contour.length];if((a.y>p.y)!==(b.y>p.y)&&p.x<(b.x-a.x)*(p.y-a.y)/(b.y-a.y)+a.x)found=!found;}return found;};
  for(let i=0;i<elements.length;i++)for(let j=0;j<elements.length;j++)if(i!==j&&enclosed(elements[i].points[0],elements[j].points))throw Error('Elements cannot overlap or be nested.');
  if (ps.length > 500) throw Error('Limit the combined geometry to 500 panels.');
  // Reject intersecting elements before solving the coupled system.
  const cross = (a: Point, b: Point, c: Point) => (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
  for (let i=0;i<ps.length;i++) for(let j=i+1;j<ps.length;j++) if(ps[i].element!==ps[j].element) {
    const p=ps[i],q=ps[j];
    if(cross(p.a,p.b,q.a)*cross(p.a,p.b,q.b)<=0 && cross(q.a,q.b,p.a)*cross(q.a,q.b,p.b)<=0 && Math.max(p.a.x,p.b.x)>=Math.min(q.a.x,q.b.x) && Math.max(q.a.x,q.b.x)>=Math.min(p.a.x,p.b.x) && Math.max(p.a.y,p.b.y)>=Math.min(q.a.y,q.b.y) && Math.max(q.a.y,q.b.y)>=Math.min(p.a.y,p.b.y)) throw Error('Airfoil elements intersect.');
  }
  const N=ps.length, E=elements.length, a=alpha*Math.PI/180;
  // Ground flow must be parallel to the wall. Rotate geometry to set incidence instead.
  if (ground !== undefined && Math.abs(alpha)>1e-12) throw Error('For ground effect use α = 0 and rotate the geometry to set incidence.');
  const cs=Math.cos(a),sn=Math.sin(a),normal:number[][]=[],tangent:number[][]=[];
  ps.forEach((p,i)=>{
    normal[i]=Array(N+E).fill(0);tangent[i]=Array(N+E).fill(0);
    ps.forEach((q,j)=>{const {s,v}=vectors(q,p.x,p.y,ground!==undefined,i===j);
      normal[i][j]=s.x*p.nx+s.y*p.ny;tangent[i][j]=s.x*p.tx+s.y*p.ty;
      normal[i][N+q.element]+=v.x*p.nx+v.y*p.ny;tangent[i][N+q.element]+=v.x*p.tx+v.y*p.ty;
    });
  });
  const A=normal.map(r=>r.slice()),b=ps.map(p=>-cs*p.nx-sn*p.ny);
  for(const indices of ranges){const first=indices[0],last=indices.at(-1)!;A.push(tangent[first].map((v,j)=>v+tangent[last][j]));b.push(-cs*(ps[first].tx+ps[last].tx)-sn*(ps[first].ty+ps[last].ty));}
  const strengths=linear(A,b);let fx=0,fy=0,cm=0;
  const data=ps.map((p,i)=>{const vt=cs*p.tx+sn*p.ty+tangent[i].reduce((sum,v,j)=>sum+v*strengths[j],0),cp=(1-vt*vt)/Math.sqrt(1-mach*mach),dx=-cp*p.nx*p.l,dy=-cp*p.ny*p.l;fx+=dx;fy+=dy;cm-=((p.x-.25)*dy-(p.y-(ground??0))*dx);return {x:p.x,y:p.y,cp,velocity:Math.abs(vt),signedVelocity:vt,side:p.side,element:p.element};});
  return {cl:-fx*sn+fy*cs,cd:fx*cs+fy*sn,cm,data,panels:ps,strengths,alpha,ground,flux:ranges.map(indices=>indices.reduce((s,i)=>s+strengths[i]*ps[i].l,0))};
}

export function velocityAt(flow: Flow, point: Point): Point {
  const a=flow.alpha*Math.PI/180;let x=Math.cos(a),y=Math.sin(a);
  for(const [j,p] of flow.panels.entries()){const {s,v}=vectors(p,point.x,point.y,flow.ground!==undefined);x+=s.x*flow.strengths[j]+v.x*flow.strengths[flow.panels.length+p.element];y+=s.y*flow.strengths[j]+v.y*flow.strengths[flow.panels.length+p.element];}
  return {x,y};
}
function inside(point:Point,panels:Panel[]){let yes=false;for(const p of panels){if((p.a.y>point.y)!==(p.b.y>point.y)&&point.x<(p.b.x-p.a.x)*(point.y-p.a.y)/(p.b.y-p.a.y)+p.a.x)yes=!yes;}return yes;}
export function streamlines(flow:Flow,seeds:Point[],steps=400,step=.012):Point[][]{
  if(!Number.isInteger(steps)||steps<1||steps>2000||!Number.isFinite(step)||step<=0||step>.1||seeds.length>80)throw Error('Invalid streamline limits.');
  return seeds.map(seed=>{const line:Point[]=[seed];for(let i=0;i<steps;i++){const p=line.at(-1)!;if(inside(p,flow.panels)||(flow.ground!==undefined&&p.y<=0)||Math.abs(p.x)>4||Math.abs(p.y)>3)break;const v=velocityAt(flow,p),norm=Math.hypot(v.x,v.y);if(!Number.isFinite(norm)||norm<1e-6)break;const mid={x:p.x+step*v.x/norm/2,y:p.y+step*v.y/norm/2},w=velocityAt(flow,mid),n=Math.hypot(w.x,w.y);if(!Number.isFinite(n)||n<1e-6)break;const q={x:p.x+step*w.x/n,y:p.y+step*w.y/n};if(inside(q,flow.panels))break;line.push(q);}return line;});
}

export type BLStation = { x:number; s:number; ue:number; theta:number; displacement:number; cf:number; shape:number; regime:'laminar'|'turbulent'|'separated' };
export type BLSide = { side:string; stations:BLStation[]; transition:number|null; separation:number|null; drag:number };
// One-way integral boundary layer: Thwaites laminar + Michel transition + Head entrainment.
// No XFOIL source is incorporated. This is NOT XFOIL's coupled e^N formulation.
export function boundaryLayer(flow:Flow,reynolds:number,trips={upper:1,lower:1}) {
  if(flow.flux.length!==1)throw Error('Boundary-layer estimates currently require one element.');
  if(!Number.isFinite(reynolds)||reynolds<1e4||reynolds>1e8)throw Error('Reynolds number must be 10,000–100,000,000.');
  if([trips.upper,trips.lower].some(x=>!Number.isFinite(x)||x<0||x>1))throw Error('Trips must be between 0 and 1 chord.');
  // Find the stagnation panel from the signed tangential-velocity reversal.
  let stag=-1;
  for(let i=0;i<flow.data.length-1;i++)if(flow.data[i].signedVelocity<=0&&flow.data[i+1].signedVelocity>=0){stag=i;break;}
  if(stag<0)throw Error('No leading-edge stagnation point found; boundary-layer estimate unavailable.');
  const left=flow.data[stag],right=flow.data[stag+1],fraction=Math.abs(left.signedVelocity)/(Math.abs(left.signedVelocity)+Math.abs(right.signedVelocity));
  const origin={x:left.x+fraction*(right.x-left.x),y:left.y+fraction*(right.y-left.y)};
  const outputs:BLSide[]=[];
  for(const side of ['upper','lower'] as const){
    const samples=side==='upper'?flow.data.slice(0,stag+1).reverse():flow.data.slice(stag+1);
    const stations:BLStation[]=[];let last=origin,s=0,integral=0,previousU=0,theta=0,H=2.59,turbulent=false,separated=false,transition:number|null=null,separation:number|null=null;
    for(const p of samples){let justTransitioned=false;const ds=Math.hypot(p.x-last.x,p.y-last.y),u=Math.max(p.velocity,.02),du=(u-previousU)/Math.max(ds,1e-8);s+=ds;
      if(!turbulent){integral+=ds/6*Array.from({length:6},(_,k)=>previousU**k*u**(5-k)).reduce((sum,v)=>sum+v,0);theta=Math.sqrt(.45/reynolds*integral/u**6);const lambda=reynolds*theta*theta*du;
        H=lambda<0?2.088+.0731/Math.max(lambda+.14,.015):2.61-3.75*lambda+5.24*lambda*lambda;
        const ret=reynolds*u*theta,res=reynolds*u*s,critical=1.174*(1+22400/Math.max(res,1))*Math.max(res,1)**.46;
        if((trips[side]<1&&p.x>=trips[side])||ret>critical){turbulent=true;justTransitioned=true;transition=p.x;H=1.4;}
        else if(lambda<=-.09){separated=true;separation=p.x;}
      }
      let cf=0;
      if(turbulent&&!separated&&!justTransitioned){
        // Head's entrainment relation, integrated in small Euler substeps.
        const sub=Math.max(1,Math.ceil(ds/.002)),delta=ds/sub;
        for(let k=0;k<sub;k++){const ue=Math.max(previousU+du*(k+.5)*delta,.02),ret=Math.max(reynolds*ue*theta,1);cf=.246*10**(-.678*H)*ret**(-.268);
          const h1=H<=1.6?3.3+.8234/Math.max(H-1.1,.01)**1.287:3.3+1.5501/Math.max(H-.6778,.01)**3.064;
          const entrainment=.0306/Math.max(h1-3,.01)**.6169,dtheta=cf/2-(H+2)*theta*du/ue;
          const nextTheta=Math.max(theta+delta*dtheta,1e-9),nextH1=Math.max((ue*theta*h1+delta*ue*entrainment)/(Math.max(ue+du*delta,.02)*nextTheta),3.31);
          H=nextH1>=5.3?1.1+(.8234/(nextH1-3.3))**(1/1.287):.6778+(1.5501/(nextH1-3.3))**(1/3.064);theta=nextTheta;
          if(H>=2.4){separated=true;separation=p.x;break;}
        }
      }else if(justTransitioned){cf=.246*10**(-.678*H)*Math.max(reynolds*u*theta,1)**(-.268);}else if(!separated){const lambda=reynolds*theta*theta*du,l=lambda<0?.22+1.402*lambda+.018*lambda/Math.max(lambda+.107,.001):.22+1.57*lambda-1.8*lambda*lambda;cf=2*Math.max(l,0)/Math.max(reynolds*u*theta,1);}
      stations.push({x:p.x,s,ue:u,theta,displacement:H*theta,cf:separated?0:cf,shape:H,regime:separated?'separated':turbulent?'turbulent':'laminar'});
      last=p;previousU=u;if(separated)break;
    }
    const tail=stations.at(-1)!;
    // Squire–Young extrapolation to wake infinity, only for attached flow.
    outputs.push({side,stations,transition,separation,drag:separation===null?2*tail.theta*tail.ue**((tail.shape+5)/2):NaN});
  }
  const attached=outputs.every(s=>s.separation===null);
  return {sides:outputs,attached,cd:attached?outputs.reduce((sum,s)=>sum+s.drag,0):null,warning:attached?'Uncoupled integral estimate; lift and moment remain inviscid.':'Separation detected: drag estimate withheld. No stall correction or post-stall prediction.'};
}

export function transform(points:Point[],scale:number,angle:number,x:number,y:number):Point[]{
  if(![scale,angle,x,y].every(Number.isFinite)||scale<=0)throw Error('Invalid geometry transform.');
  const a=angle*Math.PI/180;return points.map(p=>({x:x+scale*(p.x*Math.cos(a)-p.y*Math.sin(a)),y:y+scale*(p.x*Math.sin(a)+p.y*Math.cos(a))}));
}
export function flap(points:Point[],hinge:number,angle:number):Point[]{
  if(!Number.isFinite(hinge)||hinge<.5||hinge>.95||!Number.isFinite(angle)||Math.abs(angle)>40)throw Error('Use a hinge at 0.5–0.95 chord and deflection within ±40°.');
  const le=points.reduce((best,p,i)=>p.x<points[best].x?i:best,0);
  const height=(surface:Point[])=>{for(let i=0;i<surface.length-1;i++){const a=surface[i],b=surface[i+1];if((a.x-hinge)*(b.x-hinge)<=0&&a.x!==b.x)return a.y+(b.y-a.y)*(hinge-a.x)/(b.x-a.x);}throw Error('Hinge is outside the section.');};
  const hy=(height(points.slice(0,le+1))+height(points.slice(le)))/2,a=-angle*Math.PI/180;
  // Smooth hinge blending avoids opening or duplicating the contour.
  return points.map(p=>{const f=Math.min(1,Math.max(0,(p.x-hinge)/.04)),t=a*f;return {x:hinge+(p.x-hinge)*Math.cos(t)-(p.y-hy)*Math.sin(t),y:hy+(p.x-hinge)*Math.sin(t)+(p.y-hy)*Math.cos(t)};});
}
export function smooth(points:Point[],iterations=2):Point[]{
  validateGeometry(points);if(!Number.isInteger(iterations)||iterations<1||iterations>20)throw Error('Use 1–20 smoothing passes.');
  let out=points.map(p=>({...p}));const le=points.reduce((best,p,i)=>p.x<points[best].x?i:best,0);
  for(let k=0;k<iterations;k++)out=out.map((p,i)=>i===0||i===le||i===out.length-1?p:{x:p.x,y:.25*out[i-1].y+.5*p.y+.25*out[i+1].y});return out;
}
export function naca5(code:string,n=100):Point[]{
  if(!/^[1-5][1-5]0\d{2}$/.test(code))throw Error('Supported five-digit sections: standard non-reflex families, e.g. 23012.');
  const table=[[.058,361.4],[.126,51.64],[.2025,15.957],[.29,6.643],[.391,3.23]], [m,k]=table[+code[1]-1],k1=k*(+code[0]/2);
  const base=naca('00'+code.slice(3),n);return base.map(p=>{const x=p.x,yc=x<m?k1/6*(x**3-3*m*x*x+m*m*(3-m)*x):k1*m**3/6*(1-x),dy=x<m?k1/6*(3*x*x-6*m*x+m*m*(3-m)):-k1*m**3/6,angle=Math.atan(dy);return {x:x-p.y*Math.sin(angle),y:yc+p.y*Math.cos(angle)};});
}
export function cst(upper:number[],lower:number[],n=100):Point[]{
  if(!Number.isInteger(n)||n<20||n>300||n%2||upper.length!==lower.length||upper.length<2||upper.length>9||[...upper,...lower].some(v=>!Number.isFinite(v)))throw Error('CST needs matching 2–9 finite coefficients and 20–300 even panels.');
  const degree=upper.length-1,choose=(n:number,k:number)=>{let v=1;for(let i=1;i<=k;i++)v=v*(n-i+1)/i;return v;};
  const side=(co:number[])=>Array.from({length:n/2+1},(_,i)=>{const x=(1-Math.cos(Math.PI*i/(n/2)))/2,y=Math.sqrt(x)*(1-x)*co.reduce((s,c,j)=>s+c*choose(degree,j)*x**j*(1-x)**(degree-j),0);return {x,y};});
  const points=[...side(upper).reverse(),...side(lower).slice(1)];validateGeometry(points);return points;
}
export function designLift(points:Point[],target:number,alpha:number){
  if(!Number.isFinite(target)||Math.abs(target)>1.5)throw Error('Target CL must be within ±1.5.');
  let camber=0,geometry=points,value=solve(points,alpha).cl;
  for(let i=0;i<12;i++){if(Math.abs(value-target)<.001)return {points:geometry,cl:value,camber,iterations:i};const bump=(c:number)=>points.map(p=>({x:p.x,y:p.y+4*c*p.x*(1-p.x)})),eps=.0001,derivative=(solve(bump(camber+eps),alpha).cl-value)/eps;if(Math.abs(derivative)<1e-5)throw Error('Design sensitivity is too small.');camber+=(target-value)/derivative;if(Math.abs(camber)>.08)throw Error('Target requires over 8% additional camber; choose a closer target.');geometry=bump(camber);value=solve(geometry,alpha).cl;}
  throw Error('Lift-target design did not converge.');
}
export type AnalysisJob={elements:Element[];alpha:number;mach:number;ground?:number;reynolds?:number;trips?:{upper:number;lower:number};flowfield?:boolean};
export function analyze(job:AnalysisJob){const flow=solveElements(job.elements,job.alpha,job.mach,job.ground);const bl=job.reynolds!==undefined?boundaryLayer(flow,job.reynolds,job.trips):null;const lines=job.flowfield?streamlines(flow,Array.from({length:23},(_,i)=>({x:-.6,y:job.ground!==undefined?.015+i*.075:-.7+i*.065}))):[];return {job,flow,bl,lines};}
export function runBatch(jobs:AnalysisJob[]){if(!Array.isArray(jobs)||!jobs.length||jobs.length>81)throw Error('Use 1–81 JSON analysis jobs.');return jobs.map((job,index)=>{try{const r=analyze({...job,flowfield:false});return {index,alpha:job.alpha,reynolds:job.reynolds,cl:r.flow.cl,cm:r.flow.cm,pressureDrag:r.flow.cd,viscousDrag:r.bl?.cd??null,attached:r.bl?.attached??null,error:null};}catch(e){return {index,error:(e as Error).message};}});}
