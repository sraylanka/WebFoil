import { validateGeometry } from './geometry';
export type Point={x:number,y:number};
export function naca(code:string,n=100):Point[]{if(!Number.isInteger(n)||n<20||n>300||n%2)throw Error('Use 20–300 even panels.');if(!/^\d{4}$/.test(code))throw Error('Enter a four-digit NACA code.');if(+code.slice(2)===0)throw Error('Thickness must be positive.');const m=+code[0]/100,p=+code[1]/10,t=+code.slice(2)/100;if(m&& !p)throw Error('Camber location must be nonzero.');const upper:Point[]=[],lower:Point[]=[];for(let i=0;i<=n/2;i++){const x=(1-Math.cos(Math.PI*i/(n/2)))/2;const yt=x===1?0:5*t*(.2969*Math.sqrt(x)-.126*x-.3516*x*x+.2843*x**3-.1036*x**4);let yc=0,d=0;if(m){if(x<p){yc=m/p**2*(2*p*x-x*x);d=2*m/p**2*(p-x)}else{yc=m/(1-p)**2*(1-2*p+2*p*x-x*x);d=2*m/(1-p)**2*(p-x)}}const a=Math.atan(d);upper.push({x:x-yt*Math.sin(a),y:yc+yt*Math.cos(a)});lower.push({x:x+yt*Math.sin(a),y:yc-yt*Math.cos(a)})}return [...upper.reverse(),...lower.slice(1)]}
export function parseDat(text:string):Point[]{
  const lines=text.split(/\r?\n/).map(l=>l.trim()).filter(Boolean);const points:Point[]=[];
  for(const [index,line] of lines.entries()){
    const fields=line.split(/[\s,]+/),values=fields.map(Number);
    if(fields.length!==2||!values.every(Number.isFinite)){
      if(index===0&&/[a-zA-Z]/.test(line))continue;
      throw Error('Invalid coordinate row. Use an optional title then two finite numbers per row.');
    }
    if(values.some(v=>Math.abs(v)>900))throw Error('Multi-element separator files are not supported by DAT import.');
    points.push({x:values[0],y:values[1]});
  }
  if(points.length<20||points.length>301)throw Error('Provide 20–301 ordered surface coordinates.');
  const xmin=Math.min(...points.map(p=>p.x)),chord=Math.max(...points.map(p=>p.x))-xmin;
  if(chord<=0)throw Error('Chord must be positive.');
  let out=points.map(p=>({x:(p.x-xmin)/chord,y:p.y/chord}));
  if(out[0].x<.9||out.at(-1)!.x<.9)throw Error('Coordinates must start and end at the trailing edge.');
  if(Math.hypot(out[0].x-out.at(-1)!.x,out[0].y-out.at(-1)!.y)>1e-6)out.push({...out[0]});
  out[out.length-1]={...out[0]};
  if(out.length>301)throw Error('Use at most 300 panels after trailing-edge closure.');
  let area=0;for(let i=0;i<out.length-1;i++)area+=out[i].x*out[i+1].y-out[i+1].x*out[i].y;
  if(area<0)out=out.reverse();
  if(Math.abs(area)<1e-10)throw Error('Geometry has zero enclosed area.');
  for(let i=0;i<out.length-1;i++)if(Math.hypot(out[i+1].x-out[i].x,out[i+1].y-out[i].y)<1e-8)throw Error('Geometry has duplicate adjacent points.');
  validateGeometry(out);return out;
}
export function linear(A:number[][],b:number[]){if(!b.length||A.length!==b.length||A.some(r=>r.length!==b.length||r.some(v=>!Number.isFinite(v)))||b.some(v=>!Number.isFinite(v)))throw Error('Invalid linear system.');const original=A;const n=b.length;A=A.map((r,i)=>[...r,b[i]]);for(let k=0;k<n;k++){let piv=k;for(let i=k+1;i<n;i++)if(Math.abs(A[i][k])>Math.abs(A[piv][k]))piv=i;[A[k],A[piv]]=[A[piv],A[k]];if(Math.abs(A[k][k])<1e-12)throw Error('Singular geometry: check duplicate points.');for(let i=k+1;i<n;i++){const f=A[i][k]/A[k][k];for(let j=k;j<=n;j++)A[i][j]-=f*A[k][j]}}const x=Array(n).fill(0);for(let i=n-1;i>=0;i--){let s=A[i][n];for(let j=i+1;j<n;j++)s-=A[i][j]*x[j];x[i]=s/A[i][i]}if(x.some(v=>!Number.isFinite(v)))throw Error('Numerical solve produced non-finite results.');for(let i=0;i<n;i++){const residual=Math.abs(original[i].reduce((s,v,j)=>s+v*x[j],0)-b[i]),scale=1+Math.abs(b[i])+original[i].reduce((s,v,j)=>s+Math.abs(v*x[j]),0);if(residual/scale>1e-8)throw Error('Linear solve residual is too large.');}return x}
export function solve(points:Point[],alpha:number,mach=0){validateGeometry(points);if(points.length<20||!Number.isFinite(alpha)||Math.abs(alpha)>20)throw Error('Provide valid geometry and angle between −20° and 20°.');if(!Number.isFinite(mach)||mach<0||mach>.5)throw Error('Mach must be between 0 and 0.5.');const ps=points.slice(0,-1).map((a,i)=>{const b=points[i+1],l=Math.hypot(b.x-a.x,b.y-a.y);if(l<1e-8)throw Error('Geometry has duplicate adjacent points.');const tx=(b.x-a.x)/l,ty=(b.y-a.y)/l;return {a,b,l,tx,ty,nx:ty,ny:-tx,x:(a.x+b.x)/2,y:(a.y+b.y)/2}}),N=ps.length,cs=Math.cos(alpha*Math.PI/180),sn=Math.sin(alpha*Math.PI/180);const S:number[][]=[],T:number[][]=[],VN:number[][]=[],VT:number[][]=[];for(let i=0;i<N;i++){S[i]=[];T[i]=[];VN[i]=[];VT[i]=[];for(let j=0;j<N;j++){if(i===j){S[i][j]=.5;T[i][j]=0;VN[i][j]=0;VT[i][j]=.5;continue}const q=ps[j],r=ps[i],dx=r.x-q.a.x,dy=r.y-q.a.y,X=dx*q.tx+dy*q.ty,Y=-dx*q.ty+dy*q.tx;const u=Math.log((X*X+Y*Y)/((X-q.l)**2+Y*Y))/(4*Math.PI);const v=(Math.atan2(Y,X-q.l)-Math.atan2(Y,X))/(2*Math.PI);const ux=u*q.tx-v*q.ty,uy=u*q.ty+v*q.tx;S[i][j]=ux*r.nx+uy*r.ny;T[i][j]=ux*r.tx+uy*r.ty;VN[i][j]=-uy*r.nx+ux*r.ny;VT[i][j]=-uy*r.tx+ux*r.ty}}
const A=S.map((r,i)=>[...r,VN[i].reduce((a,b)=>a+b,0)]),b=ps.map(p=>-(cs*p.nx+sn*p.ny));A.push([...T[0].map((v,j)=>v+T[N-1][j]),VT[0].reduce((a,b)=>a+b,0)+VT[N-1].reduce((a,b)=>a+b,0)]);b.push(-(cs*(ps[0].tx+ps[N-1].tx)+sn*(ps[0].ty+ps[N-1].ty)));const z=linear(A,b);let fx=0,fy=0,cm=0;const data=ps.map((p,i)=>{const vt=cs*p.tx+sn*p.ty+T[i].reduce((s,v,j)=>s+v*z[j],0)+z[N]*VT[i].reduce((a,b)=>a+b,0),cp=(1-vt*vt)/Math.sqrt(1-mach*mach),dfx=-cp*p.nx*p.l,dfy=-cp*p.ny*p.l;fx+=dfx;fy+=dfy;cm-=((p.x-.25)*dfy-p.y*dfx);return {x:p.x,y:p.y,cp,velocity:Math.abs(vt),side:i<points.reduce((best,p,k)=>p.x<points[best].x?k:best,0)?'upper':'lower'}});return {cl:-fx*sn+fy*cs,cd:fx*cs+fy*sn,cm,data,flux:ps.reduce((s,p,i)=>s+z[i]*p.l,0)}}
