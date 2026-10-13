import type { Point } from './foil';
export function validateGeometry(points:Point[]){
  if(!Array.isArray(points)||points.length<21||points.length>301||points.some(p=>!Number.isFinite(p.x)||!Number.isFinite(p.y)))throw Error('Each element needs 20–300 finite panels.');
  const last=points.at(-1)!;
  if(Math.hypot(points[0].x-last.x,points[0].y-last.y)>1e-6)throw Error('Close the contour at its trailing edge.');
  let area=0;const cross=(a:Point,b:Point,c:Point)=>(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
  for(let i=0;i<points.length-1;i++){
    const a=points[i],b=points[i+1];area+=a.x*b.y-b.x*a.y;
    if(Math.hypot(b.x-a.x,b.y-a.y)<1e-8)throw Error('Geometry has duplicate adjacent points.');
    for(let j=i+2;j<points.length-1;j++){
      if(i===0&&j===points.length-2)continue;const c=points[j],d=points[j+1];
      if(cross(a,b,c)*cross(a,b,d)<=0&&cross(c,d,a)*cross(c,d,b)<=0&&Math.max(a.x,b.x)>=Math.min(c.x,d.x)&&Math.max(c.x,d.x)>=Math.min(a.x,b.x)&&Math.max(a.y,b.y)>=Math.min(c.y,d.y)&&Math.max(c.y,d.y)>=Math.min(a.y,b.y))throw Error('Geometry has intersecting panels.');
    }
  }
  if(area<=1e-10)throw Error('Use a nonzero-area counterclockwise contour from the trailing edge.');
}
