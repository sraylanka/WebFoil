import { analyze, runBatch, designLift, type AnalysisJob } from './advanced';
import type { Point } from './foil';
self.onmessage = (event: MessageEvent<{kind:'analysis';job:AnalysisJob}|{kind:'batch';jobs:AnalysisJob[]}|{kind:'design';points:Point[];target:number;alpha:number}>) => {
  try {
    const q=event.data;
    const result=q.kind==='analysis'?analyze(q.job):q.kind==='batch'?runBatch(q.jobs):designLift(q.points,q.target,q.alpha);
    self.postMessage({result});
  } catch(e) { self.postMessage({error:(e as Error).message}); }
};
