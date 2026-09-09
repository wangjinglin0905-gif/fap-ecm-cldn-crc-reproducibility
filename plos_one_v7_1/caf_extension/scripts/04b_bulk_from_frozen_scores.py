"""Recompute bulk estimates/intervals from the archived 380-patient scores."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
A=Path(__file__).resolve().parents[1];O=A/'results/adapted_v02'
S=json.loads((A/'derived/source_gene_sets_v02_detectability.json').read_text())
scores=pd.read_csv(O/'TCGA_380_program_scores.csv',index_col=0)
def cor_batch(x,y,c=None):
    rx=stats.rankdata(np.atleast_2d(x),axis=1);ry=stats.rankdata(np.atleast_2d(y),axis=1)
    rx-=rx.mean(axis=1,keepdims=True);ry-=ry.mean(axis=1,keepdims=True)
    if c is not None:
        rc=stats.rankdata(np.atleast_2d(c),axis=1);rc-=rc.mean(axis=1,keepdims=True)
        ss=(rc*rc).sum(axis=1,keepdims=True)
        rx-=np.divide((rx*rc).sum(axis=1,keepdims=True),ss,out=np.zeros_like(ss),where=ss>0)*rc
        ry-=np.divide((ry*rc).sum(axis=1,keepdims=True),ss,out=np.zeros_like(ss),where=ss>0)*rc
    return (rx*ry).sum(axis=1)/np.sqrt((rx*rx).sum(axis=1)*(ry*ry).sum(axis=1))
def bh(p):
    p=np.asarray(p);idx=np.argsort(p);out=np.empty(len(p));out[idx]=np.minimum(1,np.minimum.accumulate((p[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return out
rows=[]
for ki,key in enumerate(S['source']):
    for ti,target in enumerate(['SenMayo','FAP13','matrix4']):
        for mi,model in enumerate(['marginal','fib5_adjusted']):
            xx=scores[key].to_numpy();yy=scores[target].to_numpy();cc=scores.fib5.to_numpy() if mi else None
            r=float(cor_batch(xx,yy,cc)[0]);df=len(xx)-2-mi;t=r*np.sqrt(df/(1-r*r));p=float(2*stats.t.sf(abs(t),df))
            rng=np.random.default_rng(2026090940+ki*100+ti*10+mi);ix=rng.integers(0,len(xx),(5000,len(xx)))
            draws=cor_batch(xx[ix],yy[ix],cc[ix] if mi else None)
            rows.append(dict(cohort='TCGA-COAD/READ',program=key,target=target,model=model,patients=len(xx),rho=r,ci_low=float(np.quantile(draws,.025)),ci_high=float(np.quantile(draws,.975)),p=p,bootstrap=5000))
result=pd.DataFrame(rows);result['q']=result.groupby('model').p.transform(lambda p:bh(p.to_numpy()));result.to_csv(O/'TCGA_CAF_associations.csv',index=False)
print(result.to_string(index=False))
