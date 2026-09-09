"""Auxiliary TCGA associations after exact frozen-score reproduction."""
from pathlib import Path
import sys,json,hashlib,os
import numpy as np
import pandas as pd
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
A=Path(__file__).resolve().parents[1];P=A.parents[1];D=A/'derived';O=A/'results/adapted_v02';O.mkdir(parents=True,exist_ok=True)
B=A.parent/'results/bulk_composition'
S=json.loads((D/'source_gene_sets_v02_detectability.json').read_text())
source=Path(os.environ.get('FAP_TCGA_INPUT_DIR', str(A/'inputs')))
manifest=pd.read_csv(B/'tcga_primary_sample_manifest.csv');assert len(manifest)==380 and manifest.patient.nunique()==380
frames=[];receipts=[]
for project in ['COAD','READ']:
    path=source/f'TCGA.{project}.sampleMap_HiSeqV2.gz'
    frame=pd.read_csv(path,sep='\t',index_col=0)
    frame.columns=['-'.join(c.split('-')[:3]+[c.split('-')[3][:2]]) for c in frame.columns]
    assert not frame.columns.duplicated().any()
    wanted=manifest.loc[manifest.project.eq('TCGA-'+project),'sample'].tolist()
    assert set(wanted)<=set(frame.columns)
    frames.append(frame[wanted]);receipts.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
genes=frames[0].index.intersection(frames[1].index);x=pd.concat([f.loc[genes] for f in frames],axis=1)[manifest['sample']]
assert not x.index.duplicated().any() and np.isfinite(x.to_numpy()).all()
z=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1,ddof=1),axis=0).fillna(0)
def alias(gl):return ['IL8' if g=='CXCL8' else g for g in gl]
def score(gl):
    gl=[g for g in alias(gl) if g in z.index];assert gl
    return z.loc[gl].mean(axis=0).to_numpy(),gl
scores=manifest.copy();checks=[]
frozen=pd.read_csv(B/'tcga_primary_target_scores.csv').set_index('sample').loc[manifest['sample']]
bulk_source=(D/'senmayo_full_source_genes.txt').read_text().split()
bulk_sen=[g for g in bulk_source if g not in S['FAP13']]
assert all(not(set(v)&set(bulk_sen)) for v in S['purged'].values())
for key,gl in [('FAP13',S['FAP13']),('matrix4',S['matrix4']),('fib5',S['fib5']),('SenMayo',bulk_sen)]:
    sc,represented=score(gl);scores[key]=sc
    delta=float(np.max(abs(sc-frozen[key].to_numpy())))
    checks.append(dict(score=key,max_absolute_error=delta,n_genes=len(represented),pass_tolerance_1e_8=delta<1e-8))
(A/'qa/TCGA_score_reproduction.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
(D/'TCGA_consumed_files.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
assert all(r['pass_tolerance_1e_8'] for r in checks),checks
print('TCGA locked-score reproduction PASS',checks,flush=True)
coverage=[]
for key in S['source']:
    sc,gl=score(S['purged'][key]);assert len(gl)>=50
    scores[key]=sc;coverage.append(dict(program=key,genes=len(gl),gene_list=';'.join(gl)))
pd.DataFrame(coverage).to_csv(O/'TCGA_signature_coverage.csv',index=False)
scores.to_csv(O/'TCGA_380_program_scores.csv',index=False)

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
