"""Patient-level CAF-program associations; predefined overlap and null checks."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
A=Path(__file__).resolve().parents[1];D=A/'derived';O=A/'results/adapted_v02';O.mkdir(parents=True,exist_ok=True)
SEED=2026090911;NBOOT=5000;NNULL=1000
S=json.loads((D/'source_gene_sets_v02_detectability.json').read_text(encoding='utf-8'))
checks=json.loads((A/'qa/data_preflight.json').read_text());assert all(x['library_exact'] for x in checks)
cohorts=['GSE132465','GSE166555'];data={c:np.load(D/(f'{c}_fibroblast_featurefiltered_counts.npz' if c=='GSE132465' else f'{c}_fibroblast_counts.npz')) for c in cohorts}
common=set(data[cohorts[0]]['genes'])&set(data[cohorts[1]]['genes'])
sen=[g for g in S['SenMayo119'] if g in common]
sets={'SenMayo':sen}
ledger=[]
for key in S['source']:
    for variant in ['source','purged']:
        gl=[g for g in S[variant][key] if g in common]
        sets[f'{key}_{variant}']=gl
        ledger.append(dict(program=key,variant=variant,source_n=100,post_purge_n=len(S[variant][key]),common_n=len(gl),genes=';'.join(gl),removed=';'.join(sorted(set(S['source'][key])-set(gl))),evaluable=len(gl)>=50))
pd.DataFrame(ledger).to_csv(O/'signature_coverage_and_overlap.csv',index=False)
(O/'locked_evaluated_gene_sets.json').write_text(json.dumps(sets,indent=2),encoding='utf-8')
assert all(len(sets[k+'_purged'])>=50 for k in S['source']),'Program coverage gate failed'
assert all(not(set(sets[k+'_purged'])&set(sen)) for k in S['source'])

def bh(p):
    p=np.asarray(p,dtype=float);out=np.full(len(p),np.nan);ok=np.where(np.isfinite(p))[0];order=ok[np.argsort(p[ok])]
    out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1]);return out

def rank_residual_corr(scores,y,meta,mki,mode='adjusted'):
    """Rows=exposures, columns=cells; one coefficient per exposure per patient."""
    scores=np.atleast_2d(scores);result=[];patients=[];ranks=[]
    for patient in sorted(meta.patient.unique()):
        idx=np.where(meta.patient.to_numpy()==patient)[0]
        x=stats.rankdata(scores[:,idx],axis=1);yr=stats.rankdata(y[idx])
        design=[np.ones(len(idx))]
        if mode!='unadjusted':
            design.extend([stats.rankdata(np.log1p(meta.nCount_RNA.to_numpy()[idx])),stats.rankdata(mki[idx])])
        if mode=='subtype_adjusted':
            indicators=pd.get_dummies(meta.iloc[idx].original_subtype,drop_first=True,dtype=float).to_numpy()
            design.extend(indicators.T)
        des=np.column_stack(design);u,sv,_=np.linalg.svd(des,full_matrices=False)
        q=u[:,sv>sv.max()*1e-10];rank=q.shape[1]
        assert len(idx)-rank>=10,(patient,len(idx),rank)
        xr=x-(x@q)@q.T;yr=yr-q@(q.T@yr)
        den=np.sqrt((xr*xr).sum(axis=1)*np.dot(yr,yr))
        r=np.divide(xr@yr,den,out=np.full(len(scores),np.nan),where=den>1e-8)
        result.append(np.clip(r,-.999999,.999999));patients.append(patient);ranks.append(rank)
    return np.column_stack(result),patients,ranks

def summarize(r,seed):
    r=np.asarray(r);r=r[np.isfinite(r)];z=np.arctanh(r);n=len(z);rng=np.random.default_rng(seed)
    boot=np.tanh(z[rng.integers(0,n,(NBOOT,n))].mean(axis=1))
    signs=rng.choice([-1,1],(50000,n))
    p_sign=(1+np.sum(np.abs((signs*z).mean(axis=1))>=abs(z.mean())))/(len(signs)+1)
    loo=np.tanh((z.sum()-z)/(n-1))
    return dict(r=float(np.tanh(z.mean())),ci_low=float(np.quantile(boot,.025)),ci_high=float(np.quantile(boot,.975)),p=float(stats.ttest_1samp(z,0).pvalue),p_signflip=float(p_sign),patients=n,positive_patients=int((r>0).sum()),negative_patients=int((r<0).sum()),loo_min=float(loo.min()),loo_max=float(loo.max()))

all_summary=[];patient_rows=[];null_summary=[];diagnostics=[];gene_metrics=[]
for ci,cohort in enumerate(cohorts):
    f=data[cohort];meta=pd.read_csv(D/f'{cohort}_cells.csv');genes=f['genes'];raw=f['counts'].astype(np.float64)
    lookup={g:i for i,g in enumerate(genes)}
    log=np.log1p(raw/meta.nCount_RNA.to_numpy()[None,:]*10000)
    means=log.mean(axis=1);sd=log.std(axis=1,ddof=1);detection=(raw>0).mean(axis=1)
    z=np.divide(log-means[:,None],sd[:,None],out=np.zeros_like(log),where=sd[:,None]>0).astype(np.float32)
    mki=log[lookup['MKI67']];cell_scores=meta.copy();scores={}
    for name,gl in sets.items():
        scores[name]=z[[lookup[g] for g in gl]].mean(axis=0,dtype=np.float64);cell_scores[name]=scores[name]
        gene_metrics.append(dict(cohort=cohort,program=name,genes=len(gl),nonconstant_genes=sum(sd[lookup[g]]>0 for g in gl),median_detection=float(np.median([detection[lookup[g]] for g in gl])),median_expression=float(np.median([means[lookup[g]] for g in gl]))))
    # Source-specific legacy sensitivity, never overwrites the original manuscript score.
    legacy_gl=[g for g in S['SenMayo119'] if g in lookup]
    legacy_y=z[[lookup[g] for g in legacy_gl]].mean(axis=0,dtype=np.float64)
    cell_scores['SenMayo_source_specific']=legacy_y
    cell_scores.to_csv(O/f'{cohort}_cell_program_scores.csv',index=False)
    pd.DataFrame({'gene':genes,'mean_log1p_CP10K':means,'detection_fraction':detection,'sd_log1p_CP10K':sd}).to_csv(O/f'{cohort}_gene_metrics.csv',index=False)
    y=scores['SenMayo']
    for vi,variant in enumerate(['purged','source']):
        for mi,mode in enumerate(['adjusted','unadjusted','subtype_adjusted']):
            for ki,key in enumerate(S['source']):
                rr,pats,ranks=rank_residual_corr(scores[f'{key}_{variant}'],y,meta,mki,mode)
                result=summarize(rr[0],SEED+1000*ci+100*vi+10*mi+ki)
                all_summary.append(dict(cohort=cohort,program=key,variant=variant,model=mode,outcome='common_SenMayo',cells=len(meta),**result))
                for pat,r,rank in zip(pats,rr[0],ranks):patient_rows.append(dict(cohort=cohort,program=key,variant=variant,model=mode,patient=pat,r=r,design_rank=rank,n_cells=int(meta.patient.eq(pat).sum())))
    for ki,key in enumerate(S['source']):
        rr,_,_=rank_residual_corr(scores[f'{key}_purged'],legacy_y,meta,mki)
        all_summary.append(dict(cohort=cohort,program=key,variant='purged',model='adjusted',outcome='source_specific_SenMayo',cells=len(meta),**summarize(rr[0],SEED+ci*1000+400+ki)))
    # Expression/detection matching uses no target-outcome association.
    excluded=set(S['SenMayo119']+S['FAP13']+S['matrix4']+S['fib5']+['MKI67'])
    excluded.update(g for gl in S['source'].values() for g in gl)
    pool=np.array([j for j,g in enumerate(genes) if g not in excluded and sd[j]>0 and detection[j]>=.01 and not g.startswith(('MT-','RPL','RPS'))])
    features=np.column_stack([means,detection]);scale=features[pool].std(axis=0,ddof=1)
    assert (scale>0).all();features=features/scale
    for ki,key in enumerate(S['source']):
        target=np.array([lookup[g] for g in sets[f'{key}_purged']]);k=len(target)
        distance=((features[target,None,:]-features[pool][None,:,:])**2).sum(axis=2)
        # Enough ordered neighbours to permit unique matching across large sets.
        near=pool[np.argsort(distance,axis=1)[:,:min(500,len(pool))]]
        rng=np.random.default_rng(SEED+ci*1000+500+ki);draws=[]
        for b in range(NNULL):
            chosen=np.empty(k,dtype=int);used=set()
            for j in rng.permutation(k):
                candidates=[int(x) for x in near[j] if x not in used][:100]
                pick=int(rng.choice(candidates));chosen[j]=pick;used.add(pick)
            assert len(used)==k;draws.append(chosen)
        draws=np.asarray(draws);np.savez_compressed(D/f'{cohort}_{key}_matched_gene_indices.npz',indices=draws,genes=genes,targets=target)
        null_z=[]
        for start in range(0,NNULL,50):
            block=z[draws[start:start+50]].mean(axis=1,dtype=np.float64)
            rr,_,_=rank_residual_corr(block,y,meta,mki)
            assert np.isfinite(rr).all()
            null_z.extend(np.arctanh(rr).mean(axis=1))
        null_z=np.asarray(null_z);null_r=np.tanh(null_z)
        observed=[x for x in all_summary if x['cohort']==cohort and x['program']==key and x['variant']=='purged' and x['model']=='adjusted' and x['outcome']=='common_SenMayo'][0]
        obs_z=np.arctanh(observed['r']);center=float(np.median(null_z))
        empirical=(1+np.sum(abs(null_z-center)>=abs(obs_z-center)))/(NNULL+1)
        ns=dict(cohort=cohort,program=key,observed_r=observed['r'],null_median=float(np.median(null_r)),null_low=float(np.quantile(null_r,.025)),null_high=float(np.quantile(null_r,.975)),empirical_p=float(empirical),null_draws=NNULL,pool_genes=len(pool))
        null_summary.append(ns)
        pd.DataFrame({'draw':np.arange(NNULL),'r':null_r,'mean_fisher_z':null_z}).to_csv(O/f'{cohort}_{key}_null_distribution.csv',index=False)
        # Absolute standardized differences in the two matching features.
        for fi,feature in enumerate(['mean_expression','detection_fraction']):
            smd=features[draws,fi].mean(axis=1)-features[target,fi].mean()
            diagnostics.append(dict(cohort=cohort,program=key,feature=feature,target_mean=float(features[target,fi].mean()),null_mean=float(features[draws,fi].mean()),median_abs_difference=float(np.median(abs(smd))),max_abs_difference=float(abs(smd).max()),fraction_draws_abs_difference_gt_0_1=float((abs(smd)>.1).mean())))
        print('completed matched-null',cohort,key,flush=True)
    print('cohort complete',cohort,'patients',meta.patient.nunique(),'cells',len(meta),flush=True)

summary=pd.DataFrame(all_summary)
summary['q']=summary.groupby(['cohort','variant','model','outcome'],sort=False).p.transform(lambda x:bh(x.to_numpy()))
summary.to_csv(O/'single_cell_associations.csv',index=False)
pd.DataFrame(patient_rows).to_csv(O/'patient_associations.csv',index=False)
ns=pd.DataFrame(null_summary);ns['q']=ns.groupby('cohort',sort=False).empirical_p.transform(lambda x:bh(x.to_numpy()))
ns.to_csv(O/'matched_null_summary.csv',index=False)
pd.DataFrame(diagnostics).to_csv(O/'matched_null_diagnostics.csv',index=False)
pd.DataFrame(gene_metrics).to_csv(O/'score_gene_coverage_by_cohort.csv',index=False)
primary=summary[(summary.variant=='purged')&(summary.model=='adjusted')&(summary.outcome=='common_SenMayo')]
primary.to_csv(O/'primary_single_cell_results.csv',index=False)
print(primary.to_string(index=False));print(ns.to_string(index=False))
