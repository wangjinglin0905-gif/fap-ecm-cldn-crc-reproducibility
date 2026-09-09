"""Source-only marker adaptation after declared signature portability failure."""
from pathlib import Path
import json
import pandas as pd
from openpyxl import load_workbook
A=Path(__file__).resolve().parents[1];D=A/'derived'
s=json.loads((D/'source_gene_sets.json').read_text())
rows=list(load_workbook(A/'inputs/Nature_10344_Supplementary_Tables.xlsx',read_only=True,data_only=True)['Table 5 scRNA cell type markers'].values)
m=pd.DataFrame([r[:7] for r in rows[1:]],columns=rows[0][:7])
for c in ['p_val','p_val_adj','avg_log2FC','pct.1','pct.2']:m[c]=pd.to_numeric(m[c].astype(str).str.replace(',','.',regex=False),errors='raise')
full=(D/'senmayo_full_source_genes.txt').read_text().split()
protected=set(full+s['FAP13']+s['matrix4']+s['fib5']+['MKI67'])
audit=[]
for key,label in s['labels'].items():
    g=m[m.cluster.eq(label)&m.p_val_adj.lt(.05)&m.avg_log2FC.gt(0)&m['pct.1'].ge(.25)&(m['pct.1']-m['pct.2']).ge(.10)].sort_values(['avg_log2FC','gene'],ascending=[False,True]).drop_duplicates('gene')
    chosen=g.head(100).copy();assert len(chosen)==100
    s['source'][key]=chosen.gene.tolist();s['purged'][key]=[g for g in chosen.gene if g not in protected]
    chosen.insert(0,'program',key);chosen.to_csv(D/f'source_markers_{key}_v02_detectability.csv',index=False)
    audit.append(dict(program=key,source_eligible_markers=len(g),selected=100,purged_n=len(s['purged'][key]),purged_genes=';'.join(sorted(set(chosen.gene)&protected))))
s['adaptation']='Source pct.1>=0.25, pct.1-pct.2>=0.10; adjusted P<0.05, positive log2FC; top100 by log2FC. Amendment01; exploratory.'
(D/'source_gene_sets_v02_detectability.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
pd.DataFrame(audit).to_csv(D/'source_extraction_audit_v02_detectability.csv',index=False)
print(pd.DataFrame(audit).to_string(index=False))
