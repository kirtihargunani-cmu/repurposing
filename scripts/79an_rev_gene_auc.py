"""D203 criteria 3,4,5: does rev_gene rank better, and is any gain the proliferation confound?"""
import pandas as pd, numpy as np
R='/home/kirti/projects/repurposing/'
acc=[]
for ch in pd.read_csv(R+'results/lincs/regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz',
                      sep='\t',usecols=['cohort','pert_id','cmap_name','window','rev_contrib','rev_gene_contrib'],
                      chunksize=3_000_000,low_memory=False):
    acc.append(ch.groupby(['cohort','pert_id','cmap_name','window'],sort=False)
               [['rev_contrib','rev_gene_contrib']].sum().reset_index())
S=pd.concat(acc,ignore_index=True).groupby(['cohort','pert_id','cmap_name','window'],as_index=False)[
    ['rev_contrib','rev_gene_contrib']].sum().rename(columns={'rev_contrib':'rev','rev_gene_contrib':'rev_gene'})
S.to_csv(R+'results/lincs/interpretation/d203_rev_gene.tsv.gz',sep='\t',index=False,float_format='%.6g')
L=pd.read_csv(R+'results/lincs/labels/pair_labels_v2_tf.tsv.gz',sep='\t',low_memory=False)
L=L[L.unit_kind=="group"][['unit_id','pert_id','label_phase']]
D=S.merge(L,left_on=['cohort','pert_id'],right_on=['unit_id','pert_id'])
D['y']=(D.label_phase>=3).astype(int)
def auc(y,s):
    m=~np.isnan(s); y,s=np.asarray(y)[m],np.asarray(s)[m]
    n1,n0=int(y.sum()),int((1-y).sum())
    if n1<3 or n0<3: return None
    r=pd.Series(s).rank().to_numpy()
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
rows=[]
for (c,w),g in D.groupby(['cohort','window']):
    if len(g)<50 or g.y.sum()<3: continue
    a1,a2=auc(g.y.values,g.rev.values),auc(g.y.values,g.rev_gene.values)
    if a1 is None or a2 is None: continue
    rows.append(dict(cohort=c,window=w,n=len(g),n_pos=int(g.y.sum()),
                     auc_rev=round(a1,3),auc_rev_gene=round(a2,3),delta=round(a2-a1,3),
                     rho=round(g.rev.corr(g.rev_gene,method='spearman'),3)))
Rr=pd.DataFrame(rows)
Rr.to_csv(R+'results/lincs/interpretation/d203_auc.tsv',sep='\t',index=False)
d=Rr.delta
print(f"{len(Rr)} cells | median AUC rev {Rr.auc_rev.median():.3f} -> rev_gene {Rr.auc_rev_gene.median():.3f}")
print(f"CRITERION 4: median delta {d.median():+.3f} (bar +0.030), better in {int((d>0).sum())}/{len(d)} ({100*(d>0).mean():.0f}%, bar 60%) -> {'PASS' if (d.median()>=0.03 and (d>0).mean()>=0.60) else 'FAIL'}")
print(f"CRITERION 3: median Spearman(rev, rev_gene) {Rr.rho.median():+.3f}")
for w,g in Rr.groupby('window'):
    print(f"  {w:5s} n={len(g):2d} rev {g.auc_rev.median():.3f} -> rev_gene {g.auc_rev_gene.median():.3f} (delta {g.delta.median():+.3f}, better {int((g.delta>0).sum())}/{len(g)})")
print("\nbiggest gains / losses:")
print(Rr.nlargest(5,'delta')[['cohort','window','n_pos','auc_rev','auc_rev_gene','delta','rho']].to_string(index=False))
print(Rr.nsmallest(4,'delta')[['cohort','window','n_pos','auc_rev','auc_rev_gene','delta','rho']].to_string(index=False))
