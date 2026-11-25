from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
p=argparse.ArgumentParser();p.add_argument('reports');p.add_argument('output');args=p.parse_args()
q=Path(args.reports); out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
root=out
c=pd.read_csv(q/'CANDIDATES.csv');s=pd.read_csv(q/'STRATUM_MACRO_METRICS.csv');w=c[c.selected].set_index('family');near=c[c.near_family_min];ss=s.set_index(['partition','model','stratum'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'svg.hashsalt':'c19_v015'})
fig, axs=plt.subplots(1,2,figsize=(9.0,3.65),gridspec_kw={'width_ratios':[1.25,1]},layout='constrained')
mat=c.pivot(index='tau_s',columns='delay_s',values='select_macro_rmse'); best=c.select_macro_rmse.min(); loss_grid=(mat.values/best-1)*100
ax=axs[0]; im=ax.imshow(np.log10(1+loss_grid),aspect='auto',cmap='viridis',origin='lower');ax.set_xticks(range(len(mat.columns)),[f'{x:g}' for x in mat.columns]);ax.set_yticks(range(len(mat.index)),[f'{x:g}' for x in mat.index]);ax.set_xlabel('Delay candidates (s; categorical)' );ax.set_ylabel('Time-constant candidates (s; categorical)'); ax.set_title('(a) 2020 selection loss surface')
for _,r in near.iterrows():ax.scatter(list(mat.columns).index(r.delay_s),list(mat.index).index(r.tau_s),marker='o',s=29,facecolors='none',edgecolors='white',linewidths=1.0)
for _,r in w.iterrows():ax.scatter(list(mat.columns).index(r.delay_s),list(mat.index).index(r.tau_s),marker='*',s=75,facecolors='#fdae61',edgecolors='black',linewidths=.45)
cb=fig.colorbar(im,ax=ax,shrink=.88,pad=.03);ticks=[0,1,5,20,50,100];ticks=[t for t in ticks if np.log10(1+t)<=im.norm.vmax];cb.set_ticks(np.log10(1+np.array(ticks)));cb.set_ticklabels([f'{t:g}' for t in ticks]);cb.set_label('RMSE above grid minimum (%)\nlog-scaled color')
ax=axs[1]; st=['LOW_CHANGE','UP','DOWN'];xs=np.arange(3); colors=['#4378a0','#d8a143','#58846a'];
for i,model in enumerate(['LAG','FOP','FOPDT']):ax.bar(xs+(i-1)*.23,[ss.loc[('DEV_DIAGNOSTIC',model,k),'macro_rmse'] for k in st],.23,label=model,color=colors[i]);
ax.set_xticks(xs,[f'{label}\n{int(ss.loc[("DEV_DIAGNOSTIC","FOP",k),"rows"]):,} rows' for label,k in zip(('Small net change','Up','Down'),st)]);ax.set_ylabel('Mean-file RMSE (m³/h)');ax.set_title('(b) 2021 net-current-change strata');ax.legend(loc='upper left',frameon=False);ax.set_ylim(0,1.12);ax.grid(axis='y',alpha=.23);ax.set_axisbelow(True)
fig.savefig(root/'F_S3_timing_identifiability.png',dpi=300);fig.savefig(root/'F_S3_timing_identifiability.svg');plt.close(fig)
