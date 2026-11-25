from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
p=Path(sys.argv[1]); out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
s=pd.read_csv(p/'STRATUM_SSE_COMPOSITION.csv');s=s[s.model.eq('FOP')]; order=['FIT','SELECT','DEV_DIAGNOSTIC']; years=['2019 fit','2020 select','2021 diagnostic']; changes=s[s.stratum.isin(['UP','DOWN'])].groupby('partition')[['rows_pct','sse_pct']].sum().loc[order]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.hashsalt':'C19V015'}); fig,ax=plt.subplots(1,2,figsize=(9.2,3.45),layout='constrained')
x=np.arange(3); w=.34
for j,(key,lab,col) in enumerate([('rows_pct','Record share','#8795a2'),('sse_pct','Squared-error share','#216e77')]):
 vals=changes[key].values; bars=ax[0].bar(x+(j-.5)*w,vals,w,label=lab,color=col);ax[0].bar_label(bars,labels=[f'{v:.2f}' for v in vals],padding=3,fontsize=8)
ax[0].set(xticks=x,xticklabels=years,ylabel='Share of partition total (%)',title='(a) Large net changes',ylim=(0,26));ax[0].legend(loc='upper left',frameon=False,fontsize=8)
for i,(ss,col,lab) in enumerate([('UP','#b65b34','Increase'),('DOWN','#216e77','Decrease')]):
 r=s[s.stratum.eq(ss)].set_index('partition').loc[order];a=x+(i-.5)*w;ax[1].bar(a,r.macro_bias,w,color=col,label=lab+' (macro)');ax[1].scatter(a,r.pooled_bias,marker='D',s=26,facecolors='white',edgecolors='black',zorder=4,label='Pooled' if i==0 else None)
ax[1].axhline(0,color='black',lw=.7);ax[1].set(xticks=x,xticklabels=years,ylabel='Measured − predicted (m³/h)',title='(b) Signed FOP residual',ylim=(-2.0,1.05));ax[1].legend(frameon=False,fontsize=8,loc='lower right')
for a in ax:a.grid(axis='y',alpha=.20);a.set_axisbelow(True);a.tick_params(axis='x',labelsize=8)
fig.savefig(out/'F_S4_timing_composition.png',dpi=300);fig.savefig(out/'F_S4_timing_composition.svg');plt.close(fig)
