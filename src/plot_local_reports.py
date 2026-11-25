"""Plot core diagnostic quantities from local reproduced CSVs. Not historical pixel-exact templates."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

def make_figures(stage,source,out):
    out=Path(out);source=Path(source)
    if out.exists() or out.resolve().is_relative_to(source.resolve()):raise ValueError('use a new output folder outside input')
    out.mkdir(parents=True)
    def save(fig,name):fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight',pad_inches=.15);plt.close(fig)
    if stage=='nlr-lag':
        d=pd.read_csv(source/'LAG_CURVES_21.csv')
        for kind in ('level','increment'):
            fig,ax=plt.subplots(figsize=(7.5,4.6));q=d[d.kind==kind]
            for name,g in q.groupby('source'):
                s=g.groupby('lag_s').r.median();ax.plot(s.index,s,marker='o',ms=3,label=name)
            ax.axvline(0,ls=':',color='0.5');ax.axhline(0,ls=':',color='0.5');ax.set(xlabel='Lag (s); positive: command leads actual',ylabel='Median correlation',title=kind.capitalize()+' correlation');ax.legend(fontsize=7);save(fig,'lag_'+kind)
    elif stage=='nlr-replay':
        d=pd.read_csv(source/'DOMAIN_MACRO.csv');d=d[d.subset=='ALL'];fig,ax=plt.subplots(figsize=(7.5,4.6))
        for name,g in d.groupby('source'):
            g=g.sort_values('h_s');ax.plot(g.h_s,g.rmse_C-g.rmse_A,marker='o',label=name)
        ax.axhline(0,ls=':',color='0.5');ax.set(xlabel='Horizon (s)',ylabel='Mean trajectory RMSE(C) − RMSE(A) (kg/h)');ax.legend(fontsize=7);save(fig,'nlr_path_effects')
        g=d[d.h_s==30];fig,ax=plt.subplots(figsize=(8,4.5));g.set_index('source')[['mse_e_T','cross_R_T']].plot.bar(ax=ax);ax.set(ylabel='Mean trajectory squared terms ((kg/h)²)');ax.tick_params(axis='x',labelrotation=25);save(fig,'nlr_signed_terms')
    elif stage=='juelich-replay':
        d=pd.read_csv(source/'REPLAY_METRICS.csv');d=d[d.h==30];fig,ax=plt.subplots(figsize=(8,4.5));ax.bar(d.file_id.astype(str),d.rmse_H-d.rmse_A);ax.axhline(0,ls=':',color='0.5');ax.set(ylabel='RMSE(H) − RMSE(A) (m³/h)',xlabel='Measurement file identifier');ax.tick_params(axis='x',labelrotation=65);save(fig,'juelich_file_effects')
    else:raise ValueError('unsupported plot stage')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['nlr-lag','nlr-replay','juelich-replay']);p.add_argument('source',type=Path);p.add_argument('output',type=Path);a=p.parse_args();make_figures(a.stage,a.source,a.output)
