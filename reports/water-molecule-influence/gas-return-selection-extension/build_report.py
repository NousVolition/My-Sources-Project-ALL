from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'
METHODS=['position','direction','velocity','hold','linear']
NAMES={'position':'Position only','direction':'Position + direction','velocity':'Position + full velocity','hold':'Hold position','linear':'Continue current velocity'}


def main():
    r=json.loads((DATA/'results.json').read_text());rows=[x for x in r['rows'] if x['role']=='test']
    figdir=HERE/'figures';figdir.mkdir(exist_ok=True)
    colors=['#98594c','#127f8b','#9a8b50','#7563a5']
    fig,axs=plt.subplots(1,2,figsize=(13,5.4),layout='constrained')
    labels=[f'Start {x["start"]}\n'+{.2:'weaker',1:'original',10:'stronger'}[x['friction_per_ps']] for x in rows]
    xs=np.arange(6);width=.18
    for i,m in enumerate(METHODS[1:]):
        ys=[100*x['methods'][m]['mse'][0]/x['methods']['position']['mse'][0] for x in rows]
        axs[0].bar(xs+(i-1.5)*width,ys,width,color=colors[i],label=NAMES[m])
    axs[0].axhline(100,color='#666666',ls=':',lw=1.5,label='Position-only reference')
    axs[0].set(xticks=xs,xticklabels=labels,ylabel='Forecast error / position-only error (%)',title='0.02 ps ahead: full velocity helps select a passage\nContinuing the current arrow still wins three cases')
    axs[0].legend(fontsize=8.5,ncol=2,loc='upper left')
    axs[0].set_ylim(0,285)
    for i,h in enumerate([.02,.05,.1]):
        ys=[100*(1-x['methods']['velocity']['mse'][i]/x['methods']['position']['mse'][i]) for x in rows]
        axs[1].bar(xs+(i-1)*.25,ys,.25,color=['#127f8b','#be6937','#7563a5'][i],label=f'{h:.2f} ps ahead')
    axs[1].axhline(0,color='#666666',lw=1)
    axs[1].set(xticks=xs,xticklabels=labels,ylabel='Reduction in forecast MSE (%)',title='The improvement fades with longer forecasts\nPositive: better than position-only selection',ylim=(-5,29))
    axs[1].legend(fontsize=9)
    for ax in axs:
        ax.set_axisbelow(True);ax.grid(axis='y',alpha=.15);ax.spines[['top','right']].set_visible(False)
        ax.tick_params(axis='x',labelsize=9)
    fig.suptitle('Choosing an earlier passage with information available at the time',fontsize=14)
    for ext in ['png','svg']:fig.savefig(figdir/f'return_selection.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    primary='';longer='';paired='';fine='';training=''
    for x in rows:
        cells=''.join(f'<td>{x["methods"][m]["mse"][0]:.6f}</td>' for m in METHODS)
        primary+=f'<tr><td>{x["start"]}</td><td>{x["friction_per_ps"]:g}</td><td>{x["eligible"]}/990</td>{cells}</tr>'
        gain=[100*(1-x['methods']['velocity']['mse'][i]/x['methods']['position']['mse'][i]) for i in range(3)]
        longer+=f'<tr><td>{x["start"]}</td><td>{x["friction_per_ps"]:g}</td>'+''.join(f'<td>{v:+.1f}%</td>' for v in gain)+'</tr>'
        paired+=f'<tr><td>{x["start"]}</td><td>{x["friction_per_ps"]:g}</td>'+''.join(f'<td>{100*x["methods"][m]["near_same_fraction"][0]:.1f}%</td>' for m in METHODS[:3])+'</tr>'
    for x in r['fine_recording_check']:
        g=x['friction_per_ps'];n=x['common_anchors']
        for resolution,key in [('0.010 ps','coarse_common'),('0.001 ps','fine_common')]:
            m=x[key];gain=100*(1-m['velocity']['mse'][0]/m['position']['mse'][0])
            fine+=f'<tr><td>{g:g}</td><td>{resolution}</td><td>{n}</td><td>{m["position"]["mse"][0]:.6f}</td><td>{m["velocity"]["mse"][0]:.6f}</td><td>{m["linear"]["mse"][0]:.6f}</td><td>{gain:.1f}%</td></tr>'
    for x in r['training']['grid']:
        training+=f'<tr><td>{x["method"]}</td><td>{x["weight"]:g}</td><td>{x["mean_mse"]:.7f}</td></tr>'
    page=(HERE/'report-template.html').read_text(encoding='utf-8')
    for key,val in {'__PRIMARY__':primary,'__LONGER__':longer,'__PAIRED__':paired,'__FINE__':fine,'__TRAINING__':training}.items():page=page.replace(key,val)
    (HERE/'report.html').write_text(page,encoding='utf-8')
    print(json.dumps(dict(report=str(HERE/'report.html'),test_runs=len(rows))))


if __name__=='__main__':main()
