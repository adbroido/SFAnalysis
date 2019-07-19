import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D

ALMOST_BLACK = '0.125'
DARK = '0.4'
GRAY_COLOR = '0.5'
ALL_COLOR = (0.156, 0.254, 0.466)
SUPERWEAK_COLOR = (0.647, 0.772, 0.972)
WEAKEST_COLOR = (0.419, 0.627, 0.960)
WEAK_COLOR = (0.258, 0.525, 0.85)
STRONG_COLOR = (0.223, 0.368, 0.674)
STRONGEST_COLOR = (0.156, 0.254, 0.466)
LABEL_SIZE = 13
TICK_SIZE = 13

def make_domain_ploth(ax, df, xlab=False, any=False, cutoff=False, allrats=[]):
    tot = float(len(df))
    if any:
        strongest = len(df.query('Strong_Any==True'))
        strong = len(df.query('Strong_Any==True'))
        weak = len(df.query('Weak_Any==True'))
        weakest = len(df.query('Weakest_Any==True'))
        superweak = len(df.query('Super_Weak_Any==True'))
        fail = (len(df.query("Weakest_Any==False").query("Super_Weak_Any==False")))
    elif cutoff:
        strongest = len(df.query('Strongest_No_PLwC==True'))
        strong = len(df.query('Strong_No_PLwC==True'))
        weak = len(df.query('Weak_PLwC==True'))
        weakest = len(df.query('Weakest_PLwC==True'))
        superweak = len(df.query('Super_Weak_PLwC==True'))
        fail = (len(df.query("Weakest_Any==False").query("Super_Weak_Any==False")))
    else:
        strongest = len(df.query('Strongest==True'))
        strong = len(df.query('Strong==True'))
        weak = len(df.query('Weak==True'))
        weakest = len(df.query('Weakest==True'))
        superweak = len(df.query('Super_Weak==True'))
        fail = (len(df.query("Weakest==False").query("Super_Weak==False")))

    counts = [strongest, strong, weak, weakest, superweak, fail]
    barheights = [count/tot for count in counts]
    if allrats !=[]:
        allticks = [barheights[i]-allrats[i] for i in range(len(counts))]
        signs = ['+' if allticks[i] >= 0 else ' -' for i in range(len(allticks))]

    width = 4
    totwid = 4
    xlocs = np.arange(0,4*totwid,totwid)
    xlocs = np.append(xlocs, xlocs[-1]+totwid+width/2.)
    xlocs = np.append(xlocs, xlocs[-1]+totwid+width/2.)

    colors=[STRONGEST_COLOR, STRONG_COLOR, WEAK_COLOR, WEAKEST_COLOR, SUPERWEAK_COLOR, GRAY_COLOR]
    ylabels = ['Strongest', 'Strong', 'Weak', 'Weakest', 'Super-Weak', 'Not\n Scale-Free' ]
    yticks = xlocs+width/2.

    for i in range(6):
        ax.barh(y=yticks[i],width=barheights[i], height=width, color = colors[i],edgecolor='w')
        if barheights[i]>0.65:
            ax.text(barheights[i]-0.28,yticks[i],'%d (%.2f)' %(int(counts[i]) ,barheights[i]) , fontsize=LABEL_SIZE, ha='left', va='center')
        else:
            ax.text(barheights[i]+0.02,yticks[i],'%d (%.2f)' %(int(counts[i]) ,barheights[i]) , fontsize=LABEL_SIZE, ha='left', va='center')
        if allrats!=[]:
            ax.text(1.0, yticks[i], '%s %.2f' %(signs[i], np.abs(allticks[i])), fontsize=LABEL_SIZE, ha='left', va='center')
            if np.abs(allticks[i]) <=0.025:
                ax.scatter(0.97, yticks[i], marker='>', s=60, color='0.4')
            elif signs[i]=='+':
                ax.scatter(0.97, yticks[i], marker='^', s=60, color='green')
            else:
                ax.scatter(0.97, yticks[i], marker='v', s=60, color='red')
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
#     ax.set_yticklabels([])
    ax.set_xlim(0,1.0)
    ax.set_ylim(0,30)
    if xlab:
        ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
    else:
        ax.set_xticks([])
    ax.axhline((xlocs[3]+width + xlocs[4])/2., color=ALMOST_BLACK, linewidth=2)
    ax.axhline((xlocs[4]+width + xlocs[5])/2., color=ALMOST_BLACK, linewidth=2)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=TICK_SIZE)
