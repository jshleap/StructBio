#!/usr/bin/python
''' 
Plotting the principal coordinate analysis eigen vectors in 2D and 3D, for MD simulation
'''
# importing bit##########################################################################################
import matplotlib , sys , os , time
import matplotlib.pyplot as plt
from PlotPCoA import parse_eigenvectors
from numpy import array,std
# ############################################################################################

# ############################################################################################
# Some python functions#######################################################################
# ############################################################################################

def parse_fasta(prefix):
    t=[]
    fil=open(prefix+'.fasta').read().split('\n>')
    for l in fil:
        if l == '':
            continue
        else:
            name=l[:l.find('\n')]
            t.append(float(name[name.find(':')+1:]))
    a=array(t)
    a=a-t[0]
    t=list(a)
    a=a/max(a)
    return t,a

def PlotMD(ax1,ax2,t,a):
    #markerfacecolor=None
    #cm = plt.get_cmap("RdYlGn")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines['top'].set_color('none')
    ax.xaxis.tick_bottom()
    ax.spines['right'].set_color('none')
    ax.yaxis.tick_left()
    
    for i in range(len(t)):
        ax.scatter(ax1[i],ax2[i],c='%f'%(a[i]))
    
    ax.annotate('Initial conformation', xy=(ax1[0], ax2[0]), 
                xytext= (float(ax1[0])+10, float(ax2[0])+20),
                arrowprops=dict(facecolor='blue', shrink=0.05, frac=0.15))
        
    ax.annotate('Final conformation', xy=(ax1[len(t)-1], ax2[len(t)-1]),
                xytext=(float(ax1[len(t)-1])+10, float(ax2[len(t)-1])+40),
                arrowprops=dict(facecolor='blue', shrink=0.05, frac=0.15))
    plt.xlabel('Principal coordinate 1')
    plt.ylabel('Principal coordinate 2')
    #plt.show()
    fig.savefig(prefix+'_PCoA.png',dpi=300)
'''axes=[0,1]
for arg in sys.argv[1:]:
    if arg == '-ax=':
        axes=arg[4:].split(',')
        axes=[int(x)-1 for x in axes]'''
        
prefix=sys.argv[1]
t,a= parse_fasta(prefix)
ax1,ax2,ax3,nvec,allax = parse_eigenvectors(prefix)
PlotMD(ax1,ax2,t,a)
