import sys
from igraph import compare_communities
from random import sample
from scipy.stats.kde import gaussian_kde
from pylab import plot, axvline,text
from numpy import linspace

def compare_clustering(prefix1,prefix2):
    # compare two communities based on landmark and membership files
    comm1 = read_memfiles(prefix1 + '.graphcluster')
    comm2 = read_memfiles(prefix2 + '.graphcluster')
    ind1, ind2 = get_assoc_landmarks(prefix1,prefix2)
    comm1 = [comm1[x] for x in ind1]
    comm2 = [comm2[x] for x in ind2]
    nmi = compare_communities(comm1, comm2, method='nmi', remove_none=False)
    return comm1,comm2, nmi


def read_memfiles(fil):
    # read the membership vector files
    f = open(fil).read().strip().split()
    if all([x.isdigit() for x in f]):
        f = [int(x) for x in f]
    else:
        f = [ord(x)-65 for x in f]
    return f

def read_landmarks(landmfile):
    # read individual lanmarks file, return a dictionary of lists of tuples
    # the tuples contain the resid and resname
    landmarks = {}
    name = None
    indl=None
    with open(landmfile) as L:
        for line in L:
            if line.startswith('>'):
                if indl and name:
                    landmarks[name] = indl
                indl=[]
                name = line.strip().strip('>')
            else:
                bl = line.strip().split()
                indl.append((bl[1],bl[2]))
                
    return landmarks
    
    
def get_assoc_landmarks(prefix1,prefix2):
    # get the apropriate combination of landmarks
    D1 = read_landmarks(prefix1+'.landmarks')
    D2 = read_landmarks(prefix2+'.landmarks')
    match = False
    while not match:
        for k1 in D2.keys():
            len1 = len(D[k1])
            s1 = set(D1[k1])
            for k2 in D2.keys():
                len2 = len(D2[k2])
                s2 = set(D2[k2])
                if len1 < len2:
                    exist =  any(s1.difference(s2))
                    if not exist:
                        inter = s1.intersection(s2)
                        indexes1 = range(len(D1[k1]))
                        indexes2 = [D2[k2].index(x) for x in inter]
                        break
                else:
                    exist =  any(s2.difference(s1))
                    if not exist:
                        inter = s1.intersection(s1)
                        indexes2 = range(len(D2[k2]))
                        indexes1 = [D1[k1].index(x) for x in inter]
                        break
    return indexes1,indexes2
    
def nmi_null(comm1,comm2,nmi,reps=1000):
    # create a null distribution of the nmi value by shuffling communities
    null = []
    #shuff1 = sample(comm1, len(comm1))
    for i in range(reps):
        shuff2 = sample(comm2, len(comm2))
        null.append(compare_communities(comm1, shuff2, method='nmi', remove_none=False))
    
    pdf = gaussian_kde(null)
    pval = pdf(nmi)
    x = linspace(0,1,100)
    plot(x,pdf,'k')
    axvline(x=nmi,linewidth=4, color='r')
    text(nmi+0.1, min(pdf)+0.2, 'p-val = %.3f'%(pval), fontdict=None, withdash=False)
    return pval

## Application of the code
prefix1 = sys.argv[1] # this should be preferably the MD
prefix2 = sys.argv[2]

comm1,comm2, nmi = compare_clustering(prefix1,prefix2)
pval = nmi_null(comm1,comm2,nmi)

print 'the NMI is %f with a p-value of %f'%(float(nmi),float(pval)