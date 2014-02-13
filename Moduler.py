#!/usr/bin/python
'''
Graph based Modularity.
This script will evaluate the data for modules. Such modules are defined as correlating modules, so the clustering 
is performed in the correlation space. It has an optional statistical significance test for the clustering and power
analysis of the result.

ModulerV2 Copyright (C) 2012  Jose Sergio Hleap, Kyle Nguyen, Alex Safatli and Christian Blouin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@dal.ca


Requirements:
1) Python:
   a) numpy module
   b) scipy module
   c) rpy2 module
   d) matplotlib module
   ######################################
   #To install python modules in UBUNTU:#
   #sudo apt-get python-<module>        #
   ######################################

2) The igraph Python extension module:
   ######################################
   #To install python-igraph in UBUNTU: #
   #go to http://igraph.sourceforge.net/#
   #follow their instructions           #
   ######################################

3) The GNU R package:
   a) Library DAAG 
   b) Library pwr
   c) Library corpcor
   d) Library FactoMineR
   e) Library gplots
   ######################################
   #To install The R packages in UBUNTU:#
   #Get into R by typing "R" in your cmd#
   #Once in R, type:                    #
   #install.packages('<package>')       #
   ######################################   

4) set python path to:
   a) contactmapper.py (if using contacts)


'''
#importing bit####################################################################################################
import sys, pickle, os, datetime, decimal, optparse
from glob import glob
from numpy import mean
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from igraph import Graph, VertexClustering
from rpy2.robjects import r
from numpy import cov, matrix, corrcoef, array, zeros, mean, std, sqrt, var, savetxt, log
from scipy.stats import pearsonr, kendalltau, spearmanr, norm
from copy import deepcopy
from collections import Counter
from contactList.contactmapper import Contacts
from random import randint
from shutil import copyfile
from copy import deepcopy
from subprocess import Popen,call,PIPE
#from subprocess import call
#import apropriate R libraries####################################
r('library(DAAG)')
r('library(stats)')
r('library(pwr)')
r('library(corpcor)')
#r('library(gplots)')
r('library(FactoMineR)')
# End importing#################################################################################

#Some definitions###############################################################################
def Read_GM(prefix, dim, dist=False, f='gm'):
    '''
    Load data from a gm (coordinates file) file
    '''
    # Load Data ##########################################################
    data = []
    sample_size = 0
    # line for each dimension and create slush temp lists
    temps = []
    for i in range(dim):
        temps.append([])
        data.append([])
    # Set the filename with the proper suffix
    if f == 'csv':
        fname = prefix+'.csv'
    else:
        fname = prefix+'.gm'
    # Open file and read each line
    with open(fname) as F:
        for line in F:
            # Split into list
            if f == 'gm':
                if line.strip()[-1] == ';':
                    line = line.strip().split(';')[1:-1]
                else:
                    line = line.strip().split(';')[1:]
            elif f == 'csv':
                if line.strip()[-1] == ';':
                    line = line.strip().split(';')[:-1]
                else:
                    line = line.strip().split(';')			
            # Dispatch values in the right dimension list using i%dim to index
            for i in range(len(line)):
                temps[i%dim].append(float(line[i]))
            # Grow each matrix by 1 row
            for i in range(len(temps)):
                data[i].append(temps[i])
            sample_size += 1 
            # Flush temps 
            temps = []
            for i in range(dim):
                temps.append([])
    #if dist:
    #	data = as_distance(data)

    return data, sample_size

###############################################################################################
def as_distance(data):
    ''' Transform data into distances '''
    dim = len(data)
    matr=[]
    bm=[]
    for e in range(len(data[0])):
        for f in range(len(data[0][e])):
            temp=[]
            for g in range(len(data[0][e])):
                if dim == 3:
                    distance=((data[0][e][f] - data[0][e][g])**2+(data[1][e][f] - data[1][e][g])**2+
                              (data[2][e][f] - data[2][e][g])**2)
                else:
                    distance=((data[0][e][f] - data[0][e][g])**2+(data[1][e][f] - data[1][e][g])**2)
                temp.append(sqrt(distance))
        matr.append(temp)
    bm.append(matr)
    options.dim=1
    return bm

###############################################################################################
def Sigcorr(prefix,data, sample_size, method='fisher', confval=0.95, threshold='Auto', contacts=None, dim=3, absolutecov=False):
    '''
    Test if the correlation is significantly different than 0 with the method specified
    '''
    thr=[]
    matrices = []
    # Perform np magic
    for i in range(len(data)):
        mtx = SigCorrOneMatrix(prefix,data[i], sample_size, method, confval, absolutecov)
        matrices.append(mtx)

    #if threshold == 'Auto':
    #	thr = GetThreshold(matrices)

    return matrices#, thr

###############################################################################################

def SigCorrOneMatrix(prefix,data, sample_size, method='fisher', confval=0.95, absolutecov=False):
    '''
        Performs the significance of correlation test according to one of many possible tests:
        inputs:
        data : list of matrices from read_GM
        method : pearson, kendall, spearman or fisher	    
        dim : number of dimensions (integer)
    '''
    # ###############################################################
    # Convert to np matrix
    data = matrix(data)
    co = corrcoef(data.T)
    fi = open(prefix+'.bin','w')
    pickle.dump(co,fi)
    fi.close()
    # Get a matrix of zeroes.
    zero=zeros(shape=(data.shape[1], data.shape[1]))

    for e in range(len(data.T)):
        for f in range(e, len(data.T)):
            if method == 'pearson' or method == 'fisher':
                p=pearsonr(array(data.T[e]).ravel(),array(data.T[f]).ravel())
            if method == 'kendall':
                p=kendalltau(array(data.T[e]).ravel(),array(data.T[f]).ravel())
            if method == 'spearman':
                p=spearmanr(array(data.T[e]).ravel(),array(data.T[f]).ravel())					

            # Symmetrize	
            if method == 'fisher':
                if p[0] == 1.0:
                    p = (0.999,p[1])
                if absolutecov:
                    if abs(F_transf(p[0])) > Z_fisher(1-((1-confval)/2), sample_size)\
                       and Power_r(p[0],confval,options.power)[0] <= sample_size:
                        zero[e][f] = zero[f][e] = abs(p[0]) 
                    else:
                        zero[e][f] = zero[f][e] = 0.0

                elif F_transf(p[0]) > Z_fisher(confval, sample_size)\
                     and Power_r(p[0],confval,options.power)[0] <= sample_size:
                    zero[e][f] = zero[f][e] = p[0]
                else:
                    zero[e][f] = zero[f][e] = 0.0
            else:
                if p[1]<= (1-confval)\
                   and Power_r(p[0],confval,options.power)[0] <= sample_size:
                    zero[e][f] = zero[f][e] = p[0] 
                else:
                    zero[e][f] = zero[f][e] = 0.0


    return zero	


###############################################################################################

def GetThreshold(matrices):
    ''' Automatic threshold picking'''
    # Call if you need to get the vector of threshold
    thr = []
    for i in range(len(matrices)):
        for e in range(len(matrices[i])):
            thr.append(std(matrices[i][e]))
    return thr

###############################################################################################

def GetAglomeratedThreshold(lms):
    ''' Get threshold for the n dimensions'''
    # Call if you need to get the vector of threshold
    thr = []
    for i in range(len(lms)):
        for j in range(len(lms)):
            if i !=j and i < j and lms[i][j] != 0.00000:
                thr.append(lms[i][j])
        #thr.append(std(lms[i]))
    threshold = std(thr)#min(thr)
    f=open(prefix+'.parameters', 'a')
    t='Threshold mode: Auto = %f\n'%(threshold)
    print t
    f.write(t)
    f.close()
    return thr, threshold

###############################################################################################
def UseCov(matrices, data, threshold):
    '''
    Create a variance-covariance matrix
    '''
    thr=[]
    # Perform np magic
    for i in range(len(data)):
        l=88
        # Convert to np matrix
        matrices.append(matrix(data[i]))
        l+=i
        matrices[i] = cov(matrices[i].T)
        if threshold == 'Auto':
            #get the threshold from correlation matrix#
            for e in range(len(matrices[i])):
                thr.append(std(matrices[i][e]))

    return matrices, thr


###############################################################################################

def UseCorr(matrices, data, threshold):
    '''
    Use pearson correlation without a significant test
    '''
    thr=[]
    # Perform np magic
    for i in range(len(data)):
        l=88
        # Convert to np matrix
        matrices.append(matrix(data[i]))
        l+=i
        matrices[i] = corrcoef(matrices[i].T)
        if threshold == 'Auto':
            #get the threshold from correlation matrix#
            for e in range(len(matrices[i])):
                thr.append(std(matrices[i][e]))	

    return matrices, thr	


###############################################################################################

def Histocorr(prefix, method, matrices):
    '''
    Draw a histogram of the correlation/covariance matrix
    '''
    for e in range(len(matrices)):
        plt.subplot(len(matrices)-1,len(matrices)-1,e+1)
        mu, sigma = mean(matrices[e]), var(matrices[e])
        n, bins, patches = plt.hist(matrices[e], 50, normed=1, facecolor='green', alpha=0.75)
        y = mlab.normpdf( bins, mu, sigma)
        normline = plt.plot(bins, y, 'r--', linewidth=1)
        plt.xlabel('Correlation')
        plt.ylabel('Frequency')
        plt.title('Histogram of %s Coordinates'%(chr(88+e)))
    plt.savefig('%s_%s_Hist'%(prefix, method))	


###############################################################################################

def agglomerare_additive(matrices, absolutecov):
    '''
    Agglomerate landmark dimensions using euclidean distance
    '''
    lms = []
    for i in range(0,matrices[0].shape[0]):
        temp = []
        for j in range(0, matrices[0].shape[0]):
            sm = 0.0
            # Sum over all dimensions /home/lombard/LabBlouin/code/contactList/
            for m in matrices:
                d = m[i][j]
                if absolutecov:
                    d = abs(d)
                sm += (d**2)
            sq = sqrt(sm)
            temp.append(sq)
        lms.append(temp)	

    return lms


###############################################################################################

def agglomerare_mean(matrices, dim, absolutecov):
    '''Agglomerate landmark dimensions using average of correlation	'''
    lms = []
    for i in range(0,matrices[0].shape[0]):
        temp = []
        for j in range(0, matrices[0].shape[0]):
            sm = 0.0

            # Sum over all dimensions
            for m in matrices:
                d = m[i][j]
                if absolutecov:
                    d = abs(d)
                sm += (d/dim)

            temp.append(sm)
        lms.append(temp)	

    return lms


###############################################################################################

def Build_igraph(lms, options, gfilter):
    '''
    Build a graph, using igraph library
    '''
    g = Graph(len(lms))
    for i in range(len(g.vs)):
        g.vs[i]['label'] = i

    for i in range(len(lms)):
        for j in range(i+1,len(lms)): #for j in range(i,len(lms)):# just testing if redundant edges help modularity
            if lms[i][j] and lms[i][j] >= options.threshold:
                if options.contacts:
                    if (i,j) in gfilter:
                        g.add_edge(i,j)
                        g.es[len(g.es)-1]['wts'] = int(lms[i][j]*100000)
                else:
                    g.add_edge(i,j)
                    g.es[len(g.es)-1]['wts'] = int(lms[i][j]*100000)
                    #print i,j, int(lms[i][j]*100000)
    return g


###############################################################################################

def Graph_Components(g, lms):
    '''
    Clustering by components comstGreed, using igraph library
    Returns the memberships by cluster and a lis of Eigenvector centrality per cluster
    '''
    memgreed = [None]*len(lms)
    index = 0
    comp = g.components()
    #if len(comp) > 1:
    print "There are %d components in the network."%(len(comp))
    if len(comp) == len(lms):
        print 'Only singletons in the dataset'
        mem = range(len(comp))
        for e in range(len(memgreed)):
            memgreed[e] = e
    else:
        for cp in range(len(comp)):
            cpn = comp.subgraph(cp)
            if len(cpn.vs) > 1:
                mem = cpn.community_fastgreedy(cpn.es['wts']).as_clustering()
                #modmap = {}
                for i in mem:
                    for j in i:
                        memgreed[cpn.vs[j]['label']] = index
                    index += 1
            else:
                memgreed[comp[cp][0]] = index
                index += 1
    return memgreed, mem


###############################################################################################

def write_cluster(prefix, g, memgreed):
    ''' Write cluster to a file'''
    ou = open(prefix+'.parameters', 'a')
    chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y','Z',
              'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','y','z',
              '0','1','2','3','4','5','6','7','8','9']#, '~', '!', '@', '#', '$', '%', '^', '&', '(', ')', '-', '+', 
            #'=', '_', '<', '>', '/', "'", "\\", '"', '[', ']', '{', '}']	
    memgreed = VertexClustering(g,memgreed)
    m = memgreed.membership
    try:
        mod = g.modularity(m,g.es['wts'])
    except:
        mod = g.modularity(m)
    print "%d cluster(s) found."%(len(set(m)))
    print "Modularity score: %f"%(mod)
    ou.write("\nModularity score of the full graph: %f\n"%(mod))
    fout = open(prefix+'.graphcluster','w')
    if not len(set(m)) > len(chains):
        for i in m:
            fout.write(chains[i]+' ')
            print i,
        print '\n'
        fout.write('\n')
        fout.close()
    else:
        for i in m:
            fout.write(str(i)+' ')
            print i,
        print '\n'
        fout.write('\n')
        fout.close()	
    ou.close()


###############################################################################################

def clus2csv(prefix, data, sample_size, memvec = False):
    '''
    Split the data into csv files
    '''
    if not memvec:
        memvec = open(prefix+'.graphcluster').read()
    # Create the dictionaries where the output will be stored.
    fouts = {}
    #load the membership vector
    m=memvec
    mem=m.split()
    smem=set(mem)
    for i in smem:
        if not i == '?':
            fouts[i]=open(prefix+'.'+i+'.csv','w')	
    # Sort all landmarks into respective csv files ##################################
    for line in range(sample_size):
        for o in range(len(mem)):
            # Never include singletons
            if mem[o] == '?':
                continue	

            f = fouts[mem[o]]
            x = ''
            for j in range(len(data)):
                x += str(data[j][line][o]) + ';'
            f.write(x)
        for f in fouts:
            fouts[f].write('\n')
    for f in fouts:
        fouts[f].close()	

    return m, mem


###############################################################################################

def cluster_files():
    '''
    Create a list of files with individual clusters by moodularity.py
    '''
    clsfiles={}

    allfiles = os.listdir(os.getcwd())
    for name in allfiles:
        if not name.endswith('.csv'):
            continue
        if 'summary' in name:
            continue
        clsfiles[name[name.rfind('.',0,-4)+1:name.find('.csv')]]=name

    return clsfiles


###############################################################################################

def single_filter(memvec, clsfiles):
    '''
    Filter out the singletons and return a dictionary with the singletons, a list with the 
    singleton files, and clean cluster list.
    '''
    #filter singletons
    s=[x for x,v in Counter(memvec.split()).iteritems() if v == 1]
    singles={}
    singlesfiles=[]
    fls=clsfiles.keys()

    for items in fls:
        if items in s:
            singles[items]='?'
            singlesfiles.append(clsfiles[items])
            del clsfiles[items]
    return clsfiles, singles, singlesfiles


###############################################################################################

def fn2chain(clsfiles, chains, singlesfiles, itpermt):
    '''
    Map the file name with the respective chain (main use in protein data)
    '''	
    clskeys=[]
    clusters=len(clsfiles)
    for n in clsfiles:
        clskeys.append(n)
    if 0 or 1 in clskeys:
        for nu in range(len(clskeys)):
            newfilename= clsfiles[clskeys[nu]].replace(clskeys[nu],chains[nu])
            os.system('mv %s %s'%(clsfiles[clskeys[nu]],newfilename))
            clsfiles[chains[nu]]= newfilename
            del clsfiles[clskeys[nu]]
    # Exit the program when only singletons are found
    if clusters == 0 and not itpermt:
        print('Only singletons in the dataset')
        print('Bye Bye')
        for f in glob('*.csv'):
            if not 'summary' in f:
                os.remove(f)
        pass
    elif clusters == 0 and itpermt:
        print('Only singletons in the dataset in this iteration')

    # if only one cluster, test robustest against the grouped singletons
    if clusters == 1:
        lin=''
        print 'Only one cluster. Bye, Bye'
        pass

    return clsfiles, clusters


###############################################################################################

def print_memvec(m, singles):
    '''
    Read the membership vector, prints it to screen
    '''
    mem=m.split()
    me=''
    for item in mem:
        if singles.has_key(item):
            mem[mem.index(item)]=singles[item]
        me+=item + ' '
    print(me)

    return mem


###############################################################################################
def lms2Rlms(prefix, lms):
    row , col = len(lms), len(lms[0])
    r('lms <- matrix(0.0,nrow=%d,ncol=%d)'%(row,col))
    for i in range(row):
        for j in range(col):
            r('lms[%d,%d]<- %f'%(i+1,j+1,float(lms[i][j])))
            #print 'lms[%d,%d]<- %f'%(i+1,j+1,float(A[i][j]))
    r('write.table(lms,file="%s.lms",sep=";",row.names=FALSE,col.names=FALSE)'%(prefix))

def get_cluster_indices(m):
    M = m.split()
    s = set(M)
    d={}
    # get cluster indices

    for l in s:
        d[l]=[]
        for e in range(len(M)):
            if M[e] == l:
                d[l].append(e)
    return d

def filter_singletons(clusterdict,m):
    singles={}
    C = deepcopy(clusterdict)
    memb = deepcopy(m.strip())
    memb = ' '+memb+' '
    for key, value in C.iteritems():
        if len(value) == 1:
            singles[key]='?'
            del clusterdict[key]
            memb = memb.replace(' %s '%(str(key)),' ? ')

    return singles, memb[1:-1]

def AreModneighbours(A,indA,indB):
    ans = False
    for i in indA:
        for j in indB:
            if A[i][j] == 1:
                ans = True
    return ans

def Rbit_splitIntraInter(k,ke,C1,C2,options,pvals,cl,perm):
    # Split in intra and inter correlation
    r('a<-lms[c(%s),c(%s)]'%(C1,C1))
    r('a<-a[upper.tri(a)] ; a <- a[which(a > %f)]'%(options.threshold))
    if r('length(a)')[0] == 0:
        r('a<-as.vector(lms[c(%s),c(%s)])'%(C1,C1))
    r('b<-lms[c(%s),c(%s)]'%(C2,C2))
    r('b<-b[upper.tri(b)]; b <- b[which(b > %f)]'%(options.threshold))
    if r('length(b)')[0] == 0:
        r('b<-as.vector(lms[c(%s),c(%s)])'%(C2,C2))
    r('ab<-lms[c(%s),c(%s)]'%(C1,C2))
    r('ab <- ab[which(ab > %f)]'%(options.threshold))
    if r('length(ab)')[0] == 0:
        r('ab<-as.vector(lms[c(%s),c(%s)])'%(C1,C2))
    #Test A vs AB    
    pvalA=r('twotPermutation(a, ab,  nsim=%d, plotit=FALSE)'%(int(perm)))[0]
    cl['%svs%s_%s'%(k,k,ke)]=pvalA
    pvals.append(pvalA)
    print('%svs%s_%s'%(k,k,ke)+'.... Done')
    #Test B vs AB
    pvalB=r('twotPermutation(b, ab,  nsim=%d, plotit=FALSE)'%(int(perm)))[0]
    cl['%svs%s_%s'%(ke,ke,k)]= pvalB
    pvals.append(pvalB)
    print('%svs%s_%s'%(ke,ke,k)+'.... Done')
    return pvals,cl

def pairwise_lms_permt(D,perm,A,options):
    E = deepcopy(D)
    if '?' in E.keys():
        del E['?']
    cl={}
    keys=[]
    pvals=[]
    neighbours=[]
    for k,v in E.iteritems():
        keys.append(k)
        C1=','.join([str(x+1) for x in E[k]])
        indA=[int(x) for x in E[k]]
        for ke, va in E.iteritems():
            if not ke in keys:
                C2=','.join([str(x+1) for x in E[ke]])
                indB=[int(x) for x in E[ke]]
                if A:
                    if AreModneighbours(A,indA,indB):
                        pvals, cl = Rbit_splitIntraInter(k,ke,C1,C2,options,pvals,cl,perm)
                        neighbours.append((k,ke))
                    else:
                        continue
                else:
                    pvals, cl = Rbit_splitIntraInter(k,ke,C1,C2,options,pvals,cl,perm)

    return keys, pvals, cl, neighbours

def lms_permt(m,perm,A,options):
    print 'WARNING: USING LMS_PERMT... this is in BETA'
    D = get_cluster_indices(m)
    singles, m = filter_singletons(D,m)
    if len(D) == 0 and not itpermt:
        print('Only singletons in the dataset')
        print('Bye Bye')
        pass
    elif len(D) == 1:
        print 'Only one cluster. Bye, Bye'
        pass

    # print progress to screen and read in the membership vector
    print(str(len(D)*(len(D)-1))+' comparisons')
    print('Current membership vector:')
    mem = m
    print mem
    print('Progress:')
    keys, pvals, cl, neighbours = pairwise_lms_permt(D,perm,A,options)
    # FDR correction
    FDR = FDR_correction(pvals, options)
    # Perform the logical comparisons using FDR-corrected critical value
    scl = FDRc_sigtest(cl, FDR)
    popkeys={}
    for k in D.iterkeys():
        if not k == '?':
            popkeys[k]=k
    # Merge non-significant and reciprocal clusters
    mem = mem.split()
    newm = cluster_merger(D, scl, mem, popkeys, neighbours,A)

    #open Outfile and write result
    newvec = write_permt(prefix, FDR, cl, scl, newm)

    return scl, newm, newvec, FDR


def pairwise_permt(clsfiles, perm,options):
    '''
    Iterates over cluster files and perform a permutational t-test for each cluster pair 
    (A vs AB and B vs AB). Returns the keys of the comparisons, the dictionary of the numerical 
    comparisons (cl), and the list of p-values for 
    the test
    '''
    dim = options.dim
    #perm =  options.perm
    cl={}
    #Iterate over cluster files
    keys=[]
    pvals=[]
    for ke in clsfiles.iterkeys():
        if ke not in keys:
            keys.append(ke)       
        for key in clsfiles.iterkeys():
            if ke != key:
                # Read the csv files, split dimensions 
                D, ss = Read_GM(clsfiles[ke][:-4], dim, f='csv')
                anames = py2r_dim(D)
                for an in anames:
                    #change the name of dimensions for the first cluster
                    r('a_%s <- %s'%(an,an))
                D1, ss1 = Read_GM(clsfiles[key][:-4], dim, f='csv')
                bnames = py2r_dim(D1)
                for bn in bnames:
                    #change the name of dimensions for the second cluster
                    r('b_%s <- %s'%(bn,bn))	
                # R operator string : will store the call for the euclidean distance of correlation
                o = ''				
                # Iterate over dimensions, merge clusters, and estimate the correlation
                for d in range(dim):
                    r('merge_%s<-cbind(a_%s,b_%s)'%(anames[d], anames[d], bnames[d]))
                    r('C_%s <- cor(merge_%s)'%(anames[d], anames[d]))
                    o += '(C_%s^2)+'%(anames[d])
                o = o[:-1]
                # Estimate determination coefficient per dimension and add them and sqroot them
                r('R<-sqrt(%s)'%(o))
                #Ignore the diagonal
                r('diag(R) <- NA')
                # Split in intra and inter correlation
                r('a<-R[1:(dim(a_%s)[2]),1:(dim(a_%s)[2])]'%(anames[d], anames[d]))
                #r('a<-na.omit(as.vector(a))')
                r('a<-a[upper.tri(a)]')
                r('b<-R[((dim(a_%s)[2])+1):dim(R)[2],((dim(a_%s)[2])+1):dim(R)[2]]'%(anames[d], anames[d]))
                #r('b<-na.omit(as.vector(b))')
                r('b<-b[upper.tri(b)]')
                r('ab<-R[((dim(a_%s)[2])+1):dim(R)[2],1:(dim(a_%s)[2])]'%(anames[d], anames[d]))
                r('ab<-na.omit(as.vector(ab))')
                if options.histocorr:
                    #print matrix histograms
                    r('png("%s_%svs%s_permtHist.png",  width = 1200, height = 1200, units = "px")'%(prefix,key,ke))
                    r('par(mfrow=c(2,2), ps=16)')
                    r('hist(R)')
                    r('hist(a)')
                    r('hist(b)')
                    r('hist(ab)')
                    r('dev.off()')
                #Test A vs AB    
                pvalA=r('twotPermutation(a, ab,  nsim=%d, plotit=FALSE)'%(int(perm)))[0]
                cl['%svs%s_%s'%(ke,ke,key)]=pvalA
                pvals.append(pvalA)
                #print('%svs%s_%s'%(ke,ke,key)+'.... Done')
                #Test B vs AB
                pvalB=r('twotPermutation(b, ab,  nsim=%d, plotit=FALSE)'%(int(perm)))[0]
                cl['%svs%s_%s'%(key,key,ke)]= pvalB
                pvals.append(pvalB)
                print('%svs%s_%s'%(key,key,ke)+'.... Done')

    return keys, pvals, cl

###############################################################################################

def FDR_correction(pvals, options):
    '''
    Compute the False Discovery Rate correction for the critical value
    '''
    confval = options.confval
    #sort and reverse the list for FDR correction
    pvals.sort()
    pvals.reverse()
    if len(pvals)<=2:
        FDR=1-confval
    else:
        for el in pvals:
            if float(el) <= ((pvals.index(el)+1)/len(pvals))*(1-confval):
                FDR=((pvals.index(el)+1.0)/len(pvals))*(1-confval)
                break
            else:
                continue
        try:
            FDR
        except:
            FDR=((1-confval)*(len(pvals)-1))/(2*len(pvals))	

    return FDR


###############################################################################################

def FDRc_sigtest(cl, FDR):
    '''
    Perform the logical significance test using FDR corrected critical value. Returns a binary
    dictionary of the comparisons
    '''
    scl={}
    # loop over cl dictionary to write binary permt output into a dictionary
    for k, v in cl.iteritems():
        if float(v) <= FDR:
            scl[k] = True
        else:
            scl[k] = False	

    return scl

###############################################################################################
def relabel_clusters_new(singles, mem, clsfiles, chains):
    ''' Create a dictionary of new cluster names and update the 
        file map accordingly.
    '''
    # Output 1 the new map
    newmap = {'?':'?'}
    ncls={}
    for i in range(len(mem)):
        # Ignore Singletons 
        if mem[i] in singles or mem[i] == '?':
            newmap[mem[i]] = '?'
            continue

        # Not a singleton, already mapped?
        if not mem[i] in newmap.values():
            # use the next available chains, starting at chains[0]
            #newchain = chains[len(newmap)-1]
            newmap[mem[i]] = mem[i]
            # Update the map of files
            #clsfiles[newchain] = clsfiles[mem[i]]
            #ncls[newchain] = clsfiles[mem[i]]
            '''
	# Delete extra keys in the file map
	extras = list(set(clsfiles.keys()).difference(set(newmap.keys())))
	for i in extras:
		del clsfiles[i]'''

    return newmap, clsfiles



def relabel_clusters(singles, mem, clsfiles, chains):
    '''
    Label all singletons with ? and re-label the clusters with letters 
    '''
    popkeys={}
    # get the files keys
    clskeys = []
    for k in clsfiles.iterkeys():
        clskeys.append(k)
    for item in mem:
        if singles.has_key(item):
            mem[mem.index(item)]=singles[item]
        else:
            if item == '?':
                continue
            mem[mem.index(item)]=chains[clskeys.index(item)]
            if item in popkeys.itervalues():
                continue
            else:
                clsfiles[chains[clskeys.index(item)]]= clsfiles.pop(item)
                popkeys[chains[clskeys.index(item)]] = item	

    return popkeys, clsfiles

###############################################################################################

def merge(clsfiles, scl, popkeys, neighbours,A):
    '''
    Express merging events. Return a list of tuples of merging events
    '''
    merge=[]
    if not A:
        ncls=[]
        for j in clsfiles.iterkeys():
            for k in clsfiles.iterkeys():
                if not j == '?' and not k == '?':
                    try:
                        if ord(j) > ord(k):
                            ncls.append((j,k))
                    except:
                        if j > k:
                            ncls.append((j,k))
    else:
        ncls = neighbours
    for l in ncls:
        i,e = l
        if scl['%svs%s_%s'%(popkeys[i],popkeys[i],popkeys[e])] == scl['%svs%s_%s'%(popkeys[e],popkeys[e],
                                                                                   popkeys[i])] == False:
            merge.append((i,e))			

    return merge

###############################################################################################
def equiclus(merge):
    '''
    Look for equivalen clusters, and return a list with sets of equivalent clusters
    '''

    # List of sets of equivalent clusters
    newcl = []	
    for m in merge:
        m = set(m)
        merged = False
        # look for a newcl item which overlaps then merge
        for nc in range(len(newcl)):
            xnc = newcl[nc]
            if m.intersection(xnc):
                # Merge into xnc
                newcl[nc] = xnc.union(m)
                merged = True
                break
        if not merged:
            # New cluster to add to newcl
            newcl.append(m)

    return newcl


###############################################################################################

def rename_clusters(newcl, mem):
    '''
    Rename clusters using the smallest label in each set. Returns the new
    membership vector
    '''
    babel = {}
    for cluster in newcl:
        newname = min(cluster)
        for item in cluster:
            babel[item] = newname

    for i in range(len(mem)):
        if mem[i] in babel:
            mem[i] = babel[mem[i]]		

    # Create the new membership vector
    newm = ''
    for ite in mem:
        newm += ite + ' '	

    return newm


###############################################################################################

def cluster_merger(D, scl, mem, popkeys, neighbours,A):
    '''
    Merge non-significant and reciprocal clusters
    '''
    # Get the merging events
    merg = merge(D, scl, popkeys, neighbours,A)
    # List of sets of equivalent clusters
    newcl = equiclus(merg)
    # Rename clusters using the smallest label in each set
    newm = rename_clusters(newcl, mem)

    return newm

###############################################################################################

def write_permt(prefix, FDR, cl, scl, newm,count=''):
    '''
    Write and print the output of the permutation test, the new membership vector, and do some 
    cleanup
    '''
    # Rename the previous graphcluster (the original of modularity)
    os.system('mv %s %s'%(prefix+'.graphcluster',prefix+'.%sold_graphcluster'%(str(count))))
    # Open the new graphcluster file and the permt output file
    clus=open(prefix+'.graphcluster','w')
    fname=prefix+'.permt'
    fout=open(fname,'w')
    # Some lines of the output
    line0='######## Significance test of Clustering  ########\n'
    line1='#'*8+' Test of Clustering (per_t-test)  ' + '#'*8 + '\n'
    n=('#'*50)+'\n'
    print(n+line1+n)
    fout.write(n+line1+n)
    FDRc='FDR-corrected Critical value = '
    print(FDRc+str(FDR)+'\n')
    fout.write(FDRc+str(FDR)+'\n')
    theader='Comparison\tp-value\tSignificant?\n'
    print(theader)
    fout.write(theader)
    # Write the comparisons
    if cl:
        for k in cl:
            print k + '\t' + str(cl[k]) + '\t' + str(scl[k])
            fout.write(k + '\t' + str(cl[k]) + '\t' + str(scl[k])+'\n')
    else:
        cline='No adjacent clusters in the graph. Pvalue < 10^-4 or singletons.'
        print cline
        fout.write(cline)
    print(n)
    fout.write(n)
    # Print and write the new membership vector
    clustn = len(set(newm.split()))
    if '?' in newm:
        clustn = clustn-1
    print str(clustn) + ' significant clusters'
    print 'New membership vector:'
    o=newm.split()
    newvec=''
    for y in o:
        newvec+= str(y)+' '
    print newvec
    clus.write(newvec)
    clus.close()

    #remove previous CSV files... Clean up
    os.system('rm ./*.csv')	

    return newvec
###############################################################################################

def permt(prefix, data, matrices, sample_size, m, perm, options, count=''):
    '''
    test the significance of the clustering using the intra vs intercorrelation. The test is
    a permutational t-test.
    '''	
    dim, confval = options.dim, options.confval

    # Labels for chains according to PDB standards. the X,and ? are removed since they represent non-homologous
    # residues and singletons, respectively.
    chains = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y','Z',
              'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','y','z',
              '0','1','2','3','4','5','6','7','8','9', '~', '!', '@', '#', '$', '%', '^', '&', '(', ')', '-', '+', 
              '=', '_', '<', '>', '/', "'", "\\", '"', '[', ']', '{', '}']

    ########################## Statistical test of clustering ##########################
    # Create a list of files with individual clusters
    clsfiles = cluster_files()

    #filter singletons
    clsfiles, singles, singlesfiles = single_filter(m,clsfiles)

    #Change the file name to chain name
    clsfiles, clusters = fn2chain(clsfiles, chains, singlesfiles, options.itperm)

    if not len(singles) == 0:
        print 'There are %d singletons in the dataset. They will be grouped with the label "?".\n'%(len(singles))

    # print progress to screen and read in the membership vector
    print(str(clusters*(clusters-1))+' comparisons')
    print('Current membership vector:')
    mem = print_memvec(m, singles)
    print('Progress:')

    # perform pair-wise permutation test
    keys, pvals, cl = pairwise_permt(clsfiles, perm, options)

    # FDR correction
    FDR = FDR_correction(pvals, options)

    # Perform the logical comparisons using FDR-corrected critical value
    scl = FDRc_sigtest(cl, FDR)

    # Label all singletons with ? and re-label the clusters with letters 
    popkeys, clsfiles = relabel_clusters_new(singles, mem, clsfiles, chains)

    # Merge non-significant and reciprocal clusters
    newm = cluster_merger(D, scl, mem, popkeys, neighbours,A)

    #open Outfile and write result
    newvec = write_permt(prefix, FDR, cl, scl, newm)

    return scl, newm, newvec, FDR


###############################################################################################
def iterative_permt(prefix, data, matrices, m, sample_size,A, options):
    '''
    perform a permutation test iteratively, until a stable membership vector is reached
    '''
    memvec=[]
    count=0
    curr_m=deepcopy(m)
    while count <= 20:
        count+=1
        it= 'Performing iteration number %d'%(count)
        print '\n\n%s\n%s\n%s\n'%('*'*(len(it)),it,'*'*(len(it)))
        scl, newm, newvec, FDR = lms_permt(curr_m,options.itperm,A,options)
        '''scl, newm, newvec, FDR = permt(prefix, data, matrices, sample_size, curr_m, options.itperm, options, count=count)'''
        memvec.append(newm)
        if count >= 2:
            old = VectorToEdgeList(memvec[-2])
            new = VectorToEdgeList(memvec[-1])
            fsc = Fscore(Specificity(old,new), Sensitivity(old,new))
            if fsc >= 0.9999:
                break
        curr_m, mem = clus2csv(prefix,data, sample_size,newm)

    return scl, newm, newvec, FDR


###############################################################################################

def Dot_graph(prefix, scl, newm):
    '''
    Draw a graph with the significance difference among intra and intercorrelation as edges
    '''
    clus=[]
    cl = set(newm.split())
    for i in cl:
        clus.append(i)
    # Create a Dot file for the graphical representation of cluster's relationships
    dot=open(prefix+'.dot','w')
    #write the first line of the dot file
    dot.write('digraph finite_state_machine {\n')
    for e in clus:
        if not e == '?':
            dot.write('\t%s\n'%(e))
    #Iterate over cluster files	
    for ke, va in scl.iteritems():
        if va == False:
            dot.write('\t"%s" -> "%s";\n'%(ke.split('vs')[0],ke.split('vs')[-1].split('_')[-1]))
    dot.write('}')
    dot.close()
    #execute dot to create the graph and output a svg image
    os.system('dot -Tpng -o%s %s'%((prefix+'_clusgraph'+'.png'),(prefix+'.dot')))


###############################################################################################

def GetCoordinates(dim,index):
    '''
    Fetch the appropriate coordinates in the index
    '''	
    dist=0.0
    coords=[]
    #dims=[]
    for d in dim:
        M = (matrix(data[d])).T
        coor_i = mean(M[index])
        coords.append(coor_i)

    '''if dim == 2:
		xmatrix=(matrix(data[0])).T
		ymatrix=(matrix(data[1])).T
		coorX=mean(xmatrix[index])
		coorY=mean(ymatrix[index])
		coords.append(coorX)
		coords.append(coorY)
	elif dim == 3:
		xmatrix=(matrix(data[0])).T
		ymatrix=(matrix(data[1])).T
		zmatrix=(matrix(data[2])).T
		coorX=mean(xmatrix[index])
		coorY=mean(ymatrix[index])
		coorZ=mean(ymatrix[index])
		coords.append(coorX)
		coords.append(coorY)
		coords.append(coorZ)
	else:
		print 'Are you sure your data is in cartesian coordinates?'
'''
    return coords


def VectorToEdgeList(v):
    ''' Convert a membership vector to a list of edges.
        The membership vector must be a string.
    '''
    # Discretize list
    v = v.split()

    # Output vector
    out = []

    # Process edges only once (no symmetry)
    for i in range(len(v)):
        for j in range(i,len(v)):
            if v[i] == v[j]:
                out.append('%s.%s'%(i,j))

    # return edge list
    return set(out)


###############################################################################################

def Specificity(ref, test):
    '''
       Compute specificity from two edgelists
    '''
    # Size of intersection
    A = len(ref.intersection(test))

    return float(A) / len(test)


###############################################################################################

def Sensitivity(ref, test):
    '''
       Compute sensitivity from two edgelists
    '''
    # Size of intersection
    A = len(ref.intersection(test))

    return float(A) / len(ref)


###############################################################################################

def Fscore(sp, sn):
    ''' compute the F-score '''
    return 2*sp*sn/(sp+sn)


###############################################################################################

def F_transf(r):
    '''
    Compute the Fisher transformation of correlation
    '''	
    z = 0.5 * (log((1+r)/(1-r)))
    return z


###############################################################################################

def Z_fisher(confval, sample_size):
    '''
    	Compute the sample - corrected Z_alpha for hypotesis testing 
    	of Fisher transformation of correlation
    '''		
    z_critic = norm.isf(1-confval) / sqrt(sample_size-3)
    return z_critic


###############################################################################################

def Power_r(corr, confval, power):
    '''
    Compute the power of the correlation using the Z' trasnformation of correlation coefficient:
    Z'=arctang(r)+r/(2*(n-1)) (see Cohen (1988) p.546). Require the R package PWR written by
    Stephane Champely <champely@univ-lyon1.fr>. It will return the required n fo the power and
    significance chosen.
    '''
    if corr  == 0.0 or (corr < 0.01 and corr > -0.001) :
        corr = 0.0012
    elif corr < -0.001:
        corr =  corr
    elif corr >= 0.99:
        corr =  0.99			
    else:
        corr = corr	
    r('p=pwr.r.test(power=%f,r=%f, sig.level=%f)'%(float(power), float(corr), float(1-confval)))
    n = r('p[1]')[0]

    return n


###############################################################################################
'''
def Clus_Power(prefix, data, confval, power, sample_size, dim, save_r=False):

	#Make a power analysis of the clusters by igraph and outputs a table with the proportion 
	#of elements above the n required according to significance, power and correlation

	# alternative dictionary
	proportion={}
	# Open the output file
	ouf=open('%s_Power_analysis.asc'%(prefix),'w')
	ouf.write('#################################################################\n')
	ouf.write('#Power analysis and description of the intracluster correlation:#\n')
	ouf.write('#################################################################\n\n')
	ouf.write('1 - Significance(alpha)= %f;\t1 - Power(beta):%f\n'%(1-confval,1-power))
	ouf.write('Sample Size : %d\n\n'%(sample_size))
	# Create a list of files with individual clusters by moodularity.py
	cls = []
	m,mem = clus2csv(prefix,data, sample_size)
	allfiles = os.listdir(os.getcwd())
	for name in allfiles:
		if not name.endswith('.csv'):
			continue
		cls.append(name)

	for clusters in cls:
		# get name of the cluster
		clname = clusters[clusters.find('.')+1:clusters.find('.csv')]
		ouf.write('Cluster %s:\n'%(clname))
		ouf.write('%s'%('*'*len(clusters))+'\n')
		# Read the csv files split dimension and turn them into R matrices
		D, ss = Read_GM(clusters[:-4], dim, f='csv')
		mnames = py2r_dim(D)
		b=''
		# estimate the correlation for each dimension
		for mn in mnames:
			r('c%s<-cor(%s)'%(mn,mn))
			# get the lowertriangle of the correlation
			r('c%s<-c%s[lower.tri(c%s)]'%(mn,mn,mn))
			b += 'c%s'%(mn)+','
		# Bind the dimensions
		b = b[:-1]
		c = r('c <- cbind(%s)'%(b))
		# estimate quantiles
		r('q<-quantile(c)')
		q = r('q')		

		if save_r: # Save some objects in an R file for R manupulation (outside python)
			r('save(file="power",c,q,%s)'%(b))

		# estimate the proportion of variables with enough power
		enough=0.0
		for x in c:
			nsam = Power_r(x, confval, power)
			#print x, nsam[0]
			if int(round(nsam[0])) <= sample_size:
				enough += 1.0
		prop = enough/len(c)
		# write to file
		ouf.write(' 0% \t\t 25% \t\t 50% \t\t 75% \t\t 100% \t\t PVP \n')
		resline = ' %f \t %f \t %f \t %f \t %f \t %f \n\n'%(round(q[0],3),round(q[1],3),round(q[2],3),round(q[3],3),round(q[4],3),round(prop,3))
		ouf.write(resline)
		proportion[clname]=prop
	ouf.write('_'*len(resline)+'\n')
	ouf.write('Percentages are the quantiles of the lower triangle of the correlation matrix\n')
	ouf.write('PVP: Proportion of variables with enough statistical power \n')
	for f in glob('*.csv'):
		if not 'summary' in f:
			os.remove(f)

	return proportion
'''
def Clus_Power(prefix, data, confval, power, sample_size, dim, m, save_r=False):
    '''
    Make a power analysis of the clusters by igraph and outputs a table with the proportion 
    of elements above the n required according to significance, power and correlation
    '''
    # alternative dictionary
    proportion={}
    # Create a list of files with individual clusters by moodularity.py
    cls = []
    m,mem = clus2csv(prefix,data, sample_size,m)
    allfiles = os.listdir(os.getcwd())
    for name in allfiles:
        if not name.endswith('.csv'):
            continue
        cls.append(name)
    if cls:
        # Open the output file
        ouf=open('%s_Power_analysis.asc'%(prefix),'w')
        resline = '#################################################################\n'
        resline+= '#Power analysis and description of the intracluster correlation:#\n'
        resline+= '#################################################################\n\n'
        resline+= '1 - Significance(alpha)= %f;\t1 - Power(beta):%f\n'%(1-confval,1-power)
        resline+= 'Sample Size : %d\n\n'%(sample_size)
        for clusters in cls:
            # get name of the cluster
            clname = clusters[clusters.find('.')+1:clusters.find('.csv')]
            resline+= 'Cluster %s:\n'%(clname)
            resline+='%s'%('*'*len(clusters))+'\n'
            # Read the csv files split dimension and turn them into R matrices
            D, ss = Read_GM(clusters[:-4], dim, f='csv')
            mnames = py2r_dim(D)
            b=''
            # estimate the correlation for each dimension
            for mn in mnames:
                r('c%s<-cor(%s)'%(mn,mn))
                # get the lowertriangle of the correlation
                r('c%s<-c%s[lower.tri(c%s)]'%(mn,mn,mn))
                b += 'c%s'%(mn)+','
            # Bind the dimensions
            b = b[:-1]
            c = r('c <- cbind(%s)'%(b))
            # estimate quantiles
            r('q<-quantile(c)')
            q = r('q')		

            if save_r: # Save some objects in an R file for R manupulation (outside python)
                r('save(file="power",c,q,%s)'%(b))

            # estimate the proportion of variables with enough power
            enough=0.0
            for x in c:
                nsam = Power_r(x, confval, power)
                #print x, nsam[0]
                if int(round(nsam[0])) <= sample_size:
                    enough += 1.0
            prop = enough/len(c)
            # write to file
            tr = ' 0% \t\t 25% \t\t 50% \t\t 75% \t\t 100% \t\t PVP \t\t nvar \n'
            resline+= tr +' %f \t %f \t %f \t %f \t %f \t %f \t %d/%d \n\n'%(round(q[0],3),round(q[1],3),
                                                                             round(q[2],3),round(q[3],3),round(q[4],3),
                                                                             round(prop,3), enough, len(c))
            #ouf.write(resline)
            proportion[clname]=prop
        resline+= '_'*len(tr)+'\n'
        resline+= 'Percentages are the quantiles of the lower triangle of the correlation matrix\n'
        resline+= 'PVP: Proportion of variables with enough statistical power \n'
        resline+= 'nvar: Number of variables with enough power within the cluster \n'
        for f in glob('*.csv'):
            if not 'summary' in f:
                os.remove(f)
        ouf.write(resline)
        print resline
    else:
        print 'Only singletons in the dataset. Either correlation zero or not enough power to resolve'\
              'a negligible correlation.'

    return proportion

###############################################################################################

def py2r_dim(data):
    '''
    Turn the list of list of data (as in per-dimension list) into dim R matrices. Each entry in
    data should be the dimension. e.g. data[0] dimension 1, data[1]. It will return a list of 
    matrices names in R.
    '''
    # set a list of possible matrices names
    l= map(chr, range(65, 91))
    #l.extend(map(chr, range(97, 123)))
    l.reverse()
    mnames=[]
    # iterate over dimensions
    for i in range(len(data)):
        # Start an R vector
        v='c('
        # Get the dimensions of the matrix
        row = len(data[i])
        col = len(data[i][0])
        for j in range(len(data[i])):
            ro = str(data[i][j])[1:-1]
            v += ro + ','
        v = v[:-1]+')'
        r('%s <- matrix(%s,%d,%d, byrow=TRUE)'%(l[i],v,row,col))
        mnames.append(l[i])

    return mnames 


###############################################################################################

class landmark:
    def __init__(self, ID=-1):
        # The vector of data items. In the case of morphometric, this will 
        # then be a list of lists of triplets of floats [[x,y,z],[x,y,z]]
        # Where each element is a different instance of the landmark in a distinct sample.
        self.data = []

        # An identifier so the instance knows its ID in the application
        self.ID = ID

        # A list of pointers to neighbours
        self.neighbours = []

        # Private variable to save on compute time during clustering
        self.__cluster = None

    def __str__(self):
        return 'lm%d'%(self.ID)

    def __repr__(self):
        return str(self)

    def AddSample(self, v):
        ''' Add a xyz triplet to the internal data. Ensures that the data is in float format.
        '''
        self.data.append(v)

    def __len__(self):
        # returns the number of samples in the landmark
        return len(self.data)

    def GetCluster(self):
        return self.__cluster

    def SetCluster(self, c = None):
        self.__cluster = c

    def Centroid(self):
        # Return the average coordinate for this landmark
        out = [0.0]*len(self.data[0])

        n = float(len(self))

        for i in self.data:
            for j in range(len(i)):
                out[j] += i[j]/n

        return out

    def GetNeighbours(self):
        # Accessor method for the neighbours.
        return self.neighbours

    def AddNeighbour(self, N):
        # Add a neighbour to the list
        if not N in self.neighbours:
            self.neighbours.append(N)

    def OutOfClusterNeighbours(self):
        ''' Returns a lis of neighbours that are not in the same cluster '''
        out = []
        for i in self.neighbours:
            if self.GetCluster() != i.GetCluster():
                out.append(i)
        return out


###############################################################################################


class GMdata:
    def __init__(self, filename, dimension = 3):
        ''' Create a GM data object to manipulate entire datasets '''
        self.D = dimension
        # assume a .gm file extension
        self.filenameprefix = filename[:-3]

        # Landmark structure as a list. Landmarks will be indexed as integer.
        self.landmarks = []
        self.samplenames = []

        # clusters of landmarks (no particular indexing here)
        self.clusters = []
        self.singletons = []

        # Read in the data
        self.ReadFile(filename)

        # Build the connectivity between nodes
        #self.BuildContactGraph()


    def ReadFile(self, filename):
        ''' Read in the file and map landmarks to their names in the files.'''
        # Original cluster
        self.clusters.append(landmark_cluster())

        # Figure out how many landmark to create
        new = True
        with open(filename) as fin:
            for line in fin:
                line = line.split(';')

                # Last character is the name
                self.samplenames.append(line[0][-1])

                line = line[1:-1]            

                if new:
                    # Create the landmarks
                    for i in range((len(line))/self.D):
                        landm = landmark(len(self.landmarks))
                        self.landmarks.append(landm)
                        self.clusters[0].AddMember(landm)
                    new = False

                # Read line and add sample to each landmark
                for lm in range(len(self.landmarks)):
                    v = [float(line[lm*self.D])]
                    for i in range(1,self.D):
                        v.append(float(line[lm*self.D+i]))

                    # Get landmark
                    landm = self.landmarks[lm]

                    # Add data
                    landm.AddSample(v)

    def ReadGraph(self, filename):
        ''' Read in the graph and connect the landmarks
        '''
        # Opens the file
        with open(filename) as fin:
            for line in fin:
                line = line.split()
                lm = self.landmarks[int(line[0])]
                for n in line[1:]:
                    n = self.landmarks[int(n)]
                    lm.AddNeighbour(n)
                    n.AddNeighbour(lm)

    def BuildContactGraph(self, cutoff = 10):
        ''' Reciprocally add to neighbour list if two centroid's distance is
            equal or less than cutoff.
        ''' 
        ct = []
        for i in self.landmarks:
            ct.append(i.Centroid())

        # Upper triangle traversal
        for i in range(len(ct)):
            a = ct[i]
            for j in range(i+1,len(ct)):
                b = ct[j]
                d = 0.0
                for dim in range(len(a)):
                    d += (a[dim]-b[dim])**2
                d = d**0.5
                if d <= cutoff:
                    self.landmarks[i].AddNeighbour(self.landmarks[j])
                    self.landmarks[j].AddNeighbour(self.landmarks[i])

        # Ensures that all landmarks have at least 1 connection to the graph t
        for i in range(len(self.landmarks)):
            lm = self.landmarks[i]
            a = ct[i]
            if not lm.neighbours:
                # finds the nearest neighbour (possibly buggy)
                nn = None
                nnd = 1000000
                for j in range(i+1,len(self.landmarks)):
                    b = ct[j]
                    d = 0.0
                    for dim in range(len(a)):
                        d += (a[dim]-b[dim])**2
                    d = d**0.5

                    if d < nnd:
                        nn = self.landmarks[j]
                        nnd = d
                lm.AddNeighbour(nn)
                nn.AddNeighbour(lm)


    def RandomMerge(self):
        ''' Pick two clusters and merge them. 
        '''
        # Do not merge if there is only one cluster
        if len(self.clusters) == 1:
            return

        # Pick on random cluster
        random.shuffle(self.clusters)
        A = self.clusters[0]

        try:
            B = A.RandomAdjacentCluster()
        except:
            return False

        # Perform the merger
        A.MergeWith(B)

        # Clean up
        self.clusters.remove(B)

    def RandomSplit(self, N=2):
        ''' Get a random cluster, split it.
        '''
        # Pick on random cluster
        random.shuffle(self.clusters)
        A = self.clusters[0]

        # Split routine
        newclusters = A.RandomPartition(N)

        # replace
        self.clusters.remove(A)
        self.clusters.extend(newclusters)

    def RandomPointMutation(self):
        ''' Swap one landmark from one cluster to another
        '''
        # Pick on random clusterd
        random.shuffle(self.clusters)
        A = self.clusters[0]

        # Random Absorption
        return A.RandomAbsorption() 

    def ToVector(self):
        ''' returns a vector encoding of the partitonning
        ''' 
        # create a vector to map landmark into cluster
        out = []
        for i in self.landmarks:
            out.append(self.clusters.index(i.GetCluster()))
        return out

    def AssertToLandmarks(self, cluster_map=None):
        ''' Rebuild the partions using clusters as a guide.
        '''
        # Create the clusters
        self.clusters = [] 
        if cluster_map == None:
            # Reset to one big cluster
            cluster_map = [0]*len(self.landmarks)
        for i in range(len(set(cluster_map))):
            self.clusters.append(landmark_cluster())

        for i in range(len(cluster_map)):
            lm = self.landmarks[i]
            cl = self.clusters[cluster_map[i]]
            cl.AddMember(lm)

        # Delete empty clusters
        for i in range(len(self.clusters)-1,-1,-1):
            if len(self.clusters[i].members) == 0:
                self.clusters.remove(self.clusters[i])


def write_singleMatrix(prefix,matrix):
    '''
    Write to file a list of list type of matrix
    '''
    inf = open(prefix+'.mat','w')
    for e in matrix:
        t = ''
        for i in range(len(e)):
            t += str(e[i])
            t += ';'
        inf.write(t[:-1]+'\n')

    inf.close()

def pairwise_RV(clusfilename1,clusfilename2):
    '''
    Compute the RV for two given clusters files. The cluster files are cvs semicolon
    separated as created in modulerV2
    '''
    decimal.getcontext().prec = 3
    rv=()
    r('C1 <- read.table("%s", sep=";") ; C1 <- as.matrix(C1[1:dim(C1)[2]-1])'%(clusfilename1))
    r('C2 <- read.table("%s", sep=";") ; C2 <- as.matrix(C2[1:dim(C2)[2]-1])'%(clusfilename2))
    r('rv <- coeffRV(C1,C2)')
    rv = rv + (round(r('rv$rv')[0],3),)
    rv = rv + (round(r('rv$rvstd')[0],3),)
    rv = rv + (decimal.Decimal(str(r('rv$p.value')[0])).normalize(),)
    return rv

def multiple_RVs(prefix,dim=1):
    '''
    Calculate the RVs betwen all the clusters in a membership vector
    '''
    # Load Data ############################################################
    data, sample_size = Read_GM(prefix, dim)

    # Get the individual files per cluster
    m,mem = clus2csv(prefix,data, sample_size)
    cls=glob('*.csv')

    #iterate over clusters, calculate rvs and store them
    RVs={}
    for c1 in range(len(cls)):
        for c2 in range(c1+1,len(cls)):
            RVs['%svs%s'%(cls[c1][cls[c1].find('.')+1:cls[c1].rfind('.')],cls[c2][cls[c2].find('.')+1:cls[c2].rfind('.')])] = pairwise_RV(cls[c1],cls[c2])
    for f in glob('*.csv'):
        if not 'summary' in f:
            os.remove(f)
    return RVs

def write_RVfile(prefix,RVs):
    '''
    writes the RV output to file
    '''
    L= 'Calculating RVs:'
    outf = open(prefix+'.RVs','w')
    L1 = '#'*56+'\n'+'#'*10+' RV coefficient test for modularity '+'#'*10+'\n' + '#'*56+'\n'
    L2 = 'Ho = The configurations are independent.\n'
    L3 = 'Comparison\t  RV\t  Standardized RV\t p-value\n'
    outf.write(L1+L2+L3)
    print L1, L2, L3
    for k,v in RVs.iteritems():
        L4 = k+' \t'+str(v[0])+'\t\t'+str(v[1])+'\t\t '+str(v[2])+'\n'
        outf.write(L4)
        print L4
    outf.write('#'*56)
    print '#'*56
    outf.close()


def centralities(prefix, memgreed,overall,g):
    '''
    will read a landmark file, copy it, including extra columns with the centrality information
    '''
    fout = open(prefix+'.centrality','w')
    #index,names,rest = parse_landmark(prefix)
    names,rest = parse_landmark(prefix)
    evcen,btcen,clcen,degree = get_centralities(memgreed,overall,g)
    for n in range(len(names)):
        evcent,btcent,clcent,cdegree = deepcopy(evcen), deepcopy(btcen), deepcopy(clcen), deepcopy(degree)
        fout.write('>'+names[n]+'\n')
        if not overall:
            for e in range(len(memgreed)):
                fout.write(rest[n][e]+'\t'+str(evcent[memgreed[e]].pop(0)) + '\t' +\
                           str(btcent[memgreed[e]].pop(0)) + '\t' + str(clcent[memgreed[e]].pop(0)) + \
                           '\t' + str(cdegree[memgreed[e]].pop(0))+'\n')
        else:
            for e in range(len(memgreed)):
                fout.write(rest[n][e]+'\t'+str(evcent[0][e]) + '\t' + str(btcent[0][e]) + '\t' +\
                           str(clcent[0][e]) + '\t' + str(cdegree[0][e])+'\n')
    fout.write('The organization of the variables is: \nIndex\tResidue index\tAminoAcid\tEigenvalue Centrality\tBetweenness centrality\tCloseness centrality \t Degree')
    fout.close()

def parse_landmark(prefix):
    ''' 
    Parse the landmark file and return two lists with the name and the residues in that shape
    '''

    inf = open(prefix+'.landmarks').read().split('>')
    #ind = []
    nam = []
    res =[]
    for ins in inf:
        if ins == '':
            continue
        #indices = []
        names=ins[:ins.find('\n')]
        rest=[]
        inst = ins.split('\n')
        for line in inst:
            if line == '' or not '\t' in line:
                continue
            else:
                #indices.append(line[0])
                rest.append(line)
        #ind.append(indices)
        nam.append(names)
        res.append(rest)
    #return ind,nam,res
    return nam,res

def get_centralities(memgreed,overall,g):
    ''' store the centralities measures in lists '''
    mem = []
    for e in set(memgreed):
        mem.append([x for x in range(len(memgreed)) if memgreed[x] == e])
    evcen,btcen,clcen,degree=[],[],[],[]
    if not overall:
        for e in mem:
            sg = g.subgraph(e)
            if len(e) == 1:
                evcen.append(sg.evcent())
                btcen.append(sg.betweenness())
                clcen.append(sg.closeness())
                degree.append(sg.degree())
            else:
                evcen.append(sg.evcent(weights=sg.es['wts']))
                btcen.append(sg.betweenness(weights=sg.es['wts']))
                clcen.append(sg.closeness(weights=sg.es['wts']))
                degree.append(sg.degree())				
    else:
        if e != 1:
            evcen.append(g.evcent(weights=g.es['wts']))
            btcen.append(g.betweenness(weights=g.es['wts']))
            clcen.append(g.closeness(weights=g.es['wts']))
            degree.append(g.degree())
        else:
            evcen.append(g.evcent())
            btcen.append(g.betweenness())
            clcen.append(g.closeness())
            degree.append(g.degree())			
    return evcen,btcen,clcen,degree
def are_corrected():
    ''' check if the pdbs have been corrected by modeller'''
    pdbs=glob('*.pdb')
    corrected = any(['-c' in p for p in pdbs])
    return corrected

def Load_Contacts(prefix,options):
    ''' load the contacts file into the g-filter list '''
    corrected = are_corrected()
    #print corrected # for debug purposes
    if options.contacts:
        gfilter = []
        # Check if there is aprorpiate .contact file
        if not os.path.isfile('%s.contacts'%(prefix)):
            #determine if pdbs have been corrected
            print 'Need to run contactmapper.py to provide the .contacts file.'
            print 'Running contactmapper.py: Make sure the PDB file(s) is(are) in the current '\
                  'working directory and that contact mapper is in your path'
            if not options.multiple:
                if not os.path.isfile('%s.pdb'%(prefix)):
                    print 'PDB file not in directory'
                    raise
            data = Contacts(prefix,options.multiple,corrected)
            data.WriteToFile(prefix)
        with open(prefix+'.contacts') as fin:
            for line in fin:
                if line:
                    line = line.strip().strip('(').strip(')').split(',')
                    gfilter.append( (int(line[0]), int(line[1])) )
    else:
        gfilter = None	
    return gfilter

def Master_of_ceremony(options):
    # Say hi!!!#######################################################################
    print('\nWelcome to MODULER')
    print('A python script to explore modularity on coordinates data\n\n')
    l = ''
    # Print the chosen parameters and write it to a file #########################
    l+='Chosen parameters:\nPrefix=%s\nMultiple PDBs = %s\n'%(args[0],options.multiple)
    l+='Test of significance for correlation = %s\n'%(options.method)
    l+='Power analysis = %s\n'%(options.power)
    l+='RV test = %s\nAglomerate dimensions as Euclidean distance = %s\n'%(options.rv,options.additive)
    l+='Confidence Value = %f (to be use in both correlation and significance test if true)\n'%(options.confval)
    l+='Absolute value of correlation = %s\n'%(options.absolutecov)
    l+='Dimensions = %d\n'%(options.dim)
    if options.perm > 0:
        l+='Single test for significance of clusters = True with %d permutations\n'%(options.perm)
    else:
        l+='Single test for significance of clusters = False\n'
    if options.itperm > 0:
        l+='Iterative test for significance of clusters = True with %d permutations\n'%(options.itperm)
    else:
        l+='Iterative test for significance of clusters = False\n'
    l+= 'Filtering out non contact interactions = %s\n'%(options.contacts)
    l+= 'Linear discriminants prefiltering = %s\n'%(options.lda)
    l+= 'Bootstrap : %d (if 0 no bootstrap is performed)'%(options.boot)
    #l+= 'Analyse  interlandmark distances (if False will use raw coordinates) = %s'%(options.dist)
    print l
    ou = open(prefix+'.parameters', 'w')
    now = datetime.datetime.now()
    ou.write('file created on %s\n\n'%(str(now))+l)
    ou.close()

def landmarks4morph(prefix,dim):
    ''' create a dummy landmarks file if morphological data is used'''
    names=[]
    nl=open(prefix+'.landmarks','w')
    with open(prefix+'.gm') as l:
        for line in l:
            bline=line.split(';')
            if bline[-1] == '':
                bline=bline[:-1]
            names.append(bline[0].strip().strip('>'))
            lndm=len(bline[1:])/dim
    for e in names:
        nl.write('>%s\n'%(e))
        for la in range(lndm):
            nl.write(str(la)+'\t'+str(la+1)+'\tmorph\n')
    nl.close()

def main(prefix, options):
    ''' execute the code '''
    # Load Data ############################################################
    data, sample_size = Read_GM(prefix, options.dim)#, options.dist)

    # Load Contacts ########################################################
    gfilter = Load_Contacts(prefix,options)

    # Create a correlation matrix testing for significance #################
    if options.method:
        if options.absolutecov:
            matrices = Sigcorr(prefix,data, sample_size, options.method, options.confval, 
                               options.threshold, options.contacts, options.dim, 
                               options.absolutecov)
        else:
            matrices = Sigcorr(prefix,data, sample_size, options.method, options.confval,
                               options.threshold, options.contacts, options.dim)

    # Else, use covariance matrix instead or correlation without statistical
    # test #################################################################
    elif options.usecov:
        matrices, thr = UseCov(matrices, data, options.threshold)
    else:
        matrices, thr = UseCorr(matrices, data, options.threshold)

    if options.mat == 'cor':
        for m in range(len(matrices)):
            write_singleMatrix(prefix,matrices[m])
            os.system('mv %s.mat %s.%d.mat'%(prefix,prefix,m))

    # Agglomerate landmark dimensions ######################################
    if options.additive:
        lms = agglomerare_additive(matrices, options.absolutecov)
        if options.mat == 'agg':
            write_singleMatrix(prefix,lms)
    else:
        lms = agglomerare_mean(matrices, options.dim, options.absolutecov)

    # get Edge assignment threshold and print it to screen #################
    if options.threshold == 'Auto':
        thr, options.threshold = GetAglomeratedThreshold(lms)
    else:
        f=open(prefix+'.parameters', 'a')
        t = 'Threshold mode: Custom = %f'%(options.threshold)
        f.write(t)
        f.close()		
        print t
    ''' if fudge set threshold #########################################################
	if fudge:#deprecated!! is not an option in the parser... migth be usefull down the road
		print 'fudge =%s'%(fudge)
		threshold = Z_fisher(options.confval, sample_size) + F_transf(threshold)
		f= 'Edge assignment threshold = %s'%(str(threshold))
		print f
		#ou.write('\n'+f)'''

    # Build igraph ###################################################################
    g = Build_igraph(lms, options, gfilter)
    if options.contacts:
        A = g.get_adjacency() # get the adjacency matrix to filter permt
    else:
        A=[]
    g.write_dot('graph.dot')
    # Clustering by components comstGreed ############################################
    memgreed, mem = Graph_Components(g,lms)

    # write to cluster file ##########################################################
    write_cluster(prefix, g, memgreed)
    centralities(prefix, memgreed, options.overall,g)

    # create csv files per cluster
    #m, mem = clus2csv(prefix,data, sample_size)
    lms2Rlms(prefix, lms)
    m = open(prefix+'.graphcluster').read()
    c = deepcopy(m)
    nc=''
    c = c.strip().split()
    d = Counter(c)
    for e in c:
        if d[e] == 1:
            nc+= '? '
        else:
            nc+= '%s '%(str(e))	
    #copyfile(prefix+'.graphcluster',prefix+'.community')
    # run LDA to premerge clusters
    if options.lda:
        if not options.contacts:
            print 'WARNING!! Using Linear discriminant pre-filtering without contacts '\
                  'will underestimate the number of modules!! Use the non lda filtered option instead.'
        fo = open(prefix+'.community','w')
        fo.write(nc.strip())
        fo.close()
        print 'Membership vector before LDA (singletons are labelled as ?):'
        print nc
        m = LDAMerge(prefix,nc,options)
        print 'Membership vector after LDA:'
        print m
    # significance test of clustering: Outputs a new graphcluster
    if options.perm != 0:
        scl, newm, newvec, FDR = lms_permt(m,options.perm,A,options)
        '''
		scl, newm, newvec, FDR = permt(prefix, data, matrices, sample_size, m, options.perm, options)
		# Write the DOT graph
		Dot_graph(prefix, scl, newm)'''

    elif options.itperm != 0:
        scl, newm, newvec, FDR = iterative_permt(prefix, data, matrices, m, sample_size, A, options)
        # Write the DOT graph
        Dot_graph(prefix, scl, newm)

    else:
        newm = m # bootrap work around
    if options.graph:
        plot_graph(g,newm,prefix)
    # Perform options.power analyses
    if options.power:
        propor = Clus_Power(prefix, data, options.confval, options.power, sample_size,
                            options.dim,newm)
    for f in glob('*.csv'):
        if not 'summary' in f:
            os.remove(f)

    if options.rv:
        RVs = multiple_RVs(prefix,options.dim)
        write_RVfile(prefix,RVs)

    return newm,data

def plot_graph(g,m,prefix):
    ''' plot one labled and one unlabled graph '''
    from igraph import plot#,save
    col = ['red','blue', "cyan", "magenta", "yellow", 'Black','Lavender Blush', 'Medium Orchid', 
              'Midnight Blue', 'White Smoke', 'Peru', 'Seashell', 'Dark Goldenrod', 'Coral', 
              'Firebrick', 'Cornsilk 4', 'Honeydew 3', 'Dark Green', 'Violet', 'Maroon', 
              'Forest Green', 'Khaki', 'Chocolate', 'Slate Blue', 'Dark Olive Green', 'Navy', 
              'Floral White', 'Light Cyan', 'Peach Puff', 'Medium Spring Green', 'Snow 2', 
              'Gainsboro', 'Light Coral', 'Light Goldenrod Yellow', 'Bisque', 'Blanched Almond',
              'Thistle', 'Lime Green', 'Cornsilk 2', 'Honeydew', 'Green Yellow', 'Deep Sky Blue', 
              'Pale Green', 'Turquoise', 'Cornflower Blue', 'Indian Red', 'Spring Green', 
              'Dark Orange', 'Light Goldenrod', 'Peach Puff 4', 'Ivory 2', 'Steel Blue', 
              'Misty Rose', 'Pale Goldenrod', 'Bisque 3', 'Light Sea Green', 'Olive Drab', 
              'Mint Cream', 'Pale Violet Red', 'Linen', 'Sea Green', 'Dodger Blue', 'Seashell 2',
              'Plum', 'Moccasin', 'Dark Khaki', 'Blue Violet', 'Sienna', 'Dark Violet', 'Bisque 2',
              'Peach Puff 2', 'Dark Slate Gray', 'Hot Pink', 'Light Pink', 'Wheat', 'Lemon Chiffon',
              'Burlywood', 'Antique White', 'Antique White 3', 'Medium Slate Blue', 'Dark Sea Green',
              'Snow', 'Light Sky Blue', 'Old Lace', 'Pale Turquoise', 'Purple', 'Chartreuse', 
              'Medium Violet Red', 'Snow 3', 'Peach Puff 3', 'Alice Blue', 'Ghost White', 'Pink', 
              'Violet Red', 'Light Blue', 'Bisque 4', 'Dark Salmon', 'Medium Turquoise', 'Beige', 
              'Ivory 4', 'Medium Sea Green', 'Yellow Green', 'Honeydew 4', 'Deep Pink', 'Azure', 
              'Seashell 4', 'Tomato', 'Light Salmon', 'Seashell 3', 'Honeydew 2', 'Cornsilk 3', 
              'Sandy Brown', 'Orange Red', 'Papaya Whip', 'Ivory 3', 'Dark Slate Blue', 'Lawn Green',
              'Antique White 2', 'Goldenrod', 'Tan', 'Medium Purple', 'Light Yellow', 'Ivory', 
              'Light Gray', 'Medium Aquamarine', 'Saddle Brown', 'Snow 4', 'Navajo White', 'Gray', 
              'Cornsilk', 'Salmon', 'Dark Turquoise', 'Antique White 4', 'Cadet Blue', 'Rosy Brown', 
              'Gold', 'Orange', 'Light Slate Gray', 'Powder Blue', 'Lavender', 'Dark Orchid', 
              'Royal Blue', 'Medium Blue', 'Sky Blue', 'Light Slate Blue', 'Light Steel Blue', 
              'Slate Gray', 'Cyan', 'Dim Gray', 'Orchid', 'Aquamarine']
    s = list(set(m.replace(' ','')))
    colors={}
    for n in range(len(s)):
        colors[s[n]]=col[n]
    #this is just for the LDA paper, comment out when finish
    shape=['square', 'square', 'square', 'square', 'square', 'square', 'circle', 'circle', 'circle', 
           'circle', 'square', 'circle', 'triangle-up', 'circle', 'circle', 'circle', 'square', 'circle', 
           'square', 'square', 'circle', 'circle', 'square', 'square', 'circle', 'circle', 'circle', 
           'square', 'square', 'circle', 'square', 'square', 'square', 'square', 'square', 'circle', 
           'square', 'square', 'circle', 'circle', 'circle', 'circle', 'circle', 'circle', 'circle', 
           'square', 'square', 'square', 'square', 'circle', 'square', 'circle', 'circle', 'square', 
           'circle', 'square', 'circle', 'square', 'square', 'circle', 'circle', 'circle', 'square', 
           'circle', 'circle', 'square', 'square', 'circle', 'circle', 'circle', 'square', 'square', 
           'square', 'square', 'circle', 'circle', 'square', 'circle', 'circle', 'square', 'square', 
           'square', 'circle', 'square', 'square', 'square', 'square', 'square', 'square', 'triangle-up',
           'circle', 'circle', 'square', 'circle', 'circle', 'circle', 'circle', 'circle', 'circle', 
           'square', 'square', 'square']
    M = m.split(' ')
    for e in range(len(g.vs)):
        g.vs[e]['color']=colors[M[e]]
        g.vs[e]['shape']=shape[e]#uncomment this too!!
    
    #plot with lables
    plot(g, '%s_LabledGraph.png'%(prefix),layout="fr")
   # fig1.save('%s_LabledGraph.png'%(prefix))
    # plot with no lables
    plot(g, '%s_unlabledGraph.png'%(prefix),layout="fr", vertex_label=None)
    #fig2.save('%s_unlabledGraph.png'%(prefix))
################################################################################################
## Definitions from the boostrapmoduler.py file, used to create replicates of the input GM file
################################################################################################
def BootstrapReplicate(data):
    ''' returns a sample as a list of lines'''
    temp = []
    for i in range(len(data)):
        temp.append(randint(0,len(data)-1))

    for i in range(len(temp)):
        temp[i] = data[temp[i]]

    return temp

def WriteListToFile(data, filename):
    ''' write a data list into a file '''
    fout = open(filename, 'w')
    for da in data:
        if da[-1] != '\n':
            fout.write(da + '\n')
        else:
            fout.write(da)
    fout.close()

######################## vectorfscore Functions to find distance #################################
'''
This bit was mostly written by Khan Nguyen
'''
def Duplicate(pairs):
    ''' Check if there are any tuple duplicates '''
    temp_top = []
    temp_bottom = []
    for i in range(len(pairs)):
        if pairs[i][0] not in temp_top:
            temp_top.append(pairs[i][0])
        if pairs[i][1] not in temp_bottom:
            temp_bottom.append(pairs[i][1])
    return (len(pairs) != len(temp_top) and len(pairs) != len(temp_bottom))


def VectorDistance(a, b):
    '''
    Return the non-weighted distance between 2 vectors. a is the reference vector, while b is 
    the test vector
    '''
    a = a.split()
    b = b.split()
    if len(a) != len(b):
        print ''.join(a) + ' and ' + ''.join(b)
        print "These vectors are not of the same length."
        sys.exit(-1)
    else:
        pairs = [] ## this list stores the pairs of different label
        differentChars = [] ## stores the number of different characters used in test string
        dist = 0
        for i in range(len(a)):
            if (a[i], b[i]) not in pairs:
                pairs.append((a[i], b[i]))
            if b[i] not in differentChars:
                differentChars.append(b[i])
        num = 0 ## keep track of the letters used
        # Calculate the distance
        while Duplicate(pairs): ## check if there are duplicates
            found = False
            for i in range(len(pairs)):
                for j in range(i+1, len(pairs)):
                    # Check for merging possibility
                    if pairs[i][0] == pairs[j][0]:
                        num += 1
                        newChar = chr(64+len(differentChars)+num)
                        while pairs[i][1] in b or pairs[j][1] in b:
                            if pairs[i][1] in b:
                                b[b.index(pairs[i][1])] = newChar
                            if pairs[j][1] in b:
                                b[b.index(pairs[j][1])] = newChar
                        found = True
                        break
                    # Check for splitting posibility
                    elif pairs[i][1] == pairs[j][1]:
                        num += 2
                        newChar1 = chr(64+len(differentChars)+num-1)
                        newChar2 = chr(64+len(differentChars)+num)
                        for k in range(len(a)):
                            if a[k] == pairs[i][0] and b[k] == pairs[i][1]:
                                b[k] = newChar1
                            elif b[k] == pairs[i][1]:
                                b[k] = newChar2
                        found = True
                        break
                if found:
                    dist += 1
                    break
            pairs = []
            for i in range(len(a)):
                if (a[i], b[i]) not in pairs:
                    pairs.append((a[i], b[i]))

        return dist



def WeightedDistance(a, b):
    '''
    Return the weighted distance between 2 vectors. a is the reference vector, while b is the 
    test vector
    '''
    a = a.split()
    b = b.split()
    if len(a) != len(b):
        print ''.join(a) + ' and ' + ''.join(b)
        print "These vectors are not of the same length."
        sys.exit(-1)
    else:
        pairs = [] ## this list stores the pairs of different label
        differentChars = [] ## stores the number of different characters used in test string
        dist = 0
        for i in range(len(a)):
            if (a[i], b[i]) not in pairs:
                pairs.append((a[i], b[i]))
            if b[i] not in differentChars:
                differentChars.append(b[i])
        num = 0 ## keep track of the letters used

        # Calculate the distance
        while Duplicate(pairs): ## check if there are duplicates

            costs = {} ## dictionary storing the costs of all possible moves in one step
            # Find the least expensive path
            for i in range(len(pairs)):
                for j in range(i+1, len(pairs)):
                    # Check for merging possibility
                    if pairs[i][0] == pairs[j][0]:
                        cost = min([b.count(pairs[i][1]), b.count(pairs[j][1])])
                        costs[cost] = [i, j, 'merge'] ## store the indices and the type of operation along with the cost
                    # Check for splitting possibility
                    elif pairs[i][1] == pairs[j][1]:
                        cost = [0, 0]
                        for k in range(len(a)):
                            if a[k] == pairs[i][0] and b[k] == pairs[i][1]:
                                cost[0] += 1
                            elif b[k] == pairs[i][1]:
                                cost[1] += 1
                        costs[min(cost)] = [i, j, 'split'] ## store the indices and the type of operation along with the cost

            # Perform operation
            minCost = min(costs.keys())
            i = costs[minCost][0]
            j = costs[minCost][1]
            if costs[minCost][2] == 'merge':
                num += 1
                newChar = chr(64+len(differentChars)+num)
                while pairs[i][1] in b or pairs[j][1] in b:
                    if pairs[i][1] in b:
                        b[b.index(pairs[i][1])] = newChar
                    if pairs[j][1] in b:
                        b[b.index(pairs[j][1])] = newChar
            else: ## else it's a splitting operation
                num += 2
                newChar1 = chr(64+len(differentChars)+num-1)
                newChar2 = chr(64+len(differentChars)+num)
                for k in range(len(a)):
                    if a[k] == pairs[i][0] and b[k] == pairs[i][1]:
                        b[k] = newChar1
                    elif b[k] == pairs[i][1]:
                        b[k] = newChar2

            dist += minCost
            pairs = []
            for i in range(len(a)):
                if (a[i], b[i]) not in pairs:
                    pairs.append((a[i], b[i]))

        return dist

def BipartitionAgree(a, b):
    '''
    Return whether 2 strings of bipartitions agree or conflict. The strings must consist of
    1 and 0 only
    '''
    if len(a) != len(b):
        print a + ' and ' + b + ' are not of the same length to compare bipartitions.'
        sys.exit(-1)
    else:
        pairs = []
        for i in range(len(a)):
            if (a[i], b[i]) not in pairs:
                pairs.append((a[i], b[i]))
        return len(pairs) < 4


#####################################################################################
def matrix2file(prefix, matrix, typ):
    ''' write matrix to a file '''
    fout = open(prefix+'%s.score'%(typ),'w')
    for row in matrix:
        for item in row:
            fout.write(str(item)+';')
        fout.write('\n')

def scores_matrices(prefix,vectors,edges):
    '''
    Compute the scoring matrices, write them to files and return the bipartitions dictionary
    '''
    # Generate a F-score matrix
    Fmatrix=[]
    for i in range(len(edges)):
        row = []	
        for j in range(len(edges)):
            row.append(Fscore(Specificity(edges[i], edges[j]), Sensitivity(edges[i], edges[j])))
        Fmatrix.append(row)
    matrix2file(prefix,Fmatrix,'Fscore')

    i = None
    j = None
    UWmatrix = []
    # Generate a non-weighted distance matrix
    for i in range(len(vectors)):
        row = []
        for j in range(len(vectors)):
            row.append(VectorDistance(vectors[i], vectors[j]))
        UWmatrix.append(row)
    matrix2file(prefix,UWmatrix,'uwDist')

    i = None
    j = None	
    Wmatrix = []	
    # Generate a weighted distance matrix
    for i in range(len(vectors)):
        row = []
        for j in range(len(vectors)):
            row.append(WeightedDistance(vectors[i], vectors[j]))
        Wmatrix.append(row)
    matrix2file(prefix,Wmatrix,'wDist')

def bipartition_agreement(prefix, vectors, m):
    ''' Calculate the local bipartition agreement scores '''
    refvec = m.split()
    clusters = sorted(list(set(refvec)))
    refbipartitions = {} ## dictionary containing each cluster's bipartition
    for clus in clusters:
        refbipartitions[clus] = ''
        for i in range(len(refvec)):
            if refvec[i] == clus:
                refbipartitions[clus] += '1' ## 1 for dot
            else:
                refbipartitions[clus] += '0' ## 0 for dash

    # Enumerate all the bipartitions in the replicate set
    repbipartitions = {}
    for i in range(len(vectors)):
        repvec = vectors[i].split()
        cluses = sorted(list(set(repvec)))
        for clus in cluses:
            partition = ''
            for i in range(len(repvec)):
                if repvec[i] == clus:
                    partition += '1'
                else:
                    partition += '0'
            if partition not in repbipartitions:
                repbipartitions[partition] = 1
            else:
                repbipartitions[partition] += 1

    # Calculate the score
    FOUT = open(prefix+'.scores.bootstrap','w')
    bipartitions = {}
    total = sum(repbipartitions.values()) ## total number of bipartitions in the replicate set
    for clus in refbipartitions:
        agreement = 0
        for p in repbipartitions:
            if BipartitionAgree(refbipartitions[clus], p):
                agreement += repbipartitions[p]
        bipartitions[clus] = round(float(agreement)/total, 3)

    for k,v in bipartitions.iteritems():
        FOUT.write(k+'\t%f\n'%(v))

    FOUT.close()
    return bipartitions

def if_bootfile(prefix, options):
    '''
    if another bootrap intance has been called and crashed, this function will finished
    the remaining and / or compute the agreement
    '''
    vectors = []
    data = open(prefix+'.gm').read().split('\n')
    if data[-1] == '' or data[-1] =='\n':
        data = data[:-1]	
    with open(prefix+'.bootstrap') as f:
        for e in f:
            if e == '' or e == '\n':
                continue
            elif '[[' in e:
                m = e[e.find("#\t('")+4:e.find(" ',")]
                vectors.append(m)
            else:
                m = e[e.find("#\t('")+4:e.find("\n")]
                vectors.append(m)                
    dif = int(options.boot) - len(vectors)
    if dif > 0:
        # backup the original files
        F = glob(prefix+'*')
        for f in F:
            # get the extension and filename
            ext = f[f.rfind('.'):]
            name = f[:f.find(ext)]
            nf=name+'_original'+ext
            # rename the file
            os.rename(f,nf)
        # reopen boot file and keep going
        F = open(prefix+'.bootstrap','a')
        for i in range(dif):
            i = i + len(vectors)
            optionsrep = deepcopy(options)
            optionsrep.boot = 0
            optionsrep.power = False
            optionsrep.rv = False			
            l='\n#  BOOTSTRAP REPLICATE %d  #\n'%(i)
            print '#'*(len(l)-2),l,'#'*(len(l)-2)
            fn = prefix+'.%s.gm'%(i)
            D = BootstrapReplicate(data)
            WriteListToFile(D,fn)
            copyfile(prefix+'_original.landmarks', prefix+'.%d.landmarks'%(i))
            copyfile(prefix+'_original.contacts', prefix+'.%d.contacts'%(i))
            m = main(fn[:-3],optionsrep)
            F.write('Boot %d#\t%s\n'%(i,m))
            vectors.append(m[0])
            files = glob(prefix+'.%s.*'%(i))
            for fi in files:
                os.remove(fi)			
    return vectors

def bootstrap(prefix, m, options):
    '''	execute the bootstrap inferfence'''
    if os.path.isfile(prefix+'.bootstrap'):
        print 'Previous bootstrap run found. Using %s file'%(prefix+'.bootstrap')
        vectors = if_bootfile(prefix,options)
        #print vectors#debugging
    else:
        #simple data gathering
        data = open(prefix+'.gm').read().split('\n')
        if data[-1] == '' or data[-1] =='\n':
            data = data[:-1]
        # backup the original files
        F = glob(prefix+'*')
        for f in F:
            # get the extension and filename
            ext = f[f.rfind('.'):]
            name = f[:f.find(ext)]
            nf=name+'_original'+ext
            # rename the file
            os.rename(f,nf)
        # Create the out file
        bout = open(prefix+'.bootstrap','w')
        vectors=[]		
        for i in range(options.boot):
            optionsrep = deepcopy(options)
            optionsrep.boot = 0
            optionsrep.power = False
            optionsrep.rv = False
            l='\n#  BOOTSTRAP REPLICATE %d  #\n'%(i)
            print '#'*(len(l)-2),l,'#'*(len(l)-2)
            fn = prefix+'.%s.gm'%(i)
            D = BootstrapReplicate(data)
            WriteListToFile(D,fn)
            copyfile(prefix+'_original.landmarks', prefix+'.%d.landmarks'%(i))
            copyfile(prefix+'_original.contacts', prefix+'.%d.contacts'%(i))
            M = main(fn[:-3],optionsrep)
            bout.write('Boot %d#\t%s\n'%(i,M))
            vectors.append(M[0])
            files = glob(prefix+'.%s.*'%(i))
            files.extend(glob('*_clusgraph'))
            for fi in files:
                os.remove(fi)
    edges = [VectorToEdgeList(v) for v in vectors]
    #scores_matrices(prefix,vectors,edges)#very slow... comment out if really need it!
    bipartition = bipartition_agreement(prefix, vectors, m)


def LDAMerge(prefix,m,options):
    '''
    Will execute lda4mod, return the merges created by a 95% confidence ellipse on a
    LDA.
    '''
    if os.path.isfile('merges.txt'):
        os.remove('merges.txt')
    tomerge=[]
    scriptpath = '/'.join(os.path.abspath(__file__).split('/')[:-1])
    R = Popen('nohup Rscript %s/%s %s %d'%(scriptpath,'ldaellipse4mod.R', prefix+'.gm',options.dim),
              shell=True)#,stderr=PIPE,stdout=PIPE)

    #R = Popen('R CMD BATCH "--args %s %d" %s/%s'%(prefix+'.gm',options.dim,scriptpath,'ldaellipe4mod.R'),
    #          shell=True)#,stderr=PIPE,stdout=PIPE)
    R.wait()
    #call(['Rscript','%s/%s'%(scriptpath,'ldaellipe4mod.R'),'%s'%(prefix+'.gm'),'%d'%(options.dim)])
    if os.path.isfile('merges.txt'):
        with open('merges.txt') as F:
            for line in F:
                l = line.strip().split(',')
                if 'NA' in l:
                    continue
                else:
                    t = tuple(l)
                    if not t in tomerge:
                        tomerge.append(t)
        newcl = equiclus(tomerge)
        m = m.strip().split()
        memb = rename_clusters(newcl, m)
        #os.remove('merges.txt')
    else:
        memb=m
        print 'LDA cannot be performed. Perhaps variables without any variance or only singletons'\
              'in the dataset. Check your results carefully.'
    return memb

# End of definitions###########################################################################

# Aplication of the code ######################################################################
if __name__ == "__main__":
    # Command line input #############################################################
    opts = optparse.OptionParser(usage='%prog <prefix> [options]')
    opts.add_option('-o','--covariance',dest='usecov', action="store_true",
                    default=False, help='Use covariance instead of correlation to build'\
                    ' the graph. If this option is not provided the program will use cor'\
                    'relation by default.')
    opts.add_option('-a','--mean', dest='additive', action="store_false", default=True,
                    help='Use the mean of the correlation in each dimension instead of'\
                    'the euclidean distance to aglomerate the dimensions. The default '\
                    'behaviour is additconfvalive.')
    opts.add_option('-t', '--threshold', dest='threshold', action='store',# type=float, 
                    default=0.0, help= 'Set the threshold for edge assingment to the '\
                    'value provided. Otherwise is set to 0. An extra option is Auto '\
                    'which will calculate the treshold based on the mean standard '\
                    'deviation of the agglomerated dimensions matrix.')
    opts.add_option('-p','--power', dest='power', action='store', type='float', default=0.8,
                    help='Perform a power analysis with this value as the desired Power'\
                    ' (1-type II error). Default: 0.8')
    opts.add_option('-m','--method', dest='method', default='fisher', help='Test the '\
                    'significance of the correlation by the method provided(pearson, '\
                    'spearman, kendall or fisher). The latter is the fisher transform'\
                    'ation of the pearson correlation. If no test needed False should'\
                    ' be providedls.')
    opts.add_option('-i', '--confidencelevel', dest='confval', action='store', type='float',
                    default=0.95, help='Define the confidence level (1-alpha) for the'\
                    ' correlation test')
    opts.add_option('-b','--histogram', dest='histocorr', action="store_true", default=False,
                    help='Draw an histogram of the correlation matrix.')
    opts.add_option('-d','--dimensions', dest='dim',action='store', type='int', default=3,
                    help='Set the dimensions of the shape.')
    opts.add_option('-c','--contact',dest='contacts', action="store_true", default=False,
                    help='Use the all-atom contact matrix when assigning edges')
    opts.add_option('-u','--absolute',dest='absolutecov',action="store_true",default=False,
                    help= 'Use absolute values of the correlation matrix')
    opts.add_option('-e', '--permtest',dest='perm', action='store', type='int',default=0,
                    help='Test the significance of clustering with N permutations,'\
                    ' the default behaviour (0) means no permutation test')
    opts.add_option('-n','--itpermtest',dest='itperm', action='store', type='int',default=0,
                    help= 'Iterative (until convergence) test the significance of'\
                    ' clustering with N permutations. As with the previous option'\
                    ' the default is no test.')
    opts.add_option('-w','--matrix',dest='mat', default=False, help='Write to file'\
                    ' the matrix of the option chosen. It can be either cor (which '\
                    'will print the correlation matrix for each dimension), or agg '\
                    'will write to a file the agglomerated dimensions matrix.')
    opts.add_option('-r','--RV',dest='rv',action="store_false",default=True,
                    help='Do not Write the the Excouffier RV for the modules inferred.'\
                    ' The p-values will be given using the Pearson type III approximation.')
    opts.add_option('-g','--multiplepdbs',dest='multiple',action="store_true",default=False,
                    help= 'Use this flag if you have more than 50 structures and you align '\
                    'them using pairwise2multiple method')
    opts.add_option('-f','--fullcentrality',dest='overall',action="store_true",default=False,
                    help='Use this if you want to calculate the centralities on the '\
                    'full graph, as opposed as by modules ( Default behaviour ).')
    opts.add_option('-M','--morphological',dest='morph',action="store_true",default=False,
                    help='Use this if you want to estimate the modularity to morphological'\
                    ' data.')
    opts.add_option('-B','--bootstrap',dest='boot',action='store', type='int',default=0,
                    help='Compute bootstrap for sample reliability estimation.')	
    opts.add_option('-l','--lda',dest='lda',action="store_true",default=False,
                    help='Use LDA to premerge the community detection clusters')
    opts.add_option('-G','--graph',dest='graph',action="store_true",default=False,
                    help='Plot the graph colored by modules')
    options, args = opts.parse_args()

    # some workaround with the CL ####################################################
    prefix = args[0]
    try:
        options.threshold = int(options.threshold)
    except:
        options.threshold = 'Auto'
    if options.method == 'False':
        options.method = False
    fudge=False#deprecated!! is not an option in the parser... migth be usefull down the road
    # Introduce the program and write parameters to a file and screen ################
    Master_of_ceremony(options)

    # create a bogus landmarks file for morphological data
    if options.morph: landmarks4morph(prefix,options.dim)

    m,data = main(prefix, options)

    if options.boot > 0:
        bootstrap(prefix, m, options)
# End of the code ####################################################################
