#!/usr/bin/python
"""
Format GM for classification

This script will take an aligned GM file, and splitted into training and testset depending on the folders provided. It will also create
the sequence and sequence properties in conjuction with the structure and alone.
"""

blosum62 = {'*': {'*': 1, 'A':-4, 'C':-4, 'B':-4, 'E':-4, 'D':-4, 'G':-4, 'F':-4, 'I':-4, 'H':-4, 'K':-4, 'M':-4, 'L':-4, 'N':-4, 'Q': -4, 'P': -4, 'S': -4, 'R': -4, 'T': -4, 'W': -4, 'V': -4, 'Y': -4, 'X': -4, 'Z': -4},
            'A': {'*':-4, 'A': 4, 'C': 0, 'B':-2, 'E':-1, 'D':-2, 'G': 0, 'F':-2, 'I':-1, 'H':-2, 'K':-1, 'M':-1, 'L':-1, 'N':-2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2, 'X': 0, 'Z': -1}, 
            'C': {'*':-4, 'A': 0, 'C': 9, 'B':-3, 'E':-4, 'D':-3, 'G':-3, 'F':-2, 'I':-1, 'H':-3, 'K':-3, 'M':-1, 'L':-1, 'N':-3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2, 'X': -2, 'Z': -3},
            'B': {'*':-4, 'A':-2, 'C':-3, 'B': 4, 'E': 1, 'D': 4, 'G':-1, 'F':-3, 'I':-3, 'H': 0, 'K': 0, 'M':-3, 'L':-4, 'N': 3, 'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1},
            'E': {'*':-4, 'A':-1, 'C':-4, 'B': 1, 'E': 5, 'D': 2, 'G':-2, 'F':-3, 'I':-3, 'H': 0, 'K': 1, 'M':-2, 'L':-3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 4}, 
            'D': {'*':-4, 'A':-2, 'C':-3, 'B': 4, 'E': 2, 'D': 6, 'G':-1, 'F':-3, 'I':-3, 'H':-1, 'K':-1, 'M':-3, 'L':-4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1}, 
            'G': {'*':-4, 'A': 0, 'C':-3, 'B':-1, 'E':-2, 'D':-1, 'G': 6, 'F':-3, 'I':-4, 'H':-2, 'K':-2, 'M':-3, 'L':-4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3, 'X': -1, 'Z': -2},
            'F': {'*':-4, 'A':-2, 'C':-2, 'B':-3, 'E':-3, 'D':-3, 'G':-3, 'F': 6, 'I': 0, 'H':-1, 'K':-3, 'M': 0, 'L': 0, 'N':-3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3, 'X': -1, 'Z': -3},
            'I': {'*':-4, 'A':-1, 'C':-1, 'B':-3, 'E':-3, 'D':-3, 'G':-4, 'F': 0, 'I': 4, 'H':-3, 'K':-3, 'M': 1, 'L': 2, 'N':-3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1, 'X': -1, 'Z': -3}, 
            'H': {'*':-4, 'A':-2, 'C':-3, 'B': 0, 'E': 0, 'D':-1, 'G':-2, 'F':-1, 'I':-3, 'H': 8, 'K':-1, 'M':-2, 'L':-3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2, 'X': -1, 'Z': 0}, 
            'K': {'*':-4, 'A':-1, 'C':-3, 'B': 0, 'E': 1, 'D':-1, 'G':-2, 'F':-3, 'I':-3, 'H':-1, 'K': 5, 'M':-1, 'L':-2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 1}, 
            'M': {'*':-4, 'A':-1, 'C':-1, 'B':-3, 'E':-2, 'D':-3, 'G':-3, 'F': 0, 'I': 1, 'H':-2, 'K':-1, 'M': 5, 'L': 2, 'N':-2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1, 'X': -1, 'Z': -1},
            'L': {'*':-4, 'A':-1, 'C':-1, 'B':-4, 'E':-3, 'D':-4, 'G':-4, 'F': 0, 'I': 2, 'H':-3, 'K':-2, 'M': 2, 'L': 4, 'N':-3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3}, 
            'N': {'*':-4, 'A':-2, 'C':-3, 'B': 3, 'E': 0, 'D': 1, 'G': 0, 'F':-3, 'I':-3, 'H': 1, 'K': 0, 'M':-2, 'L':-3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0}, 
            'Q': {'*':-4, 'A':-1, 'C':-3, 'B': 0, 'E': 2, 'D': 0, 'G':-2, 'F':-3, 'I':-3, 'H': 0, 'K': 1, 'M': 0, 'L':-2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1, 'X': -1, 'Z': 3}, 
            'P': {'*':-4, 'A':-1, 'C':-3, 'B':-2, 'E':-1, 'D':-1, 'G':-2, 'F':-4, 'I':-3, 'H':-2, 'K':-1, 'M':-2, 'L':-3, 'N':-2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3, 'X': -2, 'Z': -1}, 
            'S': {'*':-4, 'A': 1, 'C':-1, 'B': 0, 'E': 0, 'D': 0, 'G': 0, 'F':-2, 'I':-2, 'H':-1, 'K': 0, 'M':-1, 'L':-2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2, 'X': 0, 'Z': 0}, 
            'R': {'*':-4, 'A':-1, 'C':-3, 'B':-1, 'E': 0, 'D':-2, 'G':-2, 'F':-3, 'I':-3, 'H': 0, 'K': 2, 'M':-1, 'L':-2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0}, 
            'T': {'*':-4, 'A': 0, 'C':-1, 'B':-1, 'E':-1, 'D':-1, 'G':-2, 'F':-2, 'I':-1, 'H':-2, 'K':-1, 'M':-1, 'L':-1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2, 'X': 0, 'Z': -1}, 
            'W': {'*':-4, 'A':-3, 'C':-2, 'B':-4, 'E':-3, 'D':-4, 'G':-2, 'F': 1, 'I':-3, 'H':-2, 'K':-3, 'M':-1, 'L':-2, 'N':-4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2, 'X': -2, 'Z': -3}, 
            'V': {'*':-4, 'A': 0, 'C':-1, 'B':-3, 'E':-2, 'D':-3, 'G':-3, 'F':-1, 'I': 3, 'H':-3, 'K':-2, 'M': 1, 'L': 1, 'N':-3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -2}, 
            'Y': {'*':-4, 'A':-2, 'C':-2, 'B':-3, 'E':-2, 'D':-3, 'G':-3, 'F': 3, 'I':-1, 'H': 2, 'K':-2, 'M':-1, 'L':-1, 'N':-2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7, 'X': -1, 'Z': -2}, 
            'X': {'*':-4, 'A': 0, 'C':-2, 'B':-1, 'E':-1, 'D':-1, 'G':-1, 'F':-1, 'I':-1, 'H':-1, 'K':-1, 'M':-1, 'L':-1, 'N':-1, 'Q': -1, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -2, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1}, 
            'Z': {'*':-4, 'A':-1, 'C':-3, 'B': 1, 'E': 4, 'D': 1, 'G':-2, 'F':-3, 'I':-3, 'H': 0, 'K': 1, 'M':-1, 'L':-3, 'N': 0, 'Q': 3, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 4}}

#constant dictionary with each aminoacid and its 1) hydrophobicity index at ph 7 (as proposed in
#Monera et.al. J. Prot. Sci. 1: 319:329 (1995)); 2) Van der Waals volume (as calculated by Richards
#J.Mol.Biol. 82:1-14,(1974)), and 3) charge @ ph 7
# the Unknowns are average of all others.
aa = {'A':[41,67,0],'C':[49,86,0],'D':[-55,91,-1],'E':[31,109,-1],'F':[100,135,0],'G':[0,48,0],\
      'H':[8,118,1],'I':[99,124,0],'K':[-23,135,1],'L':[97,124,0],'M':[74,124,0],'N':[-28,96,0],\
      'P':[-46,90,0],'Q':[-10,114,0],'R':[-14,148,1],'S':[-5,73,0],'T':[13,93,0],'V':[76,105,0],\
      'W':[97,163,0],'Y':[63,141,0],'X':[28,109,0.05]}

from labblouin.PDBnet import PDBstructure
from glob import glob as G
import optparse
from shutil import copy
from os.path import join,isfile
from os import getcwd,remove
from GM.classifyGM import readLandmarks
from GM.GP120classifier import gm2dict,get_list
from random import sample
from copy import deepcopy

def create_cluster(names,path):
    cls = {}
    lists = get_list(path=path)
    for i in set(lists.keys()):
        cls[i] =[]
    #gm = gm2dict(gmfile)
    for e in names:
        for k,v in lists.iteritems():
            if ':' in e:
                e = e[:e.find(':')]
            if join(path,e) in v:
                cls[k].append(e)
    return cls

def select(cls):
    l=[]
    selected=[]
    ncls=[]
    for i in cls.keys():
        if cls[i]:
            l.append(len(cls[i]))
    mini = min(l)
    for k,v in cls.iteritems():
        if v:
            S = sample(v,mini)
            selected.extend(S)
            ncls.extend([k]*mini)
    return selected,ncls
    
def get_data_gm(prefix,inputset):
    data={}
    names = []
    comb ={}
    seq ={}
    land = readLandmarks(prefix)
    with open(prefix+'.gm') as GM:
        for line in GM:
            #temp=[]
            bl=line.strip().strip(';').split(';')
            n=bl[0].strip('"').strip('>')
            stname= n[:n.find(':')]
            datum=[]
            try:
                datum.extend([float(x.strip('"')) for x in bl[1:]])
            except:
                datum.extend([float(x.strip('"')) for x in bl[2:]])
            if isfile(join(inputset,stname+'.pdb')):
                names.append(stname)
                data[stname]=deepcopy(datum)
                st = PDBstructure(join(inputset,'%s.pdb'%(stname)))
                fasta = list(st.ChainAsFASTA(st.chains.keys()[0]))
                afas = []
                for x in land[stname]:
                    T=[]
                    #T.extend([float(ord(fasta[int(x)]))]) # naive sequence
                    #T.extend([float(blosum62[fasta[int(x)]][fasta[int(x)]])])#blosum sequenc
                    [T.extend([float(blosum62[fasta[int(x)]][y])]) for y in blosum62.keys()]#dist to everything
                    T.extend([float(y) for y in aa[fasta[int(x)]]]) # include all properties
                    afas.extend(T)
                seq[stname]=afas
                datum.extend(afas)
                comb[stname]=datum
                datum=list()
    return data , seq, comb, names

def sampleSt(prefix,inputs):
    tempD={}
    tempS={}
    tempC={}
    ncls=[]
    data , seq, comb, names = get_data_gm(prefix,inputs)
    cls = create_cluster(names,inputs)
    S,ncls = select(cls)
    for e in S:
        for i in names:
            if e in i:
                tempD[i]=data[i]
                tempS[i]=seq[i]
                tempC[i]=comb[i]
    return tempD,tempS,tempC,S, ncls

def writeFiles(prefix,data,seq,comb,names,ncls=None,test=False):
    if test: prefix = 'Test'
    stout=open(prefix+'.st.gm','w')
    sqout=open(prefix+'.seq.gm','w')
    cbout=open(prefix+'.comb.gm','w')
    for e in names:
        stout.write('>'+e+';'+';'.join([str(x) for x in data[e]])+'\n')
        sqout.write('>'+e+';'+';'.join([str(y) for y in seq[e]]) +'\n')    
        cbout.write('>'+e+';'+';'.join([str(z) for z in comb[e]])+'\n')    
    cbout.close()
    sqout.close()
    stout.close()
    if not test:
        for i in ['st','seq','comb']:
            FOUT=open(prefix+'.'+i+'.cls','w')
            FOUT.write(';'.join(ncls))
            FOUT.close()

def main(options,args):
    prefix = args[0]
    inputs = args[1]
    cw = getcwd()
    #I = G(join(inputs,'*.pdb'))
    data , seq, comb, names, ncls = sampleSt(prefix,inputs)
    writeFiles(prefix,data,seq,comb,names,ncls,test=False)
    
    if options.predict: 
        problem = options.predict
        #P = G(join(problem,'*.pdb'))
        Tdata , Tseq, Tcomb, Tnames = get_data_gm(prefix,problem)
        writeFiles(prefix,Tdata , Tseq, Tcomb, Tnames,test=True)
        
if __name__ == '__main__':
    opts = optparse.OptionParser(usage='%prog [options] prefix path2inputset')
    opts.add_option('--predict','-p', action = "store", default='',
                    help="Provide the path to a problemset")
    
    options, arg = opts.parse_args()

    main(options,arg)