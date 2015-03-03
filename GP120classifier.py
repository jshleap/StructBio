#!/usr/bin/env python
'''
create classifiers for connor's GP120 datasets
'''

import sys,os
import string as S
from glob import glob as G
from random import sample
from shutil import move
from protGM import find_key


letters=list(S.uppercase)
numbers=[0,1,2,3,4,5,6,7,8,9]

def create_artificial_PDBcodes(npdb,letters):
    aPDBs=[]
    indexes=[]
    for i in range(len(numbers)):
        for j in range(len(letters)):
            for k in range(len(letters)):
                for l in range(len(letters)):
                    indexes.append((i,j,k,l))
                
    for n in range(npdb+1):
        aPDBs.append('%s%s%s%s_'%(numbers[indexes[n][0]],letters[indexes[n][1]],letters[indexes[n][2]],letters[indexes[n][3]]))
        
    return aPDBs
        
def format_gm(prefix,letters,npdb):
    PDBs=create_artificial_PDBcodes(npdb,letters)
    fi = open(prefix+'.gm')
    equi=open(prefix+'.pdbequi','w')
    lines=''
    c=-1
    eq={}
    for line in fi:
        c+=1
        oldname=line[1:line.find(':')]
        PDBname=PDBs[c]
        rest=line[line.find(';'):]
        lines+='>'+PDBname+':'+PDBname[:-1]+rest
        eq[oldname]=PDBname
        os.system('mv "%s.pdb" "%s_.pdb"'%(oldname,PDBname[:-1]))
        equi.write(oldname+'\t'+PDBname+'\n')
    fi.close()
    os.system('mv %s.gm %s_original.gm'%(prefix,prefix))
    fout=open(prefix+'.gm','w')
    fout.write(lines)
    fout.close()
    equi.close()
    return eq
    
def format_landmarks(prefix,eq):
    fil = open(prefix+'.landmarks')
    L=''
    for linea in fil:
        if linea.startswith('>'):
            na=linea.strip()[1:]
            L+='>'+eq[na]+'\n'
        else:
            L+=linea
    fil.close()
    os.system('mv %s.landmarks %s_original.landmarks'%(prefix,prefix))
    ouf = open(prefix+'.landmarks','w')
    ouf.write(L)
    ouf.close()

def get_list(path='.'):
    # open the list files and turn them into a dictionary
    lists=G(os.path.join(path,'*.list'))
    groups={}
    if lists:
        for l in lists:
            group=l[:-5]
            groups[group]=[]
            with open(l) as L:
                for line in L:
                    if line.strip() == '':
                        continue
                    bl=line.strip().strip('.pdb')
                    if not bl in groups[group]:
                        groups[group].append(bl)
    else:
        pdbs=G(os.path.join(path,'*.pdb'))
        for p in pdbs:
            q = os.path.split(p)[-1].strip().strip('.pdb')
            q = os.path.join(path,q)
            if '.CXCR4.' in p:
                if not 'X4' in groups:
                    groups['X4']=[]
                groups['X4'].append(q)
            elif '.CCR5_CXCR4.' in p:
                if not 'DUAL' in groups:
                    groups['DUAL']=[]
                groups['DUAL'].append(q)            
            elif '.CCR5.' in p:
                if not 'R5' in groups:
                    groups['R5']=[]
                groups['R5'].append(q)            
                        
    return groups

def checknsample(prefix, groups,eq,gmfile):
    linespg={}
    lenght=[]
    namesingm=[]
    for k, va in groups.iteritems():
        lenght.append(len(va))
        linespg[k]=[]
        for l in open(prefix+'.gm'):
            for sn in va:
                if sn == '':
                    continue
                if eq[sn] in l:
                    linespg[k].append(l)
    minimum=min(lenght)
    #minimum=40 #workaround to get contacts... uncomment previous line to restore
    move(prefix+'.gm',prefix+'_beforesampling.gm')
    landm = land2dict(prefix+'.landmarks')
    move(prefix+'.landmarks',prefix+'_beforesampling.landmarks')
    fout=open(prefix+'.gm','w')
    outf=open(prefix+'.landmarks','w')
    for ke, val in linespg.iteritems():
        gout = open(ke+'.gm','w')
        S = sample(val,minimum)
        for i in S:
            fout.write(i)
            gout.write(i)
            namesingm.append(find_key(eq,i[1:i.find(':')]))
            for key, value in landm.iteritems():
                if key in i:
                    outf.write('>'+key+'\n'+value+'\n')
        gout.close()
    return namesingm

def gm2dict(fname):
    d={}
    with open(fname) as F:
        for line in F:
            bl=line.split(';')
            d[bl[0][1:bl[0].find(':')]]=';'.join(bl[1:])
    return d

def land2dict(fname):
    d={}
    F = open(fname).read().split('\n>')
    for e in F:
        if e == '':
            continue
        else:
            be=e.split('\n')
            d[be[0].strip('>')]='\n'.join(be[1:])
    return d

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def ifnotGM(prefix):
    count=0
    if not os.path.isfile(prefix+'.landmarks'):
        la = land2dict('gp120_conor.landmarks')
        land = open(prefix+'.landmarks','w')
        L=True
    if not os.path.isfile(prefix+'.gm'):
        f = open(prefix+'.gm','w')
        l = open(prefix+'.list').read().split('.pdb\n')
        gm = gm2dict('gp120_conor.gm')
        for e in l:
            if e == '':
                continue
            else:
                f.write('>%s:%s;%s'%(e,e[:5],gm[e]))
                count+=1
                if L:
                    land.write('>%s\n%s\n'%(e,la[e]))
    else:
        count = file_len(prefix+'.gm')
    
    return count

def main(prefix,samp):
    npdb = ifnotGM(prefix)
    #format gm
    eq=format_gm(prefix,letters,npdb)
    format_landmarks(prefix,eq)
    groups = get_list()
    # open general file
    gmfile=open(prefix+'_original.gm')
    if samp:
        namesingm = checknsample(prefix, groups,eq,gmfile)
    else:
        #store names of the GM file in order
        namesingm=[]
        for lin in gmfile:
            n = lin[1:lin.find(':')]
            namesingm.append(n)
    
    # write the file
    fout=open(prefix+'.cls','w')
    lines=''
    for name in namesingm:
        for k, v in groups.iteritems():
            if name in v:
                lines+=k+';'
    
    fout.write(lines[:-1]+'\n')
    fout.close()
    
if __name__ == "__main__":
    prefix = sys.argv[1]
    #npdb=int(sys.argv[2])
    if '-s' in sys.argv: samp = True
    else: samp = False
    main(prefix,samp)
    
