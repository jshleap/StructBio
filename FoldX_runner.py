#!/usr/bin/python
"""
FoldX_runner Copyright (C) 2012 Jose Sergio Hleap

This script will run FoldX software for repair and mutations to alanine.
It will output a pickled dictionary of energies and the mutated files.
Check FoldX manual for more info.

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

E-mail: jshleap@squalus.org

Requires:
1) FoldX, available at http://foldx.crg.es/
2) Numpy
"""
# importing bit###################################################################################
import os,sys
import cPickle as p
from numpy import mean
from glob import glob as g
from subprocess import Popen, PIPE

# #################################################################################################
# Some python functions############################################################################
def InferSingleLettterCode(threlettercode):
    '''
    Convert the amino acid three letter code, into a single letter code
    '''
    from tl.rename.case import transform_sentence_case
    aa = {'Gly':'G' , 'Ala':'A' , 'Val':'V' , 'Leu':'L' , 'Ile':'I' , 'Met':'M' , 'Phe':'F' ,
          'Trp':'W' , 'Pro':'P' , 'Ser':'S' , 'Thr':'T' , 'Cys':'C' , 'Tyr':'Y' , 'Asn':'N' , 
          'Gln':'Q' , 'Asp':'D' , 'Glu':'E' , 'Lys':'K' , 'Arg':'R' , 'His':'H'}
    singleletter = aa[transform_sentence_case([threlettercode])[0]]
    return singleletter

def parse_landmarks(prefix,pdb):
    eq={}
    fi = open(prefix+'.landmarks').read().split('\n>')
    for st in fi:
        if st == '':
            continue
        elif st.startswith('>'):
            st=st[1:]
        else:
            if st.startswith(pdb):
                bst=st.split('\n')
                name=bst[0].replace('-c','').replace('_','').replace('-','')
                chain=name[4]
                for e in bst[1:]:
                    e=e.strip().split('\t')
                    if e == '':
                        continue
                    eq[e[1]]=InferSingleLettterCode(e[2])
    return eq,chain

def FoldXrunner(prefix,pdb):
    # get the PDB into a txt file
    cat = Popen('ls %s.pdb > PDB.txt'%(pdb),shell=True)
    cat.wait()
    # write the repair runner script
    repair = '''<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>#;
<BATCH>PDB.txt;
<COMMANDS>FOLDX_commandfile;
<RepairPDB>#;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>-CRYSTAL;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>true;
<pdb_hydrogens>false;
<END>#;
<JOBEND>#;
<ENDFILE>#;'''
    if not os.path.isfile('run_repair_%s.txt'%(pdb)):
        fname = 'run_repair_%s.txt'%(pdb)
        fout =  open(fname,'w')
        fout.write(repair)
        fout.close()
    if not os.path.isfile('RepairPDB_%s.pdb'%(pdb)):
        #run repair
        runrep = Popen('FoldX -runfile %s'%(fname),shell=True,stderr=PIPE,stdout=PIPE)
        o,e = runrep.communicate()
        print o, e
    if not os.path.isfile('nPDB.txt'):
        # get the result in a new pdb list file
        np=Popen('ls RepairPDB_%s.pdb > nPDB.txt'%(pdb),shell=True)
        np.wait()
    # get the homologous residues to be mutated
    eq,chain= parse_landmarks(prefix,pdb)
    
    #check if a partial run have been executed
    enf = g('energies_*')
    ms=[]
    if enf:
        #check if all are done
        if len(enf) == len(eq):
            print 'It seems that its been completely done.'
            #continue
        else:
            print 'It seems that it is partially donde... completing.'
            done=[]
            for f in enf:
                res = f[f.find('_')+1:]
                res = res[:res.find('_')]
                done.append(res)
            for k, v in eq.iteritems():
                if not k in done:
                    m= v + chain + k + 'a,'
                    ms.append(m)
                else:
                    continue
    else:
        for k,v in eq.iteritems():
            m= v + chain + k + 'a,'
            ms.append(m)
    m=None
    if ms:
        for m in ms:
            mutate = '''<TITLE>FOLDX_runscript;
<JOBSTART>#;
<PDBS>#;
<BATCH>nPDB.txt;
<COMMANDS>FOLDX_commandfile;
<PositionScan>#,%s;
<END>#;
<OPTIONS>FOLDX_optionfile;
<Temperature>298;
<R>#;
<pH>7;
<IonStrength>0.050;
<water>-CRYSTAL;
<metal>-CRYSTAL;
<VdWDesign>2;
<OutPDB>false;
<pdb_hydrogens>false;
<complex_with_DNA> true;
<END>#;
<JOBEND>#;
<ENDFILE>#;'''%(m[:-1])
            finame='run_mutate%s.txt'%(pdb)
            outf=open(finame,'w')
            outf.write(mutate)
            outf.close()
            mut = Popen('FoldX -runfile %s'%(finame),shell=True,stdout=PIPE,stderr=PIPE)
            o,e=mut.communicate()
            print o,e

def parse_energies(filename):
    f = open(filename).read().split('\n')
    #get the absolute values for every mutation
    l = [abs(float(x.strip().split('\t')[1])) for x in f if not x == '']
    # calculate the difference between the wt (l[0]) and the mutation (y)
    men= [abs(y-l[0]) for y in l[1:]]
    avdiff = mean(men)
    return avdiff

def energies2dict(prefix):
    nrg={}
    files=g('energies_*.txt')
    for f in files:
        res = f[f.find('_')+1:]
        res = res[:res.find('_')]
        nrg[res]=parse_energies(f)
    p.dump(nrg,open(prefix+'_energies.pckl','wb'))
# ##################################################################################################
# Aplication of the code ###########################################################################

if __name__ == "__main__":
    print 'Usage: python FoldX_runner.py <prefix> <pdb>'
    prefix = sys.argv[1]
    pdb = sys.argv[2]
    FoldXrunner(prefix,pdb)
    nrg=energies2dict(prefix)