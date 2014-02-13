#!/usr/bin/python
"""
MapModules Copyright (C) 2012 Jose Sergio Hleap

Map clustering partition into the chain field of a PDB structure, the centrality or form difference into beta and 
ocupancy. This also works for a multiple pdb(>50) with the -m option.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

Requires:
1) Utils folder, available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools
"""

#importing bit####################################################################################################
import sys,os
from utils import PDBnet
from shutil import move
from copy import deepcopy
# End importing###################################################################################################

#Some constants###################################################################################################
chainn = 'ABCDEFGHIJKLMNOPQRSTUVWYZ'
chainn += chainn.lower()
chainn += '0123456789`-=[]\;,./~!@#$%^&*()_+{}|:"<>'
## The X chain is for singletons
chainn = chainn.replace('X','')
#end constants####################################################################################################

#Some definitions#################################################################################################

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

def round_str(x,typ='oc'):
    sig = 4
    if typ == 'oc':
        if len(str(x)) > 6:
            sig = 2
    n = round(float(x),sig)
    return str(n)[:8]

def parse_centralities(prefix,centrality):
    '''
    Parse the centrality file wich also contains the landmarks file info will return a dictionary with the chain
    as key. In the case of multiple PDBs (>50), the chains are the PDB codes '''
    land = {}
    # Chains equivalencies
    chains = {}
    if centrality:
        fin = open(prefix+'.centrality')
    else:
        fin = open(prefix+'.landmarks')
        
    for line in fin:
        if line.startswith('The') or line.startswith('Index'):
            continue
        # New chain
        elif line.startswith('>'):
            ch = line[1:].strip()
            chains[ch] = []
            land[ch]={}
        else:
            # Add to chain
            line = line.split()
            chains[ch].append(int(line[1]))
            res = int(line[1])
            if centrality:
                land[ch][res] = line[3:]
    return land,chains

def map_consistency(prefix):
    ''' will create a new graphcluster consistent with the mapped chains'''
    # Read in cluster
    membership = open(prefix+'.graphcluster').read().split()

    # renumber clusters
    temp = list(set(membership))
    tempd={}
    if '?' in temp:
        temp.remove('?')
    for t in range(len(temp)):
        tempd[temp[t]]= chainn[t]
    for i in range(len(membership)):
        if not membership[i] == '?':
            membership[i] = tempd[membership[i]]
    # write a new graphcluster that is consistent with the mapping
    outf=open(prefix+'.mapped.graphcluster','w')
    for e in membership:
        outf.write(e + ' ')
    outf.close()
    return membership

def rewrite_PDB(chain,membership,chains,land,centrality,multiple):
    ''' 
    re-write the PDB with the cluster and centrality. Multiple is if more than 50 pdbs are used
    and no MStA PDB is availlable
    '''
    fout = open(prefix + '.map%s.pdb'%(chain),'w')    
    
    if not multiple:
        fin = open(prefix + ".pdb")
        #write the PDB file
        for line in fin:
            if line.startswith('ATOM') and line[21] == chain:
                chainmap = chains[chain]
                index = int(line[22:26])
                if index in chainmap:
                    # Find the cluster
                    x = chainmap.index(index)
                    x = membership[x]
                    if not x == '?':
                        newchain = x
                    else:
                        newchain = '?'
                else:
                    newchain = 'X'
                if not centrality:
                    line = line[:21] + newchain + line[22:]#55]+'-1.00'.rjust(5)+'-1.00'.rjust(6)+line[67:]
                else:
                    if newchain == 'X':
                        line = line[:21] + newchain + line[22:55] +'0.00'.rjust(5)+'0.00'.rjust(6)+line[67:]
                    else:
                        line = line[:21] + newchain + line[22:55] + round_str(land[chain][int(line[22:26])][centrality[0]])+ \
                        ' '+round_str(land[chain][int(line[22:26])][centrality[1]])+line[65:]
                fout.write(line)
        fin.close()
    else:
        fin = open(chain + ".pdb")
        st = PDBnet.PDBstructure(chain + ".pdb")
        chainmap = chains[chain]
        for line in fin:
            if not line.startswith('ATOM'):
                continue
            res = line[22:26].strip()
            if int(res) in chainmap:
                # Find the cluster
                x = chainmap.index(int(res))
                newchain = membership[x]                
            else:
                newchain = 'X'
            if centrality:
                if newchain == 'X':
                    l = line[:21]+newchain+line[22:54]+'0.00'.rjust(5)+'0.00'.rjust(6)+line[67:]
                else:
                    l = line[:21]+newchain+line[22:54]+round_str(land[chain][int(r)][centrality[0]]).rjust(5) \
                        + round_str(land[chain][int(r)][centrality[1]]).rjust(6)+line[67:]
            else:
                l = line[:21]+newchain+line[22:]
            fout.write(l)
        '''
        for r in st.chainsOrder[chain[-1]]:
            chainmap = chains[chain]
            if int(r) in chainmap:
                # Find the cluster
                x = chainmap.index(int(r))
                newchain = membership[x]
            else:
                newchain = 'X'
            for i in sorted(st.chains[chain[-1]][r].atoms, key=st.chains[chain[-1]][r].atoms.get):
                a = st.chains[chain[-1]][r].atoms[i]
                if not centrality:
                    line = 'ATOM  ' + str(a.serial).rjust(5) + '  '+a.name.ljust(3)+' '+ st.chains[chain[-1]][r].name.rjust(3)\
                        + ' ' + newchain + str(r).rjust(4) + ' '+ round_str(a.x,'coordinates').rjust(11) + \
                        round_str(a.y,'coordinates').rjust(8) + round_str(a.z,'coordinates').rjust(8) + ' ' + \
                        str(a.occupancy).rjust(5)+str(a.tempFactor).rjust(6) + ' '*11 + a.symbol.rjust(2)+a.charge+'\n'
                elif newchain == 'X':
                    line = 'ATOM  ' + str(a.serial).rjust(5) + '  '+a.name.ljust(3)+' '+ st.chains[chain[-1]][r].name.rjust(3)\
                    + ' ' + newchain + str(r).rjust(4) + ' '+ round_str(a.x,'coordinates').rjust(11) + \
                    round_str(a.y,'coordinates').rjust(8) + round_str(a.z,'coordinates').rjust(8) + ' ' + \
                    '-1.00'.rjust(5)+'-1.00'.rjust(6) + ' '*11 + a.symbol.rjust(2)+a.charge+'\n'
                else:
                    line = 'ATOM  ' + str(a.serial).rjust(5) + '  '+a.name.ljust(3)+' '+ st.chains[chain[-1]][r].name.rjust(3)\
                        + ' ' + newchain + str(r).rjust(4) + ' ' + round_str(a.x,'coordinates').rjust(11) + \
                        round_str(a.y,'coordinates').rjust(8) + round_str(a.z,'coordinates').rjust(8) + ' ' + \
                        round_str(land[chain][int(r)][centrality[0]]).rjust(5) + round_str(land[chain][int(r)][centrality[1]]).rjust(6)\
                        + ' '*10 + a.symbol.rjust(2)+a.charge+'\n'
                fout.write(line)
                '''
    fout.close()

def FDMDictByRes(prefix, chain):
    D={}
    L = land2dict(prefix+'.landmarks')
    kl = L[chain].strip().split('\n')
    for e in kl:
        bl = e.strip().split('\t')
        D[bl[1]]=int(bl[0])
    return D

def dummychain(st):
    DUMMY = deepcopy(st.chains)
    DUMMYO= deepcopy(st.chainsOrder)
    st.chains['dummy']={}
    st.chainsOrder['dummy']=[]
    for ch in DUMMY:
        st.chains['dummy'].update(st.chains[ch])
        del st.chains[ch]
    for cha in DUMMYO:
        st.chainsOrder['dummy'].extend(st.chainsOrder[cha])
        del st.chainsOrder[cha]    
    
def mapModules(st,membership,D):
    dummychain(st)
    sm = set(membership)
    sm.add('X')
    for m in sm:
        st.chains[m]={}
        st.chainsOrder[m]=[]
    for k,v in D.iteritems():
        st.AddResidueToChain(membership[D[k]], st.chains['dummy'].pop(str(k)))
    DUMMY = deepcopy(st.chains)
    for re in DUMMY['dummy']:
        st.AddResidueToChain('X', st.chains['dummy'].pop(re))
    if not st.chains['dummy']:
        del st.chainsOrder['dummy']
        del st.chains['dummy']
    else:
        print 'Variable dummy not completely remove!!! check output!'

    
def rewritePDBifFDM(prefix,chain, multiple,fdm,membership):
    D = FDMDictByRes(prefix, chain)
    if not multiple:
        print 'not coded yet'
        sys.exit(-1)
    else:
        st = PDBnet.PDBstructure(chain+'.pdb')
        for ch in st.chains:
            for res in st.chains[ch]:
                if str(res) in D.keys():
                    for a in st.chains[ch][res].atoms:
                        st.chains[ch][res].atoms[a].tempFactor=float(fdm[D[str(res)]])
                else:
                    for a in st.chains[ch][res].atoms:
                        st.chains[ch][res].atoms[a].tempFactor=-1.0
        mapModules(st,membership,D)
        st.orderofchains = st.chainsOrder.keys()
        organizeChains(st)
    st.WriteFile(prefix + '.map%s.FD.pdb'%(chain))

def organizeChains(st):
    lt=[]
    for ch in st.chains:
        st.chainsOrder[ch] = sorted(st.chainsOrder[ch], key=lambda x: int(x))
    
def tFDM2land(prefix):
    if not os.path.isfile(prefix+'.tradFDM') and os.path.isfile(prefix+'.gm') and tfd:
        os.system('Rscript ~/LabBlouin/code/GM/tradFDM.R %s %d'%(prefix+'.gm',3))
    elif not os.path.isfile(prefix+'.tradFDM') and os.path.isfile(prefix+'.gm') and tFD:
        os.system('Rscript ~/LabBlouin/code/GM/tradFDM.R %s %d %s'%(prefix+'.gm',3,prefix+'.graphcluster'))
    fdm = open(prefix+'.tradFDM').read().strip().split('\n')
    bta = []
    if os.path.isfile(prefix+'.betaeq'):
        D = land2dict(prefix+'.betaeq')
        for v in D.itervalues():
            bl = v.strip().split('\n')
            for e in bl:
                bta.append(e.strip().split('t')[-1])
    return fdm, bta
            
# End of definitions############################################################################

# Aplication of the code #######################################################################
if __name__ == "__main__": 
    if len(sys.argv) == 1 or '-help' in sys.argv:
        print 'usage MapModules.py [prefix] [chain] [option]'
        print '\tOptions:'
        print '\tprefix: Not really an option, is required. Have to be the same prefix for landmarkfile'\
              +'or centrality file, PDB file (if not using multiple PDB option), and the graphcluster file.'
        print '\tchain : Not really an option, is required. If used without the multiple option, should '\
              +'match the chains in the landmarks file. If using the multiple option, should be the prefix'\
              +'for the particular PDB file that you want to map.'
        print '\t-centrality=XXX,YYY: Two measurement of centrality are allowed to be map. You can '\
              +'choose between Eigenvector centrality (evcen), betweenness centrality (btc), closeness '\
              +'centrality (clc), and/or degree.'
        print '\t-m : choose this option if you dont have a traditional MStA but more than 50 structures '\
              +'aligned with the pairwise-multiple alignment strategy.'
        print '\t-fd: it would map the Form distortion into the beta factor for the overall structure.'
        print '\t-FD: it would map the Form distortion into the beta factor for each module in the structure.'

    # Default Parameters ###################################################
    prefix = sys.argv[1]
    chain = sys.argv[2]
    centrality=False
    multiple=False
    tfd = False
    tFD = False
    # Command line input ###################################################
    for arg in sys.argv[1:]:
        if arg.startswith('-centrality='):
            bline = arg[12:].split(',')
            if bline[0] == 'evcen':
                a=0
            elif bline[0] == 'btc':
                a=1
            elif bline[0] == 'clc':
                a=2
            elif bline[0] == 'degree':
                a=3
            elif bline[1] == 'evcen':
                b=0
            if bline[1] == 'btc':
                b=1
            elif bline[1] == 'clc':
                b=2
            elif bline[1] == 'degree':
                b=3
            centrality = (a,b)
        elif arg == '-m':
            multiple=True
        elif arg == '-fd':
            tfd = True
        elif arg == '-FD':
            tFD = True
    membership = map_consistency(prefix)
    if tfd or tFD:
        fdm, bta = tFDM2land(prefix)
        rewritePDBifFDM(prefix,chain, multiple,fdm,membership)
    #if centrality:
    land,chains = parse_centralities(prefix,centrality)
    rewrite_PDB(chain,membership,chains,land,centrality,multiple)

