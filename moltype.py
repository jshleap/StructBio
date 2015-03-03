from subprocess import Popen, PIPE
import sys
from tl.rename.case import transform_lowercase

def get_molecule_type(prefix, ext='.ent'):
    ''' 
    will get the download .extension file (either ent or pdb) and parse the molecule type out into
    a file
    '''
    d = parse_equivalencesFile(prefix)
    outf=open(prefix+'.moltype','w')
    l=[]
    for k,v in d.iteritems():
        q = Popen("cat pdb%s.ent | grep 'COMPND   2 MOLECULE'"%(transform_lowercase([k])[0]),shell=True, stdout=PIPE)
        out=q.communicate()[0].split('COMPND   2 MOLECULE:')[1].strip().strip(';')
        l.append(out+'\t'+k+'\t'+v)
    l.sort()
    for e in l:
        outf.write(e)

def parse_equivalencesFile(prefix):
    ''' 
    Parse equivalences file and return a dictionary with PDB code as key and 
    plot equivalence as value
    '''
    equ=open(prefix+'_plots.equivalences')
    d={}
    for line in equ:
        b=line.split('\t')
        d[b[1]]=b[0]+'\t'+b[2]
    return d
get_molecule_type(sys.argv[1])