import sys

def get_equivalences(prefix):
    eq = {}
    f = open(prefix+'.pdbequi')
    for line in f:
        if line == '' or line == '\n':
            continue
        else:
            bl = line.strip().split('\t')
            eq[bl[0]]=bl[1]
    return eq

def parse_landmarks(prefix):
    L = {}
    f = open(prefix + '.landmarks').read().split('\n>')
    for e in f:
        if e == '' or e == '\n':
            continue
        else:
            bl = e.split('\n')
            L[bl[0]]=bl[1:]
    return L

def parse_mFDM(prefix):
    mfdm=[]
    f = open(prefix+'_PCoAprocFdif.txt')
    for line in f:
        if line == '' or line == '\n':
            continue
        elif line.startswith('"Procrustes'):
            continue
        else:
            bl=line.strip().split('\t')
            mfdm.append('\t'.join([bl[0].strip('"'),bl[1].strip('"')]))
    return mfdm

def parse_FD(prefix):
    fds={}
    f = open(prefix+'.betaeq')
    for line in f:
        if line == '' or line == '\n':
            continue
        else:
            bl = line.split('\t')
            fds[bl[0]]=bl[-1]
    FDS = sorted(fds.items(), key=lambda x: int(x[0]))
    return FDS


prefix = sys.argv[1]

#FD =  parse_FD(prefix)
mFD = parse_mFDM(prefix)
#L = parse_landmarks(prefix)
eq = get_equivalences(prefix)
fout = open(prefix+'_summary.csv','w')

headers = 'Structure name\tIndex\tResidue index\tAminoAcid\tEigenvalue Centrality\t'\
    'Betweenness centrality\tCloseness centrality\tDegree\tmFDM (distance)\tmFDM (RMSD)\t'\
    'Normalized mFDM distance\tModule membership\n'
body = ''
C=open(prefix+'.centrality')
mv = open(prefix+'.graphcluster').read().strip().split()
for line in C:
    if line == '' or line == '\n':
        continue
    elif line.startswith('The') or line.startswith('Index'):
        continue
    elif line.startswith('>'):
        name = line[1:].strip()
	try:
		body += '>'+name+':'+eq[name]+'\t'
	except:
		for k,v in eq.iteritems():
			if name == v: n=k
        	body += '>'+name+':'+eq[n]+'\t'
        c = 0
    else:
        if c == 0:
            body += line.strip() + '\t' + mFD[c]+'\t'
        else:
            body += '\t'+ line.strip() + '\t' + mFD[c]+'\t'
        body += str(float(mFD[c].split('\t')[0])/max([float(x.split('\t')[0]) for x in mFD]))
        body += '\t'+mv[c]+'\n'
        c+=1

fout.write(headers+body)
                
