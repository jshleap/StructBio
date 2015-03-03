#!/usr/bin/python
''' 
Plotting the principal coordinate analysis eigen vectors in 2D and 3D, with taxonomic levels
'''
# importing bit##########################################################################################
import matplotlib , sys , os , time
import matplotlib.pyplot as plt
from protGM import fasta, find_key
from mpl_toolkits.mplot3d import Axes3D
from Bio import Entrez
from tl.rename.case import sentence_case
import cPickle as pickle
from subprocess import Popen, PIPE
from glob import glob as G
# ########################################################################################################

# ########################################################################################################
# Some python functions####################################################################################
# #########################################################################################################

def get_id(database,query):
	if '_' in query:
		query = query[:query.find('_')]
	handle = Entrez.esearch(db=database, term=query)
	record = Entrez.read(handle)
	ID = record['IdList']
	return ID

def fetch_info(database,ID,level):
	h = Entrez.efetch(db=database, id=ID)
	r = Entrez.read(h)
	t=r[0]['Lineage'].split(';')
	group = t[level].strip()
	return group,t    

def entrez_query(db,term,level=1):
	'''
	make a entrez query using biopython. Level is the taxonomic level we want to plot. 
	It should be an integer
	from 1 to ntaxlevels which actually varies from organism to organism
	'''
	try:
		group, lin = fetch_info(db,get_id(db,term),level)
	except:
		term = term.split()[0]
		group = fetch_info(db,get_id(db,term),level)
	return group

def number_sps(prefix,db,level,remote):
	'''
	read a fasta and gm file file, and return a dictiorary the taxonomic level 
	(as in domain, class, etc...) as keys and the sp as values. 
	Also will return a set with the different groups.
	'''
	print 'You choose taxonomic level %d:'%(level)
	print '\t1 normally is superkingdom or domain'
	print '\t>1 depends on the lineage reported in NCBI'
	if not remote:
		print "\nGoing locally:"
		# import dictionary
		P = Popen('echo $BLASTDB',shell=True, stderr=PIPE, stdout=PIPE)
		PATH = P.communicate()[0].strip()
		print '\tImporting dictionary with species and lineage from %s.'%(PATH)
		print '\tIf the path above is blank, set the $BLASTDB environment.'\
		      ' Make sure that the TaxDB.bin is in such path'
		current=time.time()
		d=pickle.load(open('%s/TaxDB.bin'%(PATH)))
		print '\tTime elapsed to load the dictionary: %f'%(time.time() - current)
	dgroup = {}
	temp = []
	#spslist=[]
	F = fasta(prefix+'.fasta')
	ax1 , ax2 , ax3 , nvec , allax = parse_eigenvectors(prefix)
	#if len(F.n) != nvec:
	#print 'Some individuals are missing.'
	print 'Looking for the final gm file.'
	pdbs, ind = get_gmfile(prefix)
	if pdbs and nvec == ind:
		count = 1
		for e in range(len(pdbs)):
			if remote:
				#print 'Sit thight, this will take a while, since it will go remotely to the NCBI taxnonomy database.'\ 
				#      ' Remember to run this in a weekend to avoid be blacklisted, or just because is going to take very'\
				#      ' long.'
				while 1:
					try:
						group = entrez_query(db,F.chains[pdbs[e][:-1]],level)
						break
					except:
						time.sleep(10)
			else:
				group = get_sp_locally(pdbs[e],d,level)

			group = sentence_case(group)

			if not group in dgroup:
				dgroup[group]=[]
			if dpc:
				dgroup[group].append((F.chains[pdbs[e][:4]],allax[e][X-1],allax[e][Y-1],ax3[e]))
			else:
				dgroup[group].append((F.chains[pdbs[e][:4]],ax1[e],ax2[e],ax3[e]))
			temp.append(group)
			#dgroup[group].sort()
			#spslist.append(F.chains[e])
			count += 1
		s= set(temp)
		return dgroup, s
	else:
		print 'Something is wrong with the files. Please check them and try again'
		sys.exit()

def get_sp_locally(pdbnchain,dictionary,level):
	''' will get the species name if locally'''
	lin=''
	#correct name
	if '_' in pdbnchain:
		pdb=pdbnchain[:pdbnchain.find('_')]
		chain=pdbnchain[pdbnchain.find('_')+1]
	else:
		pdb=pdbnchain[:4]
		chain=pdbnchain[-1]
	#get the sps name
	bl=Popen('fastacmd -d pdbaa -p T -T T -s "pdb|%s|%s"'%(pdb,chain),shell=True, stdout=PIPE, stderr=PIPE)
	o,e = bl.communicate()
	sname = o.strip().split('\nScientific name:')[-1]
	try:
		lin = dictionary[sname].split(';')
	except:
		for k, v in dictionary.iteritems():
			if sname in k:
				lin = v.split(';')
			else:
				if sname.split()[0] in k:
					lin = v.split(';')
	group = lin[level-1]
	return group

def get_gmfile(prefix):
	''' Get and parse the final GM file, and return a list of PDBs and a total number '''
	nind=0
	gmPDBs=[]
	if os.path.isfile(prefix+'.gm'):
		fil = open(prefix+'.gm')
		for l in fil:
			bl = l.split(';')
			bl = bl[0].split(':')
			pdb = bl[1]#[:-1]
			gmPDBs.append(pdb)
			nind+=1
		return gmPDBs, nind
	else:
		print 'No file named %s.gm. Check prefix or gm name.'%(prefix)
		sys.exit()

def get_list_file(prefix):
	'''
	get and parse the list filter file if MSA_pairwise_fester.py was used to align. will return a list
	of pdb codes
	'''
	if 'final' in prefix:
		nprefix = prefix[:prefix.rfind('_')]
		fname = 'list_filter_'+nprefix + '.txt'
	else:
		fname = 'list_filter_'+prefix+'.txt'
	if os.path.isfile(fname):
		fi = open(fname)
		pdbs=[]
		for line in fi:
			pdbs.append(line.strip()[:line.find('.pdb')])

		return pdbs
	else:
		print 'No filter file provided!!! Provide the filter file, check your fasta, eigenvector file and/or prefix'
		sys.exit()

def parse_eigenvectors(prefix):
	'''
	will take the file of the form <prefix>_vecs.temp, parse it and return 3 lists with the first
	three principal coordinate axes.
	'''
	ax1 = []
	ax2 = []
	ax3 = []
	allax=[]
	nvec=0
	f = open(prefix+'_vecs.temp')
	for line in f:
		allax.append([])
		nvec+=1
		if line == '':
			continue
		else:
			bl = line.split()
			ax1.append(bl[0])
			ax2.append(bl[1])
			ax3.append(bl[2])
			for i in range(len(bl)):
				allax[-1].append(bl[i])
	return ax1,ax2,ax3,nvec,allax

def read_dgroup(dgroup,groupset):
	''' will get the dgroup dictionary with the tuples and create separate lists for plotting'''
	#sort dgroup
	masterd={}
	for e in sorted(groupset):
		mastertup = ([],[],[],[])
		for el in sorted(dgroup[e]):
			mastertup[0].append(el[0])
			mastertup[1].append(float(el[1]))
			mastertup[2].append(float(el[2]))
			mastertup[3].append(float(el[3]))
		masterd[e]=mastertup
	return masterd

def plot(prefix,masterd,fontsize,symlog,threeD):
	''' will plot two axes of the PCoA, colored by taxonomic level '''
	fil = open(prefix+'_PCoA_axes%d_%d.equivalences'%(X,Y),'w')
	F = fasta(prefix+'.fasta')
	markers = ['k.','b+','g*','r.','c+','m*','y.','k+','b*','g.','r+','c*','m.','y+','k*','b.',
	           'g+','r*','c.','m+','y*']
	count = 0
	c=0
	fig = plt.figure()#figsize=(6.83 , 9.19), dpi=300)
	ax = fig.add_subplot(111)
	if threeD:
		fig3D = plt.figure()
		ax3D = fig3D.gca(projection='3d')
	ax.spines['top'].set_color('none')
	ax.xaxis.tick_bottom()
	ax.spines['right'].set_color('none')
	ax.yaxis.tick_left()
	for k in masterd:
		x = masterd[k][1]
		y = masterd[k][2]
		z = masterd[k][3]
		ax.plot(x,y, markers[c],label=k)
		if threeD:
			ax3D.plot(x, y, z, markers[c],label=k)
		c+=1
		f=0
		for e in range(len(masterd[k][1])):
			ax.annotate(count,(x[f]+0.1,y[f]+0.1),fontsize=fontsize)
			if threeD:
				ax3D.text(x[f]+0.1,y[f]+0.1,z[f]+0.1, str(count),fontsize=fontsize)
			fil.write(masterd[k][0][e]+'\t'+find_key(F.chains,masterd[k][0][e])+'\t'+str(count)+'\n')
			count+=1
			f+=1
	if symlog:
		ax.set_xscale("symlog")
		ax.set_yscale("symlog")
		ax.set_xlabel('Axis %d (symmetrical log)'%(X), fontsize=fontsize)
		ax.set_ylabel('Axis %d (symmetrical log)'%(Y), fontsize=fontsize)
	else:
		ax.set_xlabel('Axis %d'%(X), fontsize=fontsize)
		ax.set_ylabel('Axis %d'%(Y), fontsize=fontsize)
	ax.legend(loc=0, #bbox_to_anchor=(0.5, 1.1), ncol=4,
	          fancybox=True, shadow=True)
	fig.tight_layout()
	if threeD:
		ax3D.set_xlabel('Axis 1', fontsize=fontsize)
		ax3D.set_ylabel('Axis 2', fontsize=fontsize)
		ax3D.set_zlabel('Axis 3', fontsize=fontsize)
		ax3D.view_init(30, 45)
		fig3D.tight_layout()
		ax3D.legend(loc=0, #bbox_to_anchor=(0.5, -0.075),ncol=4,
		            fancybox=True, shadow=True)
	plt.show()
	fig.savefig(prefix+'_Axis%d_%dPCoA.png'%(X,Y), dpi=300)
	if threeD:
		fig3D.savefig(prefix+'_3AxesPCoA.png', dpi=300)
	fil.close()

def plot_3D(prefix, masterd):
	''' will plot three axes of the PCoA, colored by taxonomic level '''
	markers = ['k.','b.','g.','r.','c.','m.','y.','k*','b*','g*','r*','c*','m*','y*','k+','b+','g+','r+','c+','m+','y+']
	count = 0
	c=0
	fig = plt.figure()#figsize=(6.2,9.2))
	#ax = fig.add_subplot(111, projection='3d')
	ax = fig.gca(projection='3d')
	#t = ([],[])
	for k in masterd:
		count+=1
		x = masterd[k][1]
		y = masterd[k][2]
		z = masterd[k][3]
		ax.plot(x, y, z, markers[c],label=k)
		#ax.scatter(x,y,z,c=markers[c][0],s=15)
		#p = plt.Rectangle((0, 0), 0, 0)
		#t[0].append(markers[c])
		#t[1].append(k)
		c+=1
		f=0
		for e in range(len(masterd[k][1])):
			ax.text(x[f]+0.1,y[f]+0.1,z[f]+0.1, str(count),fontsize=10)
			count+=1
			f+=1
	ax.set_xlabel('Axis 1', fontsize=12)
	ax.set_ylabel('Axis 2', fontsize=12)
	ax.set_zlabel('Axis 3', fontsize=12)
	ax.view_init(30, 45)
	fig.tight_layout()
	ax.legend(loc=0, #bbox_to_anchor=(0.5, -0.075),ncol=4,
	          fancybox=True, shadow=True)
	plt.show()
	plt.savefig(prefix+'_3AxesPCoA.png', dpi=300)

def list2dict():
	ld={}
	lists=G('*.list')
	for f in lists:
		typ = f[:-5]
		inf = open(f)
		for line in inf:
			if not line == '\n' or not line == '':
				ld[line.strip()[:-4]] = typ
	return ld

def equi2dict(prefix):
	equi={}
	fi = open(prefix+'.pdbequi')
	for l in fi:
		if l == '\n' or l == '':
			continue
		else:
			bl=l.strip().split('\t')
			equi[bl[1][:4]]=bl[0]
	return equi

def masterd4nonSp(prefix):
	F = fasta(prefix+'.fasta')
	ax1,ax2,ax3,nvec,allax = parse_eigenvectors(prefix)
	pdbs, ind = get_gmfile(prefix)
	equi=equi2dict(prefix)
	ld = list2dict()
	groupset = set(ld.values())
	dgroup={}
	if pdbs and nvec == ind:
		for e in range(len(pdbs)):
			group = ld[equi[pdbs[e]]]
			if not group in dgroup:
				dgroup[group]=[]
			if dpc:
				dgroup[group].append((F.chains[pdbs[e][:4]],allax[e][X-1],allax[e][Y-1],ax3[e]))
			else:
				dgroup[group].append((F.chains[pdbs[e][:4]],ax1[e],ax2[e],ax3[e]))
		return dgroup, groupset
	else:
		print 'Something is wrong with the files. Please check them and try again'
		sys.exit()



# ######################################################################################################
# Aplication of the code ###############################################################################

if __name__ == "__main__":
	if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'Usage: python PlotPCoA.py <prefix> <email> [options]'
		print 'A <prefix>.fasta must be provided, and the headers should be ">species name:PDB"'
		print 'The level should be an integer that matches the lineage in NCBI taxonomy browser'\
		      ' index (e.g. 1 normally is the domain)'
		print 'If you want to provide a different level, you have to provide the database also'\
		      ' (even if you want taxonomy)'
		print 'Options:'
		print '\t-level=x : Choose the taxonomic level according to NCBI. Normally 1=domain,'\
		      ' 2=order, but varies according to NCBI lineage. (1 or domain by default)'
		print '\t-symlog : Use symmetrical log transformation of the axes. In only positive'\
		      ' values behaves like log. ( Default: No )'
		print '\t-fontsize=XX : Will use XX as the font size for axes titles and numbering'\
		      ' in the plot. ( Default : 12 )'
		print '\t-3D : Also will plot a 3D plot, using the first 3 axes. ( Default: No )'
		print '\t-remote : Will use the E-utilities of NCBI for the taxonomy search.'\
		      ' Otherwise will use a local copy of the databases. The latter require'\
		      ' BLAST+ tools, the pre-formatted databases, and that the TaxDB.bin is'\
		      ' in teh BLASTDB environment/folder ( Default: Local )'
		print '\t-dpc=X,Y : Will use the principal components corresponding to x and y'\
		      ' ( Default: 1,2)'


	# Default Parameters ###############################################################################
	prefix = sys.argv[1]
	email = None
	database = 'taxonomy'
	level = 1
	fontsize= 12
	symlog=False
	threeD = False
	remote=False
	fontsize=12
	dpc=False
	X=1
	Y=2
	# Get user input ###################################################################################

	for arg in sys.argv[1:]:
		if arg == '-symlog':
			symlog=True
		elif arg.startswith('-fontsize='):
			fontsize= int(arg[10:])
		elif arg == '-3D':
			threeD = True
		elif arg.startswith('-level='):
			level = int(arg[7:])
		elif arg.startswith('-dpc='):
			dpc=True
			r=arg[5:].split(',')
			X=int(r[0])
			Y=int(r[1])

	if os.path.isfile(prefix+'.pdbequi'):
		print 'Using the pdbequi and the lists files to color.'
		dgroup, groupset = masterd4nonSp(prefix)
		d = read_dgroup(dgroup, groupset)

	else:
		if os.path.isfile(prefix+'.pickle'):
			d = pickle.load( open( prefix+".pickle" , "rb" ) )
		else:
			print 'Working on %s'%(os.getcwd())
			#Tell ncbi who you are
			if remote:
				Entrez.email = sys.argv[2]
			dgroup , groupset = number_sps(prefix,database,level,remote)
			d = read_dgroup(dgroup,groupset)
			pickle.dump( d, open( prefix+".pickle", "wb" ) ) 
			#plot(prefix,d,fontsize,symlog,threeD)
	plot(prefix,d,fontsize,symlog,threeD)