#!/usr/bin/python
'''
This script will go over the SABmark DATABASE, get te PDBs and a MATT alignment to be feed to
MATT2GM and then to Raw_centroids and GMtoRMSD
'''
#importing bit####################################################################################################
import os, sys, fnmatch
from RawCentroids import get_MATTassignedChains
from subprocess import Popen, PIPE
# End importing####################################################################################################

#Some general variables ###########################################################################################
#list of characters in order for PDB
chainnames=[chr(i) for i in range(65,91)]
chainnames.extend([chr(i) for i in range(49,58)])
chainnames.extend('0')
chainnames.extend(['~','!','@','#','$','%','^','&','(',')','-','+','=','_','>','/',"'",'\\','?','"','[',']','{','}'])
chainnames.extend([chr(i) for i in range(97,123)])
# Database path
db='/home/jshleap/research_sergio/projects/modularity/Datasets/SABmark/sup/'
#SABmark results path in MATT website
smr='http://groups.csail.mit.edu/cb/matt/sabmark/'

#Some definitions##################################################################################################

def check_matt_fasta(prefix):
	'''
	will re write the fasta so is consistent
	'''
	count=0
	f = open('%s.fasta'%(prefix)).read().split('>')
	mv = Popen('mv %s.fasta %s.old.fasta'%(prefix,prefix),shell=True)
	mv.wait()
	newfasta=open('%s.fasta'%(prefix),'w')
	for e in f:
		if e == '':
			continue
		elif ':' in e:
			pdb='>'+e[:e.find(':')]
			chain=':'+chainnames[count]
			seq=e[e.find('\n'):]
			if '.' in seq:
				seq=seq.replace('.','-')
		else:
			pdb='>'+e[:e.find('\n')]
			chain=':'+chainnames[count]
			seq=e[e.find('\n'):]
			if '.' in seq:
				seq=seq.replace('.','-')			
		newfasta.write(pdb+chain+seq)
		count += 1
	newfasta.close()
			
			
	

def get_matt_results_per_group(group,smr,alphac):
	'''
	will get the PDB and the fasta from MATT, will re write the fasta so is consistent
	'''
	if not os.path.isfile(db+group+'.pdb'):
		pdb = Popen('wget %s.pdb'%(smr+group), shell=True)
		pdb.wait()	
	if not os.path.isfile(db+group+'.fasta'):
		fasta=Popen('wget %s.fasta'%(smr+group), shell=True)
		fasta.wait()
	check_matt_fasta(group)
	
	if alphac:
		try:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s -alphaC'%(group), shell=True)
			matt2gm.wait()
		except:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s -alphaC'%(group), shell=True, stderr=PIPE)
			fout = open(group+'.stderr','a')
			fout.write(matt2gm.communicate()[1])
			fout.close()
	else:
		try:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s'%(group), shell=True)
			matt2gm.wait()
		except:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s'%(group), shell=True, stderr=PIPE)
			fout = open(group+'.stderr','a')
			fout.write(matt2gm.communicate()[1])
			fout.close()	
		
def parse_fasta(group_PATH):
	'''
	will parse the fasta and return a list with the PDBs...
	'''
	PDBs=[]
	f = open('%s/group.fasta'%(group_PATH)).read().split('>')
	for e in f:
		if e == '':
			continue
		else:
			bline=e.split('\n')
			PDBs.append(bline[0].strip())
	return PDBs
	
def get_ind_PDBs(group,db, PDBs):
	'''
	wil get over fasta files and retrieve the individual PDBs
	'''
	pdbfolder = '/home/jshleap/research_sergio/projects/modularity/Datasets/SABmark/pdbs_install/'
	for p in PDBs:
		mv = Popen('cp %s.ent ./%s.pdb'%(pdbfolder+p,p),shell=True)
	try:
		rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s -sup'%(group), shell=True)
		rawcentroid.wait()
	except:
		rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s -sup'%(group), shell=True, stderr=PIPE)
		fout = open(group+'.stderr','a')
		fout.write(rawcentroid.communicate()[1])			
		fout.close()	


def set_data_n_execute(db,smr,alphac):
	'''
	will get over every folder of the superfamily dataset of SABmark and leave the individuals PDBs
	as well as the alignment from matt
	'''
	directories = fnmatch.filter(os.listdir(db),'group*')
	for d in directories:
		get_matt_results_per_group(d,smr,alphac)
		get_ind_PDBs(d,db,parse_fasta(db+d))
		try:
			gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(d), shell=True)
			gm2rmsd.wait()
		except:
			gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(d), shell=True, stderr=PIPE)
			fout = open(d+'.stderr','a')
			fout.write(gm2rmsd.communicate()[1])			
			fout.close()		
		mv = Popen('mv %s.* *.pdb *.gm *.bin %s/'%(d,db+d),shell=True)
		mv.wait()
	
	walkerextended = Popen('walker.py -extend',shell=True)
	walkerextended.wait()
	walkerav = Popen('walker.py -average',shell=True)
	walkerav.wait()
	walkeravmatt = Popen('walker.py -avmatt',shell=True)
	walkeravmatt.wait()


# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	if '-alphaC' in sys.argv:
		alphac = True
	else:
		alphac = False
	set_data_n_execute(db,smr,alphac)