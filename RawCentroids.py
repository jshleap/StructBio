#!/usr/bin/python
''' 
This script will create a GM file with the raw coordinates (not aligned) of a given set of proteins
that have been aligned with MATT. This is to get the homologous landmarks on proteins.
will require the original PDBs, the txt ouput of matt where each structure have been named, and the
landmark file created by MATT2GM
'''
#importing bit#####################################################################################################
import sys, os
from labblouin.PDBnet import PDBstructure
from rpy2.robjects import r
from subprocess import Popen
from glob import glob
# End importing####################################################################################################
#Some definitions##################################################################################################
def parse_fasta(prefix):
	'''
	will parse the fasta file returned by MATT (or MATT-like), and return a dictionary of PDB codes and chains
	'''
	PDBs={}
	fil = open(prefix+'.fasta').read().split('>')
	for el in fil:
		if el == '':
			continue
		else:
			bline = el.split('\n')
			for e in bline:
				if e == '':
					continue
				elif not ':' in e:
					continue
				else:
					bbline=bline[0].split(':')
					pdb=bbline[0][:5]
					chain=bbline[1]
					PDBs[pdb]=chain
	return PDBs
	
def get_MATTassignedChains(prefix,txt):
	'''
	Parser of the TXT file outputted by MATT after an alignment. Will return a dictionary of the
	PDB and the corresponding chain name
	'''
	PDBs={}
	if txt:
		fil = open(prefix+'.txt').read().split('\n\n')
		block = fil[1].split('\n')
		
		for e in block:
			pdb = e.split(':')[0]
			pdb = pdb.replace(' ','')
			chain = e[e.find('(')+1:e.find(')\n')]
			PDBs[pdb]=chain
	else:
		PDBs=parse_fasta(prefix)
	return PDBs

def get_homolog_res(prefix,chainname):
	'''
	Will gather the information of the residues positions for a given chainname. Requires the landmark file and
	return a list of residues numbers
	'''
	lfile = open(prefix+'.landmarks').read()
	chainplus = lfile[lfile.find('>%s'%(chainname))+2:]
	chain = chainplus[2:chainplus.find('\n>')]
	res = []
	bline = chain.split('\n')
	for e in bline:
		res.append(e.split('\t')[1])
	
	return res

def parse_landmarkfile(prefix):
	'''
	parser or the landmark file. will return a dictionary with the chain as key and the list of tuples (resindex,res) as value.
	'''
	landict={}
	lchain=[]
	f = open(prefix+'.landmarks').read().split('>')
	for e in f:
		if e == '':
			continue
		else:
			esplit = e.split('\n')
			chain = esplit[0]
			temp=[]
			if not esplit[-1] == '':
				lchain.append(int(esplit[-1].split()[0])+1)
			else:
				lchain.append(int(esplit[-2].split()[0])+1)
			for el in esplit[1:]:
				if el == '':
					continue
				else:
					elsplit = el.split()
					temp.append((elsplit[1].rjust(4),elsplit[2]))
		landict[chain]=temp
		
	return landict, lchain

def parse_pdbs(prefix):
	'''
	will parse the PDB and return a dictionary of list tuples with residue index and residue
	'''
	pdbtuples={}
	pdbs = glob('*.pdb')
	pdbs.pop(pdbs.index(prefix+'.pdb'))
	for p in pdbs:
		temp=[]
		f = open(p).read().split('\n')
		for e in f:
			if not e.startswith('ATOM'):
				continue
			else:
				if e[22+4] != ' ':
					temp.append((e[23:22+5],e[17:16+4]))
				else:
					temp.append((e[22:22+4],e[17:16+4]))
		pdbtuples[p]=temp
	return pdbtuples, pdbs

def get_matt_chain(prefix):
	'''
	Parser of the pdbs in the CWD. Will return a dictionary with the PDBcodes and the apropriate MATT chain
	'''
	PDBs={}
	landict, lchain = parse_landmarkfile(prefix)
	pdbtuples,pdbs = parse_pdbs(prefix)
	for p in pdbs:
		lpdb = pdbtuples[p]
		for k,v in landict.iteritems():
			if set(v).intersection(lpdb) and (len(set(v).intersection(lpdb)) in lchain):
				PDBs[p[:p.rfind('.')]]=k
	
	return PDBs

def write_homologous_coordinates(prefix,PDBs):
	'''
	Will get the homologous coordinates and write a gm file
	'''
	gmlines=[]
	fout=open(prefix+'_GM.gm','w')
	multiple=None
	for k,v in PDBs.iteritems():
		line=''
		line += '>'+k+';'
		if sup:
			pdbfile = k  + '.pdb'
		else:
			pdbfile = k + '.pdb'
		PDB = PDBstructure(pdbfile)
		res = get_homolog_res(prefix,v)
		usedr=[]
		for r in res:
			try:
				k[4]
				notchain=False
			except:
				notchain=True
			if not notchain:
				if k[4] == '_':
					multiple = True
					for c in PDB.chains:
						if r in PDB.chains[c]:
							if r not in usedr:
								PDB.chains[c][r].Centroid()
								line += str(PDB.chains[c][r].centroid.x)+';'+str(PDB.chains[c][r].centroid.y)+';'+str(PDB.chains[c][r].centroid.z)+';'
								usedr.append(r)
				else:
					try:
						PDB.chains[k[4]][r].Centroid()
						line += str(PDB.chains[k[4]][r].centroid.x)+';'+str(PDB.chains[k[4]][r].centroid.y)+';'+str(PDB.chains[k[4]][r].centroid.z)+';'
					except:
						multiple = True		
						for ch in PDB.chains:
							if ch == k[4]:
								continue
							else:
								if r in PDB.chains[ch]:
									if r not in usedr:
										PDB.chains[ch][r].Centroid()
										line += str(PDB.chains[ch][r].centroid.x)+';'+str(PDB.chains[ch][r].centroid.y)+';'+str(PDB.chains[ch][r].centroid.z)+';'
										usedr.append(r)
			elif len(PDB.chains) == 1:
				multiple=False
				for c in PDB.chains:
					PDB.chains[c][r].Centroid()
					line += str(PDB.chains[c][r].centroid.x)+';'+str(PDB.chains[c][r].centroid.y)+';'+str(PDB.chains[c][r].centroid.z)+';'
			
			else:
				multiple = True		
				for ch in PDB.chains:
					if ch == k[4]:
						continue
					else:
						if r in PDB.chains[ch]:
							if r not in usedr:
								PDB.chains[ch][r].Centroid()
								line += str(PDB.chains[ch][r].centroid.x)+';'+str(PDB.chains[ch][r].centroid.y)+';'+str(PDB.chains[ch][r].centroid.z)+';'
								usedr.append(r)				
		line += '\n'
		fout.write(line)
		if multiple:
			print 'Multiple chains within PDB %s will be used. Only %s is reported to be used by Homstrad'%(k[:4],k[4])
			
	fout.close()
# End of definitions###############################################################################################
# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	## Command line input ###################################################
	prefix = sys.argv[1]
	for arg in sys.argv[1:]:
		if arg == '-txt':
			txt=True
		else:
			txt=False
		if arg == '-sup':
			sup = True
		else:
			sup = False
	#####################################################
	if sup:
		PDBs = get_matt_chain(prefix)
	else:
		PDBs = get_MATTassignedChains(prefix,txt)

	write_homologous_coordinates(prefix,PDBs)
	
	os.system('cp %s_GM.gm %s_scaled.gm'%(prefix,prefix))
	os.system('/home/jshleap/LabBlouin/code/GM/procrustes.py %s 3'%(prefix+'_scaled'))
	os.system('/home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(prefix+'_scaled'))
	os.system('cp %s_GM.gm %s_noscaled.gm '%(prefix,prefix))
	os.system('/home/jshleap/LabBlouin/code/GM/procrustes.py %s 3 -scale'%(prefix+'_noscaled'))
	os.system('/home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(prefix+'_noscaled'))