#!/usr/bin/python
'''
This script will go over the HOMSTRAD DATABASE, get te PDBs and a MATT alignment to be feed to
MATT2GM and then to Raw_centroids and GMtoRMSD
'''
#importing bit####################################################################################################
import os
from RawCentroids import get_MATTassignedChains
from subprocess import Popen, PIPE
# End importing####################################################################################################

#Some definitions##################################################################################################

def get_memfile(PATH_TO_LIST=os.getcwd()):
	'''
	Parse the List file from HOMSTRAD FTP (must be provided) and return a list of FTP addresses with the memfiles
	'''
	ftps=[]
	direct = []
	f = open('LIST')
	for line in f:
		if not line.endswith('.mem\n'):
			continue
		else:
			mem=line[line.find('./')+2:].strip()
			direct.append(mem[:mem.rfind('/')])
			ftps.append('ftp://dalek.nibio.go.jp/homstrad/data/'+mem)
	return ftps, direct

def parse_memfile(memfile):
	'''
	Parser of the memfile and will return a dictionary relating the HOMSTRAD code with the actual PBD code
	'''
	#homstradcodes=[]
	pdbdict={}
	f = open(memfile).read().split('----')
	if len(f) >= 2:
		#for e in f[0].split('\n'):
		#	if not e == '':
		#		homstradcodes.append(e)
		for el in f[1][f[1].find('TRAD\n')+5:].split('\n'):
			if not el == '':
				if '-' in el:
					continue
				else:
					bline = el.split()
					if '_' in bline[0]:
						bline[0] = bline[0].replace('_','A')
					pdbdict[bline[1]]=bline[0]
	else:
		for el in f[0].split('\n'):
			if not el == '':
				bline=el.strip()
				pdbcode=bline[:4]
				try:
					chain=bline[4]
				except:
					chain='A'
				pdbdict[pdbcode]=pdbcode+chain
			
	return pdbdict	

def get_pdb_FTP(PDBcode_n_chain):
	'''
	Will use wget to get the PDB, unzip it and write a new file with the apropriate chain
	'''
	os.system('wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz'%(PDBcode_n_chain[:-1]))
	os.system('gzip -d pdb%s.ent.gz'%(PDBcode_n_chain[:-1]))
	#os.system('gzip -d %s.pdb.gz > %s.pdb'%(PDBcode_n_chain[:-1],PDBcode_n_chain[:-1]))
	os.system('grep " %s " pdb%s.ent | grep "ATOM   " > %s.pdb'%(PDBcode_n_chain[-1].upper(),PDBcode_n_chain[:-1],PDBcode_n_chain))
	os.system('rm pdb%s.ent'%(PDBcode_n_chain[:-1]))
	return PDBcode_n_chain+'.pdb'
	
def edit_MATTfasta(prefix):
	PDBs = get_MATTassignedChains(prefix,True)
	path=os.getcwd()
	for k,v in PDBs.iteritems():
		if os.path.isfile(path+'/temp.fasta'):
			os.system('cp %s.fasta %s_old.fasta'%(prefix, prefix))
			os.system('rm %s.fasta'%(prefix))
			os.system('mv temp.fasta %s.fasta'%(prefix))
		os.system('sed "/>%s/c \>%s:%s" %s.fasta > temp.fasta'%(k,k,v,prefix))
		os.system('mv temp.fasta %s.fasta'%(prefix))
	os.system('rm temp.fasta')
	
# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	ftps, dirs = get_memfile()
	for ftp in range(len(ftps)):
		os.system('wget %s'%(ftps[ftp]))
		pdbdict = parse_memfile(ftps[ftp][ftps[ftp].rfind('/')+1:])
		fout=open('list.txt','w')
		for k,v in pdbdict.iteritems():
			fout.write(get_pdb_FTP(v)+'\n')
		fout.close()
		try:
			# align using matt
			matt=Popen('Matt -o %s -L list.txt -t 16'%(dirs[ftp]), shell=True)
			matt.wait()
		except:
			matt=Popen('Matt -o %s -L list.txt -t 16'%(dirs[ftp]),shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(matt.communicate()[1])			
			fout.close()
			
		edit_MATTfasta(dirs[ftp])
		try:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s'%(dirs[ftp]), shell=True)
			matt2gm.wait()
		except:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s'%(dirs[ftp]), shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(matt2gm.communicate()[1])
			fout.close()
		try:
			rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s'%(dirs[ftp]), shell=True)
			rawcentroid.wait()
		except:
			rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s'%(dirs[ftp]), shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(rawcentroid.communicate()[1])			
			fout.close()
		try:
			gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(dirs[ftp]), shell=True)
			gm2rmsd.wait()
		except:
			gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(dirs[ftp]), shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(gm2rmsd.communicate()[1])			
			fout.close()			
		os.mkdir(dirs[ftp])
		Popen('mv %s* *.pdb ./%s'%(dirs[ftp],dirs[ftp]),shell=True)
		
