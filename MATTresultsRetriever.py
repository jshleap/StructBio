#!/usr/bin/python
'''
This script will go over the MATT reults of the homstrad DATABASE, get te PDBs and a MATT alignment to be feed to
MATT2GM and then to Raw_centroids and GMtoRMSD
'''

#importing bit####################################################################################################
import os,sys
from RawCentroids import get_MATTassignedChains
from subprocess import Popen, PIPE
from walker import *
#from bs4 import BeautifulSoup
# End importing####################################################################################################
#Some definitions##################################################################################################
def check_LIST(PATH_TO_LIST=os.getcwd()):
	'''
	will check if the LIST file is in the directory
	'''
	if not os.path.isfile(PATH_TO_LIST+'/LIST'):
		p1 = Popen('wget ftp://dalek.nibio.go.jp/homstrad/data/LIST.gz',shell=True)
		p1.wait()
		p2=Popen('gzip -d LIST.gz',shell=True)
		p2.wait()
		
def parse_list(PATH=os.getcwd()):
	check_LIST()
	f = open(PATH+'/LIST')
	files=[]
	for line in f:
		if not line.endswith('.mem\n'):
			continue
		else:
			filename = line[line.rfind('/')+1:].strip()[:-4]
			files.append(filename)
	return files

def get_n_parse_memfile(fileprefix):
	'''
	Parser of the memfile and will return a dictionary relating the HOMSTRAD code with the actual PBD code
	'''
	#get memfile
	memfile = fileprefix+'.mem'
	pathtomem=os.getcwd()+'/'+memfile
	if not os.path.isfile(pathtomem):
		os.system('wget %s'%('ftp://dalek.nibio.go.jp/homstrad/data/'+fileprefix+'/'+memfile))	
	pdbdict={}
	f = open(memfile).read().split('----')
	if len(f) >= 2:
		bpdb=[]
		for e in f[0].split('\n'):
			if e == '':
				continue
			else:
				homstrad=e.strip()
				try:
					chain = e.strip()[4].upper()
				except:
					chain= '_'	
				PDB=homstrad[:4]+chain
				bpdb.append((homstrad,PDB))
			
		for el in f[1][f[1].find('TRAD\n')+5:].split('\n'):
			apdb=[]
			if not el == '':
				if '-' in el:
					continue
				else:
					bline = el.split()
					if len(bline[1].split(',')) != 1:
						ps=bline[1].split(',')
						for p in ps:
							apdb.append((p,bline[0]))
					else:
						apdb.append((bline[0],bline[1]))
		if len(bpdb) > len(apdb):
			for e in bpdb:
				pdbdict[e[0]]=e[1]
		else:
			for el in apdb:
				pdbdict[el[0]]=el[1]
				
	else:
		for el in f[0].split('\n'):
			if not el == '':
				bline=el.strip()
				pdbcode=bline[:4]
				homstradcode=bline
				try:
					chain=bline[4]
				except:
					chain='_'
				pdbdict[homstradcode]=pdbcode+chain
				#fastaname[pdbcode]=pdbcode+chain+':'+chainnames[count]
			
	return pdbdict


def get_matt_results(fileprefix):
	present = False
	if not os.path.isfile(os.getcwd()+'/'+fileprefix+'.pdb'):
		pdb = Popen('wget http://groups.csail.mit.edu/cb/matt/homstrad/%s.pdb'%(fileprefix), shell=True)
		pdb.wait()
	
	if not os.path.isfile(os.getcwd()+'/'+fileprefix+'.fasta'):
		fasta = Popen('wget http://groups.csail.mit.edu/cb/matt/homstrad/%s.fasta'%(fileprefix), shell=True)
		fasta.wait()
	
	if os.path.isfile(os.getcwd()+'/'+fileprefix+'.pdb'):
		present = True


	return present

def get_pdb_FTP(pdbdict,dataset):
	'''
	Will use wget to get the PDB, unzip it and write a new file with the apropriate chain
	'''
	for k in pdbdict.iterkeys():
		if not os.path.isfile(os.getcwd()+'/%s.pdb'%(k)):
			#wget = Popen('wget http://www.rcsb.org/pdb/files/%s.pdb'%(v[:4]),shell=True)
			wget = Popen('wget ftp://dalek.nibio.go.jp/homstrad/data/%s.atm'%(dataset+'/'+k),shell=True)
			wget.wait()
			#gzip = Popen('gzip -d pdb%s.ent.gz'%(v[:4]),shell=True)
			#gzip.wait()
			#if not '_' in k:
			#	os.system('grep " %s " pdb%s.ent | grep "ATOM   " > %s.pdb'%(v[4].upper(),v[:4],v))
			#else:
			#Popen('grep "ATOM   " > %s.pdb'%(v[:4]),shell=True)
			#
			#Popen('rm pdb%s.ent'%(v[:4]),shell=True)
		else:
			if open(k+'.pdb').read() == '':
				Popen('rm %s.atm %s.pdb'%(k,k),shell=True)
				wget = Popen('wget http://www.rcsb.org/pdb/files/%s.pdb'%(k[:4]),shell=True)
				wget.wait()
				#gzip = Popen('gzip -d pdb%s.ent.gz'%(v[:4]),shell=True)
				#gzip.wait()				
				#Popen('grep "ATOM   " > %s.pdb'%(v[4].upper(),v[:4],v),shell=True)
				#Popen('rm pdb%s.ent'%(v[:4]),shell=True)
		Popen('mv %s.atm %s.pdb'%(k,k),shell=True)
		
def try_pipeline(fileprefix, pdbdict, single):
	try:
		matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s'%(fileprefix), shell=True)
		matt2gm.wait()
	except:
		matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s'%(fileprefix), shell=True, stderr=PIPE)
		fout = open(fileprefix+'.stderr','a')
		fout.write(matt2gm.communicate()[1])
		fout.close()
	get_pdb_FTP(pdbdict,fileprefix)
	try:
		rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s -sup'%(fileprefix), shell=True)
		rawcentroid.wait()
	except:
		rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s -sup'%(fileprefix), shell=True, stderr=PIPE)
		fout = open(fileprefix+'.stderr','a')
		fout.write(rawcentroid.communicate()[1])			
		fout.close()
	try:
		gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(fileprefix), shell=True)
		gm2rmsd.wait()
	except:
		gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s -nograph'%(fileprefix), shell=True, stderr=PIPE)
		fout = open(fileprefix+'.stderr','a')
		fout.write(gm2rmsd.communicate()[1])			
		fout.close()
	if not single:
		os.mkdir(fileprefix)
		Popen('mv %s* *.pdb ./%s'%(fileprefix,fileprefix),shell=True)	

# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	single=False
	plot=False
	PATH=os.getcwd()
	for arg in sys.argv[1:]:
		if arg.startswith('-path='):
			PATH=arg[6:]
		elif '-single' in arg:
			single=True
		elif '-plot' in arg:
			plot = True
			
	files = parse_list(PATH)
	for f in files:
		present = get_matt_results(f)
		if present:
			pdbdict = get_n_parse_memfile(f)
			try_pipeline(f, pdbdict, single)
			
	if plot:
		gmns,gms,matt = walker()
		scatter_plots(gmns,gms,matt)