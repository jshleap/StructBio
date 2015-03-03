#!/usr/bin/python
'''
This script will go over the HOMSTRAD DATABASE, get te PDBs and a MATT alignment to be feed to
MATT2GM and then to Raw_centroids and GMtoRMSD
'''
#importing bit####################################################################################################
import os, sys
from RawCentroids import get_MATTassignedChains
from subprocess import Popen, PIPE
# End importing####################################################################################################

#Some definitions##################################################################################################
#list of characters in order for PDB
chainnames=[chr(i) for i in range(65,91)]
chainnames.extend([chr(i) for i in range(49,58)])
chainnames.extend('0')
chainnames.extend(['~','!','@','#','$','%','^','&','(',')','-','+','=','_','>','/',"'",'\\','?','"','[',']','{','}'])
chainnames.extend([chr(i) for i in range(97,123)])

def check_LIST(PATH_TO_LIST=os.getcwd()):
	'''
	will check if the LIST file is in the directory
	'''
	if not os.path.isfile(PATH_TO_LIST+'/LIST'):
		p1 = Popen('wget ftp://dalek.nibio.go.jp/homstrad/data/LIST.gz',shell=True)
		p1.wait()
		p2=Popen('gzip -d LIST.gz',shell=True)
		p2.wait()
		
def get_memfile_n_alignments(PATH_TO_LIST=os.getcwd()):
	'''
	Parse the List file from HOMSTRAD FTP (must be provided) and return a list of FTP addresses with the memfiles
	'''
	alnFTPs=[]
	supFTPs=[]
	ftps=[]
	direct = []
	f = open('LIST')
	for line in f:
		if line.endswith('.mem\n'):
			mem=line[line.find('./')+2:].strip()
			direct.append(mem[:mem.rfind('/')])
			ftps.append('ftp://dalek.nibio.go.jp/homstrad/data/'+mem)
		#elif line.endswith('sup.pdb\n'):
			sup = mem[:-4]+'-sup.pdb'
			supFTPs.append('ftp://dalek.nibio.go.jp/homstrad/data/'+sup)
		#elif line.endswith('.ali\n'):
			ali = mem[:-3]+'ali'#line[line.find('./')+2:].strip()
			alnFTPs.append('ftp://dalek.nibio.go.jp/homstrad/data/'+ ali)
		else:
			continue
			
	return ftps, supFTPs, alnFTPs, direct

def get_n_parse_memfile(ftp):
	'''
	Parser of the memfile and will return a dictionary relating the HOMSTRAD code with the actual PBD code
	'''
	#get memfile
	memfile = ftp[ftp.rfind('/')+1:]
	pathtomem=os.getcwd()+'/'+memfile
	if not os.path.isfile(pathtomem):
		os.system('wget %s'%(ftp))	
	pdbdict={}
	fastaname={}
	f = open(memfile).read().split('----')
	count = 0
	if len(f) >= 2:
		for e in f[0].split('\n'):
			if e == '':
				continue
			else:
				pdbname=e.strip()[:4]
				try:
					chain = e.strip()[4].upper()
				except:
					chain= 'A'
				fastaname[pdbname]=pdbname+chain+':'+chainnames[count]
				count+=1
			
		for el in f[1][f[1].find('TRAD\n')+5:].split('\n'):
			if not el == '':
				if '-' in el:
					continue
				else:
					bline = el.split()
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
				fastaname[pdbcode]=pdbcode+chain+':'+chainnames[count]
				count+= 1
			
	return pdbdict, fastaname	


def get_n_parse_fasta(ali,chainnames,fastanames):
	'''
	parser and getter of the sequence alignment file
	'''
	alignfile = ali[ali.rfind('/')+1:]
	pathtoali=os.getcwd()+'/'+alignfile
	if not os.path.isfile(pathtoali):
		os.system('wget %s'%(ali))
		
	f = open(alignfile.strip()).read().split('>')
	fout = open(ali[ali.rfind('/')+1:].strip()[:-4]+'-sup.fasta','w')
	for e in f[1:]:
		if e == '':
			continue
		else:
			lines = e.split('\n')
			pdb = lines[0].split(';')[1][:4]
			seq = "".join(lines[2:])[:-1]
			fout.write('>'+fastanames[pdb]+'\n'+seq+'\n')
	fout.close()
					

def get_pdb_FTP(pdbdict):
	'''
	Will use wget to get the PDB, unzip it and write a new file with the apropriate chain
	'''
	for k,v in pdbdict.iteritems():
		if not os.path.isfile(os.getcwd()+'/%s.pdb'%(v)):
			wget = Popen('wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz'%(v[:4]),shell=True)
			wget.wait()
			gzip = Popen('gzip -d pdb%s.ent.gz'%(v[:4]),shell=True)
			gzip.wait()
			if not '_' in k:
				os.system('grep " %s " pdb%s.ent | grep "ATOM   " > %s.pdb'%(v[4].upper(),v[:4],v))
			else:
				Popen('grep "ATOM   " > %s.pdb'%(v[4].upper(),v[:4],v),shell=True)
			Popen('rm pdb%s.ent'%(v[:4]),shell=True)
		else:
			if open(v+'.pdb').read() == '':
				Popen('rm %s.pdb'%(v),shell=True)
				wget = Popen('wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz'%(v[:4]),shell=True)
				wget.wait()
				gzip = Popen('gzip -d pdb%s.ent.gz'%(v[:4]),shell=True)
				gzip.wait()				
				Popen('grep "ATOM   " > %s.pdb'%(v[4].upper(),v[:4],v),shell=True)
				Popen('rm pdb%s.ent'%(v[:4]),shell=True)
				
def edit_MATTfasta(prefix):
	PDBs = get_MATTassignedChains(prefix)
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
	check_LIST()
	ftps, supFTPs, alnFTPs, dirs = get_memfile_n_alignments()
	for ftp in range(len(ftps)):
		pdbdict, fastanames = get_n_parse_memfile(ftps[ftp])
		if not os.path.isfile(os.getcwd()+supFTPs[ftp][supFTPs[ftp].rfind('/'):]):
			pdb = Popen('wget %s'%(supFTPs[ftp]), shell=True)
			pdb.wait()
		get_n_parse_fasta(alnFTPs[ftp],chainnames, fastanames)
		try:
			#edit_MATTfasta(dirs[ftp])
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s-sup'%(dirs[ftp]), shell=True)
			matt2gm.wait()
		except:
			matt2gm = Popen('python /home/jshleap/LabBlouin/code/MATT2GM/MATT2GM.py %s-sup'%(dirs[ftp]), shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(matt2gm.communicate()[1])
			fout.close()
		get_pdb_FTP(pdbdict)
		try:
			rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s-sup'%(dirs[ftp]), shell=True)
			rawcentroid.wait()
		except:
			rawcentroid = Popen('python /home/jshleap/LabBlouin/code/GM/RawCentroids.py %s-sup'%(dirs[ftp]), shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(rawcentroid.communicate()[1])			
			fout.close()
		try:
			gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s-sup -nograph'%(dirs[ftp]), shell=True)
			gm2rmsd.wait()
		except:
			gm2rmsd = Popen('python /home/jshleap/LabBlouin/code/alignment/GMtoRMSD.py %s-sup -nograph'%(dirs[ftp]), shell=True, stderr=PIPE)
			fout = open(dirs[ftp]+'.stderr','a')
			fout.write(gm2rmsd.communicate()[1])			
			fout.close()			
		os.mkdir(dirs[ftp])
		Popen('mv %s* *.pdb ./%s'%(dirs[ftp],dirs[ftp]),shell=True)
			
