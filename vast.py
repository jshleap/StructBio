#!/usr/bin/python
''' 
This script will get a fasta file from a VAST search (http://www.ncbi.nlm.nih.gov/Structure/VAST/vast.shtml),
will get the pdb codes, store the sequences and fetch the PDB using the FTP server of RCSB
'''
# importing bit###############################################################################################
print 'Importing required packages:\n'
from subprocess import Popen, PIPE
import sys, os
# ############################################################################################################

# ############################################################################################################
# Some python functions#######################################################################################
# ############################################################################################################
def fasta_parser(prefix,initialpdb):
	'''
	parse a fasta file as downloaded from blast
	'''
	chains={}
	seqs={}
	fasta = open(prefix+".fasta").read()#get file to buffer
	#check if is fasta
	if not '>' in fasta:
		print 'Not a fasta file. Please provide a fasta file.'
		sys.exit()
	else:
		bfile = fasta.split('\n>gi|') #parse the file by the >gi

	#check if the first line have the pdb code in it
	if bfile[1].count('|') < 1 and not initialpdb:
		print 'There is no PDB code for the first entrance, and you did not provide it.',
		print 'Please, either edit the fasta incluing this PDB code, or provide the query PDB code in the arguments'
		sys.exit()
	else:
		for e in bfile:
			if e == '':
				continue
			else:
				bline = e.split('\n')
				if bfile.index(e) == 0:
					pdb, chain = bline[0].split('|')[3], bline[0].split('|')[4][0] #get PDB code and chain
				else:
					pdb, chain = bline[0].split('|')[2], bline[0].split('|')[3][0] #get PDB code and chain
				chains[pdb] = chain
				seqs[pdb] = ''.join(bline[1:])
	return chains,seqs

def get_PDB(chains):
	'''
	get the PDB files corresponding to the alignment and return a dictionary with the species
	'''
	sp = {}
	for k,v in chains.iteritems():
		if not os.path.isfile(os.getcwd()+'/%s%s.pdb'%(k,v)):
			wget = Popen('wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz'%(k.lower()),shell=True) #get the PDB
			wget.wait()
			gzip = Popen('gzip -d pdb%s.ent.gz'%(k.lower()),shell=True) # unzip it
			gzip.wait()
			sps = species(k.lower(),v.upper())
			sp[k]=sps
			parse_ent(k,v)
			#f=Popen('grep " %s " pdb%s.ent | grep "ATOM   " > %s%s.pdb'%(v.upper(),k.lower(),k,v),shell=True)
			#f.wait()
			if os.lstat('%s%s.pdb'%(k,v)).st_size == 0:
				print 'There is a problem with the file %s%s.pdb that you have to fix before align.'%(k,v)
				print 'File is empty. Go back to the corresponding ent file and correct the error.'
				fileproblem = True
			#Popen('rm pdb%s.ent'%(k.lower()),shell=True)			
		elif os.path.isfile(os.getcwd()+'/%s%s.pdb'%(k,v)) and open(k+v+'.pdb').read() == '':
			Popen('rm %s.pdb'%(k+v),shell=True)
			wget = Popen('wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz'%(k.lower()),shell=True)
			wget.wait()
			gzip = Popen('gzip -d pdb%s.ent.gz'%(k.lower()),shell=True) # unzip it
			gzip.wait()
			sps = species(k.lower(),v.upper())
			sp[k]=sps
			parse_ent(k,v)
			#f=Popen('grep " %s " pdb%s.ent | grep "ATOM   " > %s%s.pdb'%(v.upper(),k.lower(),k,v),shell=True)
			#f.wait()
			if os.lstat('%s%s.pdb'%(k,v)).st_size == 0:
				print 'There is a problem with the file %s%s.pdb that you have to fix before align.'%(k,v)
				print 'File is empty. Go back to the corresponding ent file and correct the error.'
				fileproblem = True
			#Popen('rm pdb%s.ent'%(k.lower()),shell=True)
		else:
			sps = species(k.lower(),v.upper())
			sp[k]=sps
			continue
	return sp

def species(pdbname,ch):
	'''
	retrieve the species name from an original pdb
	'''
	s = Popen("cat pdb%s.ent | grep 'ORGANISM_SCIENTIFIC'"%(pdbname), shell=True, stdout=PIPE, stderr=PIPE)
	spn = s.communicate()[0]
	spn = spn[spn.find(': ')+2:spn.find(';')]
	return spn


def get_list(directory):
	'''
	creates a list with the PDBs filenames
	'''
	lis = Popen('ls |grep ".pdb" > pdblist.txt',shell=True)
	lis.wait()

def create_fasta(prefix, seqs, sp, chains):
	'''
	create a fasta file with the PDBcodes and sps name and another one like matt format (is not an alignment)
	'''
	usedsp=[]
	mv = Popen('mv %s.fasta %s_original.fasta'%(prefix,prefix),shell=True) #back up the original alignment
	mv.wait()
	aln = open(prefix+'_temp_final.fasta','w')
	spec = open('splist.txt','w')
	fas = open(prefix+'.fasta','w')
	for ke,va in seqs.iteritems():
		fas.write('>'+ke+':'+sp[ke]+'\n'+va+'\n')
		spec.write(ke+'\t'+sp[ke]+'\n')
		if not sp[ke] in usedsp:
			aln.write('>'+sp[ke]+':'+ke+'\n'+va+'\n')
			usedsp.append(sp[ke])
		else:
			aln.write('>'+sp[ke]+'_'+str(sp.keys().index(ke))+':'+ke+'\n'+va+'\n')
	spec.close()
	fas.close()
	aln.close()
	
def Align(prefix, jobname, email, path2code=os.getcwd(), filelistname='pdblist.txt', quick=False):
	'''
	will run in fester the alignment
	'''
	print filelistname
	a = Popen('hostname',shell=True, stdout=PIPE, stderr=PIPE)
	host = a.communicate()[0].strip()
	print 'Host = %s'%(host)
	if host != 'fester.cs.dal.ca':
		print 'You are not on fester!!! run this piece of the code in fester'
		sys.exit()
	else:
		if quick:
			msa = Popen('python %s/MSA_pairwise_fester.py %s -name %s -L %s -email %s -quick'%(path2code,prefix,jobname,filelistname,email),shell=True)
			msa.wait()
		else:
			msa = Popen('python %s/MSA_pairwise_fester.py %s -name %s -L %s -email %s'%(path2code,prefix,jobname,filelistname,email),shell=True)

def parse_ent(pdbcode, chain):
	'''
	will parse an .ent file, extract the ATOM line of the specified chain
	'''
	o = open(pdbcode.upper()+chain+'.pdb','w')
	filename = 'pdb'+pdbcode.lower()+'.ent'
	f = open(filename)
	for line in f:
		if not line.startswith('ATOM'):
			continue
		elif not line[21:22] == chain:
			continue
		else:
			o.write(line)
	o.close()

#def get_pfam_structures(pfid):
	
# ###################################################################################################################
# Aplication of the code ############################################################################################

if __name__ == "__main__":
	if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'Usage: vast.py [prefix] [options]'
		print '\t-initialPDB=XXX : give the PDB code of the reference structure used to do the vast/blast search'\
		      'if not provided, then it is assumed that the fasta have been edited to contain that info (Default: No)'
		print '\t-dir=XXX/XXX : Give the directory where your files are (location) at if is not the current working directory (Default : <CWD>)'
		print '\t-align=W,X,Y,Z : if you are in fester, and you want to align the structures, use -align=quick,jobname,email,PATHtoCODE.'\
		      'the last argument is the path to the code MSA_pairwise_fester.py, without the last "/". The first argument is either'\
		      'True or False (case sensitive) for either to use the quick version of the alignment algorithm (Default : No)'
		print '\t-alphaC'
		#print '\t-fasta : Outputs a fasta file with the PDBcodes and sps name and another one like matt format (is not an alignment) (Default: No)'
		#print '\t-pfam=XXX : if the pfam id is known and there are structures available, provide the pfam id (XXX) (Default : No)'
		sys.exit()

	# Default Parameters ##########################################################################################
	prefix = sys.argv[1]
	initialpdb = None
	directory = os.getcwd()
	align = False
	quick=False
	fileproblem = False
	filelistname='pdblist.txt'
	#fasta = False
	#pfam = False
	# Get user input ##############################################################################################
	for arg in sys.argv[2:]:
		if arg.startswith('-initialPDB='):
			initialpdb = arg[12:]
		elif arg.startswith('-dir='):
			directory = arg[5:]
		elif arg.startswith('-align='):
			bline=arg[7:].split(',')
			q=bline[0]
			jobname = bline[1]
			email = bline[2]
			path2code = bline[3]
			align = True
			if q == 'True':
				quick=True
			'''
		elif arg.startswith('-pfam='):
			pfam = True
			pfid = arg[6:]
		elif arg == '-fasta':
			if pfam:
				print 'The pfam option does not allow you to output a fasta. Option disable.'
			else:
				fasta=True
			'''

	# Start #####################################################################################################
	#if pfam:
	#	'do'
	#else:
	chains , seqs = fasta_parser(prefix,initialpdb) ## Parse Fasta
	sp = get_PDB(chains) ## get the PDBs
	get_list(directory) ## list the PDB in directory
	#if fasta:
	create_fasta(prefix, seqs, sp, chains)
	if align:
		print 'Arguments for qsub the alignment:'
		print prefix,jobname,email,path2code
		Align(prefix, jobname, email, path2code, filelistname,quick)
	