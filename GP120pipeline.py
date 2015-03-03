#!/usr/bin/env python
'''
GP120pipeline.py
-------------------------
2014; Jose Sergio Hleap
-------------------------
Giving a set of gp120 structures (problem structures):
1) Align them to the best reference in the test set.
2) Classify them using SVM (classify GM)

This script requires PDBnet, Abermusa, GP120classifier and classifyGM to be importable and in 
the pythonpath.
'''
import optparse
from os.path import isfile, split,join
from os import getcwd, chdir, remove, mkdir
from shutil import copy, move
from subprocess import Popen
from glob import glob as G
from classifyGM import main as classifyGM
from GP120classifier import main as GP120classifier
from ABeRMuSA.ABeRMuSA import main as ABeRMuSA
from labblouin.logfile import XMLfile

class opti:
	''' bogus class for passing option to classify GM'''
	def __init__(self,options):
		self.predict        = 'Test.gm'
		self.featureselect  = True
		if options.sequence:		
			self.usesequence    = '.'
		else:
			self.usesequence    = None
		self.c_val          = options.c_val
		self.corr           = options.corr

class aberop:
	''' bogus class for passing ption to ABeRMuSA'''
	def __init__(self,options,args,ref=''):
		self.log       = True
		self.aligner   = 'matt'
		if not options.quick:
			self.reference = '%s/'%(args[1])+ref
		else:
			self.reference = None
			self.quick     = 6
		self.alphaC    = False
		self.executable= None
		self.prefix    = None
		self.scop      = None
		self.scopcache = None
		self.multi     = 0
		self.clean     = False
		self.split     = False
		self.tar       = False
		if options.optimize:
			self.optimize = True
		else:
			self.optimize = False
			
	
def align(pathtopdbs, pathtorefal,options, args):
	''' Use abermusa to align the problem set'''
	if isfile('alignment.xml'):
		xml = XMLfile('alignment.xml')
		xml.read()
		r = xml.tree.find('reference')
		ref = r.get("folder") + '.pdb'
		if isfile(r + '/ref.pickl'):
			remove(r + '/ref.pickl')		
	else:
		ref = None
	op = aberop(options,args,ref = ref)
	ABeRMuSA(op,[args[0],args[1]])
	'''
	exeline = 'ABeRMuSA.py -l'
	if options.optimize:
		exeline += ' -o '
	exeline += '-r %s/KC312379.3WARO_D16.2007.US.B.CCR5@.B99990001.pdb'%(pathtorefal)
	exeline+= ' %s %s'%(pathtorefal,pathtopdbs)
	abe=Popen(exeline,shell=True)
	abe.wait()'''

def readAlGM():
	data={}
	with open('alignment_final.gm') as A:
		for line in A:
			if line == '':
				continue
			else:
				bline = line.strip().strip(';').split(';')
				name = bline[0].strip('"').strip('>')
				name=name[:name.find(':')]
				data[name] = bline[1:]
	return data

def get_equi(filename):
	equi={}
	with open(filename) as F:
		for line in F:
			if line == '':
				continue
			else:
				bl = line.strip().split()
				equi[bl[0]]=bl[1]
	return equi

def split_aln(pathtorefal, pathtopdbs, datadict,equi):
	tout = open('Test.gm','w')
	dout = open('Alignment.gm','w')
	old = [x.strip('.pdb') for x in G(pathtorefal+'/*.pdb')]
	new = [y.strip('.pdb') for y in G(pathtopdbs+'/*.pdb')]
	for i in old:
		i = i.split('/')[-1]
		dout.write('>'+equi[i]+':'+i+';'+';'.join(datadict[i])+'\n')
	for j in new:
		j = j.split('/')[-1]
		tout.write('>'+equi[j]+':'+j+';'+';'.join(datadict[j])+'\n')
	tout.close()
	dout.close()

def translate_pred(predfile,equi):
	nfout = open(predfile.strip('prediction')+'orinames.prediction','w')
	with open(predfile) as P:
		for line in P:
			if line == '':
				continue
			else:
				bl = line.split()
				name = bl[0].strip().strip('"').strip('>')
				name = name[:name.find(':')]
				nfout('>'+equi[name]+'\t'+bl[1]+'\n')
	nfout.close()

def get_test_pdbs(testsetfile):
	'''
	Test only function to create a folder with the test pdbs. The testfile 
	is a tab delimited file with label, and pdbname. It assumes that the
	test PDBs are in the same folder as the path provided with the testfile
	'''
	path,fil = split(testsetfile)
	npath = join('.','ProblemSet')
	mkdir(npath)
	with open(testsetfile) as F:
		for line in F:
			bl=line.strip().split('\t')
			move(join(path,bl[1]),npath)
	return npath
	
def main(pathtopdbs, pathtorefal,options):
	''' execute code '''
	if options.knownTest:
		pathtopdbs = get_test_pdbs(options.knownTest)
	if not isfile('alignment_final.gm'):
		align(pathtopdbs, pathtorefal,options, [pathtopdbs, pathtorefal])
	data = readAlGM()
	for f in G('./_input/*.pdb'):copy(f,getcwd())
	GP120classifier('alignment_final',True)
	equi = get_equi('alignment_final.pdbequi')
	split_aln(pathtorefal, pathtopdbs, data,equi)
	copy('alignment_final_beforesampling.landmarks','Alignment.landmarks')
	copy('alignment_final.gm','Alignment.gm')
	copy('alignment_final.cls','Alignment.cls')
	ar = 'Alignment R5.gm R5 dual.gm dual X4.gm X4'
	op = opti(options)
	classifyGM(op,ar.split())
	#classify = Popen('classifyGM.py Alignment R5.gm R5 R5X4.gm dual X4.gm X4  -f -p Test.gm -s .',shell=True)
	#classify.wait()
	#translate_pred('Alignment.prediction',equi)
	[remove(x) for x in G('*.pdb')]

if __name__ == '__main__':
	opts = optparse.OptionParser(usage='%prog <path to problem pdbs> <path to aligned gp120 ref>'
	                             '> [options]\n The current working directory has to contain '\
	                             'the list of PDBs of each of the types to train with (R5,dual,'\
	                             ' and X4) in separate files with the extension <.list>. The pre'\
	                             'fix of those files will be the labels for the SVM classification.')

	opts.add_option('--optimize','-o', action = "store_true", default=False,
	                help="Optimize the alignment (take more time). Default: False)")
	opts.add_option('--c_val','-c', action = "store", default=None,
                    help="Provide a known optimum c value for the svm. This avoid the optimization step"\
                    " Default: No (Optimize) ")  	
	opts.add_option('--quick','-q', action = "store_true", default=False,
		                help="Instead of providing a reference structure, use the quick algorithm."\
		                " Default: No (Use best reference in a previous alignment)")	
	opts.add_option('--knownTest','-k', action = "store", default=None,
	                help="If a list of PDBs is provided here, those structures will be put into a problemset"\
	                " folder and use as testset. Default: No (This is mainly for testing purposes)")
	opts.add_option('--corr','-C', action = "store", default='',
	                help="Whether to use or not correlation of the samples giving the variables."\
	                " Default: No (Raw, use T otherwise)")	
	opts.add_option('--sequence','-s', action = "store_true", default=False,
		                help="Use sequence. Default: False)")	
	options, arg = opts.parse_args()

	main(arg[0],arg[1],options)    