#!/usr/bin/python
"""
protGM Copyright (C) 2011  (Version 2 2012) Jose Sergio Hleap


Application of some Geometric morphometrics on protein structures. That applications include 
abstraction of structures as shapes, PCoA, CVA, form differnce, MANOVA, outlier landmark removal
among others.

It uses some R functions taken and modified from Morphometrics with R(Series: Use R!) of Julien Claude.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

Requires:
1) rpy2, available at http://rpy.sourceforge.net/rpy2.html
2) Utils folder, available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools
3) Biophyton
4) Numpy
5) R packages ape, vegan, shapes, outliers, and corpcor
"""
# importing bit###################################################################################
#print 'Importing required packages:\n'
import sys,os,re,time
#from copy import deepcopy
from subprocess import Popen
from rpy2.robjects import r
import cPickle as pickle
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast.Applications import NcbiblastpCommandline
from rpy2.robjects.conversion import py2ri
from numpy import mean, median, std, sqrt
#from PDBnet import PDBstructure
from labblouin.PDBnet import PDBstructure
from labblouin.SeqMask import InferSingleLettterCode
from collections import Counter
# #################################################################################################
# Some R functions#################################################################################
def R_functions():
	'''
	Some R functions to work with shapes R package...
	'''
	
	## Convert GMfile into 'shapes'-readable array
	## A is the read.table variable of a GMfile
	## s is the sample size
	## k is the number of landmarks
	## d is the number of dimensions (2D,2;3D,3)
	r('arr <- function(A,d){s<-dim(A)[1]; k<-(dim(A)[2]/d);'\
	  'Data <- as.numeric(A);arr<-array(Data,c(s, d, k)); arr<-aperm(arr,c(3,2,1)); arr}')
	###################################################################################################
	## Mean shape coordinates##
	r('mshape <- function(A){apply(A, c(1,2), mean)}')
	###################################################################################################
	## Interlandmark distances between configurations #
	## M1: First configuration matrix of k dimensions and p landmarks.
	## M2: Second configuration matrix of k dimensions and p landmarks.
	r('ild2 <- function(M1, M2){sqrt(apply((M1-M2)^2, 1, sum))}')
	###################################################################################################
	## Euclidean distance between mean shape and a given configuration##
	##returns the diagonal vector of euclidean distances between landmarks
	r('euc.dist <- function(A){ n <- dim(A)[3]; k <- dim(A)[1];	M<-matrix(NA,n,k); msh <- mshape(A);'\
	  'for (i in 1:n){ for (j in 1:k){ M[i,j] <- sqrt(sum((msh[j,1]-A[j,1,i])^2,(msh[j,2]-A[j,2,i])^2,(msh[j,3]-A[j,3,i])^2))}}; M}')
	##################################################################################################
	## full matrix of distances FM, and the code of the fm function that returns the vectorized form
	## fm of the upper ## diagonal entry of p(p-1)/2 elements.
	r('FM<-function(M){as.matrix(dist(M))}')
	r('fm<-function(M){mat<-FM(M); mat[col(mat)<row(mat)]}')
	##################################################################################################
	## Form difference matrix (FDM) between two configurations; FDM_M1/M2 = FM1/FM2
	r('FDM<- function(M1,M2){ FDM<-FM(M1)/FM(M2)}')
	##################################################################################################
	## Convert 'shapes'-readable array into GMfile 
	## A is an array Shapes-type
	## d is the dimensions
	r('A2GM<- function(A,d,rownames=NULL){m<-matrix(NA,dim(A)[3],dim(A)[1]*d) ;'\
	  ' for (i in 1:dim(A)[3]){ for (j in 1:d){ m[i,seq(j,dim(m)[2],d)]<-A[,j,i]}};'\
	  'row.names(m)<-rownames;as.data.frame(m)}')
	##################################################################################################
	## a function to process the CVA
	r('CVA <- function(A,cls,prefix){'\
		'cls<-t(cls); png(filename = paste(prefix,"CVA.png",sep=""), width = 25, height = 25, '\
	    'units = "cm",pointsize = 12, bg = "white",  res = 600);'\
		'out<-capture.output(shapes.cva(A,as.factor(cls),scale=FALSE));dev.off();'\
		'f<-which(grepl("Proportion of trace:",out))+4 ; l<-length(out); M <- out[f:l];'\
		'mat<-matrix(NA,nrow=length(M),ncol=3);'
		'for(e in seq(length(M))){s <- strsplit(M[e], " ") ; s<- unlist(s); s<-s[s != ""];'\
		    'mat[e,1]<- s[2];mat[e,2]<- s[3]};'\
		'nm<-cbind(cls,mat);'\
		'write.table(nm,file=paste(prefix,".cva",sep=""),sep=";",row.names=FALSE,'\
	    'col.names=FALSE)}')

# ###############################################fasta###################################################
# Some python functions#############################################################################

def find_key(dic, val):
	"""return the key of dictionary dic given the value"""
	return [k for k, v in dic.iteritems() if val in v][0]

class GM:
	'''
	This class manipulate the gm file created with MATT2GM and performs a 3DScatterplot, PCoA, 
	perMANOVA, CVA, Form Difference analysis (FDA), clean the outliers from the FDA, and make a PCoA 
	of the transpose FD matrix.
	The FDM theory and functions are taken and modified from Morphometrics with R(Series: Use R!) of 
	Julien Claude (2008).
	'''	
	if __name__ == "__main__":
		'''imports... this is here to b able to import fasta class without importing all this R 
		packages'''	
		r('library(ape)')
		r('library(vegan)')
		r('library(shapes)')
		r('library(outliers)')
		r('library(corpcor)')
		
	# Load Data ###################################################################################
	def __init__(self, prefix, multiple, datatype, dim, cazy, nosingles,chain):
		self.transform_chain(chain)
		self.dim=dim
		self.prefix = prefix
		self.clasifier=None
		if multiple:
			self.currentPDB=PDBstructure(chain+'.pdb')
		else:
			self.currentPDB=PDBstructure(prefix+'.pdb')
		namelis=self.format_gm(prefix)
		self.check_fasta(namelis)
		self.F = fasta(prefix+'.fasta')
		self.clas=[]
		if cazy:
			cazy2class(self.prefix,self.F)
		if datatype == 'md':
			self.fasta4md(prefix,namelis)
		self.names, self.namestr = self.get_names(prefix,namelis,datatype)
		self.Get_classifiers()
		self.setUpR_variables()
		if nosingles:
			self.no_singletons()
			
	def transform_chain(self,chain):
		'''
		just a work around for different cases of chain names... It always assume that the fist 4 character
		correspond to the pdb code and the 5th or 6th (depending if is splitted or not) to the respective
		chain
		i.e. if chain provided is 1PPI.A.S._scrofa, will take 1PPI as the pdb code and A as chain. Likewise
		if is 1PPI_A-c.
		'''
		#cleaning the -c
		if '-c' in chain:
			self.chain=chain.replace('-c','')
			if chain[4].isalpha() and len(chain) >= 5:
				self.chain=chain[:5]
			elif not chain[4].isalpha() and chain[5].isalpha():
				self.chain=chain[:4]+chain[5]
		else:
			self.chain=chain

	def setUpR_variables(self):
			self.gm = r('D <- as.matrix(read.table("%s.r.gm",sep=";",row.names=1))'%(self.prefix))
			self.tm = r('t(D)')
			# set an array of k landmarks, d dimensions and n samples ############
			self.arr = r('A<-arr(D,%d)'%(self.dim))
			#modified for gp120... please get the other back when finished!!!!
			#self.arr = r('A<-arr(D[,-c(261, 262, 263, 264, 265, 266, 267, 268, 269 ,270, 271, 272, 273, 274, 275 ,276)],%d)'%(self.dim))
			self.n = r('dim(A)[3]')[0]
			self.d = r('dim(A)[2]')[0]
			self.k = r('dim(A)[1]')[0]
			# get the mean shape
			self.msh = r('msh<-mshape(A)')
			
	def Get_classifiers(self):
		fils= os.listdir('./')
		if self.prefix+'.cls' in fils:
			cls = set((open("%s.cls"%(self.prefix))).read().strip().split(';'))
			self.clasifier = r('cls <- as.matrix(read.table("%s.cls",sep=";", header = FALSE,'\
		                       'strip.white=TRUE))'%(prefix))
			# check if last entry is a blank
			t = r('is.na(cls[1,dim(cls)[2]])')[0]
			if t:
				self.clasifier = r('cls <- as.matrix(cls[1,1:dim(cls)[2]-1])')
			for x in cls:
				if x == '' or x == '\n':
					continue
				else:
					self.clas.append(x)
			ncls=len(self.clas)
						
	def check_fasta(self,namelis):
		if not os.path.isfile(self.prefix+'.fasta'):
			f=open(self.prefix+'.fasta','w')
			for n in namelis:
				seq=make_fasta(n+'.pdb')
				f.write(seq+'\n')
			f.close()		
			
	def no_singletons(self):
		self.gm
		self.clasifier
		f=open(self.prefix+'.cls').read().strip().split(';')
		s='c('
		c=Counter(f)
		cl=[]
		#co=-1
		for e in range(len(f)):
			#co+=1
			if e == '':
				continue
			elif c[f[e]] <= 1:
				continue
			else:
				try:
					int(f[e])
					s+=str(e+1)+','
					cl.append(f[e])
				except:
					continue
		s=s[:-1]+')'
		self.gm= r('D<-D[%s,]'%(s))
		self.clasifier = r('cls<-cls[%s]'%(s))
		self.prefix = self.prefix+'NS'
		self.arr = r('A<-A[,,%s]'%(s))
		self.n = r('dim(A)[3]')[0]
		self.d = r('dim(A)[2]')[0]
		self.k = r('dim(A)[1]')[0]
		self.clas =list(set(cl))
		
	def format_gm(self,prefix):
		'''
		Check for ; at the end and remove it, also remove the names
		'''
		namelis=[]
		self.openfile = open(prefix+'.gm')
		infile=self.openfile
		outfile=open(prefix+'.r.gm','w')
		for line in infile:
			if line == '\n' or line == '':
				continue
			bname=line[1:line.find(':'):]
			aname=line[line.find(':')+1:line.find(';')]
			if len(bname) > len(aname):
				name=bname
			else:
				name = aname
			namelis.append(name)
			if line.strip()[-1] == ';':
				outfile.write(line[:-2]+'\n')
			else:
				outfile.write(line)
		infile.close()
		outfile.close()

		return namelis
	
	def fasta4md(self,prefix,namelis):
		'''
		when a md is being used, will replicate the fasta to resemble an alignment
		'''
		fas = fasta(prefix+'.fasta')
		pdb = fas.aligned.keys()[0]
		seq = fas.seqs[pdb]
		cmd = Popen('mv %s %s'%(prefix+'.fasta',prefix+'_old.fasta'),shell=True)
		cmd.wait()
		new=open(prefix+'.fasta','w')
		for i in namelis:
			new.write('>%s:%s\n%s\n'%(pdb,str(i),seq))
		new.close()

	def get_names(self,prefix,namelis,datatype):
		'''
		Get the names from the fasta file in a dictionary relating the chainname assign by matt
		and the filename
		'''
		names=self.F.chains
		if len(namelis) != len(names):
			print 'WARNING!!! Some of the sequences were taken off because where duplicates.'\
				  ' Structures will not be taken into account.'
		namestr='c('
		for n in namelis:
			if n[:4] not in names:
				continue
			if 'FDclean' in prefix:
				namestr += n +','
			elif datatype == 'md':
				namestr += '"'+names[n]+'",'
			else:
				namestr += '"'+names[n[:4]]+'",'
		namestr = namestr[:-1]+')'

		return names, namestr

	def scatter3D(self):
		# Get each dimension
		r('x <- seq(1,dim(D)[2],%d)'%(self.d))
		r('X <- D[,x]')
		r('Vx<-as.vector(X)')
		r('y <- seq(2,dim(D)[2],%d)'%(self.d))
		r('Y <- D[,y]')
		r('Vy<-as.vector(Y)')
		r('z <- seq(3,dim(D)[2],%d)'%(self.d))
		r('Z <- D[,z]')
		r('Vz<-as.vector(Z)')

		# plot raw coordinates #####################################################################
		r('png(filename = "%s.png", width = 25, height = 25, units = "cm", pointsize = 12,'\
		  ' bg = "white",  res = 600)'%(self.prefix))
		r('scatterplot3d(Vx, Vy, Vz, col.axis="blue", col.grid="lightblue",'\
		  'main="Alpha-Amylase align coordinates", pch=20)')
		r('dev.off()')

	def PCoA(self,name,matrix,is_dist,namestr):
		dist=is_dist
		st=''
		for e in matrix:
			st+=str(e)+','
		r('mat<-matrix(c(%s),nrow=%d, ncol=%d)'%(st[:-1],matrix.dim[0],matrix.dim[1]))
		r('mite.mat <- vegdist(mat, "euclidean")')#,diag=TRUE,upper=TRUE)')
		'''positive = r('is.positive.definite(mite.mat)')
		if not positive[0]:
			print "Distance matrix is not positive definite. WARNING!!! using R's"\
			"make.positive.definite function."
			r('mite.mat<-make.positive.definite(mite.mat)')
			p = r('is.positive.definite(mite.mat)')
			print "Worked? : %s"%(p[0])'''
		if dist == False:
			try:
				r('res <- pcoa(as.dist(mite.mat), correction="cailliez", rn=%s)'%(namestr))
			except:
				# Not working... exeptions in R does not work with this try except
				print 'Using correction Lingoes after Cailliez failed'
				r('res <- pcoa(as.dist(mite.mat), correction="lingoes")')
		elif dist == True:
			try:
				r('res <- pcoa(mat, correction="cailliez")')
			except:
				# Not working... exeptions in R does not work with this try except
				print 'Using correction Lingoes after Cailliez failed'
				r('res <- pcoa(mat, correction="lingoes")')
		if spc:
			name=name+'_sub'
			r('v<-res$vectors')
			r('m=v[,1:2,drop=FALSE]')
			r('w=which(m[,1] >%f & m[,1] < %f & m[,2] > %f & m[,2] < %f)'%(x1,x2,y1,y2))
			#re-do PCoA with subset
			r('nmat<- mat[w,,drop=FALSE]')
			r('nmat<-vegdist(nmat, "euclidean")')
			r('write.table(D[w,,drop=FALSE], file = "%s.gm",sep = ";",col.names = FALSE,'
			  'row.names = TRUE)'%(name))
			r('D<-as.matrix(read.table("%s.gm",sep=";",row.names=1))'%(name))
			fnames=r('lnames<-list(rownames(D))')[0]
			nfasta=open(name+'.fasta','w')
			for e in fnames:
				nfasta.write('>'+e[-5:]+':'+e[-5:-1]+'\n'+self.F.seqs[e[-5:]]+'\n')
			nfasta.close()
			try:
				r('res <- pcoa(as.dist(nmat), correction="cailliez", rn=%s)'%(namestr))
			except:
				r('res <- pcoa(as.dist(nmat), correction="cailliez")')
		
		eigval=r('eival <- res$values')
		r('write.table(eival, file = "%s_vals.temp",sep = "\t",col.names = FALSE,'\
		  'row.names = FALSE)'%(name))
		eigvec=r('eivec <- res$vectors')
		r('write.table(eivec, file = "%s_vecs.temp",sep = "\t",col.names = FALSE,'\
		  'row.names = FALSE)'%(name))
		trace=r('trace<-res$trace')
		trace=trace[0]
		foutpcoa=open(self.prefix+'.pcoa','w')
		vatemp=open('%s_vals.temp'%(name))
		vetemp=open('%s_vecs.temp'%(name))
		axnum=len(vetemp.readline().split('\t'))
		foutpcoa.write('Values\nEigenvalues\t Relative_eig \tBroken_stick \t'
		               'Cumul_eig \tCumul_br_stick\n')
		for line in vatemp:
			foutpcoa.write(line)
		foutpcoa.write('\n')
		axi=''
		for ax in range(axnum):
			axi+= '##### Axis.%s ######\t '%(str(ax+1))
		foutpcoa.write('\n\nVectors\n'+axi+'\n')		
		for line in vetemp:
			foutpcoa.write(line)
		#os.system('cat vals.temp vecs.temp trace.temp > %s'%(self.prefix+'.pcoa'))
		r('png(filename = "%sPCoA.png", width = 35, height = 25, units = "cm", pointsize = 10,'\
		  'bg = "white",  res = 600)'%(name))
		r('biplot(res, main="Principle coordinate analysis of the %s dataset", pch=20)'%(self.prefix))
		r('dev.off()')
		foutpcoa.write('\n')
		foutpcoa.write('\nTrace = '+ str(trace))
		foutpcoa.close()
		vetemp.close()
		vatemp.close()
		#os.system('rm *.temp')

	def perMANOVA(self):
		self.clasifier
		clas=self.clas
		foutmanova=open(self.prefix+'.manova','w')
		'''r('ad<-adonis(D ~ factor(cls), permutations=%d, method="euclidean",'\
		  'contr.unordered = "contr.helmert",contr.ordered = "contr.helmert")'%(per))
		man=r('ad[1]')[0]
		foutmanova.write(str(man))
		foutmanova.close()'''
		tup=[]
		for i in range(len(clas)):
			for b in range(len(clas)):
				
				if i == b:
					continue
				elif (i,b) in tup or (b,i) in tup:
					continue
				else:
					line='Permutational MANOVA of the classifier "%s" vs All'%(clas[i])
					line2='='*(len(line)+1)
					foutmanova.write(line2+'\n'+line+'\t=\n'+line2+'\n')					
					print 'Comparing %s vs %s'%(clas[i],clas[b])
					if clas[i].isalnum() or clas[b].isalnum() :
						r('i <- which(cls == "%s")'%(clas[i]))
						r('b <- which(cls == "%s")'%(clas[b]))
					else:
						r('i <- which(as.numeric(cls)== %s)'%(clas[i]))
						r('b <- which(as.numeric(cls)== %s)'%(clas[b]))
					r('u <- c(i,b)')
					r('d <- D[u,]')
					r('BetweenGroups <- factor(cls[u])')
					tup.append((i,b))
					r('ad<-adonis(d ~ BetweenGroups, permutations=%d,method="euclidean")'%(per))
					man=r('ad[1]')
					man=man[0]
					# Write MANOVA table
					foutmanova.write('Group %s vs group %s\n'%(clas[i],clas[b])+str(man)+'\n\n')
		foutmanova.close()

	def CVA(self):
		d = int(r('length(as.numeric(cls[]))')[0])
		if d != len(self.names) and not nosingles:
			foutcva=open(self.prefix+'.cva','w')
			for e in range(len(clas)):
				foutcva.write('Canonical Variate analysis of the classifier "%s"\n\n'%(clas[e]))
				#get the classifiers
				r('clas <- as.factor(cls)')
				#perform CVA
				r('png(filename = "%sCVA%s.png", width = 28, height = 22, units = "cm",'\
				  'pointsize = 12, bg = "white",  res = 600)'%(self.prefix,clas[e]))
				r('out<-capture.output(shapes.cva(A,clas,scale=FALSE))')
				out=r('out')
				r('dev.off()')
				for o in out:
					foutcva.write(o+'\n')
				foutcva.write('\n\n')
			foutcva.close()
		else:
			r('prefix<-"%s"'%(prefix))
			r('CVA(A,cls,prefix)')
			'''
			r('png(filename = "%sCVA.png", width = 25, height = 25, units = "cm", pointsize = 12,'\
			  'bg = "white",  res = 600)'%(self.prefix))
			r('c<-capture.output(out<-shapes.cva(A,cls,scale=FALSE,pcaoutput=FALSE), file = "%s.cvaoutput")'%(self.prefix))
			r('dev.off()')
			r('nm<-cbind(cls,out)')
			r('write.table(nm,file="%s",sep=";",row.names=FALSE)'%(self.prefix+'.cva'))'''

	def formD(self, multiple=False):
		'''
		Calculate the FD and return the list of FD for each landmak
		'''
		chain=self.chain
		betares=[]
		name=self.prefix
		su=''
		count=0
		if multiple:
			cha=chain[-1]
			chain=chain[:4]
		# Iterate over samples to calculate the FDM from each one to the mean 
		for a in range(self.n):
			count+=1
			r('F%s<-FDM(A[,,%d],msh)'%(str(a),(a+1)))
			r('F%s1<-abs(F%s-median(F%s, na.rm=T))'%(str(a),str(a),str(a)))
			r('rownames(F%s1)<-1:%d'%(str(a),self.k))
			r('R%s<-round(apply(F%s1,2,sum,na.rm=T),2)'%(str(a), str(a)))
			if count < self.n:
				su+= 'R'+(str(a))+'+'
			if count == self.n:
				su+= 'R'+(str(a))
		# Estimate the average FD in all conformations
		lis=r('LIS<- round(((%s)/%d),3)'%(su,self.n))
		# Read the .landmarks file and the PDB file
		newPDB=PDBstructure() #deepcopy(self.currentPDB)
		land=open(name+'.landmarks')
		# Open the new pdb file
		#nPDB=open(name+'%s.FD.pdb'%(chain),'w')
		# Chains equivalencies 
		chains = {}
		aa={}
		for line in land:
			if line == '' or line == '\n':
				continue
			else:
					# New chain
				if line.startswith('>'):
					ch = line[1:5].strip()
					chains[ch] = []
					aa[ch]=[]
				else:
					# Add to chain
					lin = line.split()
					chains[ch].append(int(lin[1]))
					aa[ch].append(line[line.find('\t')+1:])
		if not chain in chains:
			print "Chain", chain, "not in PDB."
			sys.exit(-1)

		# write the new FD PDB file
		if multiple:
			oldchain=chain
			chain=cha
		for res in self.currentPDB.chains[chain]:
			if multiple: chainmap=chains[oldchain]
			else: chainmap = chains[chain]
			index=int(res.index)
			if not index in chainmap:
				newchain = 'X'
				temp = float(median(lis))
			else:
				newchain = chain
				temp = lis[chainmap.index(index)]
			if not newchain in newPDB.chains: newPDB.NewChain(newchain)
			newPDB.AddResidueToChain(newchain,res)
			for atom in newPDB.chains[newchain][res.index].atoms:
				newPDB.chains[newchain][res.index].atoms[atom].tempFactor=temp
			'''
			for atom in res.atoms:
				res.atoms[atom].tempFactor=temp'''
			line=(str(res.index),str(res.name),str(newchain), float(temp))
			if not line in betares:
					betares.append(line)
		betares = sorted(betares, key=lambda v: v[3])
		if multiple:
			oldch=chain
			chain=oldchain
		newPDB.WriteFile(name+'_%s.FD.pdb'%(chain))
		#nPDB.close()
		#PDB.close()
		land.close()
		self.F.equivalences(prefix)
		b = open(name+'.betaeq','w')
		for e in betares:
			e = [str(x) for x in e]
			b.write('\t'.join(e)+'\n')
			'''
			for t in e:
				b.write(str(t)+'\t')
			b.write('\n')'''
		return lis, chains, aa

	def cleanFD(self,chain='A', multiple=False):
		'''
		write a new gm file with filtered landmarks (those outliers tested with Grubbs test)
		http://en.wikipedia.org/wiki/Grubbs'_test_for_outliers
		'''
		lis, chains, aa = self.formD(chain,multiple)
		lis=list(lis)
		mlis=mean(lis)
		sdlis=std(lis)
		di={}
		dic=[]
		erased=[]
		if multiple:
			cha=chain[-1]
			chain=chain[:4]
		for i in range(len(lis)):
			di[i+1] = lis[i]
		for l in lis:
			#test statistically prescence of outliers
			r('gr<-grubbs.test(LIS)')
			sig=r('gr[3]')
			sig=sig[0][0]
			#remove values if outlier
			if sig<0.05:
				r('LIS<-rm.outlier(LIS)')
				outl=r('gr[2]')
				outl=float(outl[0][0].split()[2])
				erased.append(lis.index(l))
				di.pop(find_key(di,outl))
		for ke in di.iterkeys():
			#cou=+1
			z=ke*3
			y=z-1
			x=z-2
			dic.append(x)
			dic.append(y)
			dic.append(z)
		rs='c('+ str(dic)[1:-1] + ')'
		r('s<-%s'%(rs))
		r('nD<-D[,s]')
		r('write.table(nD, file = "%s_FDclean.gm",sep = ";",col.names = FALSE)'%(self.prefix))
		# New landmarks file
		nland=open(self.prefix+'_FDclean.landmarks','w')
		nfasta=open(self.prefix+'_FDclean.fasta','w')
		newchains={}
		sfasta=''

		#create a new chain dictionary whiout the outliers
		for k,v in chains.iteritems():
			tlist=[]
			for ele in v:
				if v.index(ele) in erased:
					tlist.append(ele)
			newchains[k]=tlist
		# prepare the new fasta file
		tempseqlist=self.F.seqslist.copy()
		for st,res in tempseqlist.iteritems():
			name=st
			chain=find_key(self.F.chains,st)
			nfasta.write('>'+name+':'+chain+'\n')
			for er in newchains[find_key(self.F.chains,st)]:
				res.pop(self.F.aligned[st][er])
			for element in res:
				nfasta.write(element)
			
			nfasta.write('\n')
			
		coun=0	
		for k, v in aa.iteritems():
			nland.write('>'+k+'\n')
			for element in range(len(v)):
				try:
					if int(v[element].split('\t')[0]) in newchains[k]:
						nland.write(str(coun)+'\t'+v[element])
						coun+=1
				except:
					continue

		nland.close()

	def FDmatrix(self,name):	
		'''
		Calculate the FD for each structure and the mean and return the transposed matrix
		of all structres
		'''
		name=name
		self.gm
		self.arr
		self.msh
		n,d,k= self.n, self.d, self.k		
		#create an empty matrix of the proper dimensions
		r('M<-matrix(NA,nrow=%d,ncol=%d)'%(n,k))
		# Iterate over samples to calculate the FDM from each one to the mean 
		for a in range(n):
			r('F%s<-FDM(A[,,%d],msh)'%(str(a),(a+1)))
			r('F%s1<-abs(F%s-median(F%s, na.rm=T))'%(str(a),str(a),str(a)))
			r('rownames(F%s1)<-1:%d'%(str(a),k))
			r('R%s<-round(apply(F%s1,2,sum,na.rm=T),2)'%(str(a), str(a)))
			r('M[%d,]<-R%s'%((a+1),str(a)))
			r('write.table(M, file = "%s_FDmatrix.M",sep = ";",col.names = FALSE)'%(self.prefix))
			#r('write.table(t(M), file = "%s_FDmatrix.tM",sep = ";",col.names = FALSE)'%(self.prefix))
		tM=r('t(M)')
		M=r('M')
		return name, M, tM

class fasta:
	def __init__(self,filein=''):
		print 'processing Fasta file.'
		self.chains={}
		self.seqs={}
		self.seqslist={}
		self.gapless={}
		self.aligned={}
		self.n=[]

		if filein:
			self.ProcessFasta(filein)
		
	def ProcessFasta(self,filename):
		try:
			fin = open(filename).read().split('\n>')
		except:
			print "file not found:", filename

		for e in fin:
			if e == '':
				continue
			elif e[0] == '>':
				e= e[1:]
			e = e.strip()
			chain=e[e.find(':')+1:e.find('\n')]
			se=re.search(':.+\n',e)
			if se:
				e=e.replace(se.group(0),'\t').replace('\n','').replace('\t','\n')
				name=e[:e.find('\n')]
			else:
				name=e[:e.find('\n')]
				e=e[e.find('\n')+1:].replace('\n','').replace('\t','\n')
			self.n.append(name)
			seq=e[e.find('\n')+1:]
			seqlist=list(seq)
			self.seqs[name]=seq
			self.seqslist[name]=seqlist
			self.chains[chain]=name
			self.Gapless()
			self.aligned_position()
		#return seqs, seqslist

	def Gapless(self):
		for k,v in self.seqs.iteritems():
			self.gapless[k]=list(v.replace('-',''))
	
	def aligned_position(self):
		for k,v in self.seqslist.iteritems():
			self.aligned[k]=[]
			itera=-1
			for e in v:
				itera+=1
				if e == '-':
					continue
				else:
					self.aligned[k].append(itera)
					
	def equivalences(self,prefix):
		f = open(prefix+'.equival','w')
		for k,v in self.aligned.iteritems():
			f.write('>'+k+':'+find_key(self.chains,k)+'\n')
			for e in v:
				f.write('(%d,%d)\t'%(int(v.index(e)),int(e)))
			f.write('\n')
			
def cazy2class(prefix,F,remote=False):
	''' 
	will take the cazy database (dictionary provided) and try to fetch subfamilies and place them
	as classifiers.
	'''
	print 'You chose to use CAZY database to classify GH13 family into subfamilies'\
          ' this will take a while, since have to go over BLAST results, etc..'
	cls=open(prefix+'.cls','w')
	# import database
	db=pickle.load(open('CazyDB.bin'))
	names = get_names(prefix+'.gm')
	for n in names:
		print 'Processing %s...'%(n)
		if remote:
			Entrez.email = 'jshleap@dal.ca'
			print '\tBlasting (Running remotely)...'
			n=n[:-1]+'_'+n[-1]
			while 1:
				try:
					b=qblast('blastp','nr',n,perc_ident=90,expect=1,gapcosts='11 1')
					print '\tBlast Done... \n\t\tAnalysing entries...'
					break
				except:
					print 'Some problem in NCBI, sleeping for 10...'
					time.sleep(10)
		else:
			print '\tBlasting (Running locally)...'
			fi=open('temp.fasta','w')
			fi.write('>%s\n%s'%(n,F.seqs[F.chains[n[:4]]]))
			fi.close()
			#blastp_cline = NcbiblastpCommandline(query="temp.fasta", db="nr", evalue=0.0001,
			#                                     outfmt=5, out="temp.xml",max_target_seqs=50, 
			#                                     num_alignments=50,num_threads=4)
			bl=Popen('blastp -db nr -outfmt "5" -query temp.fasta  -evalue 0.0001 -max_target_seqs 50 '\
			         '-seg yes -num_threads 4  -gapopen 10 -gapextend 1 -matrix BLOSUM90 -out temp.xml',
			         shell=True)
			bl.wait()
			print '\tBlast Done... \n\t\tAnalysing entries...'
			b=open('temp.xml')
		blast_record = NCBIXML.read(b)
		rm = Popen('rm temp.*',shell=True)
		rm.wait()
		nohit=True
		while nohit:
			for a in blast_record.alignments:
				print '\t\t\t%s'%(a.accession)
				h=a.hsps[0]
				if float(h.identities)/float(h.align_length) >=0.9:
					ans,k = dict_lookup(a.accession,db)
					if ans:
						cls.write(str(db[k])+';')
						print '\t\t\t\tAccession number found in CAZY!, Subfamily %s'%(db[k])
						nohit=False
						break
					else:
						if blast_record.alignments.index(a)+1 == len(blast_record.alignments):
							cls.write('%s;'%(n))
							nohit=False
							print '\tNo relative found in CAZY'
							break
				elif blast_record.alignments.index(a)+1 == len(blast_record.alignments):
					cls.write('%s;'%(n))
					nohit=False
					print '\tNo relative found in CAZY'
					break
	cls.write('\n')
	cls.close()

def dict_lookup(term,d):
	''' look if a dictionary keys contains a term'''
	answ=False
	l=None
	for k in d.iterkeys():
		if term in k:
			answ=True
			l=k
	return answ, l

def get_names(filename):
	''' get names of the proteins, assuming that are only PDBs'''
	f=open(filename)
	n=[]
	for line in f:
		if not line.startswith('>'):
			continue
		else:
			bline = line.split(';')
			h=bline[0].split(':')
			n.append(h[1])
	return n
				 
def make_fasta(filename):
	fileprefix=filename[:-4]
	if '-c' in fileprefix:
		chain=fileprefix[:fileprefix.find('-c')][-1]
	else:
		chain=fileprefix[-1]
	if "_" in fileprefix:
		name=fileprefix[:fileprefix.find('_')]
	else:
		name=filename[:-5]
	s=PDBstructure(filename)
	seq='>%s:%s\n'%(fileprefix,name)
	seq+=s.ChainAsFASTA(chain)
	return seq

def indexall(lst, value):
	return [i for i, v in enumerate(lst) if v == value]


	
# ##################################################################################################
# Aplication of the code ###########################################################################

if __name__ == "__main__":
	if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'Usage: protGM.py [prefix] [options]'
		print 'Options:'
		print '\t-nosingles : if a classifier file is provided with this option, the classifier and the'\
		      'matrix will be cropped avoiding the singletons (i.e. not classified ones or groups with a '\
		      'single member). For this, the not clasified ones must not be numeric, and the clasifications'\
		      ' numeric ( Default : No)'
		print '\t-md : will change the datatype to accept molecular dynamic simulation data (Default : No)'
		print '\t-scatter3d : Plot a 3D visualization of the landmarks ( Default : No )'
		print '\t-PCoA : Perform a Principal Coordinates Analysis and plot it. . If passed with'\
		      ' multiple, a PDB should be provided (no extension) after an equal sign( Default : No )'
		print '\t-subPCoA=X1,X2,Y1,Y2 : will perform a PCoA in a subset of a preceding PCoA. X1,X2,Y1,Y2 '\
		      'correspond to the corners of a square surrounding the points to be reanalyzed.'\
		      '( Default : No )'
		print '\t-perMANOVA=N : Perform a permutational MANOVA with the classifiers as factors. '\
		      'N permutations( Default : No )'
		print '\t-CVA : Perform a Canonical variates analysis for each classifier. If passed with'\
		      ' multiple, a PDB should be provided (no extension) after an equal sign ( Default : No )'
		print '\t-FD=X : Perform Form Difference analysis to account for the landmarks affecting the'\
		      ' form the most. '\
		      'X is the name of the visualization structure (as in the fasta file)( Default : No )'
		print '\t-cleanFD=X : Creates a new gm file stripping out the "Pinocchio" landmarks, by'\
		      ' performing a Grubbs test on the FDV for outliers'\
		      'X is the name of the chain (structure) where the analysis will be vizualized( Default :'\
		      ' No )'
		print '\t-FDmatrix : Calculates the form deformation vectors for all the shapes against the'\
		      'mean, and returns a matrix'
		print '\t-cazy=XXX : will create a clasifier file giving the cazy (http://www.cazy.org/) GH13'\
		      ' subfamilies. XXX is either'\
		      ' remote (using web-blast) or local for the blast query ( Default: No )'
		print '\t-m=PDB: will assume that a multiple PDB alignment was performed (as opposed as a MStA).'\
		      ' The pdb provided has to be in the folder ( Default: False ).'
		print '\nExample:'
		print 'python protGM.py alpha-amylase_bent -PCoa -CVA=kingdom,class -FD=A'
		print 'The clasifiers file must be in csv ";"-delimited format'\
		      'and named "[prefix].cls".'
		print 'Files required: .gm, .pdb, .landmarks, CazyDB.bin (if cazy is used)'
		sys.exit()

	# Default Parameters ###############################################################################
	prefix=sys.argv[1]
	scatter=False
	pcoa=False
	spc=False
	manova=False
	CanVA=False
	FD=False
	cfd=False
	FDm=False
	datatype='homologs'
	cazy=False
	remote=False
	nosingles=False
	multiple=False
	# Get user input ###################################################################################
	for arg in sys.argv[2:]:
		if arg == '-scatter3d':
			scatter = True
			r('library(scatterplot3d)')
		elif '-PCoA' in arg:
			pcoa=True
			if '=' in arg:
				chain=arg[6:]
		elif arg.startswith('-subPCoA='):
			pcoa=True
			spc=True
			t = arg[9:].split(',')
			x1,x2,y1,y2 = float(t[0]),float(t[1]),float(t[2]),float(t[3])
		elif arg.startswith('-perMANOVA='):
			manova=True
			per=int(arg[11:])
		elif '-CVA' in arg:
			CanVA=True
			if '=' in arg:
				chain=arg[5:]
		elif arg.startswith('-FD='):
			FD=True
			chain=arg[4:]
		elif arg.startswith('-cleanFD='):
			cfd=True
			chain=arg[9:]
		elif arg == '-FDmatrix':
			FDm=True
		elif arg == '-md':
			datatype='md'
		elif arg.startswith('-cazy='):
			cazy=True
			if arg[6:] == 'remote':
				remote=True
		elif arg == '-nosingles':
			nosingles=True
		elif arg.startswith('-m='):
			multiple=True
			chain=arg[3:]
			
	# get all R funtions
	R_functions()
	# get the class instance
	G=GM(prefix,multiple, datatype, 3, cazy, nosingles,chain)
	# 3D plot of raw coordinates #######################################################################
	if scatter:	
		G.scatter3D()

	# Principal coordinates analysys PCoA ##############################################################
	if pcoa:
		G.PCoA(G.prefix,G.gm, False,G.namestr)
		#G.PCoA(G.prefix+'_transp',G.tm, False,'NULL')
	# perMANOVA ########################################################################################
	if manova:
		G.perMANOVA()

	# CVA ##############################################################################################
	if CanVA:
		G.CVA()

	# Landmarks influencing the shape (from the mean) ##################################################
	if FD:
		G.formD(multiple)

	# Strip off landmarks influencing the shape (from the mean) ########################################
	if cfd:
		G.cleanFD(chain,multiple)
	# PCoA in the transposed FD matrix #################################################################
	if FDm:
		name , M, tM = G.FDmatrix(prefix+'_FDA')
		G.PCoA(name+'_distBclean',M, False, G.namestr)
		G.PCoA(name+'_transpBclean',tM, False,'NULL')
		H=GM(prefix+'_FDclean',3,cazy)
		name , M, tM = H.FDmatrix(prefix+'_FDA')
		H.PCoA(name+'_distAclean',M, False, H.namestr)
		H.PCoA(name+'_transpAclean',tM, False,'NULL')	
