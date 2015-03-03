#!/usr/bin/python

''' 
Simulator of wing, using cholesky decomposition to generate correlated variables.
Requieres an initial shape in GM format, and a partition file of the form:

Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation within partition 1]
inter: [the ammount of correlation betwen partition 1 and the rest]

Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation within partition 2]
inter: [the ammount of correlation betwen partition 2 and the rest]
.
.
.
Partitionn
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation within partition n]
inter: [the ammount of correlation betwen partition n and the rest]

Inter has to be the same!
'''
#importing bit####################################################################################################
from rpy2.robjects import r
import sys, os
from random import normalvariate, shuffle
import matplotlib.pyplot as plt
r('library(corpcor)')
r('library(shapes)')
# End importing####################################################################################################

#Start fuctions####################################################################################################
## Convert GMfile into 'shapes'-readable array
## A is the read.table variable of a GMfile
## s is the sample size
## k is the number of landmarks
## d is the number of dimensions (2D,2;3D,3)
r('arr <- function(A,d){s<-dim(A)[1]; k<-(dim(A)[2]/d); Data <- as.numeric(A);arr<-array(Data,c(s, d, k)); arr<-aperm(arr,c(3,2,1)); arr}')
###################################################################################################
###################################################################################################
## Convert 'shapes'-readable array into GMfile 
## A is an array Shapes-type
## d is the dimensions
r('A2GM<- function(A,d){m<-matrix(NA,dim(A)[3],dim(A)[1]*d) ; for (i in 1:dim(A)[3]){ for (j in 1:d){ m[i,seq(j,dim(m)[2],2)]<-A[,j,i]}};as.data.frame(m)}')
###################################################################################################
###################################################################################################
## Transform a random matrix into a correlated matrix using cholesky decomposition and a given correlation matrix
r('Cholesky<- function(random,correlation){t<-random%*%solve(chol(var(random))) ; cor.mat<-make.positive.definite(correlation) ;  t.new<-t%*%chol(cor.mat) ; as.matrix(t.new)}')
###################################################################################################


def plot_gm(filename, col):
	''' 
	plot a scatter GM file of 2D 
	'''
	#read and split file
	x=[]
	y=[]
	f = open(filename)
	plots = []
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_axis_off()
	count=0
	for line in f:
		count+=1
		bline = line.strip().split(';')[1:]
		for i in range(0,len(bline),2):
			#x.append(bline[i])
			#y.append(bline[i+1])
			x = float(bline[i])
			y = float(bline[i+1])
			if col:
				lin, = ax.plot(x,y, ls='None', marker='o', color=colors[(i+2)/2])
			else:
				lin, = ax.plot(x,y, ls='None', marker='o')
			if count == 1:
				ax.annotate('%s'%(str((i+2)/2)), xy=(x,y), fontsize=14, style='oblique', xytext=(x+2*2,y+2*2))

	plt.savefig('%s.png'%(filename[:filename.find('.')]), transparent=True)

def ReadPartitionFile(pfilename):
	'''
	Take a partition file and extract the landmarks and correlation in each partition
	'''
	f = open(pfilename).read().strip().split('Partition')
	partitions=[]
	# get the landmarks and corfor each partition
	for el in f[1:]:
		la = el[el.find(': ')+1:el.find('\ncorr:')].strip(' ').split(',')
		co = el[el.rfind('corr: ')+5:el.find('\ninter:')].strip(' ').split()
		inter = el[el.rfind('inter: ')+6:].strip(' ').split()
		partitions.append((la,co,inter))
	landmarks =[]
	for p in partitions:
		landmarks.extend((int(x) for x in p[0]))
	landmarksN = max(landmarks)*2
	
	return partitions, landmarksN

def CorrMarPart(partitions, hist=False):
	'''
	Creates a correlation matrix per partition
	'''
	part=len(partitions)
	l=65
	mnames=[]
	r('png("CorrelationMatrices.png")')
	r('par(mfrow=c(%d,1))'%(part))
	for p in range(len(partitions)):
		tempx=[]
		tempy=[]
		coors=[]
		landinpart=len(partitions[p][0])
		r('%s <- matrix(%f,%d,%d)'%(chr(l+p), float(partitions[p][2][0]),landinpart*2,landinpart*2))
		mnames.append(chr(l+p))
		for e in partitions[p][0]:
			tempx.append((int(e)*2)-1)
			tempy.append((int(e)*2))
		coors.extend(tempx)
		coors.extend(tempy)
		coors.sort()
		
		for i in range(len(coors)):
			for j in range(len(coors)):
				if coors[i] == coors[j]:
					r('%s[%d,%d]<- %f'%(chr(l+p),int(coors.index(coors[i]))+1,int(coors.index(coors[j]))+1,1.0))
				else:
					if coors[i] in tempx and coors[j] in tempx:
						r('%s[%d,%d]<- %f'%(chr(l+p),int(coors.index(coors[i]))+1,int(coors.index(coors[j]))+1,float(partitions[p][1][0])))
					elif coors[i] in tempy and coors[j] in tempy:
						r('%s[%d,%d]<- %f'%(chr(l+p),int(coors.index(coors[i]))+1,int(coors.index(coors[j]))+1,float(partitions[p][1][0])))
					else:
						r('%s[%d,%d]<- %f'%(chr(l+p),int(coors.index(coors[i]))+1,int(coors.index(coors[j]))+1,0.00))
		
		r('diag(%s)<-1'%(chr(l+p)))
		r('hist(as.matrix(%s), breaks=50)'%(chr(l+p)))
	r('dev.off()')
	
	r('png("images.png")')
	r('par(mfrow=c(%d,1))'%(len(mnames)))	
	for m in mnames:
		r('image(as.matrix(%s))'%(m))
	r('dev.off()')
	
	return mnames

def R_table2array2GPA2gm(prefix,dim):
	r('Ar<-arr(as.matrix(r),%d)'%(dim))
	r('pc <- procGPA(Ar)')
	r('plotshapes(pc$rotated)')
	r('gm <- A2GM(pc$rotated, %d)'%(dim))
	r('res <- A2GM(Ar,%d) - gm'%(dim))
	r('save(Ar, res, pc, gm, arr,A2GM, file ="proc.R" )')
	#r('write.table(gm,file="%s",sep=";",row.names=TRUE, col.names=FALSE)'%(prefix+'.gm'))
	
	
def RandomProcrustes(prefix, landmarksN, samplen, sd, dim):
	'''
	Create a random matrix with a given sd and shape, and then perform procrustes
	'''
	r('s <- read.table("%s",sep=";",row.names=1)'%(prefix+'.txt'))# read the shape file
	r('m <- matrix(NA,%d,%d)'%(samplen,landmarksN))
	for i in range(landmarksN):
		r('m[,%d]<-rnorm(%d,sd=%f)'%(i+1,int(samplen),float(sd)))	
	r('s <- matrix(as.numeric(s), %d, dim(s)[2], byrow=TRUE)'%(samplen))#create a matrix repeating the initial shape
	r('r <- s+m')
	R_table2array2GPA2gm(prefix,dim)

def correlatedDistribution(mnames, partitions, samplen,sd):
	'''
	creates a correlated matrix per partition starting with random variables
	with mean given by the fist procrustes coordinate
	'''
	s=''
	cholm=[]
	for p in range(len(partitions)):
		linp=''
		for e in partitions[p][0]:
			linp += str((int(e)*2)-1)+','+str(int(e)*2)+','
		linps=linp[:-1].split(',')
		r('%s%s <- matrix(NA, %d,%d)'%(mnames[p],chr(ord(mnames[p])+1),samplen,len(partitions[p][0])*2))
		for i in range(samplen):
			for j in range(len(linps)):
				r('%s%s[%d,%d]<- rnorm(1,mean=gm[%d,%d],%f)'%(mnames[p],chr(ord(mnames[p])+1),i+1,j+1,i+1,int(linps[j]),float(sd)))
		r('Chol%s<-Cholesky(%s%s,%s)'%(mnames[p],mnames[p],chr(ord(mnames[p])+1),mnames[p]))
		r('colnames(Chol%s)<-c(%s)'%(mnames[p],linp[:-1]))
		s+= mnames[p]+',Chol%s,'%(mnames[p])+mnames[p]+chr(ord(mnames[p])+1)+','
		cholm.append(('Chol%s'%(mnames[p]),linps))
	r('save(%s file="randommat")'%(s))
	return cholm

def SimulateProcrustes(landmarksN, samplen, mnames, partitions, cholm):
	'''
	starting from a procrustes coordinates simulate with correlated random numbers
	'''
	r('tr<-transformations(arr(as.matrix(gm),2),Ar)')
	r('corProc <- matrix(NA,%d,%d)'%(int(samplen),int(landmarksN)))
	for c in cholm:
		r('corProc[,as.numeric(colnames(%s))] <- %s'%(c[0],c[0]))
	r('transformed <- (corProc*tr$scale)- as.vector(t(tr$translation))')
	r('new.p <- procGPA(arr(as.matrix(transformed),2))$rotated')
	r('new.gm <- A2GM(pc$rotated,2)')
	r('save(tr,corProc,transformed,new.p,new.gm,file="sp")')
	r('png("correlated.png")')
	r('par(mfrow=c(2,2))')
	r('hist(cor(as.matrix(new.p)))')
	r('hist(cor(as.matrix(new.gm)))')
	r('hist(cor(as.matrix(transformed)))')
	r('dev.off()')
	
			
		
# End functions####################################################################################################

# Aplication of the code ##########################################################################################
if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'usage ShapeSim.py [prefix] [option]'
		print 'Options:'
		print '\t-sample=XXX : Create XXX samples ( Default: 100 )'
		print '\t-mean= YYY : Create a distribution of points with mean YYY ( Default: 0.0 )'
		print '\t-sd= ZZZ: Create a distribution of points with ZZZ standard deviation ( Default: 1.00 )'
		print '\t-partition=ZZZ : Use ZZZ file to get the partitions to be simulated ( Default: [prefix].par )'
		print '\t-histograms : Will plot some histograms of the correlation ( Default : No)'
		
##Default Parameters ###################################################
sd = 1.00
samplen = 100
prefix = sys.argv[1]
pfilename = prefix+'.par'
colors = ['#8b8989','#cdb79e','#000000','#2f4f4f','#d3d3d3','#191970','#0000cd','#87cefa','#b0c4de','#b0c4de',
          '#5f9ea0','#66cdaa','#7fffd4','#006400','#556b2f','#8fbc8f','#20b2aa','#7cfc00','#bdb76b','#f0e68c', 
          '#ffff00','#ffd700','#b8860b','#b8860b','#bc8f8f','#cd5c5c','#8b4513','#a0522d','#b22222','#fa8072',
          '#ff0000', '#ff69b4', '#ff1493','#9400d3']
shuffle(colors)
col = False
procrustes = False
histogram=False
## #####################################################################

partitions, landmarksN = ReadPartitionFile(pfilename)
mnames = CorrMarPart(partitions)
RandomProcrustes(prefix, landmarksN, samplen, sd, 2)
cholm = correlatedDistribution(mnames, partitions, samplen,sd)
SimulateProcrustes(landmarksN, samplen, mnames, partitions, cholm)