#!/usr/bin/python

''' 
Simulator of wing, using cholesky decomposition to generate correlated variables.
Requieres an initial shape in GM format, and a partition file of the form:

Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation among partition 1]

Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation among partition 2]
.
.
.
Partitionn
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation among partition n]

'''
#importing bit####################################################################################################
from rpy2.robjects import r
import sys, os
from random import normalvariate, shuffle
import matplotlib.pyplot as plt
r('library(corpcor)')
#r('library(MASS)')
# End importing####################################################################################################

#Start fuctions####################################################################################################
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
		co = el[el.rfind(': ')+1:].strip(' ').split()
		partitions.append((la,co))
	landmarks =[]
	for p in partitions:
		landmarks.extend((int(x) for x in p[0]))
	landmarksN = max(landmarks)*2
	
	return partitions, landmarksN

def CreateCorrelationMatrix(partitions, landmarksN):
	'''
	Given the partition list, creates a correlation matrix in R
	'''
	r('c <- matrix(NA,%d,%d)'%(landmarksN,landmarksN))
	for p in partitions:
		tempx=[]
		tempy=[]
		for e in p[0]:
			tempx.append((int(e)*2)-1)
			tempy.append((int(e)*2))
		#p[0].extend(temp)
		for i in range(1,landmarksN+1):
			for x in tempx:#p[0]:
				if i in tempx:#(int(x) for x in p[0]):
					if i == int(x):
						r('c[%d,%d]<- %f'%(int(i),int(x),1.0))
					else:
						r('c[%d,%d]<- %f'%(int(i),int(x),float(p[1][0])))
				else:
					r('c[%d,%d]<- %f'%(int(i),int(x),0.00))
		#for j in range(2,landmarksN+1,2):
			for y in tempy:#p[0]:
				if i in tempy:#(int(x) for x in p[0]):
					if i == int(y):
						r('c[%d,%d]<- %f'%(int(i),int(y),1.0))
					else:
						r('c[%d,%d]<- %f'%(int(i),int(y),float(p[1][0])))
				else:
					r('c[%d,%d]<- %f'%(int(i),int(y),0.00))
	r('png("CorrelationMatrix.png")')
	r('hist(c)')
	r('dev.off()')
	#r('save(c,file ="test" )')
	
'''def CreateRandomMatrix(landmarksN, samplen, sd):
	
	Create a random matrix which each entry has a given mean and sd
	'''

			
def CreateCorrelatedRandomMatrix(prefix, samplen, sd, landmarksN):
	'''
	Reads the shapefile (must have only one entry of the reference shape in gm format,
	and must have a txt extension) and creates a random matrix with apropriate dimensions 
	and apply cholesky decomposition to create a correlated matrix
	'''
	r('s <- read.table("%s",sep=";",row.names=1)'%(prefix+'.txt'))# read the shape file
	r('m <- matrix(NA,%d,%d)'%(samplen,landmarksN))
	for i in range(landmarksN):
		r('m[,%d]<-rnorm(%d,rnorm(1,0,%f),%f)'%(i+1,int(samplen),float(sd),float(sd)))	
	r('s <- matrix(as.numeric(s), %d, dim(s)[2], byrow=TRUE)'%(samplen))#create a matrix repeating the initial shape
	#r('m<- matrix(rnorm(dim(s)[1]*dim(s)[2], sd=%f), dim(s)[1],dim(s)[2])'%(sd))#create a random matrix with desired sd
	#r('xp<-scale(m)') # Scale the random matrix
	r('t<-m%*%solve(chol(var(m)))')# comment out this if exact correlation is not needed
	r('cor.mat<-make.positive.definite(c)') # Transform the correlation matrix in positive definite
	r('t.new<-t%*%chol(cor.mat)') # Create the new correlated variables using Choleski decomposition or the correlation matrix
	#r('t.new[,1]<-t.new[,1]*attr(xp,"scaled:scale")+attr(xp,"scaled:center")') # Create the new dataset with correlation
	r('correlated <- t.new+s')
	r('write.table(correlated, "%s", sep = ";", row.names = TRUE, col.names = FALSE)'%(prefix+'correlated.gm'))
	r('png("CorrelatedMatrix.png")')
	r('hist(cor(correlated))')
	r('dev.off()')	
	r('save(c,m,s,t,t.new,pc,file ="test" )')
# End functions####################################################################################################

# Aplication of the code ##########################################################################################
if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'usage ShapeSim.py [prefix] [option]'
		print 'Options:'
		print '\t-sample=XXX : Create XXX samples ( Default: 100 )'
		print '\t-mean= YYY : Create a distribution of points with mean YYY ( Default: 0.0 )'
		print '\t-sd= ZZZ: Create a distribution of points with ZZZ standard deviation ( Default: 1.00 )'
		print '\t-partition=ZZZ : Use ZZZ file to get the partitions to be simulated ( Default: [prefix].par )'
		
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
## #####################################################################
for arg in sys.argv:
		if arg.startswith('-sample='):
			samplen = int(arg[8:])
		elif arg.startswith('-sd='):
			sd = float(arg[4:])
		elif arg.startswith('-partition='):
			pfilename = arg[11:]
		elif arg == '-procrustes':
			procrustes = True
			
samples = []
f=open(prefix+'.txt').read().strip().split(';')
if f[-1] == '':
	f=f[1:-1]
else:
	f=f[1:]
if len(f) <= len(colors):
	col = True
	
samples.append(f)

plot_gm(prefix+'.txt',col)

partitions, landmarkN = ReadPartitionFile(pfilename)

CreateCorrelationMatrix(partitions, landmarkN)

#CreateRandomMatrix(landmarkN, samplen, sd)

CreateCorrelatedRandomMatrix(prefix, samplen, sd, landmarkN)

plot_gm(prefix+'correlated.gm',col)

if procrustes:
	os.system('python /home/jshleap/LabBlouin/code/GM/procrustes.py %s 2'%(prefix+'correlated'))