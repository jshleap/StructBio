#!/usr/bin/python
"""
ShapeSim.py Copyright (C) 2012 Jose Sergio Hleap

Simulates shapes using cholesky decomposition to generate correlated variables (coordinates) with a shape constraint.
Requires an initial shape in GM format, and a partition file of the form:

Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation among partition 1]

Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation among partition 2]
.
.
.
Partition
Landmarks: [list of the landmark number in the order of entries in the GM file]
corr: [the ammount of correlation among partition n]

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

Requires:
1) rpy2, available at http://rpy.sourceforge.net/rpy2.html
2) Matplotlib, available at http://matplotlib.org/
3) R packages corpcor and shapes
"""

#importing bit####################################################################################################
from rpy2.robjects import r
import sys, os, glob
from random import shuffle, uniform
from math import sqrt
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
r('arr <- function(A,d){s<-dim(A)[1]; k<-(dim(A)[2]/d); Data <- as.numeric(A);arr<-array(Data,c(s, d, k));arr<-aperm(arr,c(3,2,1)); arr}')
###################################################################################################
###################################################################################################
## Convert 'shapes'-readable array into GMfile 
## A is an array Shapes-type
## d is the dimensions
r('A2GM<- function(A,d){m<-matrix(NA,dim(A)[3],dim(A)[1]*d) ; for (i in 1:dim(A)[3]){ for (j in 1:d){ m[i,seq(j,dim(m)[2],2)]<-A[,j,i]}};as.data.frame(m)}')
###################################################################################################
###################################################################################################
## Transform a random matrix into a correlated matrix using cholesky decomposition and a given correlation matrix
r('Cholesky<- function(random,correlation){t<-random ; cor.mat<-make.positive.definite(correlation) ;  t.new<-t%*%chol(cor.mat) ; as.matrix(t.new)}')
# Old with variance control: r('Cholesky<- function(random,correlation){t<-random%*%solve(chol(var(random))) ; cor.mat<-make.positive.definite(correlation) ;  t.new<-t%*%chol(cor.mat) ; as.matrix(t.new)}')

###################################################################################################
r('save(arr,A2GM,Cholesky,file="functions.R")')


def create_pars_inter(pfilename,intra):
	'''
	Given an original par file, will create the same partitions with intercorelation ranging from 0.0 to inter
	'''
	partitions, landmarksN = ReadPartitionFile(pfilename,prot)
	os.system('rm %s'%(pfilename))
	for i in range(0,(int(round(intra*10)))*2):
		inter = i*0.05
		ouf = open(pfilename[:-3]+str(intra)+'_'+str(inter)+'.par','w')
		for p in partitions:
			la=''
			for e in p[0]:
				la+=e+','
			ouf.write('Partition\nLandmarks: %s\ncorr: %f \ninter: %f \n\n'%(la[:-1],float(intra), inter))
		ouf.close()
			
def simulate_intercorrelation(pfilename, samplen,sd,moduler):
	'''
	 will simulate shapes with intercorrelations ranging from 0.0 to intracorrelation
	'''
	cwd = os.getcwd()
	print cwd
	for d in range(1,11):
		intra = d/10.0
		os.system('mkdir %.1f'%(d/10.0))
		os.system('cp ./*.txt ./*.par ./%.1f'%(d/10.0))
		os.chdir('./%.1f'%(d/10.0))
		print os.getcwd()
		create_pars_inter(pfilename, d/10.0)
		parfiles = glob.glob('*.par')
		shape=pfilename[:-3]+'txt'
		for f in parfiles:
			os.system('mv %s %s'%(shape,f[:-3]+'txt'))
			shape = f[:-3]+'txt'
			os.system('python ~/LabBlouin/code/GM/ShapeSim.py %s -sd=%s -sample=%s -partition=%s'%(f[:-4],str(sd), str(samplen),f))
			os.system('python ~/LabBlouin/code/modularity/modulerV2.py %s -ndim=2 -permt=1000 -power=0.8'%(f[:-4]))
		os.chdir(cwd)
	
		

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

def ProteinPartitionGenerator(prefix,sample):
	'''
	Giving a txt (a single line GM file), a original landmarks file and a grahcluster file (membership vector), create a
	partition file and a landmarks for protein data
	'''
	print 'Creating a partition file for a protein dataset. Be aware that the graphcluster file'\
	      ' and the txt file must match the number of landmarks'
	land=open(prefix+'.landmarks').read().split('\n>')[0].split('\n')
	nl=open(prefix+'.landmarks','w')
	for la in range(sample):
		nl.write('>"%d"\n'%(la)+'\n'.join(land[1:])+'\n')
	nl.close()
	p=open(prefix+'.part','w')
	g=open(prefix+'.graphcluster').read().strip().split()
	s=set(g)
	l=list(s)
	for i in l:
		p.write('Partition\nLandmarks: ')
		li=''
		for e in range(len(g)):
			if i == g[e]:
				li+=str(e+1)+','
		p.write(li[:-1]+'\ncorr: %f\n\n'%(uniform(0.15,0.85)))
	p.close()

def ReadPartitionFile(pfilename,prot):
	'''
	Take a partition file and extract the landmarks and correlation in each partition
	'''
	f = open(pfilename).read().strip().split('Partition')
	partitions=[]
	# get the landmarks and corfor each partition
	for el in f[1:]:
		la = el[el.find(': ')+1:el.find('\ncorr:')].strip(' ').split(',')
		co = el[el.rfind('corr: ')+5:].strip(' ').split()
		#inter = el[el.rfind('inter: ')+6:].strip(' ').split()
		partitions.append((la,co))#,inter))
	landmarks =[]
	for p in partitions:
		landmarks.extend((int(x) for x in p[0]))
	if prot:
		landmarksN = max(landmarks)*3
	else:
		landmarksN = max(landmarks)*2
	
	return partitions, landmarksN

###############################################################################################

def GetDistance(coords_i,coords_j):
	'''
	calculate the euclidean distances between coordinates i and j
	'''	
	#if dim == 2:
	d= sqrt(((float(coords_i[0])-float(coords_j[0]))**2)+((float(coords_i[1])-float(coords_j[1]))**2))
	#elif dim == 3:
	#	d= sqrt(((coords_i[0]-coords_j[0])**2)+((coords_i[1]-coords_j[1])**2)+((coords_i[2]-coords_j[2])**2))
	return d

###############################################################################################

def contact2D(prefix,sd):
	'''
	Create a contact file and a landmarks file for a 2d shape. Assumes 1 unit grid 
	'''
	outf=open(prefix+'.landmarks','w')
	outf.write('>Shape\n')
	fout = open(prefix+'.contacts','w')
	radius = 1+float(sd)
	coor = open(prefix+'.txt').read().split(';')[1:]
	x = [coor[x] for x in range(0,len(coor),2)]
	y = [coor[y] for y in range(1,len(coor),2)]
	for i in range(len(x)):
		outf.write(str(i)+'\t'+str(i)+'\t'+'Point\n')
		for j in range(len(x)):
			if not j == i:
				dist = GetDistance((x[i],y[i]),(x[j],y[j]))
				if dist <= radius:
					fout.write('(%d,%d)\n'%(i,j))
	fout.close()
	
def correlation_matrix(prefix, partitions, landmarksN,sd):
	'''
	Creates a correlation matrix with sligthly random deviation from original shape
	'''
	if prot:
		dim=3
		correction=(2,1)
	else:
		dim=2
		correction=(1,0)
	r('s <- read.table("%s",sep=";",row.names=1)'%(prefix+'.txt'))# read the shape file
	r('s <- matrix(as.numeric(s), %d, dim(s)[2], byrow=TRUE)'%(samplen))#create a matrix repeating the initial shape
	#r('x <- matrix(rnorm(dim(s)[1]*dim(s)[2], sd=%f),dim(s)[1],dim(s)[2])'%(float(sd))) # random matrix with slight variation
	#r('s.new <- s+x') # add the variation to the shape matrix
	if inter:
		r('co <-  matrix(%f,%d,%d)'%(uniform(0,maxIc),landmarksN,landmarksN)) # get the correlation matrix with the inter and intra of the partitionfile
	else:
		r('co <-  matrix(%f,%d,%d)'%(0,landmarksN,landmarksN))# get the correlation matrix with intra of the partitionfile (inter 0)
	#get the index of the partitions' landmarks
	for p in range(len(partitions)):
		tempx=[]
		tempy=[]
		tempz=[]
		coors=[]
		for e in partitions[p][0]:
			tempx.append((int(e)*dim)-correction[0])
			tempy.append((int(e)*dim-correction[1]))
			if prot:
				tempz.append((int(e)*dim))
		coors.extend(tempx)
		coors.extend(tempy)
		if prot:
			coors.extend(tempz)
		coors.sort()
		for x1 in coors:
			for x2 in coors:
				if x1 == x2:
					r('co[%d,%d]<- %f'%(int(x1),int(x2),1.0))
				else:
					r('co[%d,%d]<- %f'%(int(x1),int(x2),float(partitions[p][1][0])))

def PCM(prefix, partitions):
	'''
	Creates a correlation matrix given a GPA of a sligthly random deviation from original shape
	'''
	r('s <- read.table("%s",sep=";",row.names=1)'%(prefix+'.txt'))# read the shape file
	r('s <- matrix(as.numeric(s), %d, dim(s)[2], byrow=TRUE)'%(samplen))#create a matrix repeating the initial shape
	r('x <- matrix(rnorm(dim(s)[1]*dim(s)[2], sd=0.01),dim(s)[1],dim(s)[2])') # random matrix with slight variation
	r('s.new <- s+x') # add the variation to the shape matrix
	r('M <- procGPA(arr(s.new,2))') # perform GPA
	r('co <- cor(A2GM(M$rotated,2))') # get the correlation matrix of the procrustes superimposition
	r('d <- density(co[upper.tri(co)])')
	#get the index of the partitions' landmarks
	for p in range(len(partitions)):
		tempx=[]
		tempy=[]
		coors=[]
		for e in partitions[p][0]:
			tempx.append((int(e)*2)-1)
			tempy.append((int(e)*2))
		coors.extend(tempx)
		coors.extend(tempy)
		coors.sort()
		for x1 in coors:
			for x2 in coors:
				if x1 == x2:
					r('co[%d,%d]<- %f'%(int(x1),int(x2),1.0))
				else:
					r('co[%d,%d]<- %f'%(int(x1),int(x2),float(partitions[p][1][0])))
	r('save(s,x,s.new,M,co, file="PCM.R")')
	
def CorrelatedMatrix(prefix, sd, hist, proc):
	'''
	Creates the correlated matrix by cholesky decomposition
	'''
	r('m <- matrix(rnorm(dim(s)[1]*dim(s)[2], sd=%f),dim(s)[1],dim(s)[2])'%(sd))
	#r('m<-make.positive.definite(c)')
	r('t.new <- Cholesky(m,co)') 
	#r('plotshapes(arr(t.new,2))')
	r('correlated <- t.new+s')
	if hist and not proc:
		r('png("correlated.png"); X <- cor(correlated)')# ; dev.off()')
		r('hist(X[upper.tri(X)], breaks = 100, main = paste("Histogram of the upper triangle of the correlation matrix before GPA"),freq=FALSE)')
		#r('d <- matrix(rnorm(dim(correlated)[1]*dim(correlated)[2]), dim(correlated)[1],dim(correlated)[2])')
		#r('lines(d,col="red")')
		r('dev.off()')
	#r('plotshapes(n$rotated)')
	if hist and proc:
		r('n <- procGPA(arr(correlated,2))')
		r('png("correlatedProc.png"); X <- cor(A2GM(n$rotated,2)); dev.off()')
		r('hist(X[upper.tri(X)], breaks = 100, main = paste("Histogram of the upper triangle of the correlation matrix after GPA"),freq=FALSE)')
		r('d <- matrix(rnorm(dim(correlated)[1]*dim(correlated)[2]), dim(correlated)[1],dim(correlated)[2])')
		r('dp <- procGPA(arr(d,2)); cd <- cor(A2GM(dp$rotated,2)) ; d <- cd[upper.tri(cd)]')
		r('lines(density(d),col="red")')
		r('dev.off()')
	if proc:
		r('n <- procGPA(arr(correlated,2))')
		r('write.table(A2GM(n$rotated,2), "%s.gm", sep = ";", row.names = TRUE, col.names = FALSE)'%(prefix+'_proc'))
	r('write.table(correlated, "%s.gm", sep = ";", row.names = TRUE, col.names = FALSE)'%(prefix))
	if not prot: plot_gm(prefix+'.gm',col)
	if hist and proc:
		r('save(m,t.new,correlated,n,co,file="CorrelatedMatrix.R")')
	else:
		r('save(m,t.new,correlated,co,file="CorrelatedMatrix.R")')
		
def includeNonShapeVariables(prefix, samplen, angle, landmarksN):
	'''
	includes in the correlated matrix variables as scale, traslation, and rotation
	'''
	a=angle.split(',')
	r('odd <- seq(1,%d,2) ; even <- seq(2,%d,2)'%(landmarksN,landmarksN))
	r('scale <- runif(%d,0,1.5)'%(samplen))
	r('angle <- runif(%d,%f,%f)'%(samplen,float(a[0]),float(a[1])))
	r('angle[1]<-0')
	r('tra <- matrix(runif(%d*2,max(correlated)*0.05,max(correlated)*0.2),%d,2)'%(samplen,samplen))
	r('ts <- correlated*scale ; tt <- matrix(NA,dim(ts)[1],dim(ts)[2])')
	r('tt[,odd] <- ts[,odd]+tra[,1] ; tt[,even]<- ts[,even]+tra[,2]')
	r('tr <- matrix(NA,dim(tt)[1],dim(tt)[2]); tr[,odd] <- tt[,odd]*cos(angle)+ tt[,even]*sin(angle) ; tr[,even] <- tt[,odd]*-sin(angle)+ tt[,even]*cos(angle)')
	r('save(odd,even, scale, angle, tra, ts, tt, tr, file="nonshape.R")')
	r('write.table(tr, "%sintrinsic_rotated.gm", sep = ";", row.names = TRUE, col.names = FALSE)'%(prefix))
	r('write.table(A2GM(procGPA(arr(tr,2))$rotated,2), "%sintrinsic_rotatedProc.gm", sep = ";", row.names = TRUE, col.names = FALSE)'%(prefix))
# End functions####################################################################################################

# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	if len(sys.argv) == 1 or 'help' in sys.argv or '-h' in sys.argv:
			print 'usage ShapeSim.py [prefix] [option]'
			print 'Options:'
			print '\t-sample=XXX : Create XXX samples ( Default: 100 )'
			print '\t-mean= YYY : Create a distribution of points with mean YYY ( Default: 0.0 )'
			print '\t-sd= ZZZ: Create a distribution of points with ZZZ standard deviation ( Default: 1.00 )'
			print '\t-partition=ZZZ : Use ZZZ file to get the partitions to be simulated ( Default: [prefix].part )'
			print '\t-nonshape=VV,WW : Use the range VV-WW to randomly assign degrees to rotate the final configurations. This option also',
			print 'translate and scale the configuration ( Default: No )'
			print '\t-moduler : Apply modulerV2 to the resulting (procrustean) data ( Default: no )'
			print '\t-PCM : Correlation matrix will include the correlation imposed by procrustes ( Default: No )'
			print '\t-SI : Simulates intercorrelation from 0.0 to intracorrelation in step of 0.05 ( Default: No )'
			print '\t-histograms : Will plot some histograms of the correlation ( Default : No)'
			print '\t-procrustes : Write a gm file with the full procrustes coordinates and perform modularity inference on that (Default : No)'
			print '\t-protein: simulating a protein in 3d'
			print '\t-c: create a contact file of a 2D shape'
			print '\t-inter=MAX: will create intercorrelation taken from a uniform random distribution between'\
			      ' 0 and a max value that has to be specified ( Default: 0)'
	##Default Parameters ###################################################
	sd = 1.00
	samplen = 100
	prefix = sys.argv[1]
	pfilename = prefix+'.part'
	colors = ['#8b8989','#cdb79e','#000000','#2f4f4f','#d3d3d3','#191970','#0000cd','#87cefa','#b0c4de','#b0c4de',
		      '#5f9ea0','#66cdaa','#7fffd4','#006400','#556b2f','#8fbc8f','#20b2aa','#7cfc00','#bdb76b','#f0e68c', 
		      '#ffff00','#ffd700','#b8860b','#b8860b','#bc8f8f','#cd5c5c','#8b4513','#a0522d','#b22222','#fa8072',
		      '#ff0000', '#ff69b4', '#ff1493','#9400d3']
	shuffle(colors)
	col = False
	nonshape=False
	moduler = False
	pcm = False
	si = False
	hist = False
	proc = False
	prot = False
	contact = False
	inter = False
	## #####################################################################
	for arg in sys.argv:
			if arg.startswith('-sample='):
				samplen = int(arg[8:])
			elif arg.startswith('-sd='):
				sd = float(arg[4:])
			elif arg.startswith('-partition='):
				pfilename = arg[11:]
			elif arg.startswith('-nonshape='):
				nonshape = True
				angle = arg[10:]
			elif arg == '-moduler':
				moduler = True
			elif arg == '-PCM':
				pcm = True
			elif arg == '-SI':
				si = True		
			elif arg == '-histograms':
				hist = True
			elif arg == '-procrustes':
				proc = True
			elif arg == '-protein':
				prot = True
			elif arg == '-c':
				contact = True
			elif arg.startswith('-inter='):
				inter = True
				maxIc = float(arg[7:])
								
	samples = []
	f=open(prefix+'.txt').read().strip().split(';')
	if f[-1] == '':
		f=f[1:-1]
	else:
		f=f[1:]
	if len(f) <= len(colors):
		col = True
		
	samples.append(f)
	print 'Contact = %s'%contact
	if contact:
		contact2D(prefix,sd)
	
	if prot:
		dim = 3
	else:
		plot_gm(prefix+'.txt',col)
		dim=2
	if si:
		print 'simulating intercorrelation'
		simulate_intercorrelation(pfilename, samplen,sd, moduler)
		sys.exit(-1)
	else:
		if prot:
			if not os.path.isfile(prefix+'.part'):
				ProteinPartitionGenerator(prefix,samplen)
		partitions, landmarksN = ReadPartitionFile(pfilename,prot)
		print 'Reading partition file'
		if pcm:
			print 'Creating a correlated matrix with GPA noise'
			PCM(prefix, partitions)
		else:
			print 'Creating a correlation matrix with desired sd'
			correlation_matrix(prefix, partitions, landmarksN,sd)
		
		print 'Performing Cholesky decomposition on the matrix to correlate the variables'	
		CorrelatedMatrix(prefix, sd, hist, proc)
		
		if nonshape:
			print 'Including rotation, translation and scaling'
			includeNonShapeVariables(prefix, samplen, angle, landmarksN)
			if moduler:
				print "perform modularity on procrustean data"
				os.system('python ~/LabBlouin/code/modularity/modulerV2.py %s -d %d -p 0.8 -c -n 999 -r -g'%(prefix+'intrinsic_rotatedProc',dim))
		else:
			if moduler:
				print 'Executing modulerV3 on the correlated dataset'
				os.system('python ~/LabBlouin/code/modularity/modulerV3.py %s -d %d -n 999 -l -r'%(prefix,dim))
				if proc:
					print 'Executing modulerV2 on the correlated and "procrutean" dataset'
					os.system('python ~/LabBlouin/code/modularity/modulerV2.py %s -d %d -n 999 -p 0.8'%(prefix+'_proc',dim))
		os.system('python /home/jshleap/LabBlouin/code/ABeRMuSA/alignment/GMtoRMSD.py %s'%(prefix))
