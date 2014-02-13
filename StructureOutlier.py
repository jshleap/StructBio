#!/usr/bin/python
"""
StructureOutlier.py Copyright (C) 2012 Jose Sergio Hleap

This script will identify full structure outliers from a gm file. It is based on the Grubss test with RMSD as 
variable.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

Requires:
1) rpy2, available at http://rpy.sourceforge.net/rpy2.html
2) Numpy
3) R packages outliers
"""

#importing bit####################################################################################################
import sys, os, math
import numpy as np
from rpy2.robjects import r
r('library(outliers)')
# End importing####################################################################################################

#Some definitions##################################################################################################

def find_key(dic, val):
	"""return the key of dictionary dic given the value"""
	return [k for k, v in dic.iteritems() if v == val][0]

def bylettersplit(string):
	split = []
	for l in string:
		split.append(l)
	return split

def listelinlist(test, reference):
	boole = False
	for e in test:
		if e in reference:
			boole = True
	return boole

def fixgmfile(gmfileprefix):
	forbiden = map(chr,range(32,65))
	forbiden.pop(14)
	inf = open(gmfileprefix+'.gm')
	outfn = gmfileprefix+'.r'
	ouf=open(outfn+'.gm','w')	
	for line in inf:
		bline = line.split(';')
		if listelinlist(bylettersplit(bline[0][1:]),forbiden):
			ouf.write('"'+bline[0]+'"'+line[len(bline[0]):])
	return outfn
			

def GM2R_meanshape(gmfileprefix):
	'''
	Passes the GM file to R, computes the mean shape and write a file including the mean shape
	'''
	forbiden = map(chr,range(32,65))
	if listelinlist(bylettersplit(open(gmfileprefix+'.gm').readline().split(';')[0]),forbiden):
		gmfileprefix = fixgmfile(gmfileprefix)
	#read the gm file
	r('m<-read.table("%s.gm", sep=";", row.names=1)'%(gmfileprefix))
	r('m<-m[1:dim(m)[2]-1]')
	# Set it as matrix
	r('m<-as.matrix(m)')
	# Create an empty matrix with an extra row
	r('n<-matrix(NA, nrow=dim(m)[1]+1, ncol=dim(m)[2])')
	#fill the matrix with the previous data
	r('n[1:dim(m)[1],]<-m')
	#Set the names of rows
	r('row.names(n)<-c(row.names(m),">mean")')
	# get the average shape and put it as the last element
	r('n[dim(m)[1]+1,]<- apply(m,c(2),mean)')
	# write to file
	filename= '%s_mean.csv'%(gmfileprefix)
	r('write.table(n, file = "%s",sep = ";",col.names = FALSE)'%(filename))
	return filename


def GM2npMatrix(gmfileprefix):
	'''
	Convert the GM file into a numpy matrix and a list of names
	'''
	inf = open(gmfileprefix+'.gm')
	data = []
	
	names = []
	for line in inf:
		if line.strip().endswith(';'):
			line = line.strip()[:-1]
		temp=[]
		bline= line[line.find(';')+1:].strip().split(';')
		for el in bline:
			temp.append(float(el))
		data.append(temp)
		names.append(line[:line.find(';')])
	m = np.matrix(data)
	return m, names

def meanshape(matrix):
	'''
	Calculates the meanshape using a numpy matrix
	'''
	msh = matrix.mean(0)
	return msh

def rmsd_dict(matrix,names, msh):
	'''
	returns a dictionary with the name of the structure to be passed to
	RMSDfromMEAN using iteritems
	'''
	lines = {}
	for i in range(matrix.shape[0]):
		l = ''
		for j in range(matrix.shape[1]):
			l += str(matrix[i,j]) +';'
		lines[names[i]] = l[:-1]
	ml=''
	for k in range(msh.shape[1]):
		ml += str(msh[0,k]) + ';'
	lines['>mean']=ml[:-1]
	return lines
			
			
def storelines(csvfilename):
	'''
	Read the csv file generated by GM2R_meanshape, and store each line into an 
	element of a dictionary
	'''
	lines = {}
	f = open(csvfilename)
	token = f.read()
	fsplit = token.split('">')
	fsplit=fsplit[1:]
	for i in range(len(fsplit)):
		tok = fsplit[i].split('";')
		lines[tok[0]] = tok[1][:-1]

	return lines   


def lines2list(dictval):
	ll = []
	for d in range(3):
		temp=[]
		for i in range(d,len(dictval.split(';')),3):
			temp.append(float(dictval.split(';')[i]))
		ll.append(temp)
	return ll

def RMSDfromMEAN(p,m):
	'''
	Calculates the RMSD for the meanshape (m) and the problem shape (p)
	'''      
	sum_dist=0.0
	count = 0
	for j in range(len(m[0])):
		d =(m[0][j]-p[0][j])**2 + (m[1][j]-p[1][j])**2 + (m[2][j]-p[2][j])**2
		sum_dist += d
		count += 1
	# calculate the sqrt of the mean deviation
	RMSD = math.sqrt(sum_dist/count)

	return RMSD

def outlier_test(RMSDstr, linesdict,revdict):
	'''
	Takes a string of the RMSD in the form "f1,f2,f3 ... fn", and pass it to R to do a Grubbs test
	'''
	#popsout=[]
	r('library(outliers)')
	r('rmsd <- c(%s)'%(RMSDstr))
	for s in range(len(rmsd)):
		#test statistically prescence of outliers
		r('gr<-grubbs.test(rmsd)')
		sig=r('gr[3]')
		sig=sig[0][0]
		#remove values if outlier
		#count = 0
		if sig<0.05:
			r('rmsd<-rm.outlier(rmsd)')
			outl=r('gr[2]')
			outl=float(outl[0][0].split()[2])
			#popsout.append(l)
			try:
				linesdict.pop(revdict['%.10f'%round(outl,10)])#find_key(linesdict,outl))
				rm = revdict['%.10f'%round(outl,10)]
			except:
				try:
					linesdict.pop(revdict[str(outl)[:12]])
					rm = revdict[str(outl)[:12]]
				except:
					linesdict.pop(revdict['%.10f'%outl])
					rm = revdict['%.10f'%outl]
			print 'structures removed: %s'%(rm)
	print 'Quantiles after removal:'
	print r('quantile(rmsd)')
	return linesdict


def GM_clean(linesdict, prefix):
	'''
	Write a new GM file without the outliers
	'''
	ouf = open(prefix+'_clean.gm','w')
	for k,v in linesdict.iteritems():
		ouf.write('%s;%s\n'%(k,v))
	ouf.close()

	
#End definitions ##################################################################################################

# Application of the code #########################################################################################

prefix =  sys.argv[1]

m,names = GM2npMatrix(prefix)
lines = rmsd_dict(m,names, meanshape(m))
rmsd=''
rmsdlist=[]
revdict={}
mean = lines2list(lines.pop('>mean'))
for k, v in lines.iteritems():
	p = lines2list(v)
	rms = RMSDfromMEAN(p,mean)
	rmsdlist.append(rms)
	rmsd += str(rms) + ','
	revdict['%.10f'%rms]=k
print '0%=', min(rmsdlist), '50%=', np.median(rmsdlist), '100%=', max(rmsdlist)
GM_clean(outlier_test(rmsd[:-1],lines,revdict),prefix)