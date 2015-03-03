#!/usr/bin/python
''' 
Geometric morphometrics. Most R functions are taken and modified from Morphometrics with R(Series: Use R!) of Julien Claude.
require package shapes (R) and 
'''
# importing bit###################################################################################
import sys
import os
from rpy2.robjects import r
from protGM import R_functions
r('library(shapes)')
# #################################################################################################
# Some python functions#############################################################################
def tps2gm(filename):
	'''
	Convert TPS files into GM
	'''
	try:
		filen = filename+'.tps'
		inf=open(filen)
	except:
		filen = filename+'.TPS'
		inf=open(filen)
	foun = filename+'.gm'
	fou=open(foun,'w')
	d=inf.read()
	if '\r' in d:
		d = d.replace('\r','')
	d=d.split('LM=')
	data={}
	for l in d:
		if l == '':
			continue
		else:
			temp=[]
			coord=l[l.find('\n')+1:l.find('IMAGE=')]
			bline=coord.strip().split()
			for e in range(len(bline)):
				temp.append(bline[e])
			ID=l[l.find('ID='):].split()
			data[ID[0][3:]]=temp

	for k,v in data.iteritems():
		fou.write('>'+str(k)+';')
		for el in range(len(v)):
			fou.write(str(v[el])+';')
		fou.write('\n')
	inf.close()
	fou.close()
	return foun

def R_read_table(infile):
	line=open(infile).readline()
	r('m <- read.table("%s", sep=";", row.names=1)'%(infile)) 
	if line.strip()[-1] == ';':
		r('m<- m[1:dim(m)[2]-1]')	
	r('m<- as.matrix(m)')
	r('shnames<- row.names(m)')
	#os.system('mv %s %s'%(infile,infile[:-3]+'_before.gm'))

def R_table2array2GPA2gm(prefix,dim, scale):
	r('A<-arr(as.matrix(m),%d)'%(dim))
	r('pc <- procGPA(A, scale=%s)'%(scale))
	#r('png("%s_procrustes.png")'%(prefix))
	#r('plotshapes(pc$rotated)')
	#r('dev.off()')
	r('gm <- A2GM(pc$rotated, %d,rownames=shnames)'%(dim))
	r('mshape<-pc$mshape')
	r('dim(mshape)[3]<-1')
	r('save(A,pc, gm, arr,A2GM, file ="proc.R" )')
	#r('png("%s_procCorrMat.png")'%(prefix))
	#r('hist(cor(gm))')
	#r('dev.off()')
	r('write.table(gm,file="%s",sep=";",row.names=TRUE, col.names=FALSE)'%(prefix+'.gm'))
	r('write.table(A2GM(mshape,%d), file="mshape.gm", sep=";",row.names=TRUE, col.names=FALSE)'%(dim))

def main(prefix,dim,extension='gm',scale='FALSE',tps=False):
	R_functions()
	if tps:
		tps2gm(prefix)
	R_read_table(prefix+'.'+extension)
	if extension != 'gm' and scale == 'FALSE':
		R_table2array2GPA2gm(prefix+'_noScaled',dim, scale)
	elif extension != 'gm' and scale == 'TRUE':
		R_table2array2GPA2gm(prefix+'_Scaled',dim, scale)
	else:
		R_table2array2GPA2gm(prefix,dim, scale)
	
# End of definitions###############################################################################################

# Aplication of the code ##########################################################################################
if __name__ == '__main__':
	if len(sys.argv) == 1 or '-help' in sys.argv:
			print 'usage procrustes.py [prefix] [dim] [option]'
			print 'Options:'
			print '\t-tps : Input is a tps file ( Default: gmfile )'
			print '\t-scale: Will NOT use scale'
			sys.exit()
	
	##Default Parameters ###################################################
	prefix = sys.argv[1]
	dim = int(sys.argv[2])
	tps = False
	scale='TRUE'
	## #####################################################################
	
	for arg in sys.argv:
			if arg == '-tps':
				tps=True
			if arg == '-scale':
				scale="FALSE"
	
	main(prefix,dim,scale=scale,tps=tps)