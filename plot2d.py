#!/usr/bin/python

#importing bit####################################################################################################
import matplotlib.pyplot as plt
from Moduler import Read_GM
import sys
#End importing####################################################################################################

#Some definitions##################################################################################################
def mask_graphcluster2GM(prefix):
	try:
		graph = open(prefix+'.final.graphcluster').readline().strip().split()
	except:
		graph = open(prefix+'.graphcluster').readline().strip().split()
	modulesindex = {}
	for i in range(len(set(graph))):
		tuples = ()
		for j in range(len(graph)):
			if graph[j] == list(set(graph))[i]:
				tuples+= (j,)
		modulesindex[list(set(graph))[i]]=tuples
		
	return modulesindex

def strip_coordGM(prefix, data, moduletuple):
	x=[]
	y=[]
	for i in moduletuple:
		for linex in data[0]:
			x.append(linex[i])
		for liney in data[1]:
			y.append(liney[i])
	return x,y

# End of definitions###############################################################################################

# Aplication of the code ##########################################################################################
prefix = sys.argv[1]

modulesindex=mask_graphcluster2GM(prefix)
data, sample_size=Read_GM(prefix,2)
fig = plt.figure(1)
for k,v in modulesindex.iteritems():
	x,y = strip_coordGM(prefix,data,v)
	ax = fig.add_subplot(111)
	ax.plot(x,y, marker='o', linestyle='None', label=str(k))


plt.axis('equal')
plt.legend()
if '-show' in sys.argv:
	plt.show()
else:
	plt.savefig('%s.svg'%(prefix+'_plot'), dpi=600, format='svg')
	plt.close()
		
	
