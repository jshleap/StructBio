#!/usr/bin/python
'''
This script will walk thoughout the Homstrad dataset and verify if the three binary files are there and will 
plot GM scaled vs GM not scale, and GM not scaled vs MATT alignment
'''
#importing bit####################################################################################################
import os,sys,fnmatch
from glob import glob
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import Subplot
from pickle import load
from subprocess import Popen
from scipy import mean, std
from scipy.stats import pearsonr
# End importing####################################################################################################

#Some definitions##################################################################################################
def get_binaries(PATH=os.getcwd()):
	'''
	Given a PATH to a directory, will fetch the binaries and return them. Also will check if there is
	a landmark file and return a boolean
	'''
	landmark=False
	binaries = glob(PATH+'/*.bin')
	if len(binaries) != 0:
		landmark = os.path.isfile(PATH+'%s.landmarks'%(PATH[PATH.rfind('/'):]))
	return binaries, landmark

def solve_problem(path):
	'''
	will solve problem of len(gmns) vs len(matt)
	'''
	i=path
	rmpdbs = Popen('rm %s/*.pdb'%(i), shell=True)
	rmpdbs.wait()
	greplist = Popen('grep "%s" LIST > %s/LIST'%(i[i.rfind('/'):],i), shell=True)
	greplist.wait()
	rerun = Popen('python /home/jshleap/LabBlouin/code/GM/MATTresultsRetriever.py -path=%s -single'%(i), shell=True)
	rerun.wait()
	Popen('mv %s* *.pdb ./%s'%(i[i.rfind('/')+1:],i[i.rfind('/')+1:]),shell=True)


def walker(PATH_to_HOMSTRAD_folder= os.getcwd(), gmns=[],gms=[],matt=[],avmatt=False):
	'''
	This function will walk thoughout the HOMSTRAD folder, fetch the binaries, unpickle them and return an extended list of it
	'''
	errors=open('errors.txt','w')
	allroots=[]
	notok=[]
	ok=True
	for root, dirs, files in os.walk(PATH_to_HOMSTRAD_folder):
		if not root in allroots:
			allroots.append(root)
		if root == PATH_to_HOMSTRAD_folder:
			continue
		else:
			binaries, landmark = get_binaries(root)
			if not len(binaries)==3:
				errors.write('Not enough binaries in %s\n'%(root[root.rfind('/'):]))
				print 'Error found'
			elif not landmark:
				errors.write('Not landmarkfile in %s\n'%(root[root.rfind('/'):]))
				print 'Error found'
			else:
				for b in binaries:
					temp = load(open(b))
					if '_noscaled_' in b:
						if avmatt:
							gmns.append(mean(temp))
						else:
							gmns.extend(temp)
					elif '_scaled_' in b:
						if avmatt:
							gms.append(mean(temp))
						else:
							gms.extend(temp)
					else:
						if avmatt:
							matt.append(mean(temp))
						else:
							matt.extend(temp)
							
					for e in temp:
						if int(e) >= 5:
							errors.write('Outlier in RMSD (%f) of %s dataset\n'%(e,b))
							print 'Error found'
					temp=[]
			if len(gmns) != len(gms):
				errors.write('Not the same number of RMSDs between GMns and GMs in %s dataset\n'%(b[b.rfind('/')+1:b.rfind('_')]))
			if len(gmns) != len(matt):
				ok=False
				errors.write('Not the same number of RMSDs between GMns and MATT in %s dataset\n'%(b[b.rfind('/')+1:b.rfind('_')]))		
				notok.append(root)
				solve_problem(root)
	errors.close()
	return gmns,gms,matt,ok,notok

def parse_txt(txtfile):
	'''
	will parse the matt's TXT file and return the core residues used in the aln
	'''
	f = open(txtfile)
	for line in f:
		if not line.startswith('Core Residues: '):
			continue
		else:
			bline = line.split(': ')
			res = bline[1]
	clean = Popen('rm %s'%(txtfile),shell=True)
	clean.wait()
	return float(res)

def averagediffVScoresize(PATH_to_HOMSTRAD_folder= os.getcwd(),SABmark=False):
	'''
	Will get the diferences in average RMSD and plot them against core size of the aliged proteins
	'''
	coreRes=[]
	Avdiff=[]
	if SABmark:
		directories = fnmatch.filter(os.listdir(PATH_to_HOMSTRAD_folder),'group*')
	else:
		directories=[]
		for root, dirs, files in os.walk(PATH_to_HOMSTRAD_folder):
			if root == PATH_to_HOMSTRAD_folder:
				continue
			else:
				directories.append(root[root.rfind('/'):])
				
	for d in directories:
		if not os.path.isfile(d+d[d.rfind('/'):]):
			if not SABmark:
				txt = Popen('wget http://groups.csail.mit.edu/cb/matt/homstrad%s.txt'%(d),shell=True)
				txt.wait()
			else:
				txt = Popen('wget http://groups.csail.mit.edu/cb/matt/sabmark/%s.txt'%(d),shell=True)
				txt.wait()				
		coreRes.append(parse_txt(d+'.txt'))
		binaries, landmark = get_binaries(os.getcwd()+'/'+d)
		count=0
		for b in binaries:
			temp = load(open(b))
			if '_noscaled_' in b:
				avns=mean(temp)
				count+=1
			elif '_scaled_' in b:
				continue
			else:
				avmatt=mean(temp)
				count += 1
				
			if count == 2:
				Avdiff.append(avmatt - avns)
	gms=[]
	return Avdiff, gms , coreRes
				

	
def scatter_plots(gmns=[],gms=[],matt=[],average=False, avmatt=False):
	'''
	will plot to scatters of the RMSDs of GMnoscaled vs GMscaled and GMnoscaled vs MATT
	'''
	if gms:
		r,pval=pearsonr(gmns,gms)
		line=r'''
		$r$=%.2f, pval=%.4f
		$R^2$=%.2f'''%(r,pval,r**2)
		fig = plt.figure()
		ax = Subplot(fig, 1,1,1)
		fig.add_subplot(ax)
		ax.axis["right"].set_visible(False)
		ax.axis["top"].set_visible(False)
		ax.axis([min(gmns), max(gmns)+0.5,min(gms),max(gms)+0.5])
		plt.plot(range(int(max(gmns))+3),'r')
		plt.scatter(gmns,gms)
		ax.annotate(line,xy=(min(gmns)+0.5,max(gms)), xytext=(min(gmns)+0.5,max(gms)))
		if avmatt:
			plt.xlabel("Average Procrustes RMSD without scaling")
			plt.ylabel("Average Procrustes RMSD with scaling")
			plt.savefig('AvGMnsvsAvGMs.png')
		else:
			plt.xlabel("Procrustes RMSD without scaling")
			plt.ylabel("Procrustes RMSD with scaling")
			plt.savefig('GMnsvsGMs.png')#, transparent=True)
		
	
	r,pval=pearsonr(gmns,gms)
	line=r'''
	$r$=%.2f, pval=%.4f
	$R^2$=%.2f'''%(r,pval,r**2)	
	fig2 = plt.figure()
	ax2 = Subplot(fig2, 1,1,1)
	fig2.add_subplot(ax2)
	ax2.axis["right"].set_visible(False)
	ax2.axis["top"].set_visible(False)
	ax2.axis([min(matt),max(matt)+0.5, min(gmns), max(gmns)+0.5])
	plt.plot(range(int(max(gmns))+3),'r')
	plt.scatter(matt,gmns)
	ax2.annotate(line,xy=(min(gmns)+0.5,max(gms)), xytext=(min(gmns)+0.5,max(gms)))
	if average:
		plt.xlabel("Number of Core residues")
		plt.ylabel("Difference of average RMSD (mean(MATT) - mean(procrustes))")
		plt.savefig('diffvscore.png')
	elif avmatt:
		plt.xlabel("Average MATT RMSD")
		plt.ylabel("Average Procrustes RMSD without scaling")
		plt.savefig('AvMATTvsAvgmns.png')
	else:
		plt.xlabel("MATT alignment RMSD")
		plt.ylabel("Procrustes RMSD without scaling")
		plt.savefig('MattvsGMns.png')#, transparent=True)
	
	
# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	avmatt=False
	average=False
	SABmark=False
	if '-extend' in sys.argv:
		gmns,gms,matt,ok,notok = walker(os.getcwd(),[],[],[],avmatt)
	elif '-average' in sys.argv:
		if '-SABmark' in sys.argv:
			SABmark=True
		gmns,gms,matt = averagediffVScoresize(os.getcwd(),SABmark)
		average=True
	elif '-avmatt' in sys.argv:
		avmatt=True
		gmns,gms,matt,ok,notok = walker(os.getcwd(),[],[],[],avmatt)
	
	scatter_plots(gmns,gms,matt,average,avmatt)