import sys

print 'Usage: <prefix> <FDM filename>'
prefix = sys.argv[1]
fdmfile=sys.argv[2]
lnmk=open(prefix+'.landmarks')
fdms = open(fdmfile).read().split('\n')
fd=[]
for i in fdms:
	if i == '':
		continue
	else:
		fd.append(i)

nlnmk=open(prefix+'.fdm','w')
c=0
for line in lnmk:
	if line.startswith('>'):
		nlnmk.write(line)
	else:
		l = line.strip()
		nlnmk.write(l+'\t'+fd[int(line.split('\t')[0])]+'\n')
nlnmk.close()
