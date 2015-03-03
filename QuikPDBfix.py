import sys,os
'''
quick naive fix of pdb ordeing
'''
pdb = open(sys.argv[1]).read().strip().split('\n')
os.rename(sys.argv[1],sys.argv[1]+'.bkup')
npdb= open(sys.argv[1],'w')
ordered={}
for l in pdb:
    bl=l.split()
    ordered[int(bl[1])]=l

sortedkeys=ordered.keys()
sortedkeys.sort()
for k in sortedkeys:
    npdb.write(ordered[k]+'\n')

npdb.close()