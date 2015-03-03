#!/usr/bin/python
'''
This script will calculate the correlation of the form deformtation residuals and either entropy
or evolutionary trace scrores
'''
#importing bit##########################################################################################
import sys, optparse, os
from rpy2.robjects import r
import numpy as np
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from labblouin.FASTAnet import FASTAstructure as F
from numpy import log
from FoldX_runner import *
import cPickle as P
from matplotlib import rc

def format_fasta(prefix):
    fasta = F(prefix+'.fasta')
    for i in fasta.orderedSequences:
        on=i.name
        nn=on[on.find(':')+1:]
        i.name = nn
    fasta.writeFile(prefix+'.fasta')
    
def get_FD(prefix,args):
    #process betaeq
    if not os.path.isfile(prefix+'.betaeq'):
        pgm=Popen('python /home/jshleap/LabBlouin/code/GM/protGM.py %s'\
                  ' -FD=%s -m=%s'%(prefix,args[1],args[1]),shell=True)
    betas={}
    fil=open(prefix+'.betaeq')
    for line in fil:
        if line == '':
            continue
        else:
            bline=line.strip().split('\t')
            betas[bline[0]]=bline[-1]
    
    return betas

def process_landmarks(prefix,chain):
    #process landmarks
    eq={}
    fi = open(prefix+'.landmarks').read().split('\n>')
    for st in fi:
        if st == '':
            continue
        elif st.startswith('>'):
            st=st[1:]
        else:
            if st.startswith(chain):
                bst=st.split('\n')
                for e in bst[1:]:
                    e=e.strip().split('\t')
                    if e == '':
                        continue
                    eq[e[0]]=e[1]    
    return eq

def get_tFDM(prefix,eq):
    tfd={}
    if not os.path.isfile(prefix+'.tradFDM'):
        t = Popen('Rscript /home/jshleap/LabBlouin/code/GM/tradFDM.R %s.gm 3'%(prefix),
                  shell=True)
        t.wait()
    f=open(prefix+'.tradFDM').read().split('\n')
    for l in range(len(f)):
        if f[l] == '':
            continue
        else:
            tfd[eq[str(l)]]=f[l]
    return tfd

def get_entropies(prefix):
    # process web logo
    if not os.path.isfile(prefix+'.hom.fasta'):
        hom = Popen('/home/jshleap/LabBlouin/code/GM/HomologousFasta.py %s'%(prefix),shell=True)
        hom.wait()
    if not os.path.isfile(prefix+'..weblogo.data'):
        wl=Popen('weblogo -f %s.hom.fasta -o %s.weblogo.data -F logodata'%(prefix,prefix),shell=True)
        wl.wait()
    logo={}
    f = open(prefix+'.weblogo.data')
    for l in f:
        if l.startswith('#'):
            continue
        elif l == '':
            continue
        elif l == '\n':
            continue
        else:
            bl=l.strip().split('\t')
            logo[bl[0].strip()]=bl[21]
    return logo

def AppAli4ET(prefix,extension,pdb):
    if not os.path.isfile(prefix+'.msf'):
        c=Popen('clustalw -CONVERT -INFILE=%s -OUTPUT=GCG'%(prefix+'.'+extension),
                shell=True)
        c.wait()
    rn=Popen('cat %s.msf | grep %s| grep "Name:"'%(prefix,pdb[:4]),shell=True,
             stdout=PIPE,stderr=PIPE)
    o,e = rn.communicate()
    refname=o.split()[1]
    return refname

def parse_ET_result(prefix):
    ET={}
    f=open(prefix+'.ranks_sorted')
    for line in f:
        if line.startswith('%'):
            continue
        if line == '':
            continue
        if line == '\n':
            continue
        else:
            bline=line.strip().split()
            ET[bline[1]]=bline[6]
    return ET      

def get_ETs(prefix,extension,pdb):
    if not os.path.isfile(prefix+'.tre'):
        tree = Popen('FastTree -wag %s.fasta > %s.tre'%(prefix,prefix),shell=True)
        tree.wait()
    refname = AppAli4ET(prefix,extension,pdb)
    if not os.path.isfile(prefix+'.ranks_sorted'):
        et=Popen('wetc -p %s.msf -x %s %s -similarity -readtree %s.tre -c '\
                 '-o %s'%(prefix,refname,pdb+'.pdb',prefix,prefix),shell=True)
        et.wait()
    ET = parse_ET_result(prefix)
    return ET

def transform_small_values(xl,yl):
    if max(xl) < 0.1:
        xl = [round((x/max(xl)),2) for x in xl]
    else:
        xl = xl
    if max(yl) < 0.1:
            yl = [round((y/max(yl)),2) for y in yl]
    else:
        yl=yl
    return xl, yl

def plot_correlation(xlist,ylist,opts,string,fout,n):
    x, data_y = transform_small_values(xlist,ylist)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, data_y, 'b.')
    coefs = np.lib.polyfit(x, data_y, 1)
    fit_y = np.lib.polyval(coefs, x)
    ax.plot(x, fit_y, 'r-')
    #ax.text(max(x)-0.5,max(data_y)-1,string)
    ax.annotate(string,  xy=(max(x), min(data_y)), xycoords='data',
                #xytext=(0, max(data_y)-5) , textcoords='offset points',
                xytext=(1, 1), textcoords='axes fraction',
                bbox=dict(boxstyle="round", fc="0.9"), horizontalalignment='right', 
                verticalalignment='top')
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    if opts.entropy:
        plt.xlabel('Entropy (nats)')
        if opts.PDB:
            plt.ylabel('Evolutionary Trace')
        elif opts.FD:
            plt.ylabel('Form Difference')
        elif opts.cent:
            if n == 0:
                plt.ylabel('Eigenvalue centrality')
            elif n == 1:
                plt.ylabel('Betweeness centrality')
            elif n == 2:
                plt.ylabel('Closeness centrality')
            else:
                plt.ylabel('Degree centrality')            
        elif opts.foldx:
            #plt.ylabel('Average delta energy')
            plt.ylabel(r'Average $\Delta  \Delta$G')
        elif opts.tfd:
            plt.ylabel('Form Difference')

    elif opts.PDB:
        plt.xlabel('Evolutionary Trace')
        if opts.FD:
            plt.ylabel('Form Difference')
        elif opts.cent:
            if n == 0:
                plt.ylabel('Eigenvalue centrality')
            elif n == 1:
                plt.ylabel('Betweeness centrality')
            elif n == 2:
                plt.ylabel('Closeness centrality')
            else:
                plt.ylabel('Degree centrality')            
        elif opts.foldx:
            #plt.ylabel('Average delta energy')
            plt.ylabel(r'Average $\Delta  \Delta$G')
        elif opts.tfd:
            plt.ylabel('Form Difference')        

    elif opts.cent and not opts.FD or not opts.tfd:
        if n == 0:
            plt.ylabel('Eigenvalue centrality')
        elif n == 1:
            plt.ylabel('Betweeness centrality')
        elif n == 2:
            plt.ylabel('Closeness centrality')
        else:
            plt.ylabel('Degree centrality')
    
    elif opts.cent and (opts.FD or opts.tfd):
        plt.ylabel('Form Difference')
        if n == 0:
            plt.xlabel('Eigenvalue centrality')
        elif n == 1:
            plt.xlabel('Betweeness centrality')
        elif n == 2:
            plt.xlabel('Closeness centrality')
        else:
            plt.xlabel('Degree centrality')
    else:
        plt.ylabel('Form Difference')
    #plt.show()       
    plt.savefig(fout)
    plt.close()

def parse_centralities(prefix,pdb):
    eige={}
    betw={}
    clos={}
    degr={}
    fi = open(prefix+'.centrality').read().split('\n>')
    for l in fi:
        if not pdb in l:
            continue
        bline=l.split('\n')
        for e in bline:
            if e == '':
                continue
            elif 'The organization' in e:
                continue
            elif pdb in e:
                continue
            else:
                bl=e.strip().split('\t')
                eige[bl[1]]=bl[3]
                betw[bl[1]]=bl[4]
                clos[bl[1]]=bl[5]
                degr[bl[1]]=bl[6]
    return eige,betw,clos, degr

def get_centralities(prefix,pdb):
    if not os.path.isfile(prefix+'.centrality'):
        mod = Popen('python /home/jshleap/LabBlouin/code/modularity/modulerV2.py '\
                    '%s -p 0.8 -c -n 999 -r -g -f'%(prefix),shell=True)
        mod.wait()
    eige,betw,clos, degr = parse_centralities(prefix,pdb)
    return [eige,betw,clos,degr]

def vectorNlist(xdict,ydict,eqdict,options,fout,n):
    # create the vectors for R
    x, xl, y, yl = 'c(' , [] , 'c(' , []
    for k,v in eqdict.iteritems():
        b=ydict[v]
        yl.append(float(b))
        y+=str(b)+','
        if options.entropy:
            en=xdict[str(int(k)+1)]
            xl.append(float(en))
            x+=str(en)+','
        else:
            e=xdict[v]
            xl.append(float(e))
            x+=str(e)+','
    y , x = y[:-1]+')' , x[:-1]+')'
    r('e <- %s'%(x))
    r('f<-%s'%(y))
    r('C<-cor.test(e,f)')
    pval=r('C$p.value')[0]
    cor=r('C$estimate')[0]
    string='''
    r: %s
    p-value: %s
    '''%(round(cor,3),round(pval,5))
    print  string
    plot_correlation(xl,yl,options,string,fout,n)
    string=''

def get_energies(prefix,pdb):
    if not os.path.isfile(prefix+'_energies.pckl'):
        FoldXrunner(prefix,pdb)
        nrg = energies2dict(prefix)
    else:
        nrg = P.load(open(prefix+'_energies.pckl'))
    return nrg

def main(args,options):
    prefix=args[0]
    chain=args[1]
    format_fasta(prefix)
    eq=process_landmarks(prefix,chain)
    if options.entropy:
        logo = get_entropies(prefix)
        xdict = logo
        na='Entropy_vs_'
        if options.PDB:
            na += 'ET.png'
            ydict = get_ETs(prefix,options.extfile,options.PDB)

        elif options.FD:
            na += 'FD.png'
            ydict = get_FD(prefix,args)

        elif options.cent:
            C = get_centralities(prefix,chain)
            names=['Eigen','Betw','Clos', 'Degree']
            for i in range(len(C)):
                if names[i] == 'Eigen':
                    print na.replace('_',' '),'Eigenvalue centrality:'
                elif names[i] ==  'Betw':
                    print na.replace('_',' '),'Betweeness centrality:'
                elif names[i] ==  'Clos':
                    print na.replace('_',' '),'Closeness centrality:'
                else:
                    print na.replace('_',' '),'Degree centrality:'
                vectorNlist(xdict,C[i],eq,options,na+names[i]+'.png',i)            
        elif options.foldx:
            ydict = get_energies(prefix,chain)
            na += 'foldx.png'
        elif options.tfd:
            ydict=get_tFDM(prefix,eq)
            na += 'tFD.png'

    elif options.PDB:
        ET=get_ETs(prefix,options.extfile,options.PDB)
        xdict = ET
        na='ET_vs_'
        if options.FD:
            na += 'FD.png'
            ydict = get_FD(prefix,args)
        
        
        elif options.tfd:
            ydict=get_tFDM(prefix,eq)
            na += 'tFD.png'        

        elif options.cent:
            C = get_centralities(prefix,chain)
            names=['Eigen','Betw','Clos', 'Degree']
            for i in range(len(C)):
                if names[i] == 'Eigen':
                    print na.replace('_',' '),'Eigenvalue centrality:'
                elif names[i] ==  'Betw':
                    print na.replace('_',' '),'Betweeness centrality:'
                elif names[i] ==  'Clos':
                    print na.replace('_',' '),'Closeness centrality:'
                else:
                    print na.replace('_',' '),'Degree centrality:'
                vectorNlist(xdict,C[i],eq,options,na+names[i]+'.png',i)

        elif options.foldx:
            ydict = get_energies(prefix,chain)
            na += 'foldx.png'        
        
    elif options.cent:
        C = get_centralities(prefix,chain)
        names=['Eigen','Betw','Clos', 'Degree']
        if options.FD:
            ydict = get_FD(prefix,args)
            na = 'FD.png'

        elif options.foldx:
            ydict = get_energies(prefix,chain)
            na = 'foldx.png'

        elif options.tfd:
            ydict=get_tFDM(prefix,eq)
            na = 'tFD.png'        

        for i in range(len(C)):
            if names[i] == 'Eigen':
                print na.replace('_',' '),'Eigenvalue centrality:'
            elif names[i] ==  'Betw':
                print na.replace('_',' '),'Betweeness centrality:'
            elif names[i] ==  'Clos':
                print na.replace('_',' '),'Closeness centrality:'
            else:
                print na.replace('_',' '),'Degree centrality:'
            vectorNlist(C[i],ydict,eq,options,names[i]+'_vs_'+na,i)
    
    elif options.foldx:
        xdict = get_energies(prefix,chain)
        na = 'foldx_vs_'
        if options.FD:
            ydict = get_FD(prefix,args)
            na+='FD.png'
        elif options.tfd:
            ydict = get_tFDM(prefix,eq)    
            na+= 'tFD.png'

    elif options.tfd and options.FD:
            xdict = get_FD(prefix,args)
            na='FD_vs_tFD.png'
            ydict = get_tFDM(prefix,eq)    
    
    if not options.cent:
            vectorNlist(xdict,ydict,eq,options,na,1)    
'''
    elif options.FD:
        xdict = get_FD(prefix,args)
        na='FD_vs_foldx.png'
        ydict = get_energies(prefix,chain)
        
    elif options.tfd:
        xdict=get_tFDM(prefix,eq)
        na = 'tFD_vs_foldx.png'
        ydict = get_energies(prefix,chain)

    if not options.cent:
        vectorNlist(xdict,ydict,eq,options,na,1)'''
# End of definitions##################################################################

# Aplication of the code #############################################################
if __name__ == "__main__":
    # Command line input #############################################################
    opts = optparse.OptionParser(usage='%prog <prefix> <pdbcode> [options]')
    opts.add_option('-e', '--entropy',dest='entropy', action="store_true",
	                default=False, help='Use entropy as X for the correlation test.')
    opts.add_option('-t', '--ET',dest='PDB', action="store", default=False, 
                    help='Use evolutionary trace as X for the correlat'\
                    'ion test. A PDB code where the evolutionary trace want to be displayed'\
                    ' has to be provided and must match the pdb file in folder minus the ".PDB".'\
                    ' This will need a <prefix>.tre, and an alignent file.')
    opts.add_option('-a','--AliFileext', dest='extfile', action='store', default='fasta',
                    help='The extension of the alignment file. Is Fasta by default.')
    opts.add_option('-c', '--centrality',dest='cent', action="store_true", default=False,  
                    help='Use centrality instead of FD to test correlation with entropy '\
                    'and evolutionary trace.')
    opts.add_option('-f', '--FD',dest='FD', action="store_true", default=False,  
                        help='Use Form Difference residuals (FD) as the Y component in '\
                        'the correlation analyses or plot. This has to be used with -e, -t or -c.')   
    opts.add_option('-x', '--FoldX',dest='foldx', action="store_true", default=False,  
                            help='Use difference in energies between the wild type (pdb) and mutations'\
                            ' to alanine ')
    opts.add_option('-r', '--tFD',dest='tfd', action="store_true", default=False,  
                                help='Compute the traditional univariate FDM.')     
    options, args = opts.parse_args()
    
    main(args,options)




