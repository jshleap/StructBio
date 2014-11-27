#!/usr/bin/python
"""
classifyGM Copyright (C) 2012 Jose Sergio Hleap with greatly appreated contributions (almost all; svm module) from
Alex Safatli

Given gmfiles (each one containing elements belonging to the same group and therefore with 
an assigned label) create a Support Vector Machine model, classifying them. the gmfiles given
are treated as the training set. Use -h to see options.
A gmfile is a semicolon-delimited text file, where the first element is the name of the structure/shape.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

Requirements:
1) Utils folder, available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools
2) Numpy module
3) svmutil module

"""
from subprocess import Popen
from os.path import isfile
#from shutil import copyfileobj
from labblouin.PDBnet import PDBstructure
import numpy as np
from glob import glob as G

blosum62 = {
"C":{"C":9, "S":-1, "T":-1, "P":-3, "A":0,  "G":-3, "N":-3, "D":-3, "E":-4, "Q":-3, "H":-3, "R":-3, "K":-3, "M":-1, "I":-1, "L":-1, "V":-1, "F":-2, "Y":-2, "W":-2},
"S":{"C":-1,"S":4,  "T":1,  "P":-1, "A":1,  "G":0,  "N":1,  "D":0,  "E":0,  "Q":0,  "H":-1, "R":-1, "K":0,  "M":-1, "I":-2, "L":-2, "V":-2, "F":-2, "Y":-2, "W":-3},
"T":{"C":-1,"S":1,  "T":4,  "P":1,  "A":-1, "G":1,  "N":0,  "D":1,  "E":0,  "Q":0,  "H":0,  "R":-1, "K":0,  "M":-1, "I":-2, "L":-2, "V":-2, "F":-2, "Y":-2, "W":-3},
"P":{"C":-3,"S":-1, "T":1,  "P":7,  "A":-1, "G":-2, "N":-1, "D":-1, "E":-1, "Q":-1, "H":-2, "R":-2, "K":-1, "M":-2, "I":-3, "L":-3, "V":-2, "F":-4, "Y":-3, "W":-4},
"A":{"C":0, "S":1,  "T":-1, "P":-1, "A":4,  "G":0,  "N":-1, "D":-2, "E":-1, "Q":-1, "H":-2, "R":-1, "K":-1, "M":-1, "I":-1, "L":-1, "V":-2, "F":-2, "Y":-2, "W":-3},
"G":{"C":-3,"S":0,  "T":1,  "P":-2, "A":0,  "G":6,  "N":-2, "D":-1, "E":-2, "Q":-2, "H":-2, "R":-2, "K":-2, "M":-3, "I":-4, "L":-4, "V":0,  "F":-3, "Y":-3, "W":-2},
"N":{"C":-3,"S":1,  "T":0,  "P":-2, "A":-2, "G":0,  "N":6,  "D":1,  "E":0,  "Q":0,  "H":-1, "R":0,  "K":0,  "M":-2, "I":-3, "L":-3, "V":-3, "F":-3, "Y":-2, "W":-4},
"D":{"C":-3,"S":0,  "T":1,  "P":-1, "A":-2, "G":-1, "N":1,  "D":6,  "E":2,  "Q":0,  "H":-1, "R":-2, "K":-1, "M":-3, "I":-3, "L":-4, "V":-3, "F":-3, "Y":-3, "W":-4},
"E":{"C":-4,"S":0,  "T":0,  "P":-1, "A":-1, "G":-2, "N":0,  "D":2,  "E":5,  "Q":2,  "H":0,  "R":0,  "K":1,  "M":-2, "I":-3, "L":-3, "V":-3, "F":-3, "Y":-2, "W":-3},
"Q":{"C":-3,"S":0,  "T":0,  "P":-1, "A":-1, "G":-2, "N":0,  "D":0,  "E":2,  "Q":5,  "H":0,  "R":1,  "K":1,  "M":0,  "I":-3, "L":-2, "V":-2, "F":-3, "Y":-1, "W":-2},
"H":{"C":-3,"S":-1, "T":0,  "P":-2, "A":-2, "G":-2, "N":1,  "D":1,  "E":0,  "Q":0,  "H":8,  "R":0,  "K":-1, "M":-2, "I":-3, "L":-3, "V":-2, "F":-1, "Y":2,  "W":-2},
"R":{"C":-3,"S":-1, "T":-1, "P":-2, "A":-1, "G":-2, "N":0,  "D":-2, "E":0,  "Q":1,  "H":0,  "R":5,  "K":2,  "M":-1, "I":-3, "L":-2, "V":-3, "F":-3, "Y":-2, "W":-3},
"K":{"C":-3,"S":0,  "T":0,  "P":-1, "A":-1, "G":-2, "N":0,  "D":-1, "E":1,  "Q":1,  "H":-1, "R":2,  "K":5,  "M":-1, "I":-3, "L":-2, "V":-3, "F":-3, "Y":-2, "W":-3},
"M":{"C":-1,"S":-1, "T":-1, "P":-2, "A":-1, "G":-3, "N":-2, "D":-3, "E":-2, "Q":0,  "H":-2, "R":-1, "K":-1, "M":5,  "I":1,  "L":2,  "V":-2, "F":0,  "Y":-1, "W":-1},
"I":{"C":-1,"S":-2, "T":-2, "P":-3, "A":-1, "G":-4, "N":-3, "D":-3, "E":-3, "Q":-3, "H":-3, "R":-3, "K":-3, "M":1,  "I":4,  "L":2,  "V":1,  "F":0,  "Y":-1, "W":-3},
"L":{"C":-1,"S":-2, "T":-2, "P":-3, "A":-1, "G":-4, "N":-3, "D":-4, "E":-3, "Q":-2, "H":-3, "R":-2, "K":-2, "M":2,  "I":2,  "L":4,  "V":3,  "F":0,  "Y":-1, "W":-2},
"V":{"C":-1,"S":-2, "T":-2, "P":-2, "A":0,  "G":-3, "N":-3, "D":-3, "E":-2, "Q":-2, "H":-3, "R":-3, "K":-2, "M":1,  "I":3,  "L":1,  "V":4,  "F":-1, "Y":-1, "W":-3},
"F":{"C":-2,"S":-2, "T":-2, "P":-4, "A":-2, "G":-3, "N":-3, "D":-3, "E":-3, "Q":-3, "H":-1, "R":-3, "K":-3, "M":0,  "I":0,  "L":0,  "V":-1, "F":6,  "Y":3,  "W":1},
"Y":{"C":-2,"S":-2, "T":-2, "P":-3, "A":-2, "G":-3, "N":-2, "D":-3, "E":-2, "Q":-1, "H":2,  "R":-2, "K":-2, "M":-1, "I":-1, "L":-1, "V":-1, "F":3,  "Y":7,  "W":2},
"W":{"C":-2,"S":-3, "T":-3, "P":-4, "A":-3, "G":-2, "N":-4, "D":-4, "E":-3, "Q":-2, "H":-2, "R":-3, "K":-3, "M":-1, "I":-3, "L":-2, "V":-3, "F":1,  "Y":2,  "W":11},
"X":{"X":8,"A":-0.9}
}
#constant dictionary with each aminoacid and its 1) hydrophobicity index at ph 7 (as proposed in
#Monera et.al. J. Prot. Sci. 1: 319:329 (1995)); 2) Van der Waals volume (as calculated by Richards
#J.Mol.Biol. 82:1-14,(1974)), and 3) charge @ ph 7
# the Unknowns are average of all others.
aa = {'A':[41,67,0],'C':[49,86,0],'D':[-55,91,-1],'E':[31,109,-1],'F':[100,135,0],'G':[0,48,0],\
      'H':[8,118,1],'I':[99,124,0],'K':[-23,135,1],'L':[97,124,0],'M':[74,124,0],'N':[-28,96,0],\
      'P':[-46,90,0],'Q':[-10,114,0],'R':[-14,148,1],'S':[-5,73,0],'T':[13,93,0],'V':[76,105,0],\
      'W':[97,163,0],'Y':[63,141,0],'X':[28,109,0.05]}


#normalize props... comment out

maxi = np.array([100.0,163.0,1.0])
mini = np.array([-55,48,-1])
norm = maxi - mini
keys = aa.keys()
for k in keys:
    aa[k]=list((np.array(aa.pop(k)) - mini)/norm)


def get_data_gm(prefix, gmfile,options):
    data=[]
    names = []
    if options.usesequence:
        land = readLandmarks(prefix)
    with open(gmfile) as GM:
        for line in GM:
            bl=line.strip().strip(';').split(';')
            names.append(bl[0])
            try:
                datum = [float(x.strip('"')) for x in bl[1:]]
            except:
                datum = [float(x.strip('"')) for x in bl[2:]]
            if options.usesequence:
                name = bl[0].strip('"').strip('>')
                stname= name[:name.find(':')]
                st = PDBstructure('%s/%s.pdb'%(options.usesequence,stname))
                fasta = list(st.ChainAsFASTA(stname[-1]))
                afas = []
                for x in land[stname]:
                    T=[]
                    T.extend([float(ord(fasta[int(x)]))]) # naive sequence
                    T.extend([float(blosum62[fasta[int(x)]][fasta[int(x)]])])#blosum sequenc
                    T.extend([float(blosum62[fasta[int(x)]]["A"])])#dist to ala
                    #T = list(((np.array(T)-min(np.array(T)))+1)/((max(np.array(T))-min(np.array(T)))+1))# normalize naive seq
                    #T.extend([float(aa[fasta[int(x)]][2])]) # only charge
                    T.extend([float(y) for y in aa[fasta[int(x)]]]) # include all properties
                    afas.extend(T)
                datum.extend(afas)            
            data.append(datum)
    return data, names

def readLandmarks(prefix):
    '''read landmarks file and return a dict'''
    land={}
    l = open(prefix+'.landmarks').read().split('>')
    for e in l:
        if e == '':
            continue
        else:
            bl = e.split('\n')
            land[bl[0]]=[]
            for i in bl[1:]:
                if i == '':
                    continue
                else:
                    bline= i.split()
                    land[bl[0]].append(bline[1])
    return land
    
def setFilesLabels(prefix, arg,options):
    data=[]
    lab=[]
    labels=[]
    labdict={}
    files=[]
    for i in range(0,len(arg),2):
        dat, names = get_data_gm(prefix, arg[i],options)
        files.append(arg[i])
        label=arg[i+1]
        for d in range(len(dat)): labels.append(label)
        data.extend(dat)
    s=list(set(labels))
    for i in labels:
        lab.append(s.index(i))
        labdict[float(s.index(i))]= i 
    return files, data, lab, labdict

def train(prefix, data, labels,options):
    if options.c_val:
        model = svm.svm_model(data,labels, kernel=svmutil.RBF,c=float(options.c_val))
    else:
        model = svm.svm_model(data,labels, kernel=svmutil.RBF)
        model.optimize(low=0.1, up=45.1, steps=0.1)
        print "Optimum C: %f"%(model.c)
    if options.predict:
        cv = svm.cross_validate(model,get_data_gm(prefix, options.predict,options))
    else:
        cv = svm.cross_validate(model)
    cv.perform()
    Fscore = cv.f1stats
    Cmatrix = cv.stats
    F = cv.proportion
    print 'Percentage of correct classifications in the cross-validation = %f'%(F*100)
    return cv

def featSel_n_useSeq(options,prefix,args):
    if options.predict:
        testdata, testnames = get_data_gm(prefix,options.predict,options)
        fname ='test.gm'
        f = open(fname,'w')
        for i in range(len(testdata)):
            line = testnames[i]
            line += ';' + ';'.join([str(x) for x in testdata[i]])
            f.write(line+'\n')
        f.close()
        options.predict = fname
    if options.featureselect and not options.usesequence:
        rs = Popen('Rscript ./featuresel.R %s.gm %s.cls 3 FALSE %s'%(prefix,prefix,options.predict),shell=True)
        rs.wait()
        if options.corr:
            rs1=Popen('Rscript ./GM2MDS.R %s.gm %s.cls %s'%(prefix,prefix,options.predict),shell=True)
            rs1.wait()        
        options.predict = options.predict[:-2]+'sv.gm'
        for i in range(0,len(args[1:]),2):
            args[i+1] = args[i+1][:-2]+'sv.gm'
    if options.featureselect and options.usesequence:
        files, data, labels, labdict = setFilesLabels(prefix,args[1:],options)
        nfile = open(prefix+'.seq.gm','w')
        #prefix = prefix+'.seq'
        for i in range(len(data)):
            nfile.write('>%s;%s\n'%(str(i),';'.join([str(x) for x in data[i]])))
        nfile.close()
        rs = Popen('Rscript ./featuresel.R %s.seq.gm %s.cls 3 FALSE %s'%(prefix,prefix,options.predict),shell=True)
        rs.wait()
        if options.corr:
            rs1=Popen('Rscript ./GM2MDS.R %s.seq.gm %s.cls %s'%(prefix,prefix,options.predict),shell=True)
            rs1.wait()
        options.predict = options.predict[:-2]+'sv.gm'
        for i in range(0,len(args[1:]),2):
            if not 'sv' in args[i+1]:
                args[i+1] = args[i+1][:-2]+'sv.gm'
            options.usesequence = None    
    
def main(options,args):
    spline=''
    snline=''
    fsline=''    
    prefix=args[0]
    featSel_n_useSeq(options,prefix,args)
    files, data, labels, labdict = setFilesLabels(prefix,args[1:],options)
    if options.predict:
        files.append(options.predict)
        #if not isfile(prefix+'.svm'):
        cv = train(prefix, data, labels,options)
        #else:
         #   model = svm.svm_model()
          #  model.load(prefix+'.svm')
            #cv = svm.cross_validate(model)
            #cv = svm.cross_validate(model.model,get_data_gm(prefix, options.predict,options))
        print "Prediction:"
        pred = cv.model.predict(cv.test[0])
        fout = open(prefix+'.prediction','w')
        for i in range(len(pred)):
            line = cv.test[1][i].strip('"') + '\t' + labdict[pred[i]]
            fout.write(line+'\n')
            print line
        fout.close()
        log = open(prefix+'.log','w')
        print "Optimum C: %f"%(cv.model.c)
        for k,v in cv.f1stats.iteritems():
            spline += labdict[k] + ': %f'%(v[0]) + '\t'
            snline += labdict[k] + ': %f'%(v[1]) + '\t'
            fsline += labdict[k] + ': %f'%(v[2]) + '\t'
        whole = "Performance in the cross-validation:\n Specificity :\n"+spline+'\nSensitivity :'
        whole += "\n"+snline+'\nF-score :\n'+fsline+"\nOptimum C: %f"%(cv.model.c)
        log.write(whole)
        log.close()
        
    else:
        m = train(prefix, data, labels,options)
        print "Optimum C: %f"%(m.model.c)
        for k,v in m.f1stats.iteritems():
            spline += labdict[k] + ': %f'%(v[0]) + '\t'
            snline += labdict[k] + ': %f'%(v[1]) + '\t'
            fsline += labdict[k] + ': %f'%(v[2]) + '\t'
                
        print 'Specificity :'
        print spline
        print 'Sensitivity :'
        print snline        
        print 'F-score :'
        print fsline

if __name__ == '__main__':
    import svmutil, optparse
    from labblouin import svm
    opts = optparse.OptionParser(usage='%prog [options] prefix GMfile1 label1 [GMfilename2 label2]...\n'\
                                 ' The prefix will be use to save and load the model, which will be'\
                                 ' <prefix>.svm. It will be also used to load the gm file and the cls'\
                                 ' file if feature selection is selected, as well as the landmarks file'\
                                 ' if usesequence is selected. The Rscript featuresel.R should be in the'\
                                 ' folder.')
    opts.add_option('--predict','-p', action = "store", default='',
                    help="Use an existing trained model to test a file for predictions. With"\
                    " this flag you have to provide the name of a test file in gm format (semi"\
                    "colon-delimited csv file') of the aligned structures (aligned with the "\
                    "training set)Default: None (train first)")
    opts.add_option('--featureselect','-f', action = "store_true", default=False,
                    help=" Use univariate Kuskall-Wallis ANOVA per variable to select significant variables. Default: No ")    
    opts.add_option('--usesequence','-s', action = "store", default=None,
                    help="When a protein dataset, use the sequence information to cluster. You have to provide"\
                    "the path where the pdbs are to be found. Default: No ")  
    opts.add_option('--c_val','-c', action = "store", default=None,
                        help="Provide a known optimum c value for the svm. This avoid the optimization step"\
                        " Default: No (Optimize) ")
    opts.add_option('--corr','-C', action = "store", default='',
                            help="Wheather to use or not correlation of the samples giving the variables"\
                            " Default: No (Raw, option: T) ")    
    options, arg = opts.parse_args()

    main(options,arg)

