#!/usr/bin/python
'''
Equi2Latex Copyright (C) 2012 Jose Sergio Hleap

This script (pure gradware) will take the equivalence file from the 
PCoA plot (if you use some other of my scripts) and turn it into a latex table, 
given a number of columns. the source file should have the following architechture:

SUS SCROFA_197 \t 1VAH \t 0 \n
SUS SCROFA_146 \t 1KXQ \t 1 \n
HOMO SAPIENS_190 \t 1KBK \t 2 \n
HOMO SAPIENS_160 \t 3DHP \t 3 \n

Field one: Species name and number after an underscore
Field two: pdb code for that entry
Field three: the index in which the plot was made or any other thing important to you

To execute this script:

python Equi2Latex.py [filename] [columns number]

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

Requires: Numpy
'''
import sys
from tl.rename.case import transform_lowercase, transform_sentence_case
from numpy import median

filename = sys.argv[1]
col = int(sys.argv[2]) #number of columns
f = open(filename).read().split('\n')
#parse the equivalences file
sps = {}
for e in f:
    if e == '' or e == '\t\t':
        continue
    else:
        bline = e.split('\t')
        sp = bline[0].split('_')[0]
        sp = sp.split()        
        if 'SP.' in bline[0]:
            sp = transform_sentence_case([sp[0]])[0]+' '+transform_lowercase([sp[1]])[0]
        else:
            sp = ' '.join(sp)#sp[0][0]+'. '+transform_lowercase([sp[1]])[0]
            
        pdb=bline[1]
        ind=bline[2]
        if not sp in sps:
            sps[sp]= []
            sps[sp].append((pdb,ind))
        else:
            sps[sp].append((pdb,ind))

#create the table
t = open(filename[:filename.rfind('.')]+'.latextable','w')
header ='\\begin{table}[ht]\n\t\\tiny\n\t\\caption{}\n\t\\begin{center}\n\\begin{tabular}{ ccc'+' cc'*(col-1)
header += '}\n\t\\toprule\n\t Species'+ '& PDB & Plot'*col +'\\\\\n'
header +='&code&equivalences'*col+'\\\\\n\midrule\n'
t.write(header)
for k,v in sps.iteritems():
    bigline = ''
    p = len(v)
    r=range(p)
    if len(r) <= col:
        m = 0
    rows=(p/col)+1
    if p == col:
        rows = 1
    #m=int(median(range(rows)))
    for ro in range(rows):
        tbw=r[ro*col:(ro*col)+col]
        d=col-len(tbw)
        if tbw:
            #if ro == m:
            if ro == 0:
                bigline+= k
            for el in tbw:
                bigline+= '&' + v[el][0] + '&' + v[el][1]
                
            if d != 0:
                for i in range(d):
                    bigline+='&&'
        bigline+='\\\\\n'
    bigline+='\\\\[-0.2in]\n'
    t.write(bigline)
bottom ='\t\\bottomrule\n\t\\end{tabular}\n\t\\end{center}\n\\end{table}'
t.write(bottom)
t.close()                
