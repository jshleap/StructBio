#!/usr/bin/python
"""
landmark2latex Copyright (C) 2012 Jose Sergio Hleap

Really simple script to transform a landmark file into latex format (yes I know, stupid, 
but somre reviwers ask or it)

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

"""


import sys

def parse_matt_txt(prefix,md):
    d={}
    try:
        f=open(prefix+'.txt').read()
        f= f.split('Reference structure:')[1].split('\n\n')[1].split('\n')
        for e in f:
            if e == '':
                continue
            else:
                e=e.strip()
                code=e[:4]
                chain=e[5]
                sp=e[7:e.find(':')]
                key=e[-2]
                d[key]=(code,chain,sp)
    except:
        d['A']=('3GKH', 'A' , 'H._sapiens')
    
    return d
    
    
    
prefix = sys.argv[1]
try:
    sys.argv[2]
    md=True
    protein='NPC1'
except:
    md=False
    protein='$\\alpha$-amylase'
    
outf=open(prefix+'.supl','w')
land=open(prefix+'.landmarks').read().split('\n>')
memb=open(prefix+'.mapped.graphcluster').read().strip().split()
d = parse_matt_txt(prefix,md)

for f in land:
    if f == '':
        continue
    else:
        if not md:
            k=f[0]
        else:
            k='A'
        outf.write('\\begin{center} \n \\begin{longtable}{c c c c} \n \\caption{Residues membership'\
                   ' for the \\textit{ %s } %s (PDB code %s, chain %s)protein}\\\\ \n \\toprule \n'\
                   '\\textbf{Homology index}&\\textbf{Residue Index}&\\textbf{Residue abbreviation}'\
                   '&\\textbf{Module}\\\\ \n \\midrule \n \\endfirsthead \n \\multicolumn{4}{c} \n'\
                   '{\\tablename\\ \\thetable\\ -- \\textit{Continued from previous page}}\\\\ \n'\
                   '\\hdashline \n \\textbf{Homology index}&\\textbf{Residue Index}&\\textbf{Residue'\
                   ' abbreviation}&\\textbf{Module}\\\\ \n \\midrule \n \\endhead \n \\midrule \n'\
                   '\\multicolumn{4}{r}{\\textit{Continued on next page}}\\\\ \n \\endfoot \n'\
                   '\\bottomrule \n \\endlastfoot \n'%(d[k][2].replace('_',' '),protein,d[k][0],d[k][1]))
        bl=f[2:].split('\n')
        c=-1
        for e in bl:
            if e == '':
                continue
            else:
                outf.write(e.replace('\t','&'))
                c+=1
                outf.write('& %s \\\\ \n'%(memb[c]))
        
        outf.write('\\end{longtable} \n \\end{center} \n\n')
