"""
Given a pdb.xml output from modweb, split the pdb into individuals and align them.
This assume matt aligner on your path
"""

import xml.etree.ElementTree as ET
from subprocess import Popen
import sys, os

tree = ET.parse(sys.argv[1])
files = tree.getroot()

if not os.path.isdir(sys.argv[2]):
    os.mkdir(sys.argv[2])
for pdbfile in files:
    content = pdbfile[1].text
    name = content[content.find('Original ID:'):]
    name = name[: name.find('\n')].strip().split()[2]
    with open(os.path.join(sys.argv[2],'%s.pdb' % name), 'w') as f:
        f.write(content)

try:
    cpus = int(sys.argv[3])
except:
    cpus = 1

outfn = sys.argv[1][:sys.argv[1].find('.pdb')]
matt = Popen('Matt -t %d -o %s_aln %s/*.pdb' % (outfn, cpus, sys.argv[2]),
             shell=True)


