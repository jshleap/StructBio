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

outfn = sys.argv[1][:sys.argv[1].find('.pdb')]
matt = Popen('Matt -o %s_aln %s/*.pdb' % (outfn, sys.argv[2]),
             shell=True)


