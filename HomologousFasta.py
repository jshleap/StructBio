#!/usr/bin/python
''' 
plotting sequence identities and RMSDs from homologous landmarks
'''
from SeqMask import InferSeqs_landmarks, dict2fasta
import sys

prefix=sys.argv[1]


seqs=InferSeqs_landmarks(prefix)
dict2fasta(prefix+'.hom.fasta',seqs)
