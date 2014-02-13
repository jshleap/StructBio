StructBio
=========

Scripts for structural biology. This includes modularity testing, svm classification for structures, feature selection, among others.

classifyGM.py:
-------------
Given gmfiles (each one containing elements belonging to the same group and therefore with 
an assigned label) create a Support Vector Machine model, classifying them. the gmfiles given
are treated as the training set. Use -h to see options.
A gmfile is a semicolon-delimited text file, where the first element is the name of the structure/shape.
