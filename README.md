StructBio
=========

Scripts for structural biology. This includes modularity testing, svm classification for structures, feature selection, among others.

classifyGM.py:
-------------
  Given gmfiles (each one containing elements belonging to the same group and therefore with 
  an assigned label) create a Support Vector Machine model, classifying them. the gmfiles given
  are treated as the training set. Use -h to see options.
  A gmfile is a semicolon-delimited text file, where the first element is the name of the structure/shape.


Equi2Latex.py:
-------------
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
  
  Requires: Numpy
  
featuresel.R:
------------
  R script for feature selection for GM data. Will pick variables that have 
  a significant kuskal-wallis difference given a classification scheme.
  Usage:  Rscript featureSel.R <gmfile> <membershipfile> <dimensions> <FDM>
          Arguments:
          1) File name of the GM data
          2) File name of the membership (must include the membership for all 
             structures in the gm file)
          3) Number of dimensions of the structures (must be multiple of the 
             total of variables)
          4) a boolean (TRUE/FALSE) if the feature selection should be made on 
             the raw coordinates (FALSE) or in the form difference matrix
            
