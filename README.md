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
            
FoldX_runner.py:
---------------
  This script will run FoldX software for repair and mutations to alanine.
  It will output a pickled dictionary of energies and the mutated files.
  Check FoldX manual for more info.
  Requires:
  1) FoldX, available at http://foldx.crg.es/
  2) Numpy
  

protGM.py:
---------
  Application of some Geometric morphometrics on protein structures. That applications include 
  abstraction of structures as shapes, PCoA, CVA, form differnce, MANOVA, outlier landmark removal
  among others.
  
  It uses some R functions taken and modified from Morphometrics with R(Series: Use R!) of Julien Claude.
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  Requires:
  1) rpy2, available at http://rpy.sourceforge.net/rpy2.html
  2) Utils folder, available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools
  3) Biophyton
  4) Numpy
  5) R packages ape, vegan, shapes, outliers, and corpcor


ShapeSim.py:
-----------
  Simulates shapes using cholesky decomposition to generate correlated variables (coordinates) with a shape constraint.
  Requires an initial shape in GM format, and a partition file of the form:
  
    Partition
    Landmarks: [list of the landmark number in the order of entries in the GM file]
    corr: [the ammount of correlation among partition 1]
    
    Partition
    Landmarks: [list of the landmark number in the order of entries in the GM file]
    corr: [the ammount of correlation among partition 2]
    .
    .
    .
    Partition
    Landmarks: [list of the landmark number in the order of entries in the GM file]
    corr: [the ammount of correlation among partition n]
  
  Requires:
  1) rpy2, available at http://rpy.sourceforge.net/rpy2.html
  2) Matplotlib, available at http://matplotlib.org/
  3) R packages corpcor and shapes

