StructBio
=========

Scripts for structural biology. This includes modularity testing, svm classification for structures, feature selection, among others.

classifyGM.py:
-------------
  Given gmfiles (each one containing elements belonging to the same group and therefore with 
  an assigned label) create a Support Vector Machine model, classifying them. the gmfiles given
  are treated as the training set. Use -h to see options.
  A gmfile is a semicolon-delimited text file, where the first element is the name of the structure/shape.

RFclassGM.R:
------------
  Given a gmfile (each one containing elements belonging to the same group and therefore with 
  an assigned label) and a membership vector file (space delimited labels in a single line, of the same lenght of entries in  
  the GM file). Run it once for options:
  -LDfs XX: Where XX can be RF (random forest) or KW (kruskall-wallis). Feature selection with RF or KW followed by an linear 
    discriminant pre-classification to curate the training set based on a 95% confidence ellipses.
  -MDS : A non-metric MDS transformation of the data.
  -FS XX: XX can be also RF or KW and will perform feature selection without a posterior LD curation
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
  
  Usage:  Rscript featureSel.R gmfile membershipfile dimensions FDM
  
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
  

StructureOutlier.py:
-------------------
  This script will identify full structure outliers from a gm file. It is based on the Grubss test with RMSD as 
  variable.
  
  Requires:
  
  1) rpy2, available at http://rpy.sourceforge.net/rpy2.html
  
  2) Numpy
  
  3) R packages outliers


landmark2latex.py:
-----------------
Really simple script to transform a landmark file into latex format (yes I know, stupid, 
but somre reviwers ask or it).


ldaellipse4mod.R:
----------------
  R script for Linear discriminat analysis and 95% ellipse collision test for GM data. 
  This script will classify the variables in a GM file given a cluster file with an LDA,
  compute the 95% confidence ellipses around the groups, and compute the collisions among
  all compured ellipses
  
  Usage: Rscript featureSel.R prefix dimensions
  
  Arguments:
  
  1) Prefix: The name scheme use in the analysis (normally modularity using Moduler).
             Two files must be available an [prefix].lms with the agglomerated correlation
             vector magnitudes of a shape, and [prefix].community, with a membership vector
             as a single line text file with the membership vector as string splitted by
             space.
             
  2) Number of dimensions of the structures (must be multiple of the total of variables)


MapModules.py:
-------------
  THIS IS MAINLY DEPRECATED AND NOT MAINTAINED SINCE IT WAS INCORPORATED INTO MODULER!!!
  Map clustering partition into the chain field of a PDB structure, the centrality or form difference into beta and 
  ocupancy. This also works for a multiple pdb(>50) with the -m option. it requires [prefix].landmark, [prefix].gm, and   [prefix].graphcluster files. If used without the -m option it also needs a [prefix].pdb file.
  
  See help for more details.
  
  
  Requires:
  
  1) Utils folder, available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools
  
  
Moduler.py:
----------
  Graph based Modularity. This script is a joint effort between me and Kyle Nguyen, Alex Safatli and Christian Blouin
  
  
  This script will evaluate the data for modules. Such modules are defined as correlating variables, so the clustering 
  is performed in the correlation space. It has an optional statistical significance test for the clustering and power
  analysis of the result, as well as a bootstrap analysis. See options for more details.
  
  If you use this software, please cite:
  Hleap, J.S., Susko, E., & Blouin, C.Defining structural and evolutionary modules in proteins: a community
  detection approach to explore sub-domain architecture. BMC Structural Biology, 13, 20.
  
  Requirements:
  
  1) Python:
  
     a) numpy module
     
     b) scipy module
     
     c) rpy2 module
     
     d) matplotlib module
     
     ######################################
     #To install python modules in UBUNTU:#
     #sudo apt-get python-<module>        #
     ######################################
  
  2) The igraph Python extension module:
  
     ######################################
     #To install python-igraph in UBUNTU: #
     #go to http://igraph.sourceforge.net/#
     #follow their instructions           #
     ######################################
  
  3) The GNU R package:
     Libraries DAAG, pwr, corpcor, FactoMineR, gplots
     
     ######################################
     #To install The R packages in UBUNTU:#
     #Get into R by typing "R" in your cmd#
     #Once in R, type:                    #
     #install.packages('<package>')       #
     ######################################   
  
  4) Folder contactlist, available in https://github.com/jshleap/Collaboration
  
  5) Matplotlib, avalable at Matplotlib, available at http://matplotlib.org/
  
  6) rpy2, available at http://rpy.sourceforge.net/rpy2.html
  
  7) Scipy

comparecommunities.py:
----------
  Given two membership vectors of communities (i.e. Moduler.py output) compute the comparison of the two using      
  estandardized mutual information and a null distribution permuting the communities:

  Requirements:
  
  1) Python:
      a) igraph
      b) Scipy
      c) Matplotlib
  
  THIS STILL IN BETA!!!
  
comparecommunities.py:
----------
  Searching for optimal 1D partition for 3D correlation modularity. Essentially looking for AFUs-like modules within 
  proteins, although can be used for general purposes when a correlation graph wants to be further constrained by 
  sequentiality.
  Requirements:
  
  1) Python:
      a) igraph
      b) numpy
      c) pickle


