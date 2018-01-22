HGT Detection 3.4 - March, 2012
--

###Citation

If you use this program in a publication, please cite the following reference:	
Boc, A., Philippe, H. and Makarenkov, V. (2010), Inferring and validating horizontal gene transfer events using bipartition dissimilarity, Systematic Biology.


###Requirements
With Linux/UNIX :
Perl v5.10.0 or higher
g++ compiler (or any other C++ compiler)


###Summary
1. How to create the executable file.
2. How to run the HGT-DETECTION program.
3. Options.
4. Input file example.
5. How to specify species and gene tree root in a file (speciesRoot and geneRoot files).
6. Output file example.
7. Important notes.


####1) How to create the executable file

With Linux/UNIX, you need to compile the following files in order to create the executable file (hgt*):

 |-> hgt3.4.cpp
 |-> structures.h
 |-> utils_tree.cpp
 |-> fonctions.cpp
 |-> Makefile
 
To do it, carry out the command:
 > make hgt3.4


####2) How to call the HGT-DETECTION program

    run :> perl run_hgt3.4.pl -inputfile=inputfilename [options]
the program will provide a formated results files (output.txt - see an example below)

You can also run directly the c++ program:   

    run :> hgt -inputfile=inputfilename [options]
the program will provide a raw results file only (results.txt - see an example below) and HGT bootstrapping is not be available with this command


####3) Options

* -inputfile : input file must contain the species and the gene phylogenetic trees in the Newick format (the species tree must be reported first) the trees can be rooted; if they are unrooted, they will be rooted by midpoint (default option)

* -outputfile : default name of the output file is output.txt

* -speciesroot or -generoot :  you can also specify the root of the species (or gene) tree in a file, otherwise the tree will be rooted by midpoint
 |-> file                 :  the root position is specified in the file called speciesRootLeaves.txt (or geneRootLeaves.txt) as indicated in Point 5 below
 |-> bestbipartition      :  if the root position in the gene tree is uncertain, this option locates the gene tree root in such a way that the bipartition induced by the gene tree root is the closest possible to the bipartition induced by the root of the species tree (this option is available for the gene tree only)
 |-> midpoint             :  the tree is rooted by midpoint (default option)
 
* -mode :  
 |-> multicheck : several HGTs by iteration are detected (default option)
 |-> monocheck : one HGT by iteration is detected (the execution with this option is a little bit slower)
 
* -criterion : 
 |-> bd : Bipartiton dissimilarity (default)
 |-> rf : Robinson and Foulds distance
 |-> ls : Least-squares

* -bootstrap :
 |-> no  : default
 |-> yes : if you want to carry out HGT bootstrap, the first line of your input file should specify the species tree, the second - the gene tree, and all the following lines - the gene tree replicates. 
  HGT boostrap is computed taking into account the robustness of the gene tree and all possible minimum-cost HGT scenarios (for more details, see Boc et al., 2010). The root of the gene tree replicates may be not specified.

* -nbhgt : maximum number of detected HGTs in the computation (default value : 100)

* -path : location of your data files (input file, output file, ....) 
  
  
####4) Input file example

    ((J:1.0,A:1.0):1.0,((C:1.0,E:1.0):10.0,(G:1.0,(F:1.0,H:1.0):1.0):1.0):1.0); |-> species tree
    ((J:1.0,F:1.0):1.0,((C:1.0,G:1.0):1.0,(E:1.0,(A:1.0,H:1.0):1.0):1.0):1.0);  |-> gene tree

 
####5) How to specify species and gene tree root (speciesRoot and geneRoot files)

If the species and (or) gene tree(s) are not rooted, the root position(s) can be specified in the species root (speciesRootLeaves.txt) and (or) gene root (or geneRootLeaves.txt) file(s) as follows: 
 
    C E F H G
    <>
    J A
Example of a command line call with this option: 

    perl run_hgt3.4.pl -inputfile=Input_7_taxa.txt -outputfile=output.txt -speciesroot=file

 
####6) Output file example

An example of a FORMATTED OUTPUT file:

    =================================================================================    
    Program :   HGT Detection 3.2 - November, 2009							
    Authors   : Alix Boc and Vladimir Makarenkov (Universite du Quebec a Montreal)	
    This program computes a unique scenario of horizontal gene transfers (HGTs) for the given pair of species and gene phylogenetic trees.				
    =================================================================================					    
    Subtree constraint   :yes								
    Criterion            :bipartition dissimilarity						
    Bootstrap            :no								
    
    Species tree :
    ((J:1.0,A:1.0):1.0,((C:1.0,E:1.0):10.0,(G:1.0,(F:1.0,H:1.0):1.0):1.0):1.0);
    
    Gene Tree :	
    ((J:1.0,F:1.0):1.0,((C:1.0,G:1.0):1.0,(E:1.0,(A:1.0,H:1.0):1.0):1.0):1.0);											
    =============================================									    
    = Criteria values before the computation 									
    =============================================									    
    Robinson and Foulds distance (RF) = 10									    
    Least-squares coefficient(LS)     = 30869.737									      
    Bipartition dissimilarity         = 11.0									
    
    ================================================================						
    | Iteration #1 : 3 HGTs were found										   
    ================================================================						
    |															    
    | HGT 1 / 3 												    
    | From subtree (H) to subtree (A)											    
    | RF = 8 , LS = 30789.952 , BD = 6.0									    
    | 															    
    | HGT 2 / 3 												    
    | From subtree (C) to subtree (G)											    
    | RF = 8 , LS = 29533.007 , BD = 8.0									    
    | 															    
    | HGT 3 / 3 												    
    | From subtree (F) to subtree (J)											    
    | RF = 8 , LS = 29671.656 , BD = 7.5									    
    | 															    
    | ================================================================						
    | After this iteration the criteria values are as follows :							
    | RF = 4 , LS = 10367.302 , BD = 2.5										    
    | ================================================================						
    |															    
    | ================================================================						
    | Iteration #2 : 1 HGT was found											    
    | ================================================================						
    |															    
    | HGT 1 / 1 												    
    | From subtree (E) to subtree (A, H)										    
    | RF = 0 , LS = 6807.857 , BD = 0.0										    
    | 															    
    | ================================================================						
    | After this iteration the criteria values are as follows :							
    | RF = 0 , LS = 0.000 , BD = 0.0											    
    | ================================================================						
    | 															    
    | Total number of HGTs : 4 
    -------------------------------------------------------------------


An example of a RAW OUTPUT file (e.g., it can be used in simulations):

    ----------------------------------------------------
    | 10,30869.736516,11.000000                        |  -> RF, LS and BD values
    | mode=multicheck                                  |  -> selected mode
    | 3                                                |  -> number of HGT(s) found at the current iteration
    | HGT 1 / 3                                        |  -> first HGT found at the current iteration
    | H                                                |  -> list of species of the source subtree
    | A                                                |  -> list of species of the recepient subtree
    | From branch H--12 to branch A--9                 |  -> two species tree branches affected by this HGT
    | RF = 8 , LS = 30789.952463 , BD = 6.000000       |  -> criteria values (RF, LS and BD) after this HGT	
    | HGT 2 / 3                                        |  -> .....
    | C                                                |
    | G                                                |
    | From branch C--10 to branch G--13                |
    | RF = 8 , LS = 29533.006815 , BD = 8.000000       |
    | HGT 3 / 3                                        |
    | F                                                |
    | J                                                |
    | From branch F--12 to branch 9--J                 |
    | RF = 8 , LS = 29671.655834 , BD = 7.500000       |
    | RF = 4 , LS = 10367.301891 , BD = 2.500000       |  -> criteria values (RF, LS and BD) after the 3 HGTs above
    | 1                                                |
    | HGT 1 / 1                                        |
    | E                                                |
    | A H                                              |
    | From branch E--10 to branch H--12                |
    | RF = 0 , LS = 6807.856562 , BD = 0.000000        |
    | RF = 0 , LS = 0.000000 , BD = 0.000000           |
    | Stat= 4 0                                        | -> number of regular HGT(s) found (followed by the number of trivial HGT(s), if any)
    ----------------------------------------------------

 
####7) Important notes

* If the branch specified (as a root branch) in the speciesRoot or geneRoot file is not found, the branch closest to your selection will be used to root the tree(s)
* If a species is present in the species tree and absent in the gene tree (or vice versa), then this species will be ruled out from the computation
 
