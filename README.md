# crnpred
Prediction of 1D protein structures.

* 0. Introduction
This directory contains a computer program for predicting one-dimensional 
protein structures (secondary structures [SS], contact numbers [CN], and 
residue-wise contact orders [RWCO]) by the method of critical random networks 
described in:

Ref. 1 (Description of the software)
  "CRNPRED: Highly accurate prediction of one-dimensional protein structures
    by large-scale critical random networks."
  Kinjo AR, Nishikawa K.
  7:401 (2006) (DOI: https://doi.org/10.1186/1471-2105-7-401)

and 

Ref. 2 (Method of critical random networks)
  "Predicting secondary structures, contact numbers, and residue-wise contact 
    orders of native protein structure from amino acid sequence using critical 
    random networks."
  Kinjo AR, Nishikawa K.
  BIOPHYSICS, 1:67-74 (2005) (DOI: https://doi.org/10.2142/biophysics.1.67).

This software is in public domain. You can use, modify and/or destroy it freely,
but we do not take any responsibility for the consequences of your use.

* 1. INSTALLING CRNPRED.
 To install the CRNPRED program, you need the following:
 (0) UNIX-like operating system (Linux, MacOS X, *BSD, etc.)
 (1) bash (or zsh)
 (2) make
 (3) gcc 
 (4) PSI-BLAST and related databases (amino acid sequences and BLOSUM 
     scoring matrices).

  First, set the environment variable CRNPRED_DIR to this directory (that is,
  the directory containing this file "00README").
  If you are using sh, ksh, bash, or zsh, write

    export CRNPRED_DIR=/path/to/this/directory

  in your ~/.profile and do 
	% . ~/.profile

  If you are using csh or tcsh, write

    setenv CRNPRED_DIR /path/to/this/directory

  in your ~/.cshrc and do
	% source ~/.cshrc

  To compile the program, do 

	% ./make.sh

  Then the program named "xpredm" is installed under the directory 
  ${CRNPRED_DIR}/bin.

   After xpredm has been installed, test it by running

        % ${CRNPRED_DIR}/bin/xpredm sample/d3nul__.prof > hoge.out
  
   Compare hoge.out with sample/d3nul__.out. There are a few sample 
   inputs and outputs in the directory named "sample".

* 2. RUNNING CRNPRED.  
  Make sure you have set the environment variable CRNPRED_DIR appropriately.
  A utility shell script "run_crn.sh" is available for your convenience.
  If you have FASTA format amino acid sequence file (say, "test.seq"), do

  ${CRNPRED_DIR}/bin/run_crn.sh -d uniref100 test.seq 

  where "uniref100" is the sequence database used by PSI-BLAST. 
  Then, after some time,  you have a file named "test.seq.d.out" which contains
  the result of the prediction. If it does not work, check the content of 
  "run_crn.sh" and modify the environment variables such as BLASTDB, BLASTMAT, 
  and CRNPRED_DIR, or you may have to change the first line "#!/bin/sh" to 
  something like "#!/usr/bin/env bash" or "#!/usr/bin/env/ zsh".
  Run 
    ${CRNPRED_DIR}/bin/run_crn.sh -h
  to see other options.

  Alternatively, you can directly run the program. You first need to run 
  PSI-BLAST to make a position-specific scoring matrix:

  blastpgp -d nr -h 0.0005 -j 3 -i test.seq -Q test.prof > /dev/null 

  Then do 

  ${CRNPRED_DIR}/bin/xpredm test.prof > test.out

  The result is saved in "test.out".

* 3. INTERPRETING THE RESULTS.
  Below is an example of prediction.

    * Lines starting with "AA" show the amino acid sequence you fed.
    * Lines starting with "SS" show the predicted secondary structures 
      where "H", "E", and "C" mean "alpha-helix", "beta-strand", and "coils", 
      respectively.
    * Lines starting with "CN" show the predicted contact numbers in 2-state 
      description where "B" and "E" mean "buried" and "exposed", respectively. 
      The threshold values are the average contact number for each residue 
      type (see Appendix below for the list of the average contact numbers).
    * Lines after "># AA : SS P_H P_E P_C : CN : RWCO" are the details of the 
      prediction:
          o The column corresponding to "AA" indicates the residue numbers 
            and the amino acid residues.
          o The column corresponding to "SS" indicates the predicted secondary 
            structure followed by the ad hoc probability for each secondary 
            structure class (i.e., "P_H" for the probability for the residue to
            be in the alpha-helix class, etc.).
          o The column corresponding to "CN" indicates the predicted contact 
            numbers in 2-state description ("B" or "E") followed by the real 
            predicted contact numbers.
          o The column corresponding to "RWCO" indicates the predicted 
            residue-wise contact orders (real numbers). 

---------sample output starts here--------
>prediction for: test.prof


#                  *         *         *         *         *         *
AA:       SWQSYVDDHLMCDVEGNHLTAAAILGQDGSVWAQSAKFPQLKPQEIDGIKKDFEEPGFLA
SS:       CCHHHHHHHHHCCCCCCCCHEEEEECCCCCEEEECCCCCCCCHHHHHHHHHCCCCCCCCC
CN:       BBBBBBEBBEBBBBBBBBEEEEEEEEEEEEEEEEEBBBBBEEBBEBBBEEBEBBBBBBBB
#                  *         *         *         *         *         *
AA:       PTGLFLGGEKYMVIQGEQGAVIRGKKGPGGVTIKKTNQALVFGFYDEPMTGGQCNLVVER
SS:       CCEEEECCCEEEEEECCCCEEEEECCCCCEEEEEECCCEEEEEEECCCCCCHHHHHHHHH
CN:       EBEEBEBBEEEEEEEBBBBBBEEEEEBBBEEEEEEEEEEEEEEEBBBBBBBBEEEBEEEB
#                  *
AA:       LGDYLIESEL
SS:       HHHHHHHCCC
CN:       EEEBEBBBBB
//

>#   AA : SS P_H P_E P_C : CN     : RWCO
    1 S : C   11   7  82 : B   14 :  840
    2 W : C   23  10  67 : B   22 : 1221
    3 Q : H   59  11  30 : B   18 :  864
    4 S : H   79   8  12 : B   18 :  860
    5 Y : H   86   6   7 : B   25 : 1199
    6 V : H   89   5   6 : B   27 : 1276
    7 D : H   90   5   6 : E   21 :  855
    8 D : H   90   4   6 : B   17 :  728
    9 H : H   89   5   6 : B   22 :  954
   10 L : H   85   6   8 : E   30 : 1188
   11 M : H   72   9  18 : B   24 :  850
   12 C : C   44  11  46 : B   22 :  826
   13 D : C   18   8  73 : B   18 :  669
   14 V : C   10   7  83 : B   22 :  751
   15 E : C    8   6  86 : B   17 :  593
   16 G : C    8   7  85 : B   18 :  640
   17 N : C   10   8  82 : B   17 :  696
   18 H : C   16  11  73 : B   24 :  808
   19 L : C   30  16  54 : E   32 : 1103
   20 T : H   45  22  32 : E   26 :  962
   21 A : E   37  43  20 : E   27 : 1017
   22 A : E   15  75  10 : E   37 : 1265
   23 A : E    7  88   5 : E   38 : 1286
   24 I : E    6  89   5 : E   41 : 1341
   .
   .
   .
   .
---------sample output ends here--------

* 4. CONTACT INFORMATION

Akira Kinjo
Center for Information Biology and DNA Data Bank of Japan,
National Institute of Genetics,
Mishima, 411-8540, JAPAN
email: akinjo @ genes . nig . ac . jp

* Appendix A.
 The average contact number of each residue type is listed below:
-------------------
  25.430 , /* A */
  21.038 , /* R */
  20.093 , /* N */
  18.594 , /* D */
  29.647 , /* C */
  20.206 , /* Q */
  18.008 , /* E */
  22.505 , /* G */
  23.572 , /* H */
  29.469 , /* I */
  28.173 , /* L */
  18.452 , /* K */
  26.466 , /* M */
  28.057 , /* F */
  20.350 , /* P */
  21.420 , /* S */
  22.747 , /* T */
  26.913 , /* W */
  26.627 , /* Y */
  28.656 , /* V */
-------------------

* Appendix B.
** Faster but less accurate predictions.
The default implementation of CRNPRED uses 5000 dimensional state vectors for
critical random networks. This makes the prediction process quite slow when 
you use the program on old computers or when you predict large proteins.

If you want predictions quickly, there are two options: 
(1) linear predictor or (2) 2000 dimensional state vectors.

*** Using linear predictor
The linear predictor as described in Ref. 2 is implemented as a separate 
program named "lpredm" which is installed along with xpredm (CRNPRED).
Use it as follows:
  ${CRNPRED_DIR}/bin/lpredm test.prof > test.out

*** Using 2000 dimensional CRNPRED
To use CRNPRED with 2000 dimensional state vectors, you need to recompile the 
program. Do it as follows:
 
  cd ${CRNPRED_DIR}/src
  make realclean
  make NDIM=2000 install
  cd ..
  cp w2000/WMATS .
  cp w2000/WMAT_ENS .

This produces the executable file "xpredm" just like before, but it now 
uses 2000-dimensional state vectors.

*** Comparison of predictors
Here is a brief summary of speed and accuracy of the linear predictor (lpredm),
xpredm(2000),  and xpredm (5000). The CPU times were measured for 
the sample file "sample/d8abp__.prof" (305 AA) on Mac OS X (PPC G5, 2.5GHz).
The CPU time is (almost) linearly proportional to the protein length.


program	       speed		accuracy	note 
--------------------------------------------------------------
xpredm		very slow	SS:Q3=81	default
(5000)       	5min52s		CN:Cor=0.75
				RWCO:Cor=0.61	

xpredm		slow		SS:Q3=79	
(2000)       	1min12s		CN:Cor=0.74
				RWCO:Cor=0.61	

lpredm		fast		SS:Q3=76	
	       	0.558s		CN:Cor=0.72
				RWCO:Cor=0.59	
--------------------------------------------------------------

Note that the accuracies are the average values based on a benchmark.
The difference between Q3=81 and Q3=79 may seem insignificant on average, but 
there can be a big difference for individual predictions [e.g., an incorrectly
predicted alpha helix with xpred(2000) may be correctly predicted as a beta
strand with xpred(5000)].

Updated: 2009-06-10

# Local variables:
# mode: outline
# End:
