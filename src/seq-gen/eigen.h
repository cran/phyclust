/*  Header file for eigen.c                                */

/*  
   Sequence Generator - seq-gen, version 1.3.2
   Andrew Rambaut & Nick Grassly
   Department of Zoology, University of Oxford			
	
   The code in this file is taken from Ziheng Yang's PAML package.
   http://abacus.gene.ucl.ac.uk/

   Any feedback is very welcome.
   http://evolve.zoo.ox.ac.uk/software/Seq-Gen/
   email: andrew.rambaut@zoo.ox.ac.uk
*/


//WCC:add
#include "R_seq-gen.h"


#ifndef _EIGEN_H_
#define _EIGEN_H_

//WCC int abyx (double a, double x[], int n);
int seq_gen_abyx (double a, double x[], int n);

//WCC int xtoy (double x[], double y[], int n);
int seq_gen_xtoy (double x[], double y[], int n);

//WCC int matinv( double x[], int n, int m, double space[]);
int seq_gen_matinv( double x[], int n, int m, double space[]);

int eigen(int job, double A[], int n, double rr[], double ri[],
          double vr[], double vi[], double w[]);

#endif /* _EIGEN_H_ */
