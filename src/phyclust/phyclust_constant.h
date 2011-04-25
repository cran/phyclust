/* This file contains constants for phyclust. */

#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


#ifndef __PHYCLUST_CONSTANT_
#define __PHYCLUST_CONSTANT_

#include <float.h>
#define Inf DBL_MAX

#ifdef __GNUC__
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#else
#define VARIABLE_IS_NOT_USED
#endif


/* Nucleotides. */
#define NN 4								/* Number of nucleotides. */
#define NNG 5								/* Index of nucleotides with gap. */
enum {A, G, C, T, GAP};							/* Nucleotides. */
static const char NUCLEOTIDE_CODE[NNG]  = {'A', 'G', 'C', 'T', '-'};
static const char NUCLEOTIDE_lower[NNG]  = {'a', 'g', 'c', 't', '-'};
static const char NUCLEOTIDE_ID[NNG]    = {'0', '1', '2', '3', '4'};

/* SNPs. */
#define NSNP 2								/* Number of SNPs. */
#define NSNPG 3								/* Index of SNPs with unknown. */
enum {NL, SU, UN};							/* Normal/Substitution/Unknown. */
static const char SNP_CODE[NSNPG] = {'0', '1', '-'};
static const char SNP_ID[NSNPG] = {'1', '2', '-'};

/* CODEs. */
#define N_CODE_TYPE 2							/* Total number of code type. */
enum {NUCLEOTIDE, SNP};							/* Number of codes. */
static char VARIABLE_IS_NOT_USED *CODE_TYPE[N_CODE_TYPE] =
	{"NUCLEOTIDE", "SNP"};
static const int NCODE[N_CODE_TYPE] = {NN, NSNP};			/* Indicate the dimension. */
#define PRINT_CODE_TYPE 0						/* 0 for CODE, 1 for CODE_ID. */
static const int MISSING_INDEX[N_CODE_TYPE] = {NN, NSNP};		/* Indicate the missing index. */


/* EM procedures. */
#define N_INIT_PROCEDURE 4
enum {exhaustEM, emEM, RndEM, RndpEM};					/* Initialization procedures. */
static char VARIABLE_IS_NOT_USED *INIT_PROCEDURE[N_INIT_PROCEDURE] =
	{"exhaustEM", "emEM", "RndEM", "RndpEM"};

/* Initialization methods. */
#define N_INIT_METHOD 6
enum {randomMu, NJ, randomNJ, PAM, kMedoids, manualMu};			/* Initialization methods. */
static char VARIABLE_IS_NOT_USED *INIT_METHOD[N_INIT_METHOD] =
	{"randomMu", "NJ", "randomNJ", "PAM", "K-Medoids", "manualMu"};

/* Substitution models. */
#define N_SUB_MODEL 9							/* Number of evolution models. */
enum {JC69, K80, F81, HKY85, SNP_JC69, SNP_F81,
	E_F81, E_HKY85, E_SNP_F81};
static char VARIABLE_IS_NOT_USED *SUBSTITUTION_MODEL[N_SUB_MODEL] =
	{"JC69", "K80", "F81", "HKY85", "SNP_JC69", "SNP_F81",
		"E_F81", "E_HKY85", "E_SNP_F81"};

/* Distance models. */
#define N_EDIST 3
enum {D_JC69, D_K80, D_HAMMING};					/* Distance models. */
static char VARIABLE_IS_NOT_USED *EDISTANCE_MODEL[N_EDIST] =
	{"D_JC69", "D_K80", "D_HAMMING"};


/* Identifier for constraints. */
#define N_IDENTIFIER 4
enum {EE, EV, VE, VV};							/* First letter for Q, second for Tt. */
static char VARIABLE_IS_NOT_USED *IDENTIFIER[N_IDENTIFIER] =
	{"EE", "EV", "VE", "VV"};


/* EM methods. */
#define N_EM_METHOD 3
enum {EM, ECM, AECM};
static char VARIABLE_IS_NOT_USED *EM_METHOD[N_EM_METHOD] =
	{"EM", "ECM", "AECM"};


/* Boundary methods. */
#define N_BOUNDARY_METHOD 2
enum {ADJUST, IGNORE};
static char VARIABLE_IS_NOT_USED *BOUNDARY_METHOD[N_BOUNDARY_METHOD] =
	{"ADJUST", "IGNORE"};


/* Rprintf: print message to R console. */
#ifdef R_EXT_PRINT_H_
#undef printf
#define printf Rprintf
#undef exit
#define exit(a) error("%d\n", a)
#endif


/* EM Debugging only. */
#define EMDEBUG 0		/* 0 for no output, 1 for em steps, >1 for m step, >2 for nm, >3 for -logpL.
				 * 0 = 0000 for no output,
				 * 1 = 0001 for em steps,
				 * 2 = 0010 for m step,
				 * 4 = 0100 for nm step,
				 * 8 = 1000 for -logpL. */
#define INITDEBUG 0		/* 1 */		/* 0 for no output, >0 for tracing initialization. */
#define verbosity_em_step 0	/* 3 */		/* 0 for no output, >0 for print information. */
#define verbosity_exhaust_EM 0	/* 1 */		/* 0 for no output, >0 for print information. */
#define PRINT_ERROR 0		/* 1 */		/* 0 for no output, >0 for error messages. */


#endif	/* End of __PHYCLUST_CONSTANT_. */

