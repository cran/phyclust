/* This file contains declarations for phyclust. */

#ifndef __PHYCLSUT_STRUCT_
#define __PHYCLSUT_STRUCT_

#include "phyclust_qmatrix_array.h"

typedef struct _phyclust_struct	phyclust_struct;

/* Phyloclustering structure for storing data and results. */
struct _phyclust_struct{
/* Fixed variables, used in initial and duplicate functions. */
	/* Define code type. */
	int		code_type;		/* NUCLEOTIDE/SNP. */
	int		ncode;			/* = NN/NNG(4/5) or NSNP/NSNPG(2/3). */
	/* Storage. */
	int		n_param;		/* Number of parameters, including Mu's, and eta's. */
	int		compress_method;	/* 1 for compress, 0 otherwise.*/
	int		N_X_org;		/* Number of original sequences. */
	int		N_X_unique;		/* Number of unique sequnces, used for initialization only*/
	int		N_X;			/* Number of sequences. */
	int		N_seg_site;		/* Total segregating sites. */
	int		L;			/* Number of loci. */
	int		K;			/* Number of clusters. */
	int		**X_org;		/* Original data pointer, dim = N_X_org * L. */
	int		**X_unique;		/* Unique data pointer, dim = N_X_unique * L. */
	int		**X;			/* Data pointer, dim = N_X * L if compress_method = 1, dim = N_X_org * L otherwise. */
	int		*map_X_org_to_X_unique;	/* Map indexes from X_org to X_unique, dim = N_X_org, used for initialization only. */
	int		*map_X_unique_to_X_org;	/* Map indexes from X_unique to X_org, dim = N_X_unique, used for initialization only. */
	int		*replication_X_unique;	/* Count replications of each unique sequence, dim = N_X_unique. */
	int		*map_X_org_to_X;	/* Map indexes from X_org to X, dim = N_X_org. */
	int		*map_X_to_X_org;	/* Map indexes from X to X_org, dim = N_X. */
	int		*replication_X;		/* Count replications of each unique sequence, dim = N_X. */
	/* Summarizd information for X. */
	int		*seg_site_id;		/* Segregating site id, dim = N_seg_site. */

/* Dynamical variables, used in copy functions only. */
	/* Parameters. */
	int		**Mu;			/* Centers, dim = K * L. */
	double		*Eta;			/* Proportion, dim = K. */
	double		**Z_normalized;		/* Normalized Z, dim = N_X * K. */
	/* For summary. */
	double		logL_observed;		/* Observed logL. */
	double		logL_entropy;		/* Observed logL + entropy. */
	double		bic;			/* bic. */
	double		aic;			/* aic. */
	double		icl;			/* icl. */
	int		*class_id;		/* Class id of each sequence, dim = N_X. */
	int		*n_class;		/* Number of sequences for each class, dim = K. */
};

phyclust_struct* initialize_phyclust_struct(int code_type, int N_X_org, int L, int K);
void free_phyclust_struct(phyclust_struct *pcs);
void update_phyclust_struct(phyclust_struct *pcs, int compress_method);


/* ----- Summary tool. ----- */
void assign_class(phyclust_struct *pcs);
void update_ic(phyclust_struct *pcs, Q_matrix_array *QA);


/* ----- For debug. ----- */
void print_X(phyclust_struct *pcs);
void print_Mu(phyclust_struct *pcs);
void print_class_id(phyclust_struct *pcs);

#endif	/* End of __PHYCLSUT_STRUCT_. */
