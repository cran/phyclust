/* This file contains declarations for em control, em phyclust struct, and
 * em function points. */


#ifndef __PHYCLSUT_EM_
#define __PHYCLSUT_EM_

#include "phyclust_struct.h"
#include "phyclust_qmatrix_array.h"

typedef struct _em_control		em_control;
typedef struct _em_phyclust_struct	em_phyclust_struct;
typedef struct _em_fp			em_fp;


/* Control em iterations and convergence in "phyclust_em_control.c". */
struct _em_control{
/* Fixed variables, used in initial and duplicate functions. */
	/* For em step. */
	int	exhaust_iter;			/* Exhaust em iterations. */
	int	fixed_iter;			/* Fixed em iterations for em+.EM. */
	int	short_iter;			/* Short em iterations. */
	int	EM_iter;			/* Long EM iterations. */
	double	short_eps;			/* Short em tolerance. */
	double	EM_eps;				/* Long EM tolerance. */

	/* For cm optim. */
	double	cm_reltol;			/* Relative tolerance. */
	int	cm_maxit;			/* Maximum iterations. */

	/* For NM optim. */
	double	nm_abstol_Mu_given_QA;		/* Absolute tolerance, update_flag = 0. */
	double	nm_abstol_QA_given_Mu;		/* Absolute tolerance, update_flag = 1. */
	double	nm_reltol_Mu_given_QA;		/* Relative tolerance. */
	double	nm_reltol_QA_given_Mu;		/* Relative tolerance. */
	int	nm_maxit_Mu_given_QA;		/* Maximum iterations. */
	int	nm_maxit_QA_given_Mu;		/* Maximum iterations. */
	/* For simplified optim. */
	int	est_non_seg_site;		/* Estimate non-segregating sites. */

	/* For initialization and EM steps. */
	int	max_init_iter;			/* Max initial steps. */
	int	min_n_class;			/* Min n for each class. */
	int	init_procedure;			/* Initialization procedure. */
	int	init_method;			/* Initialization method. */
	int	substitution_model;		/* Substition model. */
	int	edist_model;			/* Distance model for initialization. */
	int	identifier;			/* EE, EV, VE, VV. */
	int	code_type;			/* NUCLEOTIDE/SNP. */
	int	em_method;			/* EM method. */
	int	boundary_method;		/* Method to deal with boundary solutions. */
	double	Eta_lower_bound;		/* Lower bound of Eta, 1/N_X_org. */
	double	Eta_upper_bound;		/* Upper bound of Eta, 1 - Eta_lower_bound. */

/* Dynamical variables, used in copy functions only. */
	/* For report, reset by reset_em_control(). */
	double	converge_eps;			/* Convergent tolerance. */
	double	converge_error;			/* Convergent error as likelihood decreasing. */
	int	converge_flag;			/* Convergent flag. */
	int	converge_iter;			/* Convergent iteration. */
	int	converge_inner_iter;		/* Convergent inner iteration. */
	int	converge_cm_iter;		/* Convergent CM iteration. */
	int	update_flag;			/* 0 for update Mu, 1 for update Q. */
};

em_control* initialize_em_control();
void free_em_control(em_control *EMC);
em_control* duplicate_em_control(em_control *org_EMC);
void update_em_control(em_control *EMC);	/* update EMC when initial settings are changed. */
void reset_em_control(em_control *EMC);		/* For long EM steps. */


/* EM phyloclustering structure for EM steps only. */
struct _em_phyclust_struct{
/* Fixed variables, used in initial and duplicate functions, copy or point to pcs. */
	/* Define code type. */
	int		code_type;		/* NUCLEOTIDE/SNP. */
	int		ncode;			/* = NN/NNG(4/5) or NSNP/NSNPG(2/3). */
	/* Constant points to phylcust_struct *pcs. */
	int		N_X_org;		/* Number of original sequences. */
	int		N_X_unique;		/* Number of unique sequnces, used for initialization only*/
	int		N_X;			/* Number of sequences. */
	int		N_seg_site;		/* Total segregating sites. */
	int		L;			/* Number of loci. */
	int		K;			/* Number of clusters. */
	int		**X_org;		/* Original data pointer, dim = N_X_org * L. */
	int		**X_unique;		/* Unique data pointer, dim = N_X_unique * L. */
	int		**X;			/* Data pointer, dim = N_X * L if compress = 1, dim = N_X_org * L otherwise. */
	int		*map_X_org_to_X_unique;	/* Map indexes from X_org to X_unique, dim = N_X_org, used for initialization only. */
	int		*map_X_unique_to_X_org;	/* Map indexes from X_unique to X_org, dim = N_X_unique, used for initialization only. */
	int		*replication_X_unique;	/* Count replications of each unique sequence, dim = N_X_unique. */
	int		*map_X_org_to_X;	/* Map indexes from X_org to X, dim = N_X_org. */
	int		*map_X_to_X_org;	/* Map indexes from X to X_org, dim = N_X. */
	int		*replication_X;		/* Count replications of each unique sequence, dim = N_X. */
	/* Summarized information for X. */
	int		*seg_site_id;		/* Segregating site id, dim = N_seg_site. */
	int		*class_id;		/* For manually initialization, dim = N_X_org. */

/* Dynamical variables, used in copy functions only, owned memory and need to be freeed. */
	/* For EM. */
	int		*n_class;		/* Number of sequences for each class, dim = K. */
	int		**Mu;			/* Centers, dim = K * L. */
	double		**Z_modified;		/* unnormalized log Z, dim = N_X * K. */
	double		**Z_normalized;		/* Normalized Z, dim = N_X * K. */
	double		*Z_total;		/* Total posterior given k, dim = K. */
	double		*Eta;			/* Proportion, dim = K. */
	double		*log_Eta;		/* Log of proportion, dim = K. */
	double		logL_observed;		/* Observed logL. */

	/* Computing storage. */
	int		****count_Mu_X;		/* Used in update_Z_modified(), logL_observed(), update_count_Mu_X(). */
};

em_phyclust_struct* initialize_em_phyclust_struct(phyclust_struct *pcs);
void free_em_phyclust_struct(em_phyclust_struct *empcs);
em_phyclust_struct* duplicate_em_phyclust_struct(em_phyclust_struct *org_empcs);

/* E_step's. */
void e_step_with_stable_exp(int *K, double *a_Z_normalized, double *total_sum, double *scale_exp, int *flag_out_range);
void E_step_and_logL_observed_AU(em_phyclust_struct *empcs, Q_matrix_array *QA);
void E_step_and_logL_observed_NU(em_phyclust_struct *empcs, Q_matrix_array *QA);
void E_step_simple(em_phyclust_struct *empcs, Q_matrix_array *QA);

/* m_step's. */
int M_step_simple(em_phyclust_struct *empcs, Q_matrix_array *new_QA, em_control *EMC, em_fp *EMFP, em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA);
/* All CM steps require to set EMC->update_flag = 1. */
int M_step_CM(em_phyclust_struct *empcs, Q_matrix_array *new_QA, em_control *EMC, em_fp *EMFP, em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA);
int M_step_ACM(em_phyclust_struct *empcs, Q_matrix_array *new_QA, em_control *EMC, em_fp *EMFP, em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA);

/* Method to update Eta.
 * AU = all unique (treating all sequences are distinct).
 * NU = not unique (deal with summarized unique sequences).
 * IGNORE = ignore the results and stop if Eta's are out of boundary.
 * ADJUST = adjust the Eta's to be fixed on the boundary. */
int Update_Eta_given_Z_ADJUST_AU(em_phyclust_struct *empcs, em_control *EMC);
int Update_Eta_given_Z_ADJUST_NU(em_phyclust_struct *empcs, em_control *EMC);
int Update_Eta_given_Z_IGNORE_AU(em_phyclust_struct *empcs, em_control *EMC);
int Update_Eta_given_Z_IGNORE_NU(em_phyclust_struct *empcs, em_control *EMC);

/* This special function updates empcs->Z_modified (log and unnormalized). */
void update_Z_modified(em_phyclust_struct *empcs, Q_matrix_array *QA);

/* Check convergence. */
int Check_convergence_em(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, em_control *EMC);
int Check_convergence_org(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, em_control *EMC);

/* Different EM methods. */
void Em_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);
void Short_em_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);


/* Control function pointers in "phyclust_em_fp.c". */
struct _em_fp{
/* Deal with different initialization methods. */
	/* In "phyclust_init_method.c". */
	void	(*Update_init)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*);	/* Initialization method. */

/* Deal with different EM algorithms. */
	/* In "phyclust_em.c". */
	int	(*M_step)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*, em_phyclust_struct*, Q_matrix_array*);		/* M_step. */
	int	(*Check_convergence)(em_phyclust_struct*, em_phyclust_struct*, em_control*);	/* Check convergence. */
	void	(*Em_step)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*);		/* Em_step. */
	void	(*Short_em_step)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*);	/* Short_em_step. */

/* Deal with unique or non-unique sequences. */
	/* In "phyclust_em.c". */
	void	(*E_step_and_logL_observed)(em_phyclust_struct*, Q_matrix_array*);			/* LogL by E_step. */
	int	(*Update_Eta_given_Z)(em_phyclust_struct*, em_control*);		/* Update Eta given Z. */
	/* In "phyclust_em_tool.h". */
	double	(*LogL_observed)(em_phyclust_struct*, Q_matrix_array*);		/* Observed logL. */
	double	(*LogL_complete)(em_phyclust_struct*, Q_matrix_array*);		/* Complete logL. */
	double	(*LogL_profile)(em_phyclust_struct*, Q_matrix_array*);		/* Complete logL. */
	void	(*Copy_empcs_to_pcs)(em_phyclust_struct*, phyclust_struct*);	/* copy empcs to pcs. */
	void	(*Copy_pcs_to_empcs)(phyclust_struct*, em_phyclust_struct*);	/* copy empcs to pcs. */
	/* In "phyclust_logpL.c". */
	void	(*Update_Mu_given_QA)(em_phyclust_struct*, Q_matrix_array*);	/* Update Mu given QA in logpL. */
	double	(*Compute_R)(em_phyclust_struct*);				/* Update function R(). */
};

em_fp* initialize_em_fp(em_control *EMC, phyclust_struct *pcs);
void free_em_fp(em_fp *EMFP);

#endif	/* End of __PHYCLSUT_EM_. */
