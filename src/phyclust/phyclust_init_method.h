/* This file contains declarations for phyclust. */


#ifndef __PHYCLSUT_INIT_METHOD_
#define __PHYCLSUT_INIT_METHOD_

#include "phyclust_em.h"
#include "phyclust_ape_nj.h"


/* Tools for initialization. */
int rdunif(int n);
void srswor(int n, int k, int *x);
int init_m_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);
int check_all_min_n_class(int k, int *n_class, int min_n_class);
void assign_Mu_by_class(int N_X_org, int K, int L, int ncode, int *class_id, int **X_org, int **Mu);


/* Initialization functions.
 * All functions will update empcs, QA.
 * Basically, these functions only provides Mu and update QA.
 * Further updates for em are implemented "phyclust_init_procedure.c". */

/* Randomly pick Mu's and assign by edist. */
void Update_init_random_Mu_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By neighbor-joining. */
void Update_init_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By random neighbor-joining. */
void Update_init_random_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By pam. */
void Update_init_pam(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By k-medoids. */
void Update_init_k_medoids(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By hand-coding or a prespecified input file. */
void Update_init_manually(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);


/* File: "phyclust_init_method_nj.c". */
void search_largest_branch(nj_struct *njs, int *largest_branch_id);
void random_branch(nj_struct *njs, int *random_branch_id);
int assign_class_by_njs_branch(int K, nj_struct *njs, int *branch_id, int *class_id);

/* File: "phyclust_init_method_kmed.c" and "phyclust_init_method_kmed_ex.c". */
void assign_class_by_k_medoids(int N_X, int K, double **EDM, int *center_id, int *class_id);
void assign_class_unique_by_k_medoids(int N_X_org, int K, double **EDM, int N_X_unique, int *map_X_unique_to_X_org, int *center_id, int *class_id);

/* File: "phyclust_init_method_pam.c". */
void assign_class_by_pam(int N_X_unique, int K, double **EDM_LT_pam, int *center_id, int *class_id);

#endif	/* End of __PHYCLSUT_INIT_METHOD_. */
