/* This file contains declarations for em steps. */

#ifndef __PHYCLSUT_EM_TOOL_
#define __PHYCLSUT_EM_TOOL_

#include "phyclust_struct.h"
#include "phyclust_em.h"
#include "phyclust_qmatrix_array.h"

/* ----- Summary tool. ----- */
double LogL_observed_AU(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_observed_NU(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_complete_AU(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_complete_NU(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_profile_AU(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_profile_NU(em_phyclust_struct *empcs, Q_matrix_array *QA);

void update_count_Mu_X(em_phyclust_struct *empcs);

/* ----- Checking tool. ----- */
int is_finite(double x);

/* ----- For copy. ----- */
void copy_EMC(em_control *EMC_from, em_control *EMC_to);
void copy_empcs(em_phyclust_struct *empcs_from, em_phyclust_struct *empcs_to);

void Copy_empcs_to_pcs_AU(em_phyclust_struct *empcs, phyclust_struct *pcs);
void Copy_empcs_to_pcs_NU(em_phyclust_struct *empcs, phyclust_struct *pcs);
void Copy_pcs_to_empcs_AU(phyclust_struct *pcs, em_phyclust_struct *empcs);
void Copy_pcs_to_empcs_NU(phyclust_struct *pcs, em_phyclust_struct *empcs);


/* ----- For debug. ----- */
void print_empcs(em_phyclust_struct *empcs);
void print_EMC(em_control *EMC);
void print_result(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC);
void print_rich_result(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC);
void print_Z_modified(em_phyclust_struct *empcs);
void print_Z_normalized(em_phyclust_struct *empcs);
void print_Eta(em_phyclust_struct *empcs);
void print_empcs_Mu(em_phyclust_struct *empcs);
void print_empcs_Mu_seg_site(em_phyclust_struct *empcs);
void print_count_Mu_X(em_phyclust_struct *empcs, int n_X, int k);

#endif	/* End of __PHYCLSUT_EM_TOOL_. */
