/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_tool.h"


/* Initial a em_phyclust structure by a phyclust structure and a QA. */
em_phyclust_struct* initialize_em_phyclust_struct(phyclust_struct *pcs){
	int i, j, N_X_org = pcs->N_X_org, N_X = pcs->N_X, L = pcs->L, K = pcs->K;
	em_phyclust_struct *empcs;

	empcs = (em_phyclust_struct*) malloc(sizeof(em_phyclust_struct));
	empcs->code_type = pcs->code_type;
	empcs->ncode = pcs->ncode;
	empcs->N_X_org = N_X_org;
	empcs->N_X_unique = pcs->N_X_unique;
	empcs->N_X = N_X;
	empcs->N_seg_site = pcs->N_seg_site;
	empcs->L = L;
	empcs->K = K;
	empcs->X_org = pcs->X_org;
	empcs->X_unique = pcs->X_unique;
	empcs->X = pcs->X;
	empcs->map_X_org_to_X_unique = pcs->map_X_org_to_X_unique;
	empcs->map_X_unique_to_X_org = pcs->map_X_unique_to_X_org;
	empcs->replication_X_unique = pcs->replication_X_unique;
	empcs->map_X_org_to_X = pcs->map_X_org_to_X;
	empcs->map_X_to_X_org = pcs->map_X_to_X_org;
	empcs->replication_X = pcs->replication_X;
	empcs->seg_site_id = pcs->seg_site_id;

	empcs->class_id = pcs->class_id;
	empcs->n_class = allocate_int_1D(K);
	empcs->Mu = allocate_int_2D_AP(K); 
	for(i = 0; i < K; i++){
		empcs->Mu[i] = allocate_int_1D(L); 
		/* This loop will be replaced by initialization functions. */
		for(j = 0; j < L; j++){
			empcs->Mu[i][j] = pcs->Mu[i][j];
		}
	}
	empcs->Z_modified = allocate_double_2D_AP(N_X); 
	for(i = 0; i < N_X; i++){
		empcs->Z_modified[i] = allocate_double_1D(K); 
	}
	empcs->Z_normalized = allocate_double_2D_AP(N_X); 
	for(i = 0; i < N_X; i++){
		empcs->Z_normalized[i] = allocate_double_1D(K); 
	}
	empcs->Z_total = allocate_double_1D(K); 
	empcs->Eta = allocate_double_1D(K); 
	empcs->log_Eta = allocate_double_1D(K); 
	empcs->count_Mu_X = allocate_int_RT_4D(N_X, pcs->K, pcs->ncode, pcs->ncode);

	for(i = 0; i < K; i++){
		empcs->Eta[i] = pcs->Eta[i];
		empcs->log_Eta[i] = log(pcs->Eta[i]);
	}
	empcs->logL_observed = 0.0;
	update_count_Mu_X(empcs);
	return(empcs);
} /* End of initialize_em_phyclust_struct(). */

void free_em_phyclust_struct(em_phyclust_struct *empcs){
	free(empcs->n_class);
	free_int_RT(empcs->K, empcs->Mu);
	free_double_RT(empcs->N_X, empcs->Z_modified);
	free_double_RT(empcs->N_X, empcs->Z_normalized);
	free(empcs->Z_total);
	free(empcs->Eta);
	free(empcs->log_Eta);
	free_int_RT_4D(empcs->N_X, empcs->K, empcs->ncode, empcs->count_Mu_X);
	free(empcs);
} /* End of free_em_phyclust_struct(). */

em_phyclust_struct* duplicate_em_phyclust_struct(em_phyclust_struct *org_empcs){
	em_phyclust_struct *new_empcs;
	int i, N_X = org_empcs->N_X, L = org_empcs->L, K = org_empcs->K;

	new_empcs = (em_phyclust_struct*) malloc(sizeof(em_phyclust_struct));
	new_empcs->code_type = org_empcs->code_type;
	new_empcs->ncode = org_empcs->ncode;

	new_empcs->N_X_org = org_empcs->N_X_org;
	new_empcs->N_X_unique = org_empcs->N_X_unique;
	new_empcs->N_X = N_X;
	new_empcs->N_seg_site = org_empcs->N_seg_site;
	new_empcs->L = L;
	new_empcs->K = K;
	new_empcs->X_org = org_empcs->X_org;
	new_empcs->X_unique = org_empcs->X_unique;
	new_empcs->X = org_empcs->X;
	new_empcs->map_X_org_to_X_unique = org_empcs->map_X_org_to_X_unique;
	new_empcs->map_X_unique_to_X_org = org_empcs->map_X_unique_to_X_org;
	new_empcs->replication_X_unique = org_empcs->replication_X_unique;
	new_empcs->map_X_org_to_X = org_empcs->map_X_org_to_X;
	new_empcs->map_X_to_X_org = org_empcs->map_X_to_X_org;
	new_empcs->replication_X = org_empcs->replication_X;
	new_empcs->seg_site_id = org_empcs->seg_site_id;

	new_empcs->class_id = org_empcs->class_id;
	new_empcs->n_class = allocate_int_1D(K);
	new_empcs->Mu = allocate_int_2D_AP(K); 
	for(i = 0; i < K; i++){
		new_empcs->Mu[i] = allocate_int_1D(L); 
	}
	new_empcs->Z_modified = allocate_double_2D_AP(N_X); 
	for(i = 0; i < N_X; i++){
		new_empcs->Z_modified[i] = allocate_double_1D(K); 
	}
	new_empcs->Z_normalized = allocate_double_2D_AP(N_X); 
	for(i = 0; i < N_X; i++){
		new_empcs->Z_normalized[i] = allocate_double_1D(K); 
	}
	new_empcs->Z_total = allocate_double_1D(K);
	new_empcs->Eta = allocate_double_1D(K);
	new_empcs->log_Eta = allocate_double_1D(K);
	new_empcs->count_Mu_X = allocate_int_RT_4D(org_empcs->N_X, org_empcs->K, org_empcs->ncode, org_empcs->ncode);

	copy_empcs(org_empcs, new_empcs);
	return(new_empcs);
} /* End of duplicate_em_phyclust_struct(). */

