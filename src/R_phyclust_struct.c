/* This file contains functions for initialization and memory deallocation.
 * Copy from "phyclust_struct.c" to avoid coping some memory. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust/phyclust.h"


/* Initial a phyclust structure without assigning data X.
 * Assign X later and call update_phyclust_struct(). */
phyclust_struct* R_initialize_phyclust_struct(int code_type, int N_X_org, int L, int K){
	phyclust_struct *pcs = NULL;

	pcs = (phyclust_struct*) malloc(sizeof(phyclust_struct));
	pcs->code_type = code_type;
	pcs->ncode = NCODE[code_type];
	pcs->missing_index = MISSING_INDEX[code_type];
	pcs->missing_flag = 0;
	pcs->n_param = K - 1 + K * L;
	pcs->N_X_org = N_X_org;
	pcs->N_X_unique = 0;
	pcs->N_X = 0;
	pcs->N_seg_site = 0;
	pcs->L = L;
	pcs->K = K;
	pcs->X_org = allocate_int_2D_AP(N_X_org);
	pcs->X_unique = NULL;
	pcs->X = NULL;
	pcs->map_X_org_to_X_unique = NULL;
	pcs->map_X_unique_to_X_org = NULL;
	pcs->replication_X_unique = NULL;
	pcs->map_X_org_to_X = NULL;
	pcs->map_X_to_X_org = NULL;
	pcs->replication_X = NULL;
	pcs->seg_site_id = NULL;

	pcs->Mu = allocate_int_2D_AP(K);			/* Assigned by R. */
	pcs->Eta = NULL;					/* Assigned by R. */
	pcs->Z_normalized = allocate_double_2D_AP(N_X_org);	/* Assigned by R. */

	pcs->logL_observed = 0.0;
	pcs->logL_entropy = 0.0;
	pcs->bic = 0.0;
	pcs->aic = 0.0;
	pcs->icl = 0.0;
	pcs->class_id = NULL; 
	pcs->n_class = NULL; 

	return(pcs);
} /* End of R_initialize_phyclust_struct(). */

void R_free_phyclust_struct(phyclust_struct *pcs){
	free(pcs->X_org);
	free(pcs->X_unique);
	free(pcs->map_X_org_to_X_unique);
	free(pcs->map_X_unique_to_X_org);
	free(pcs->replication_X_unique);
	free(pcs->seg_site_id);
	free(pcs->Mu);
	free(pcs->Z_normalized);
	free(pcs);
} /* End of free_phyclust_struct(). */

