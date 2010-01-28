/* This file contains all function pointers required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_init_method.h"
#include "phyclust_logpL.h"

/* Initial a em_fp. */
em_fp* initialize_em_fp(em_control *EMC, phyclust_struct *pcs){
	em_fp *EMFP = (em_fp*) malloc(sizeof(em_fp));
	
	/* Assign init_method FP. */
	switch(EMC->init_method){
		case randomMu:
			EMFP->Update_init = &Update_init_random_Mu_unique;
			break;
		case NJ:
			EMFP->Update_init = &Update_init_nj_unique;
			break;
		case randomNJ:
			EMFP->Update_init = &Update_init_random_nj_unique;
			break;
		case PAM:
			EMFP->Update_init = &Update_init_pam;
			break;
		case kMedoids:
			EMFP->Update_init = &Update_init_k_medoids;
			break;
		case manualMu:
			EMFP->Update_init = &Update_init_manually;
			break;
		default:
			fprintf(stderr, "PE: The initial method is not found.\n");
			exit(1);
	}

	/* Assign EM method. */
	switch(EMC->em_method){
		case EM:
			EMFP->M_step = &M_step_simple;
			EMFP->Check_convergence = &Check_convergence_em;
			EMFP->Em_step = &Em_step;
			EMFP->Short_em_step = &Short_em_step;
			break;
		case ECM:
			EMFP->M_step = &M_step_CM;
			EMFP->Check_convergence = &Check_convergence_org;
			EMFP->Em_step = &Em_step;
			EMFP->Short_em_step = &Short_em_step;
			break;
		case AECM:
			EMFP->M_step = &M_step_ACM;
			EMFP->Check_convergence = &Check_convergence_org;
			EMFP->Em_step = &Em_step;
			EMFP->Short_em_step = &Short_em_step;
			break;
		default:
			fprintf(stderr, "PE: The EM method is not found.\n");
			exit(1);
	}

	/* Update functions. */
	if(pcs->N_X < pcs->N_X_org){	/* Some are not unique. */
		/* In "phyclust_em.c". */
		EMFP->E_step_and_logL_observed = &E_step_and_logL_observed_NU;
		switch(EMC->boundary_method){
			case IGNORE:
				EMFP->Update_Eta_given_Z = &Update_Eta_given_Z_IGNORE_NU;
				break;
			case ADJUST:
				EMFP->Update_Eta_given_Z = &Update_Eta_given_Z_ADJUST_NU;
				break;
			default:
				fprintf(stderr, "PE: The boundary method is not found.\n");
				exit(1);
		}
		/* In "phyclust_em_tool.h". */
		EMFP->LogL_observed = &LogL_observed_NU;
		EMFP->LogL_complete = &LogL_complete_NU;
		EMFP->LogL_profile = &LogL_profile_NU;
		EMFP->Copy_empcs_to_pcs = &Copy_empcs_to_pcs_NU;
		EMFP->Copy_pcs_to_empcs = &Copy_pcs_to_empcs_NU;
		/* In "phyclust_logpL.c". */
		if(EMC->est_non_seg_site != 0){
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_full_NU;
		} else{
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_skip_non_seg_NU;
		}
		EMFP->Compute_R = &Compute_R_NU;
	} else{				/* All are unique. */
		/* In "phyclust_em.c". */
		EMFP->E_step_and_logL_observed = &E_step_and_logL_observed_AU;
		switch(EMC->boundary_method){
			case IGNORE:
				EMFP->Update_Eta_given_Z = &Update_Eta_given_Z_IGNORE_AU;
				break;
			case ADJUST:
				EMFP->Update_Eta_given_Z = &Update_Eta_given_Z_ADJUST_AU;
				break;
			default:
				fprintf(stderr, "PE: The boundary method is not found.\n");
				exit(1);
		}
		/* In "phyclust_em_tool.h". */
		EMFP->LogL_observed = &LogL_observed_AU;
		EMFP->LogL_complete = &LogL_complete_AU;
		EMFP->LogL_profile = &LogL_profile_AU;
		EMFP->Copy_empcs_to_pcs = &Copy_empcs_to_pcs_AU;
		EMFP->Copy_pcs_to_empcs = &Copy_pcs_to_empcs_AU;
		/* In "phyclust_logpL.c". */
		if(EMC->est_non_seg_site != 0){
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_full_AU;
		} else{
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_skip_non_seg_AU;
		}
		EMFP->Compute_R = &Compute_R_AU;
	}
	
	return(EMFP);
} /* End of initialize_em_control(). */


void free_em_fp(em_fp *EMFP){
	free(EMFP);
} /* End of free_em_fp(). */

