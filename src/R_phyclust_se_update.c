/* This file is modified from "R_phyclust_em_step.c" for updating the results
 * of phyclust() to  phyclust.se() and find.best.se() in R.
 *
 * Writen: Wei-Chen Chen on 2012/03/01. */


#include "R_phyclust_se.h"

/* Input:
 *   R_N_X_org: SEXP[1], number of sequences.
 *   R_L: SEXP[1], length of sequences.
 *   R_X: SEXP[1], sequences.
 *   R_K: SEXP[1], number of clusters.
 *   R_Eta: SEXP[1], Eta.
 *   R_Mu: SEXP[1], Mu.
 *   R_vect: SEXP[1], vect contains pi, kappa, and Tt.
 *   R_substitution_model: SEXP[1], substitution model.
 *   R_identifier: SEXP[1], identifier.
 *   R_code_type: SEXP[1], code_type.
 *   R_label: SEXP[1], labeles.
 *   R_se_model: SEX[1], se_model.
 *   R_se_constant: SEX[1], se_constant.
 * Output:
 *   ret: a list contains everythings returned from phyclust in C. */
SEXP R_phyclust_se_update(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label, SEXP R_se_model, SEXP R_se_constant){
	/* Declare variables for calling C. */
	int *C_N_X_org, *C_L, *C_K;
	double *C_vect;
	em_control *EMC;
	phyclust_struct *pcs;
	Q_matrix_array *QA;
	em_phyclust_struct *empcs;
	em_fp *EMFP;

	/* Declare variables for R's returning. */
	EMPTR_SE emptr = allocate_emptr_se();
	SEXP emobj;
	int C_protect_length;
	
	/* Declare variables for processing. */
	int i, j, *tmp_ptr;
	double *tmp_ptr_double;

	/* Set initial values. */
	C_N_X_org = INTEGER(R_N_X_org);
	C_L = INTEGER(R_L);
	C_K = INTEGER(R_K);
	C_vect = REAL(R_vect);

	/* Assign controler. */
	EMC = initialize_em_control();
	EMC->substitution_model = INTEGER(R_substitution_model)[0];
	EMC->identifier = INTEGER(R_identifier)[0];
	EMC->code_type = INTEGER(R_code_type)[0];
	EMC->em_method = EM;				/* For se model. */
	EMC->se_type = SE_YES;				/* For se model. */
	EMC->se_model = INTEGER(R_se_model)[0];		/* For se model. */
	EMC->se_constant = REAL(R_se_constant)[0];	/* For se model. */
	update_em_control(EMC);

	/* Assign data, read only. */
	pcs = R_initialize_phyclust_struct(EMC->code_type, *C_N_X_org, *C_L, *C_K);
	emobj = initialize_emptr_se(emptr, pcs);		/* !! Don't move this. */
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->X_org[i] = tmp_ptr;			/* Assign poiners. */
		tmp_ptr += *C_L;
	}

	/* Assign parameters. Updates are required, so make a copy. */
	tmp_ptr = INTEGER(R_Mu);				/* Read only. */
	for(i = 0; i < *C_K; i++){
		for(j = 0; j < *C_L; j++){
			pcs->Mu[i][j] = *tmp_ptr;		/* Copy from the original input. */
			tmp_ptr++;
		}
	}
	tmp_ptr_double = REAL(R_Eta);				/* Read only. */
	for(i = 0; i < *C_K; i++){
		pcs->Eta[i] = *tmp_ptr_double;			/* Copy from the original input. */
		tmp_ptr_double++;
	}
	update_phyclust_struct(pcs);

	/* Assign labels. */
	R_update_phyclust_label(pcs, R_label);

	/* Assign function pointers. */
	EMFP = initialize_em_fp(EMC, pcs);

	/* Assign QA. */
	QA = initialize_Q_matrix_array(EMC->code_type, *C_K, EMC->substitution_model, EMC->identifier);
	QA->Convert_vect_to_Q_matrix_array(C_vect, QA);	/* Copy from the original input. */
	QA->Update_log_Pt(QA);

	/* Compute for se. */
	if(EMC->code_type == NUCLEOTIDE){
		update_phyclust_se_struct(pcs, EMC);
		update_em_fp_se(EMFP, EMC, pcs);

		/* Initialize empcs. */
		empcs = initialize_em_phyclust_struct(pcs);
	
		/* EM steps. */
		EMFP->Em_step(empcs, QA, EMC, EMFP);
		EMFP->Copy_empcs_to_pcs(empcs, pcs);

		/* Print results. */
		assign_class(pcs);
		update_ic(pcs, QA);

		free_em_phyclust_struct(empcs);
	}

	/* For return. */
	copy_all_to_emptr_se(pcs, QA, EMC, emptr);

	/* Free memory and release protectation. */
	free_em_control(EMC);
	free_phyclust_se_struct(pcs);
	R_free_phyclust_struct(pcs);
	free_em_fp(EMFP);
	free_Q_matrix_array(QA);
	C_protect_length = emptr->C_protect_length;
	free(emptr);

	UNPROTECT(C_protect_length);
	return(emobj);
} /* End of SEXP R_phyclust_se_update(). */
