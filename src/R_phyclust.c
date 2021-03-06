/* This file contains a function R_phyclust() called by R wraps phyclust() in
 * "R/f_phyclust.r", and this function calls the relative functions in
 * "src/phyclust/" and "src/".
 *
 * Writen: Wei-Chen Chen on 2009/08/20. */

#include "R_phyclust.h"

/* For getting list elements. Copied from R document. */
SEXP getListElement(SEXP list, const char *str){
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return(elmt);
} /* End of getListElement(). */

void copy_R_EMC_to_EMC(SEXP R_EMC, em_control *EMC){
	EMC->exhaust_iter = INTEGER(getListElement(R_EMC, "exhaust.iter"))[0];
	EMC->fixed_iter = INTEGER(getListElement(R_EMC, "fixed.iter"))[0];
	EMC->short_iter = INTEGER(getListElement(R_EMC, "short.iter"))[0];
	EMC->EM_iter = INTEGER(getListElement(R_EMC, "EM.iter"))[0];
	EMC->short_eps = REAL(getListElement(R_EMC, "short.eps"))[0];
	EMC->EM_eps = REAL(getListElement(R_EMC, "EM.eps"))[0];

	EMC->cm_reltol = REAL(getListElement(R_EMC, "cm.reltol"))[0];
	EMC->cm_maxit = INTEGER(getListElement(R_EMC, "cm.maxit"))[0];
	
	EMC->nm_abstol_Mu_given_QA = REAL(getListElement(R_EMC, "nm.abstol.Mu.given.QA"))[0];
	EMC->nm_reltol_Mu_given_QA = REAL(getListElement(R_EMC, "nm.reltol.Mu.given.QA"))[0];
	EMC->nm_maxit_Mu_given_QA = INTEGER(getListElement(R_EMC, "nm.maxit.Mu.given.QA"))[0];
	EMC->nm_abstol_QA_given_Mu = REAL(getListElement(R_EMC, "nm.abstol.QA.given.Mu"))[0];
	EMC->nm_reltol_QA_given_Mu = REAL(getListElement(R_EMC, "nm.reltol.QA.given.Mu"))[0];
	EMC->nm_maxit_QA_given_Mu = INTEGER(getListElement(R_EMC, "nm.maxit.QA.given.Mu"))[0];
	EMC->est_non_seg_site = INTEGER(getListElement(R_EMC, "est.non.seg.site"))[0];

	EMC->max_init_iter = INTEGER(getListElement(R_EMC, "max.init.iter"))[0];
	EMC->init_procedure = INTEGER(getListElement(R_EMC, "init.procedure"))[0];
	EMC->init_method = INTEGER(getListElement(R_EMC, "init.method"))[0];
	EMC->substitution_model = INTEGER(getListElement(R_EMC, "substitution.model"))[0];
	EMC->edist_model = INTEGER(getListElement(R_EMC, "edist.model"))[0];
	EMC->identifier = INTEGER(getListElement(R_EMC, "identifier"))[0];
	EMC->code_type = INTEGER(getListElement(R_EMC, "code.type"))[0];
	EMC->em_method = INTEGER(getListElement(R_EMC, "em.method"))[0];
	EMC->boundary_method = INTEGER(getListElement(R_EMC, "boundary.method"))[0];

	EMC->min_n_class = INTEGER(getListElement(R_EMC, "min.n.class"))[0];
} /* End of copy_R_EMC_to_EMC(). */

EMPTR allocate_emptr(void){
	EMPTR emptr = (EMPTR) malloc(sizeof(struct _emptr));
	if(emptr == NULL){
		printf("Memory allocation fails!\n");
		exit(1);
	}
	return(emptr);
} /* End of allocate_emptr(). */

SEXP initialize_emptr(EMPTR emptr, phyclust_struct *pcs){
	SEXP emobj, emobj_names, QA_names, converge_names;
	SEXP N_X_org, N_X, L, K, Eta, Z_normalized, Mu, QA, logL, p, bic, aic, icl, N_seg_site,
	     class_id, n_class, converge, label_method;
	SEXP pi, kappa, Tt;
	SEXP converge_eps, converge_error, converge_flag, converge_iter, converge_inner_iter, converge_cm_iter, check_param;
	char *names_emobj[] = {"N.X.org", "N.X", "L", "K", "Eta", "Z.normalized", "Mu", "QA", "logL", "p",
				"bic", "aic", "icl", "N.seg.site", "class.id", "n.class", "conv", "label.method"};
	char *names_QA[] = {"pi", "kappa", "Tt"};
	char *names_converge[] = {"eps", "error", "flag", "iter", "inner.iter", "cm.iter", "check.param"};
	int emobj_length = 18, QA_length = 3, converge_length = 7;
	int i, j, *tmp_ptr_int;
	double *tmp_ptr_double;

	/* Allocate and protect storages. */
  	PROTECT(emobj = allocVector(VECSXP, emobj_length));
  	PROTECT(emobj_names = allocVector(STRSXP, emobj_length));
  	PROTECT(QA_names = allocVector(STRSXP, QA_length));
  	PROTECT(converge_names = allocVector(STRSXP, converge_length));
	
  	PROTECT(N_X_org = allocVector(INTSXP, 1));
  	PROTECT(N_X = allocVector(INTSXP, 1));
  	PROTECT(L = allocVector(INTSXP, 1));
  	PROTECT(K = allocVector(INTSXP, 1));
  	PROTECT(Eta = allocVector(REALSXP, pcs->K));
	PROTECT(Z_normalized = allocVector(REALSXP, pcs->N_X_org * pcs->K));
  	PROTECT(Mu = allocVector(INTSXP, pcs->K * pcs->L));
  	PROTECT(QA = allocVector(VECSXP, QA_length));
	  	PROTECT(pi = allocVector(REALSXP, pcs->ncode * pcs->K));
  		PROTECT(kappa = allocVector(REALSXP, 1 * pcs->K));
  		PROTECT(Tt = allocVector(REALSXP, 1 * pcs->K));
  	PROTECT(logL = allocVector(REALSXP, 1));
  	PROTECT(p = allocVector(INTSXP, 1));
  	PROTECT(bic = allocVector(REALSXP, 1));
  	PROTECT(aic = allocVector(REALSXP, 1));
  	PROTECT(icl = allocVector(REALSXP, 1));
  	PROTECT(N_seg_site = allocVector(INTSXP, 1));
  	PROTECT(class_id = allocVector(INTSXP, pcs->N_X_org));
  	PROTECT(n_class = allocVector(INTSXP, pcs->K));
  	PROTECT(converge = allocVector(VECSXP, converge_length));
		PROTECT(converge_eps = allocVector(REALSXP, 1));
		PROTECT(converge_error = allocVector(REALSXP, 1));
		PROTECT(converge_flag = allocVector(INTSXP, 1));
		PROTECT(converge_iter = allocVector(INTSXP, 1));
		PROTECT(converge_inner_iter = allocVector(INTSXP, 1));
		PROTECT(converge_cm_iter = allocVector(INTSXP, 1));
		PROTECT(check_param = allocVector(INTSXP, 1));
	PROTECT(label_method = allocVector(INTSXP, 1));

	/* Set the elements and names. */
	i = 0;
	SET_VECTOR_ELT(emobj, i++, N_X_org);
	SET_VECTOR_ELT(emobj, i++, N_X);
	SET_VECTOR_ELT(emobj, i++, L);
	SET_VECTOR_ELT(emobj, i++, K);
	SET_VECTOR_ELT(emobj, i++, Eta);
	SET_VECTOR_ELT(emobj, i++, Z_normalized);
	SET_VECTOR_ELT(emobj, i++, Mu);
	SET_VECTOR_ELT(emobj, i++, QA);
	SET_VECTOR_ELT(emobj, i++, logL);
	SET_VECTOR_ELT(emobj, i++, p);
	SET_VECTOR_ELT(emobj, i++, bic);
	SET_VECTOR_ELT(emobj, i++, aic);
	SET_VECTOR_ELT(emobj, i++, icl);
	SET_VECTOR_ELT(emobj, i++, N_seg_site);
	SET_VECTOR_ELT(emobj, i++, class_id);
	SET_VECTOR_ELT(emobj, i++, n_class);
	SET_VECTOR_ELT(emobj, i++, converge);
	SET_VECTOR_ELT(emobj, i++, label_method);

	i = 0;
	SET_VECTOR_ELT(QA, i++, pi);
	SET_VECTOR_ELT(QA, i++, kappa);
	SET_VECTOR_ELT(QA, i++, Tt);

	i = 0;
	SET_VECTOR_ELT(converge, i++, converge_eps);
	SET_VECTOR_ELT(converge, i++, converge_error);
	SET_VECTOR_ELT(converge, i++, converge_flag);
	SET_VECTOR_ELT(converge, i++, converge_iter);
	SET_VECTOR_ELT(converge, i++, converge_inner_iter);
	SET_VECTOR_ELT(converge, i++, converge_cm_iter);
	SET_VECTOR_ELT(converge, i++, check_param);
	
	for(i = 0; i < emobj_length; i++){
		SET_STRING_ELT(emobj_names, i, mkChar(names_emobj[i])); 
	}
	setAttrib(emobj, R_NamesSymbol, emobj_names);
	for(i = 0; i < QA_length; i++){
		SET_STRING_ELT(QA_names, i, mkChar(names_QA[i])); 
	}
	setAttrib(QA, R_NamesSymbol, QA_names);
	for(i = 0; i < converge_length; i++){
		SET_STRING_ELT(converge_names, i, mkChar(names_converge[i])); 
	}
	setAttrib(converge, R_NamesSymbol, converge_names);

	/* Assign emptr. */
	tmp_ptr_int = INTEGER(Mu);
	for(i = 0; i < pcs->K; i++){
		pcs->Mu[i] = tmp_ptr_int;
		tmp_ptr_int += pcs->L;
		for(j = 0; j < pcs->L; j++){
			pcs->Mu[i][j] = 0;
		}
	}
	pcs->Eta = REAL(Eta);
	for(i = 0; i < pcs->K; i++){
		pcs->Eta[i] = 1 / (double) pcs->K;
	}
	tmp_ptr_double = REAL(Z_normalized);
	for(i = 0; i < pcs->N_X_org; i++){
		pcs->Z_normalized[i] = tmp_ptr_double;
		tmp_ptr_double += pcs->K;
		for(j = 0; j < pcs->K; j++){
			pcs->Z_normalized[i][j] = 0;
		}
	}

	pcs->class_id = INTEGER(class_id);
	pcs->n_class = INTEGER(n_class);

	emptr->C_N_X_org = INTEGER(N_X_org);
	emptr->C_N_X = INTEGER(N_X);
	emptr->C_L = INTEGER(L);
	emptr->C_K = INTEGER(K);
	emptr->C_logL = REAL(logL);
	emptr->C_p = INTEGER(p);
	emptr->C_bic = REAL(bic);
	emptr->C_aic = REAL(aic);
	emptr->C_icl = REAL(icl);
	emptr->C_N_seg_site = INTEGER(N_seg_site);
	emptr->C_pi = REAL(pi);
	emptr->C_kappa = REAL(kappa);
	emptr->C_Tt = REAL(Tt);
	emptr->C_converge_eps = REAL(converge_eps);
	emptr->C_converge_error = REAL(converge_error);
	emptr->C_converge_flag = INTEGER(converge_flag);
	emptr->C_converge_iter = INTEGER(converge_iter);
	emptr->C_converge_inner_iter = INTEGER(converge_inner_iter);
	emptr->C_converge_cm_iter = INTEGER(converge_cm_iter);
	emptr->C_check_param = INTEGER(check_param);
	emptr->C_label_method = INTEGER(label_method);

	/* emptr->C_protect_length = 4 + emobj_length + QA_length + converge_length; */

	/* Do NOT call UNPROTECT() within this constructor!! */
	// UNPROTECT(emptr->C_protect_length);
	UNPROTECT(32);
	return(emobj);
} /* End of initialize_emptr(). */

void copy_all_to_emptr(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, EMPTR emptr){
	int i, k, i2;
	*emptr->C_N_X_org = pcs->N_X_org;
	*emptr->C_N_X = pcs->N_X;
	*emptr->C_L = pcs->L;
	*emptr->C_K = pcs->K;
	*emptr->C_logL = pcs->logL_observed;
	*emptr->C_p = pcs->n_param + QA->total_n_param;
	*emptr->C_bic = pcs->bic;
	*emptr->C_aic = pcs->aic;
	*emptr->C_icl = pcs->icl;
	*emptr->C_N_seg_site = pcs->N_seg_site;
	i2 = 0;
	for(k = 0; k < pcs->K; k++){
		for(i = 0; i < pcs->ncode; i++){
			emptr->C_pi[i2++] = QA->Q[k]->pi[i];
		}
	}
	for(k = 0; k < pcs->K; k++){
		emptr->C_kappa[k] = *QA->Q[k]->kappa;
		emptr->C_Tt[k] = *QA->Q[k]->Tt;
	}
	*emptr->C_converge_eps = EMC->converge_eps;
	*emptr->C_converge_error = EMC->converge_error;
	*emptr->C_converge_flag = EMC->converge_flag;
	*emptr->C_converge_iter = EMC->converge_iter;
	*emptr->C_converge_inner_iter = EMC->converge_inner_iter;
	*emptr->C_converge_cm_iter = EMC->converge_cm_iter;
	*emptr->C_check_param = QA->check_param;
	*emptr->C_label_method = pcs->label->label_method;
} /* End of copy_all_to_emptr(). */




/* This function calls init_em_step() in
 * "src/phyclust/phyclust_init_procedure.c" and is
 * called by phyclust() using .Call() in "R/f_phyclust.r".
 * Input:
 *   R_N: SEXP[1], number of sequences.
 *   R_L: SEXP[1], length of sequences.
 *   R_K: SEXP[1], number of clusters.
 *   R_X: SEXP[1], sequences.
 *   R_EMC: SEXP[1], EM controler.
 *   R_manual_id: SEXP[1], manual class id.
 *   R_label: SEXP[1], labels.
 * Output:
 *   ret: a list contains everythings returned from phyclust in C. */
SEXP R_phyclust(SEXP R_N_X_org, SEXP R_L, SEXP R_K, SEXP R_X, SEXP R_EMC, SEXP R_manual_id, SEXP R_label){
	/* Declare variables for calling phyclust. */
	int *C_N_X_org, *C_L, *C_K, *C_manual_id;
	em_control *EMC;
	phyclust_struct *pcs;
	Q_matrix_array *QA;
	em_fp *EMFP;

	/* Declare variables for R's returning. */
	EMPTR emptr = allocate_emptr();
	SEXP emobj;
	/* int C_protect_length; */

	/* Declare variables for processing. */
	int i, *tmp_ptr;

	/* Set initial values. */
	C_N_X_org = INTEGER(R_N_X_org);
	C_L = INTEGER(R_L);
	C_K = INTEGER(R_K);
	C_manual_id = INTEGER(R_manual_id);

	/* Assign controler. */
	EMC = initialize_em_control();
	copy_R_EMC_to_EMC(R_EMC, EMC);
	update_em_control(EMC);

	/* Assign data. */
	pcs = R_initialize_phyclust_struct(EMC->code_type, *C_N_X_org, *C_L, *C_K);
	PROTECT(emobj = initialize_emptr(emptr, pcs));		 /* !! Don't move this. */
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->X_org[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}
	if(EMC->init_method == manualMu){
		for(i = 0; i < *C_N_X_org; i++){
			pcs->class_id[i] = C_manual_id[i];
		}
	}
	update_phyclust_struct(pcs);

	/* Assign labels. */
	R_update_phyclust_label(pcs, R_label);

	/* Assign function pointers. */
	EMFP = initialize_em_fp(EMC, pcs);

	/* Assign QA. */
	QA = initialize_Q_matrix_array(EMC->code_type, *C_K, EMC->substitution_model, EMC->identifier);

	/* Compute. */
	init_em_step(pcs, QA, EMC, EMFP);
	assign_class(pcs);
	update_ic(pcs, QA);

	/* For return. */
	copy_all_to_emptr(pcs, QA, EMC, emptr);

	/* Free memory and release protectation. */
	free_em_control(EMC);
	R_free_phyclust_struct(pcs);
	free_em_fp(EMFP);
	free_Q_matrix_array(QA);
	/* C_protect_length = emptr->C_protect_length; */
	free(emptr);

	/* UNPROTECT(C_protect_length); */
	UNPROTECT(1);
	return(emobj);
} /* End of SEXP R_phyclust(). */

