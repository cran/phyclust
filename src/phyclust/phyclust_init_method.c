/* This file contains all fucntions for initialization. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_edist.h"
#include "phyclust_em_tool.h"
#include "phyclust_init_method.h"
#include "phyclust_logpL.h"
#include "phyclust_tool.h"
#include "phyclust_ape_nj.h"


/* Define double rdunif(). */
#ifdef RMATH_H
	#ifdef MATHLIB_STANDALONE
		/* Require set_seed(SEED1, SEED2) & get_seed(SEED1, SEED2) to call this. */
		#include <Rmath.h>
		int rdunif(int n){
			return((int) floor(n * unif_rand()));
		}
	#else
		#include <R.h>
		#include <Rmath.h>
		int rdunif(int n){
			int ret = 0;
			GetRNGstate();
			ret = (int) floor(n * runif(0, 1));
			PutRNGstate();
			return(ret);
		}
	#endif
#else
	/* Require srand(seed) to call this. */
	#include <math.h>
	int rdunif(int n){
		return((int) floor(n * (double) rand() / ((double) RAND_MAX + 1.0)));
	}
#endif


/* Simple random sample withor replacement. Choose k from 0 to n-1.
 * Modified from Dr. Maitra's code. */
void srswor(int n, int k, int *x){
	int i, j;
	int *tmp_x = allocate_int_1D(n);

	for(i = 0; i < n; i++){
		tmp_x[i] = i;
	}

	for(i = 0; i < k; i++){
		j = rdunif(n);
		x[i] = tmp_x[j];
		tmp_x[j] = tmp_x[--n];
	}

	free(tmp_x);
}


/* Copy from phyclust_em.c and modified for initialization method. */
int init_m_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	update_count_Mu_X(empcs); /* Update count_Mu_X. */
	ret_stop = EMFP->Update_Eta_given_Z(empcs, EMC); /* Find Eta. */
	if(ret_stop){
		return(ret_stop);
	}

	/* Find QA. */
	EMC->update_flag = 1;		/* For update QA, given Mu. */
	maximize_logpL(empcs, QA, EMC, EMFP);
	QA->Update_log_Pt(QA);
	EMC->update_flag = 0;		/* Reset to 0 for update Mu, given QA. */
	return(ret_stop);
} /* End of init_m_step(). */

int check_all_min_n_class(int K, int *n_class, int min_n_class){
	int k, ret = 1;

	for(k = 0; k < K; k++){
		ret &= (n_class[k] >= min_n_class);
	}

	return(ret);
} /* End of check_n_class(). */

void assign_Mu_by_class(int N_X_org, int K, int L, int ncode, int *class_id, int **X_org, int **Mu){
	int i, n_X_org, k, l;
	int count_N[ncode], tmp_count, tmp_Mu;

	for(l = 0; l < L; l++){
		/* Find the most common nucleotide in all sequences for replacing when tie in each k. */
		for(i = 0; i < ncode; i++){
			count_N[i] = 0;
		}
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			count_N[X_org[n_X_org][l]]++;
		}
		tmp_count = -1;
		tmp_Mu = 0;
		for(i = 0; i < ncode; i++){
			if(count_N[i] > tmp_count){
				tmp_count = count_N[i];
				tmp_Mu = i;
			}
		}

		/* Find the most common nucleotide for each k. */
		for(k = 0; k < K; k++){
			for(i = 0; i < ncode; i++){
				count_N[i] = 0;
			}
			for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
				if(class_id[n_X_org] == k){
					count_N[X_org[n_X_org][l]]++;
				}
			}
			tmp_count = -1;
			for(i = 0; i < ncode; i++){
				if(count_N[i] > tmp_count){
					tmp_count = count_N[i];
					Mu[k][l] = i;
				} else if(count_N[i] == tmp_count && (Mu[k][l] == tmp_Mu || i == tmp_Mu)){
					tmp_count = count_N[i];
					Mu[k][l] = tmp_Mu;	/* Return the most common nucleotide if tie. */
				}
			}
		}
	}
} /* End of assign_Mu_by_class(). */




void Update_init_manually(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int n_X, k, ret_stop = 0;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			empcs->Z_normalized[n_X][k] = 0.0;
		}
		empcs->Z_normalized[n_X][empcs->class_id[empcs->map_X_to_X_org[n_X]]] = 1.0;
	}

	reset_Q_matrix_array(QA);
	assign_Mu_by_class(empcs->N_X_org, empcs->K, empcs->L, empcs->ncode, empcs->class_id, empcs->X_org, empcs->Mu);
	ret_stop = init_m_step(empcs, QA, EMC, EMFP);
	if(ret_stop){
		fprintf(stderr, "PE: Initialization error.\n");
		exit(1);	
	}
	if(!is_finite(EMFP->LogL_observed(empcs, QA))){
		fprintf(stderr, "PE: manual initialization leads to non-finite observed log likelihood\n");
		exit(1);
	}
} /* End of Update_init_manually(). */


/* Randomly pick Mu from X. */
void Update_init_random_Mu_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X, k, l, N_X_unique = empcs->N_X_unique, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], tmp_id;
	double tmp, tmp_min, init_logL_observed;
	edist_struct *eds;

	eds = initialize_edist_struct_UT(EMC->edist_model, N_X, L, empcs->X);
	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);

		/* Randomly pick mu from X_unique. */
		srswor(N_X_unique, K, center_id);
		
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				empcs->Mu[k][l] = empcs->X_unique[center_id[k]][l];
			}
			empcs->n_class[k] = 0;
			center_id[k] = empcs->map_X_org_to_X[empcs->map_X_unique_to_X_org[center_id[k]]];
		}

		/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
		for(n_X = 0; n_X < N_X; n_X++){
			tmp_min = eds->get_pair_edist(eds, n_X, center_id[0]);
			tmp_id = 0;
			for(k = 1; k < K; k++){
				tmp = eds->get_pair_edist(eds, n_X, center_id[k]);
				if(tmp < tmp_min){
					tmp_min = tmp;
					tmp_id = k;
				}
			}
	
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][tmp_id] = 1.0;
			empcs->n_class[tmp_id] += empcs->replication_X[n_X];
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			ret_stop = init_m_step(empcs, QA, EMC, EMFP);
			if(ret_stop){
				continue;
			}
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop){
			fprintf(stderr, "PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			fprintf(stderr, "PE: Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	}

	free_edist_struct(eds);
} /* End of Update_init_random_Mu_unique(). */


/* Pick clusters by cutting the longest K internal branches of neighbor-joining tree.
 * There is no randomness for this method, so the following settings may be suggested.
 * EMC->init_procedure = exhaustEM;
 * EMC->exhaust_iter = 1; */
void Update_init_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int n_X_org, n_X_unique, n_X, k, ret_stop = 0;
	int N_X_org = empcs->N_X_org, N_X_unique = empcs->N_X_unique, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int largest_branch_id[N_X_unique - 3], class_id[N_X_org], class_id_unique[N_X_unique];
	double init_logL_observed;
	edist_struct *eds;
	nj_struct *njs;
	
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X_unique, L, empcs->X_unique);
	njs = initialize_nj_struct(N_X_unique);
	njs->D = eds->EDM[0];
	phyclust_ape_nj(njs);
	if(! check_njs(njs)){
		fprintf(stderr, "PE: NJ may be not valid!\n");
		print_njs(njs->n_edge, njs);
		exit(1);
	}
	#if INITDEBUG > 0
		print_njs(njs->n_edge, njs);
	#endif

	search_largest_branch(njs, largest_branch_id);
	ret_stop = assign_class_by_njs_branch(K, njs, largest_branch_id, class_id_unique);
	if(ret_stop != 0){
		fprintf(stderr, "PE: Class assignment fails.\n");
		exit(1);
	}

	/* Assign Mu and recreate Z_normalized. */
	for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
		class_id[n_X_org] = class_id_unique[empcs->map_X_org_to_X_unique[n_X_org]];
	}
	assign_Mu_by_class(empcs->N_X_org, empcs->K, empcs->L, empcs->ncode, empcs->class_id, empcs->X_org, empcs->Mu);
	for(k = 0; k < K; k++){
		empcs->n_class[k] = 0;
	}
	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = 0.0;
		}
		empcs->Z_normalized[n_X][class_id_unique[empcs->map_X_org_to_X_unique[empcs->map_X_to_X_org[n_X]]]] = 1.0;
	}
	for(n_X_unique = 0; n_X_unique < N_X_unique; n_X_unique++){
		empcs->n_class[class_id_unique[n_X_unique]] += empcs->replication_X_unique[n_X_unique];
	}

	if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop){
			fprintf(stderr, "PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			fprintf(stderr, "PE: Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	} else{
		fprintf(stderr, "PE: Initialization is not valid for min_n_class = %d. (%s)\n", EMC->min_n_class, INIT_METHOD[EMC->init_method]);
		exit(1);
	}

	free_edist_struct(eds);
	free_nj_struct(njs);
} /* End of Update_init_nj_unique(). */


/* Pick clusters by randomly cutting the longest K internal branches of neighbor-joining tree. */
void Update_init_random_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X_org, n_X_unique, n_X, k;
	int N_X_org = empcs->N_X_org, N_X_unique = empcs->N_X_unique, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int random_branch_id[N_X_unique - 3], class_id[N_X_org], class_id_unique[N_X_unique];
	double init_logL_observed;
	edist_struct *eds;
	nj_struct *njs;
	
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X_unique, L, empcs->X_unique);
	njs = initialize_nj_struct(N_X_unique);
	njs->D = eds->EDM[0];
	phyclust_ape_nj(njs);
	if(! check_njs(njs)){
		fprintf(stderr, "PE: NJ may be not valid!\n");
		print_njs(njs->n_edge, njs);
		exit(1);
	}
	#if INITDEBUG > 0
		print_njs(njs->n_edge, njs);
	#endif

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);

		/* Randomly pick mu from X. */
		random_branch(njs, random_branch_id);
		ret_stop = assign_class_by_njs_branch(K, njs, random_branch_id, class_id_unique);
		if(ret_stop != 0){
			continue;
		}

		/* Assign Mu and recreate Z_normalized. */
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			class_id[n_X_org] = class_id_unique[empcs->map_X_org_to_X_unique[n_X_org]];
		}
		assign_Mu_by_class(empcs->N_X_org, empcs->K, empcs->L, empcs->ncode, empcs->class_id, empcs->X_org, empcs->Mu);
		for(k = 0; k < K; k++){
			empcs->n_class[k] = 0;
		}
		for(n_X = 0; n_X < N_X; n_X++){
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][class_id_unique[empcs->map_X_org_to_X_unique[empcs->map_X_to_X_org[n_X]]]] = 1.0;
		}
		for(n_X_unique = 0; n_X_unique < N_X_unique; n_X_unique++){
			empcs->n_class[class_id_unique[n_X_unique]] += empcs->replication_X_unique[n_X_unique];
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			ret_stop = init_m_step(empcs, QA, EMC, EMFP);
			if(ret_stop){
				continue;
			}
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop){
			fprintf(stderr, "PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			fprintf(stderr, "PE: Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	}

	free_edist_struct(eds);
	free_nj_struct(njs);
} /* End of Update_init_random_nj_unique(). */


void Update_init_k_medoids(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X_org, n_X, k, l, N_X_org = empcs->N_X_org, N_X_unique = empcs->N_X_unique;
	int N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X_org];
	double init_logL_observed;
	edist_struct *eds;

	eds = initialize_edist_struct_UT(EMC->edist_model, N_X_org, L, empcs->X_org);

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);

		assign_class_unique_by_k_medoids(N_X_org, K, eds->EDM, N_X_unique, empcs->map_X_unique_to_X_org, center_id, class_id);

		/* Pick mu from X. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				empcs->Mu[k][l] = empcs->X_org[center_id[k]][l];
			}
			empcs->n_class[k] = 0;
		}

		/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
		for(n_X = 0; n_X < N_X; n_X++){
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][class_id[empcs->map_X_to_X_org[n_X]]] = 1.0;
		}
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			empcs->n_class[class_id[n_X_org]]++;
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			init_m_step(empcs, QA, EMC, EMFP);
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop){
			fprintf(stderr, "PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			fprintf(stderr, "PE: Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	}

	free_edist_struct(eds);
} /* End of Update_init_k_medoids(). */


/* Pick clusters by PAM.
 * There is no randomness for this method, so the following settings may be suggested.
 * EMC->init_procedure = exhaustEM;
 * EMC->exhaust_iter = 1; */
void Update_init_pam(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	int n_X_org, n_X, k, l, N_X_org = empcs->N_X_org, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X_org];
	double init_logL_observed;
	edist_struct *eds;
	
	eds = initialize_edist_struct_LT_pam(EMC->edist_model, N_X_org, L, empcs->X_org);
	assign_class_by_pam(N_X_org, K, eds->EDM, center_id, class_id);

	/* Pick mu from X. */
	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			empcs->Mu[k][l] = empcs->X_org[center_id[k]][l];
		}
		empcs->n_class[k] = 0;
	}

	/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = 0.0;
		}
		empcs->Z_normalized[n_X][class_id[empcs->map_X_to_X_org[n_X]]] = 1.0;
	}
	for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
		empcs->n_class[class_id[n_X_org]]++;
	}

	if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop){
			fprintf(stderr, "PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			fprintf(stderr, "PE: Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	} else{
		fprintf(stderr, "PE: Initialization is not valid for min_n_class = %d. (%s)\n", EMC->min_n_class, INIT_METHOD[EMC->init_method]);
		exit(1);
	}

	free_edist_struct(eds);
} /* End of Update_init_pam(). */
