/* This file contains functions about log profile likelihood. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_logpL.h"
#include "phyclust_optim_nmmin.h"
#include "phyclust_tool.h"


/* This function will update Mu given QA for full length sequences. */
void Update_Mu_given_QA_full_AU(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int org_s_mu, new_s_mu, s, n_X, k, l, flag = 0;
	int K = empcs->K, L = empcs->L, Ncode = empcs->ncode, N_X = empcs->N_X;
	double org_psi_k, new_psi_k, psi_k[empcs->ncode];

	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			for(s = 0; s < Ncode; s++){
				psi_k[s] = 0.0;
			}
			for(n_X = 0; n_X < N_X; n_X++){
				psi_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
			}

			org_s_mu = empcs->Mu[k][l];
			org_psi_k = 0.0;
			for(s = 0; s < Ncode; s++){
				org_psi_k += psi_k[s] * QA->Q[k]->log_Pt[org_s_mu][s];
			}
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				new_psi_k = 0.0;
				for(s = 0; s < Ncode; s++){
					new_psi_k += psi_k[s] * QA->Q[k]->log_Pt[new_s_mu][s];
				}
				if(new_psi_k > org_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					org_psi_k = new_psi_k;
					flag |= 1;
				}
			}
			/* Dynamically update empcs->count_Mu_X:
			 * If empcs->Mu[k][l] were changed from org_s_mu to empcs->Mu[k][l]=new_s_mu,
			 * then count_Mu_X[][][][] should be updated by
			 *      count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]-- and
			 *      count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++
			 * for all n_X = 1, ..., N_X. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_AU(). */

void Update_Mu_given_QA_full_NU(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int org_s_mu, new_s_mu, s, n_X, k, l, flag = 0;
	int K = empcs->K, L = empcs->L, Ncode = empcs->ncode, N_X = empcs->N_X;
	double org_psi_k, new_psi_k, psi_k[empcs->ncode];

	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			for(s = 0; s < Ncode; s++){
				psi_k[s] = 0.0;
			}
			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->replication_X[n_X] == 1){
					psi_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
				} else{
					psi_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
				}
			}

			org_s_mu = empcs->Mu[k][l];
			org_psi_k = 0.0;
			for(s = 0; s < Ncode; s++){
				org_psi_k += psi_k[s] * QA->Q[k]->log_Pt[org_s_mu][s];
			}
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				new_psi_k = 0.0;
				for(s = 0; s < Ncode; s++){
					new_psi_k += psi_k[s] * QA->Q[k]->log_Pt[new_s_mu][s];
				}
				if(new_psi_k > org_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					org_psi_k = new_psi_k;
					flag |= 1;
				}
			}
			/* Dynamically update empcs->count_Mu_X:
			 * If empcs->Mu[k][l] were changed from org_s_mu to empcs->Mu[k][l]=new_s_mu,
			 * then count_Mu_X[][][][] should be updated by
			 *      count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]-- and
			 *      count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++
			 * for all n_X = 1, ..., N_X. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_NU(). */

/* This function will update Mu given QA, but skip non-segregating sites. */
void Update_Mu_given_QA_skip_non_seg_AU(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int org_s_mu, new_s_mu, s, n_X, k, l, i, flag = 0;
	int K = empcs->K, N_seg_site = empcs->N_seg_site, Ncode = empcs->ncode, N_X = empcs->N_X;
	double org_psi_k, new_psi_k, psi_k[empcs->ncode];

	for(k = 0; k < K; k++){
		for(i = 0; i < N_seg_site; i++){
			l = empcs->seg_site_id[i];

			for(s = 0; s < Ncode; s++){
				psi_k[s] = 0.0;
			}
			for(n_X = 0; n_X < N_X; n_X++){
				psi_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
			}

			org_s_mu = empcs->Mu[k][l];
			org_psi_k = 0.0;
			for(s = 0; s < Ncode; s++){
				org_psi_k += psi_k[s] * QA->Q[k]->log_Pt[org_s_mu][s];
			}
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				new_psi_k = 0.0;
				for(s = 0; s < Ncode; s++){
					new_psi_k += psi_k[s] * QA->Q[k]->log_Pt[new_s_mu][s];
				}
				if(new_psi_k > org_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					org_psi_k = new_psi_k;
					flag |= 1;
				}
			}
			/* Dynamically update empcs->count_Mu_X:
			 * If empcs->Mu[k][l] were changed from org_s_mu to empcs->Mu[k][l]=new_s_mu,
			 * then count_Mu_X[][][][] should be updated by
			 *      count_Mu_X[n][k][org_s_mu][empcs->X[n][l]]-- and
			 *      count_Mu_X[n][k][empcs->Mu[k][l]][empcs->X[n][l]]++
			 * for all n = 1, ..., N_X and where l's should belong to segregating sites. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_skip_non_seg_AU(). */

void Update_Mu_given_QA_skip_non_seg_NU(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int org_s_mu, new_s_mu, s, n_X, k, l, i, flag = 0;
	int K = empcs->K, N_seg_site = empcs->N_seg_site, Ncode = empcs->ncode, N_X = empcs->N_X;
	double org_psi_k, new_psi_k, psi_k[empcs->ncode];

	for(k = 0; k < K; k++){
		for(i = 0; i < N_seg_site; i++){
			l = empcs->seg_site_id[i];

			for(s = 0; s < Ncode; s++){
				psi_k[s] = 0.0;
			}
			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->replication_X[n_X] == 1){
					psi_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
				} else{
					psi_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
				}
			}

			org_s_mu = empcs->Mu[k][l];
			org_psi_k = 0.0;
			for(s = 0; s < Ncode; s++){
				org_psi_k += psi_k[s] * QA->Q[k]->log_Pt[org_s_mu][s];
			}
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				new_psi_k = 0.0;
				for(s = 0; s < Ncode; s++){
					new_psi_k += psi_k[s] * QA->Q[k]->log_Pt[new_s_mu][s];
				}
				if(new_psi_k > org_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					org_psi_k = new_psi_k;
					flag |= 1;
				}
			}
			/* Dynamically update empcs->count_Mu_X:
			 * If empcs->Mu[k][l] were changed from org_s_mu to empcs->Mu[k][l]=new_s_mu,
			 * then count_Mu_X[][][][] should be updated by
			 *      count_Mu_X[n][k][org_s_mu][empcs->X[n][l]]-- and
			 *      count_Mu_X[n][k][empcs->Mu[k][l]][empcs->X[n][l]]++
			 * for all n = 1, ..., N_X and where l's should belong to segregating sites. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_skip_non_seg_NU(). */




double Compute_R_AU(em_phyclust_struct *empcs){
	int n_X, k;
	double ret = 0.0;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			/* empcs->Z_normalized[n_X][k] is not updated in m-step, and only
			 * empcs->Z_modified[n_X][k] is updated by update_Z_modified(). */
			ret += empcs->Z_normalized[n_X][k] * empcs->Z_modified[n_X][k];
		}
	}

	return(ret);
} /* End of Compute_R_AU(). */

double Compute_R_NU(em_phyclust_struct *empcs){
	int n_X, k;
	double ret = 0.0, tmp_ret;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		tmp_ret = 0.0;
		for(k = 0; k < empcs->K; k++){
			/* empcs->Z_normalized[n_X][k] is not updated in m-step, and only
			 * empcs->Z_modified[n_X][k] is updated by update_Z_modified(). */
			tmp_ret += empcs->Z_normalized[n_X][k] * empcs->Z_modified[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			ret += tmp_ret;
		} else{
			ret += tmp_ret * empcs->replication_X[n_X];
		}
	}

	return(ret);
} /* End of Compute_R_NU(). */




/* Negative of profile log-likelihood for minimizing, -R(Mu, QA, Tt).
 * parameters: eta_1, ..., eta_(K-1), QA(model)
 * model JC69: Tt.
 *       K80: kappa, Tt.
 *       HKY85: pi_A, pi_G, pi_C, kappa, Tt.
 * vect stores parameters: QA(model).
 */
double negative_logpL_Mu_given_QA(int m, double *vect, void *ex){
	ex_struct *in = (ex_struct*) ex;
	double ret = 0.0;

	in->QA->Convert_vect_to_Q_matrix_array(vect, in->QA);	/* Update QA, Tt and check_param. */

	/* empcs->Z_normalized is fixed and should not be updated in M-step.
	 * QA and Tt are updated by NM. Mu are updated given QA, Tt and empcs->Z_normalized. */
	if(in->QA->check_param){
		in->QA->Update_log_Pt(in->QA);				/* Update log(P(t)). */
		in->EMFP->Update_Mu_given_QA(in->empcs, in->QA);	/* Update Mu. */
		update_Z_modified(in->empcs, in->QA);			/* Update empcs->Z_modified (log and unnormalized). */
		ret = -in->EMFP->Compute_R(in->empcs);			/* Compute R(Mu, QA, Tt). */
	} else{
		ret = Inf;				/* NM failed, return Inf to stop NM. */
	}

	#if EMDEBUG > 3
		printf("    Update Mu given QA\n");
		printf("      vect:");
		print_vect(m, vect);
		if(is_finite(ret)){
			printf("    -logpL: %.4f\n", ret);
		} else{
			printf("    -logpL: %.4e\n", ret);
		}
	#endif

	return(ret);
} /* End of negative_logpL_Mu(). */


double negative_logpL_QA_given_Mu(int m, double *vect, void *ex){
	ex_struct *in = (ex_struct*) ex;
	double ret = 0.0;

	in->QA->Convert_vect_to_Q_matrix_array(vect, in->QA);	/* Update QA, Tt and check_param. */

	/* empcs->Z_normalized is fixed and should not be updated in M-step.
	 * QA and Tt are updated by NM. Mu are updated given QA, Tt and empcs->Z_normalized. */
	if(in->QA->check_param){
		in->QA->Update_log_Pt(in->QA);		/* Update log(P(t)). */
		update_Z_modified(in->empcs, in->QA);	/* Update empcs->Z_modified (log and unnormalized). */
		ret = -in->EMFP->Compute_R(in->empcs);	/* Compute R(Mu, QA, Tt). */
	} else{
		ret = Inf;				/* NM is fail, return Inf to stop NM. */
	}

	#if EMDEBUG > 3
		printf("    Update QA given Mu\n");
		printf("      vect:");
		print_vect(m, vect);
		if(is_finite(ret)){
			printf("    -logpL: %.8f\n", ret);
		} else{
			printf("    -logpL: %.8e\n", ret);
		}
	#endif

	return(ret);
} /* End of negative_logpL_QA_given_Mu(). */


/* Maximize profile log-likelihood R(Mu, QA, Tt). */
void maximize_logpL(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	ex_struct exs;
	nm_struct *nms;
	double *vect = allocate_double_1D(QA->total_n_param);

	QA->Convert_Q_matrix_array_to_vect(QA, vect);
	exs.empcs = empcs;
	exs.EMFP = EMFP;
	exs.QA = QA;

	nms = initialize_nm_struct(QA->total_n_param);
	nms->Bvec = vect;
	nms->ex = &exs;
	if(EMC->update_flag == 0){
		nms->fminfn = &negative_logpL_Mu_given_QA;
		nms->abstol = EMC->nm_abstol_Mu_given_QA;
		nms->reltol = EMC->nm_reltol_Mu_given_QA;
		nms->maxit = EMC->nm_maxit_Mu_given_QA;
	} else{
		nms->fminfn = &negative_logpL_QA_given_Mu;
		nms->abstol = EMC->nm_abstol_QA_given_Mu;
		nms->reltol = EMC->nm_reltol_QA_given_Mu;
		nms->maxit = EMC->nm_maxit_QA_given_Mu;
	}

	#if EMDEBUG > 2
		int i, j;
		double tmp_i, tmp_logpL, tmp_R, tmp_Q, tmp_obs;
		tmp_logpL = -nms->fminfn(nms->n_param, nms->Bvec, nms->ex);
		tmp_R = EMFP->LogL_profile(empcs, QA);
		tmp_Q = EMFP->LogL_complete(empcs, QA);
		tmp_obs = EMFP->LogL_observed(empcs, QA);
		printf("    maximize_logpL:\n");
		printf("    init. coefs:");
		print_vect(nms->n_param, vect);
		if(is_finite(tmp_logpL) && is_finite(tmp_R) && is_finite(tmp_Q) && is_finite(tmp_obs)){
			printf("        logpL: %.6f\n", tmp_logpL);
			printf("        R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
			for(i = 0 ; i < empcs->N_X; i++){
				tmp_i = 0.0;
				for(j = 0; j < empcs->K; j++){
					tmp_i += empcs->Z_normalized[i][j] * empcs->log_Eta[j];
				}
				if(empcs->replication_X[i] == 1){
					tmp_logpL += tmp_i;
				} else{
					tmp_logpL += tmp_i * empcs->replication_X[i];
				}
			}
			printf("        R+head terms = %.6f. [indep fcn]\n", tmp_logpL);
		} else{
			printf("        logpL: %.6e\n", tmp_logpL);
			printf("        R = %.6e, Q = %.6e, Obs = %.6e. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		}
	#endif

	phyclust_optim_nmmin(nms);

	#if EMDEBUG > 2
		tmp_logpL = -nms->fminfn(nms->n_param, nms->Bvec, nms->ex);
		tmp_R = EMFP->LogL_profile(empcs, QA);
		tmp_Q = EMFP->LogL_complete(empcs, QA);
		tmp_obs = EMFP->LogL_observed(empcs, QA);
		printf("    conv. coefs:");
		print_vect(nms->n_param, vect);
		if(is_finite(tmp_logpL) && is_finite(tmp_R) && is_finite(tmp_Q) && is_finite(tmp_obs)){
			printf("        logpL: %.6f\n", tmp_logpL);
			printf("        R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
			for(i = 0 ; i < empcs->N_X; i++){
				tmp_i = 0.0;
				for(j = 0; j < empcs->K; j++){
					tmp_i += empcs->Z_normalized[i][j] * empcs->log_Eta[j];
				}
				if(empcs->replication_X[i] == 1){
					tmp_logpL += tmp_i;
				} else{
					tmp_logpL += tmp_i * empcs->replication_X[i];
				}
			}
			printf("        R+head terms = %.6f.\n", tmp_logpL);
		} else{
			printf("        logpL: %.6e.\n", tmp_logpL);
			printf("        R = %.6e, Q = %.6e, Obs = %.6e. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		}
		printf("    convergence: %d\n", *nms->fail);
	#endif
	
	EMC->converge_inner_iter += *nms->fncount;
	/*
	QA->Convert_vect_to_Q_matrix_array(vect, QA);
	QA->Update_log_Pt(QA);
	*/

	free(vect);
	free_nm_struct(nms);
} /* End of maximize_logpL(). */




/* ----- For debug. ----- */
void print_vect(int m, double *vect){
	int i;

	for(i = 0; i < m; i++){
		printf(" %.4f", vect[i]);
	}
	printf("\n");
} /* End of print_vect(). */

