/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_logpL.h"
#include "phyclust_tool.h"
#include "phyclust_init_method.h"


/* Gamma[n][k]
   = pi_k * L_k(X_n) / sum_i(pi_i * L_i(X_n))
   = 1 / sum_i(pi_i * L_i(X_n) / (pi_k * L_k))
   = 1 / sum_i(exp(log(pi_i) + log(L_i(X_n)) - log(pi_k) - log(L_k(X_n)))) */
void e_step_with_stable_exp(int *K, double *a_Z_normalized, double *total_sum, double *scale_exp, int *flag_out_range){
	int k;
	double tmp_exp, max_exp;

	*total_sum = 0.0;
	*scale_exp = 0.0;
	*flag_out_range = 0;
	max_exp = a_Z_normalized[0];
	for(k = 1; k < *K; k++){
		if(a_Z_normalized[k] > max_exp){
			max_exp = a_Z_normalized[k];
		}
	}

	/* tmp_exp = HUGE_VAL for overflow, 0 for underflow.
	 * Scale it by 2 such that close to +0 or -0 until errno is not ERANGE. */
	errno = 0;
	tmp_exp = exp(max_exp);
	if(errno == ERANGE){
		*flag_out_range = 1;
		*scale_exp = (tmp_exp == HUGE_VAL) ? max_exp : -max_exp;
		do{
			errno = 0;
			*scale_exp *= 0.5;
			tmp_exp = exp(*scale_exp);
		} while(errno == ERANGE);
		*scale_exp = max_exp - *scale_exp;
	}
	
	if(*flag_out_range){
	  for(k = 0; k < *K; k++){
	  	a_Z_normalized[k] -= *scale_exp;
	  }
	}

	*total_sum = 0.0;
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = exp(a_Z_normalized[k]);
		*total_sum += a_Z_normalized[k];
	}
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = a_Z_normalized[k] / *total_sum;
	}
} /* End of e_step_with_stable_exp(); */

void E_step_and_logL_observed_AU(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int n_X, k, K = empcs->K, flag_out_range;
	double total_sum, scale_exp;

	update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	empcs->logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized and logL_observed. */
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_modified[n_X][k] + empcs->log_Eta[k];
		}

		e_step_with_stable_exp(&K, empcs->Z_normalized[n_X], &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		empcs->logL_observed += log(total_sum);
		if(flag_out_range){
			empcs->logL_observed += scale_exp;
		}
	}
} /* End of E_step_and_logL_observed_AU(). */

void E_step_and_logL_observed_NU(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int n_X, k, K = empcs->K, flag_out_range;
	double total_sum, scale_exp;

	update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	empcs->logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized and logL_observed. */
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_modified[n_X][k] + empcs->log_Eta[k];
		}

		e_step_with_stable_exp(&K, empcs->Z_normalized[n_X], &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		if(empcs->replication_X[n_X] == 1){
			empcs->logL_observed += log(total_sum);
		} else{
			empcs->logL_observed += log(total_sum) * empcs->replication_X[n_X];
		}
		if(flag_out_range){
			total_sum += scale_exp;
		}
	}
} /* End of E_step_and_logL_observed_NU(). */

void E_step_simple(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int n_X, k, K = empcs->K, flag_out_range;
	double total_sum, scale_exp;

	update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized. */
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_modified[n_X][k] + empcs->log_Eta[k];
		}

		e_step_with_stable_exp(&K, empcs->Z_normalized[n_X], &total_sum, &scale_exp, &flag_out_range);
	}
} /* End of E_step_simple(). */

/* This is a original version of E-step, and has some numerical problems such as overflow or underflow. */
void E_step_original(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int n_X, k, K = empcs->K;
	double total_sum;

	update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized. */
		total_sum = 0.0;
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = exp(empcs->Z_modified[n_X][k] + empcs->log_Eta[k]);
			total_sum += empcs->Z_normalized[n_X][k];
		}

		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_normalized[n_X][k] / total_sum;
		}
	}
} /* End of E_step_original(). */


int M_step_simple(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP, em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA){
	int ret_stop = 0;

	ret_stop = EMFP->Update_Eta_given_Z(empcs, EMC);	/* Find Eta. */
	if(ret_stop){
		return(ret_stop);
	}
	maximize_logpL(empcs, QA, EMC, EMFP);		/* Find QA, Tt, and Mu. */
	return(ret_stop);
} /* End of M_step_simple(). */

/* All of these CM/ACM steps require to set EMC->update_flag = 1. */
int M_step_CM(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP, em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA){
	int ret_stop = 0, cm_iter = 1, flag = 1;
	double cm_reltol = 0.0, R = 0.0, tmp_R = 0.0;

	copy_empcs(empcs, tmp_empcs);
	tmp_QA->Copy_Q_matrix_array(QA, tmp_QA);

	ret_stop = EMFP->Update_Eta_given_Z(tmp_empcs, EMC);	/* Update Eta given Mu and QA. */
	if(ret_stop){
		return(ret_stop);
	}

	EMFP->Update_Mu_given_QA(tmp_empcs, tmp_QA);	/* Update Mu given Eta and QA. */ 
	update_Z_modified(tmp_empcs, tmp_QA);		/* Update empcs->Z_modified (log and unnormalized). */
	maximize_logpL(tmp_empcs, QA, EMC, EMFP);	/* Update QA given Eta and Mu. */ 
	tmp_QA->Update_log_Pt(tmp_QA);
	update_Z_modified(tmp_empcs, tmp_QA);		/* Update empcs->Z_modified (log and unnormalized). */
	tmp_R = EMFP->Compute_R(tmp_empcs);		/* Compute R(Mu, QA, Tt). */
	do{
		copy_empcs(tmp_empcs, empcs);
		tmp_QA->Copy_Q_matrix_array(tmp_QA, QA);
		R = tmp_R;

		EMFP->Update_Mu_given_QA(tmp_empcs, tmp_QA);	/* Update Mu given Eta and QA. */ 
		update_Z_modified(tmp_empcs, tmp_QA);		/* Update empcs->Z_modified (log and unnormalized). */
		maximize_logpL(tmp_empcs, tmp_QA, EMC, EMFP);	/* Update QA given Eta and Mu. */ 
		tmp_QA->Update_log_Pt(tmp_QA);
		update_Z_modified(tmp_empcs, tmp_QA);		/* Update empcs->Z_modified (log and unnormalized). */
		tmp_R = EMFP->Compute_R(tmp_empcs);		/* Compute R(Mu, QA, Tt). */

		if(tmp_R < R){
			flag = 0;
			break;
		}

		cm_reltol = fabs(tmp_R / R - 1.0);
		cm_iter++;
	} while((cm_reltol > EMC->cm_reltol) && (cm_iter < EMC->cm_maxit));

	if(flag){
		copy_empcs(tmp_empcs, empcs);
		tmp_QA->Copy_Q_matrix_array(tmp_QA, QA);
	}
	EMC->converge_cm_iter += cm_iter;
	return(ret_stop);
} /* End of M_step_CM(). */

int M_step_ACM(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP, em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA){
	int ret_stop = 0;

	ret_stop = EMFP->Update_Eta_given_Z(empcs, EMC);	/* Update Eta given Mu and QA. */
	if(ret_stop){
		return(ret_stop);
	}
	E_step_simple(empcs, QA);
	EMFP->Update_Mu_given_QA(empcs, QA);	/* Update Mu given Eta and QA. */ 
	E_step_simple(empcs, QA);
	maximize_logpL(empcs, QA, EMC, EMFP);	/* Update QA given Eta and Mu. */ 
	QA->Update_log_Pt(QA);
	EMC->converge_cm_iter++;
	return(ret_stop);
} /* End of M_step_ACM(). */


int Update_Eta_given_Z_ADJUST_AU(em_phyclust_struct *empcs, em_control *EMC){
	int n_X, k;
	double total_sum = 0.0;

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] = 0.0;
		for(n_X = 0; n_X < empcs->N_X; n_X++){
			empcs->Eta[k] += empcs->Z_normalized[n_X][k];
		}
		total_sum += empcs->Eta[k];
	}
	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] /= total_sum;
		empcs->log_Eta[k] = log(empcs->Eta[k]);
	}
	for(k = 0; k < empcs->K; k++){
		if(empcs->Eta[k] < EMC->Eta_lower_bound){
			empcs->Eta[k] = EMC->Eta_lower_bound;
		} else if(empcs->Eta[k] > EMC->Eta_upper_bound){
			empcs->Eta[k] = EMC->Eta_upper_bound;
		}
	}

	#if EMDEBUG > 2
		printf("    Eta:");
		for(k = 0; k < empcs->K; k++){
			printf(" %f", empcs->Eta[k]);
		}
		printf("\n");
	#endif

	return(0);
} /* End of Update_Eta_given_Z_ADJUST_AU(). */

int Update_Eta_given_Z_ADJUST_NU(em_phyclust_struct *empcs, em_control *EMC){
	int n_X, k;
	double total_sum = 0.0;

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] = 0.0;
		for(n_X = 0; n_X < empcs->N_X; n_X++){
			if(empcs->replication_X[n_X] == 1){
				empcs->Eta[k] += empcs->Z_normalized[n_X][k];
			} else{
				empcs->Eta[k] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
			}
		}
		total_sum += empcs->Eta[k];
	}
	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] /= total_sum;
		empcs->log_Eta[k] = log(empcs->Eta[k]);
	}
	for(k = 0; k < empcs->K; k++){
		if(empcs->Eta[k] < EMC->Eta_lower_bound){
			empcs->Eta[k] = EMC->Eta_lower_bound;
		} else if(empcs->Eta[k] > EMC->Eta_upper_bound){
			empcs->Eta[k] = EMC->Eta_upper_bound;
		}
	}

	#if EMDEBUG > 2
		printf("    Eta:");
		for(k = 0; k < empcs->K; k++){
			printf(" %f", empcs->Eta[k]);
		}
		printf("\n");
	#endif

	return(0);
} /* End of Update_Eta_given_Z_ADJUST_NU(). */


int Update_Eta_given_Z_IGNORE_AU(em_phyclust_struct *empcs, em_control *EMC){
	int n_X, k;
	double total_sum = 0.0;

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] = 0.0;
		for(n_X = 0; n_X < empcs->N_X; n_X++){
			empcs->Eta[k] += empcs->Z_normalized[n_X][k];
		}
		total_sum += empcs->Eta[k];
	}
	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] /= total_sum;
		empcs->log_Eta[k] = log(empcs->Eta[k]);
	}
	for(k = 0; k < empcs->K; k++){
		if(empcs->Eta[k] < EMC->Eta_lower_bound || empcs->Eta[k] > EMC->Eta_upper_bound){
			#if EMDEBUG > 2
				printf("  Eta[%d]=%f is out of (%f, %f)\n", k, empcs->Eta[k], EMC->Eta_lower_bound, EMC->Eta_upper_bound);
			#endif
			return(1);
		}
	}

	#if EMDEBUG > 2
		printf("    Eta:");
		for(k = 0; k < empcs->K; k++){
			printf(" %f", empcs->Eta[k]);
		}
		printf("\n");
	#endif

	return(0);
} /* End of Update_Eta_given_Z_IGNORE_AU(). */

int Update_Eta_given_Z_IGNORE_NU(em_phyclust_struct *empcs, em_control *EMC){
	int n_X, k;
	double total_sum = 0.0;

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] = 0.0;
		for(n_X = 0; n_X < empcs->N_X; n_X++){
			if(empcs->replication_X[n_X] == 1){
				empcs->Eta[k] += empcs->Z_normalized[n_X][k];
			} else{
				empcs->Eta[k] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
			}
		}
		total_sum += empcs->Eta[k];
	}
	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] /= total_sum;
		empcs->log_Eta[k] = log(empcs->Eta[k]);
	}
	for(k = 0; k < empcs->K; k++){
		if(empcs->Eta[k] < EMC->Eta_lower_bound || empcs->Eta[k] > EMC->Eta_upper_bound){
			#if EMDEBUG > 2
				printf("  Eta[%d]=%f is out of (%f, %f)\n", k, empcs->Eta[k], EMC->Eta_lower_bound, EMC->Eta_upper_bound);
			#endif
			return(1);
		}
	}

	#if EMDEBUG > 2
		printf("    Eta:");
		for(k = 0; k < empcs->K; k++){
			printf(" %f", empcs->Eta[k]);
		}
		printf("\n");
	#endif

	return(0);
} /* End of Update_Eta_given_Z_IGNORE_NU(). */


/* For each sequence, compute log_Pt for all (n, k) cells and save in empcs->Z_modified (log and unnormalized). */
/* ToDo: extend count_Mu_X[empcs->ncode, empcs->ncode] to count_Mu_X[empcs->ncode, empcs->ncode_gap] to allow GAP. */
void update_Z_modified(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			empcs->Z_modified[n_X][k] = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					empcs->Z_modified[n_X][k] += QA->Q[k]->log_Pt[s_from][s_to] * empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
		}
	}
} /* End of update_Z_modified(). */


void print_status(em_phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, int verbosity){
	int k;
	if(!verbosity) return;
	if(verbosity==1){
		printf(".");
		/* fflush(stdout); */
		return;
	}
	if(verbosity==2){
		printf("%5d %12.3f\n", EMC->converge_iter, pcs->logL_observed);
		return;
	}
	if(verbosity==3) {
		printf("%5d eta", EMC->converge_iter);
		for(k = 0; k<pcs->K; k++){
			printf(" %6.4f", pcs->Eta[k]);
		}
		print_QA(QA);
		printf(" %12.3f\n", pcs->logL_observed);
		return;
	}
} /* End of print_status(). */


/* converge_flag = 0, successed to converge.
 *               = 1, out of iterations.
 *               = 2, maybe converged where logL_observed decreasing or Eta's go to 0.
 *               = 3, fail to converge.
 * update_flag = 0 for update Mu given QA. (for profile logL)
 *             = 1 for update QA given Mu. (for profile logL and ECM/AECM) */
int Check_convergence_em(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, em_control *EMC){
	int ret_stop = 0;
	if(new_empcs->logL_observed < org_empcs->logL_observed){	/* logL decreasing. */
		if(EMC->update_flag == 0){				/* change to update Q, Tt given Mu. */
			copy_empcs(org_empcs, new_empcs);		/* Repalce object with higher logL. */
			EMC->update_flag = 1;
		} else{							/* otherwise stop. */
			EMC->converge_flag = 9;
			EMC->converge_error = new_empcs->logL_observed - org_empcs->logL_observed;
			ret_stop = 1;
		}
	} else{
		if(EMC->update_flag == 1){				/* logL increasing and current is update QA, Tt given Mu*/
			EMC->update_flag = 0;				/* change to original plan. */
		}
	}
	return(ret_stop);
} /* End of Check_convergence_em(). */

int Check_convergence_org(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, em_control *EMC){
	int ret_stop = 0;
	if(new_empcs->logL_observed < org_empcs->logL_observed){	/* logL decreasing. */
		EMC->converge_flag = 9;
		EMC->converge_error = new_empcs->logL_observed - org_empcs->logL_observed;
		ret_stop = 1;
	}
	return(ret_stop);
} /* End of Check_convergence_org(). */


void Em_step(em_phyclust_struct *org_empcs, Q_matrix_array *org_QA, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	em_phyclust_struct *new_empcs, *tmp_empcs;
	Q_matrix_array *new_QA, *tmp_QA;

	reset_em_control(EMC);
	new_QA = duplicate_Q_matrix_array(org_QA);
	new_empcs = duplicate_em_phyclust_struct(org_empcs);
	tmp_QA = duplicate_Q_matrix_array(org_QA);
	tmp_empcs = duplicate_em_phyclust_struct(org_empcs);

	#if EMDEBUG > 0
		int k;
		printf("Start: EM\n");
		printf("  Eta:");
		for(k = 0; k < new_empcs->K; k++){
			printf(" %.4f", new_empcs->Eta[k]);
		}
		printf("\n");
		print_QA(new_QA);
	#endif

	EMC->update_flag = (EMC->em_method == EM) ? 0 : 1;
	EMFP->E_step_and_logL_observed(new_empcs, new_QA);
	#if EMDEBUG > 0
		double tmp_logL;
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter - 1, EMC->update_flag);
		tmp_logL = EMFP->LogL_observed(new_empcs, new_QA);
		if(is_finite(new_empcs->logL_observed)){
			printf("  logL: %.8f, %.8f [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
		} else{
			printf("  logL: %.8e, %.8e [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
		}
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
	#endif
	do{
		#if verbosity_em_step > 0
			print_status(new_empcs, new_QA, EMC, verbosity_em_step);
		#endif
		copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(new_QA, org_QA);

		#if EMDEBUG > 1
			double tmp_R, tmp_Q, tmp_obs;
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("  M-step:\n");
			printf("    init.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif
		ret_stop = EMFP->M_step(new_empcs, new_QA, EMC, EMFP, tmp_empcs, tmp_QA);
		if(ret_stop){
			EMC->converge_flag = 2;	/* Eta < 1/N or > 1 - 1/N. */
			break;
		}
		#if EMDEBUG > 1
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("    conv.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif

		EMFP->E_step_and_logL_observed(new_empcs, new_QA);

		EMC->converge_eps = fabs(new_empcs->logL_observed / org_empcs->logL_observed - 1.0);
		EMC->converge_iter++;
		#if EMDEBUG > 0
			tmp_logL = EMFP->LogL_observed(new_empcs, new_QA);
			if(is_finite(new_empcs->logL_observed)){
				printf("  logL: %.8f, %.8f [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			} else{
				printf("  logL: %.8e, %.8e [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			}
			printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
		#endif

		ret_stop = EMFP->Check_convergence(new_empcs, org_empcs, EMC);
		if(ret_stop){
			break;	/* logL is decreasing and fails after switch. */
		}
	} while((EMC->converge_eps > EMC->EM_eps) && (EMC->converge_iter < EMC->EM_iter));

	#if verbosity_em_step > 1
		printf("\n");
	#endif

	if(EMC->converge_iter > EMC->EM_iter){
		EMC->converge_flag = 1;
	}

	if(EMC->converge_flag < 2){
		copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(new_QA, org_QA);
	}

	free_Q_matrix_array(new_QA);
	free_em_phyclust_struct(new_empcs);
	free_Q_matrix_array(tmp_QA);
	free_em_phyclust_struct(tmp_empcs);
} /* End of Em_step(). */

void Short_em_step(em_phyclust_struct *org_empcs, Q_matrix_array *org_QA, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	em_phyclust_struct *new_empcs, *tmp_empcs;
	Q_matrix_array *new_QA, *tmp_QA;
	double logL_0;

	reset_em_control(EMC);
	new_QA = duplicate_Q_matrix_array(org_QA);
	new_empcs = duplicate_em_phyclust_struct(org_empcs);
	tmp_QA = duplicate_Q_matrix_array(org_QA);
	tmp_empcs = duplicate_em_phyclust_struct(org_empcs);

	#if EMDEBUG > 0
		int k;
		double tmp_logL;
		printf("Start: EM\n");
		printf("  Eta:");
		for(k = 0; k < new_empcs->K; k++){
			printf(" %.4f", new_empcs->Eta[k]);
		}
		printf("\n");
		print_QA(new_QA);
	#endif

	EMC->update_flag = (EMC->em_method == EM) ? 0 : 1;
	EMFP->E_step_and_logL_observed(new_empcs, new_QA);
	logL_0 = new_empcs->logL_observed;
	#if EMDEBUG > 0
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter - 1, EMC->update_flag);
		update_Z_modified(org_empcs, org_QA);
		tmp_logL = EMFP->LogL_observed(org_empcs, org_QA);
		if(is_finite(org_empcs->logL_observed)){
			printf("  logL: %.8f, %.8f [indep fcn]\n", org_empcs->logL_observed, tmp_logL);
		} else{
			printf("  logL: %.8e, %.8e [indep fcn]\n", org_empcs->logL_observed, tmp_logL);
		}
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
	#endif
	do{
		copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(new_QA, org_QA);

		#if EMDEBUG > 1
			double tmp_R, tmp_Q, tmp_obs;
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("  M-step:\n");
			printf("    init.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif
		ret_stop = EMFP->M_step(new_empcs, new_QA, EMC, EMFP, tmp_empcs, tmp_QA);
		if(ret_stop){
			EMC->converge_flag = 2;	/* Eta < 1/N or > 1 - 1/N. */
			break;
		}
		#if EMDEBUG > 1
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("    conv.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif

		EMFP->E_step_and_logL_observed(new_empcs, new_QA);

		EMC->converge_eps = (org_empcs->logL_observed - new_empcs->logL_observed) / (logL_0 - new_empcs->logL_observed);
		EMC->converge_iter++;
		#if EMDEBUG > 0
			tmp_logL = EMFP->LogL_observed(new_empcs, new_QA);
			if(is_finite(new_empcs->logL_observed)){
				printf("  logL: %.8f, %.8f [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			} else{
				printf("  logL: %.8e, %.8e [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			}
			printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
		#endif

		ret_stop = EMFP->Check_convergence(new_empcs, org_empcs, EMC);
		if(ret_stop){
			break;	/* logL is decreasing and fails after switch. */
		}
	} while((EMC->converge_eps > EMC->short_eps) && (EMC->converge_iter < EMC->short_iter));

	if(EMC->converge_iter > EMC->short_iter){
		EMC->converge_flag = 1;
	}

	if(EMC->converge_flag < 2){
		copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(new_QA, org_QA);
	}

	free_Q_matrix_array(new_QA);
	free_em_phyclust_struct(new_empcs);
	free_Q_matrix_array(tmp_QA);
	free_em_phyclust_struct(tmp_empcs);
} /* End of Short_em_step(). */

