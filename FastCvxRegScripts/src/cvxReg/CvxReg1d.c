/**
 * @file CvxReg1d.c
 * @author Ang Li
 * @date 2017-11-6
 * @brief Convex regression with the ADMM algorithm.
 * The main function reads simulated data from csv, run convex regression
 * and write the fitted curve, duality_gap series into csv files
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#include "tf.h"
#include "utils.h"

#include "isoReg.h"
#include "csvIO.h"

/**
* @brief              1-d convex regression
*
* @param n            the length of input vector
* @param y            a vector of responses
* @param x            a vector of x
* @param w            a vector of weights
* @param output       allocated space for output
* @param duality_gap  allocated space for the time (iteration) series of duality gap
* @param max_iter     max number of ADMM iterations
* @param rho          tuning parameter for the ADMM
* @param obj_tol      stopping criteria tolerance
*/
void cr_1d_admm(int n, double *y, double *x, double *w, double *output, double *duality_gap, int max_iter, double rho, double obj_tol)
{
    int it;
    int i;

    double *beta;
    double *z;
    double *u;

    double *tmp1;
    double *tmp2;
    double *Dtu;

    cs *D1;
    cs *D1t;
    cs *D1tD1;
    cs *kernmat;
    gqr *kernmat_qr;

    beta = (double *) malloc(n * sizeof(double));
    z = (double *) malloc((n - 1) * sizeof(double));
    u = (double *) malloc((n - 1) * sizeof(double));

    tmp1 = (double *) malloc((n-1) * sizeof(double));
    tmp2 = (double *) malloc(n * sizeof(double));
    Dtu = (double *) malloc((n-1) * sizeof(double));

    D1 	= tf_calc_dktil(n, 1, x);
    D1t = cs_transpose(D1, 1);
    D1tD1 = cs_multiply(D1t, D1);
    kernmat = scalar_plus_diag(D1tD1, rho, w);
    kernmat_qr = glmgen_qr(kernmat);

    /* Initialize z, u */
    for (i = 0; i < n - 1; ++i) {
    	z[i] = 0;
    	u[i] = 0;
    }

    for (it = 0; it < max_iter; ++it){
        /* Update beta */
        for (i = 0; i < n-1; ++i) tmp1[i] = z[i] - u[i];
        tf_dtxtil(x, n, 1, tmp1, tmp2);
        for (i = 0; i < n; ++i) beta[i] = w[i] * y[i] + rho * tmp2[i];
        /* Solve the least squares problem with sparse QR */
        glmgen_qrsol(kernmat_qr, beta);

        /* Update z */
        tf_dxtil(x, n, 1, beta, tmp1);
        for (i = 0; i < n-1; i++) tmp1[i] += u[i]; /* tmp1 = D(1)beta + u */
        /* Use PAVA */
        iso_reg(n, tmp1, z);

        /* Update u: dual update */
        for (i = 0; i < n-1; i++) u[i] = tmp1[i] - z[i];

        /* Stop if duality_gap < obj_tol */
        /* Note that weight is not being considered here (refine later) */
        for (i = 0; i < n-1; i++) Dtu[i] = u[i];
        tf_dtxtil(x, n, 1, Dtu, Dtu);
        for (i = 0; i < n; i++) duality_gap[it] += (y[i] - beta[i])*(y[i] - beta[i])/2 + (rho - rho*rho/2)*Dtu[i]*Dtu[i] - Dtu[i]*y[i];
        /* Rescale the duality gap with ((x_n - x_1)/n)^2 */
        for (i = 0; i < 2; i++) duality_gap[it] *= (x[n-1] - x[0])/n;
        printf("Duality_gap = %f\n", duality_gap[it]);
        if (it > 0 && fabs(duality_gap[it]) < obj_tol) break;

    }
    printf("Actual num of iterations = %d\n", it + 1);

    memcpy(output, beta, n * sizeof(double));

    free(beta);
    free(z);
    free(u);
    free(tmp1);
    free(tmp2);
    free(Dtu);

    cs_spfree(D1);
    cs_spfree(D1t);
    cs_spfree(D1tD1);
    cs_spfree(kernmat);
    glmgen_gqr_free(kernmat_qr);
}


/**
* @brief              retrieve data from input file, run convex regression, and write results into output files
*
* @param n            the length of input vector
* @param max_iter     max number of ADMM iterations
* @param rho          tuning parameter for the ADMM
* @param obj_tol      stopping criteria tolerance
* @param input_filename         input file contains a table of 3 columns (response, x, weight)
* @param output_filename_1      output fitted curve
* @param output_filename_2      output duality gap history
* @param verbose                if true, print the input data and output curve
*/
void convex_regression(int n, int max_iter, double rho, double obj_tol, char *input_filename, char *output_filename_1, char *output_filename_2, bool verbose)
{
	double *y = (double *) malloc(n * sizeof(double));
	double *x = (double *) malloc(n * sizeof(double));
	double *w = (double *) malloc(n * sizeof(double));

	double *output = (double *) malloc(n * sizeof(double));
	double *duality_gap = (double *) malloc(max_iter * sizeof(double));

	read_csv(n, y, x, w, input_filename);

	if(verbose){
		printf("input:\n");
		for(int i = 0; i < n; ++i) {
			printf("%f %f %f\n", y[i], x[i], w[i]);
		}
		printf("\n");
	}

	cr_1d_admm(n, y, x, w, output, duality_gap, max_iter, rho, obj_tol);
	create_csv(output_filename_1, output, n);
	create_csv(output_filename_2, duality_gap, max_iter);

	if(verbose){
		printf("output:\n");
		for(int i = 0; i < n; ++i) {
			printf("%f\n", output[i]);
		}
	}

	free(y);
	free(x);
	free(w);
	free(output);
	free(duality_gap);
}



int main(){

	int n = 1e6;
	int max_iter = 1e4;
	double rho = 5e-8;
	double obj_tol = 1e-6;
	bool verbose = false;

	char *input_filename = "cvxReg_input.csv";
	char *output_filename_1 = "cvxReg_output_1.csv";
	char *output_filename_2 = "cvxReg_output_2.csv";

	clock_t tic = clock();
	convex_regression(n, max_iter, rho, obj_tol, input_filename, output_filename_1, output_filename_2, verbose);
	clock_t toc = clock();
	printf("\nElapsed: %f seconds", (double)(toc - tic) / CLOCKS_PER_SEC);

	return 0;
}




