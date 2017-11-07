/**
 * @file CvxReg1d.h
 * @author Ang Li
 * @date 2017-11-6
 * @brief Convex regression with the ADMM algorithm.
 * The main function reads simulated data from csv, run convex regression
 * and write the fitted curve, duality_gap series into csv files
 */

#ifndef CVXREG1D_H
#define CVXREG1D_H

#include <stdbool.h>

/* 1-d convex regression */
void cr_1d_admm(int n, double *y, double *x, double *w, double *output, double *duality_gap, int max_iter, double rho, double obj_tol);

/* retrieve data from input file, run convex regression, and write results into output files */
void convex_regression(int n, int max_iter, double rho, double obj_tol, char *input_filename, char *output_filename_1, char *output_filename_2, bool verbose);
#endif
