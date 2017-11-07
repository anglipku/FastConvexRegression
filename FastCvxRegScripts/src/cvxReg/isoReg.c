/**
 * @file isoReg.c
 * @author Ang Li
 * @date 2017-11-6
 * @brief isotonic regression function, based on the pool adjacent violators algorithm (PAVA)
 */

#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "isoReg.h"
#include "csvIO.h"


/**
* @brief         Solving: min{||x-y||^2_2: x_1<=...<=x_n}
*                (Implementation of PAVA on linear graph)
*                (this function scales linearly)
* @param n       number of observations
* @param y       response vector
* @param output  allocated space for the output
* @return void
*/

void iso_reg(int n, double *y, double *output){
	double *yprime;
	double *w;
	int *S;
	int i, j;

	yprime = (double *) malloc((n+1) * sizeof(double));
	w = (double *) malloc((n+1) * sizeof(double));
	S = (int *) malloc((n+1) * sizeof(int));

	yprime[1] = y[0];
	w[1] = 1;
	j = 1;
	S[0] = 0;
	S[1] = 1;
	for(i = 2; i <= n; ++i){
		j += 1;
		yprime[j] = y[i-1];
		w[j] = 1;
		while((j > 1) && (yprime[j] < yprime[j-1])){
			yprime[j-1] = (w[j] * yprime[j] + w[j-1] * yprime[j-1]) / (w[j] + w[j-1]);
			w[j-1] = w[j] + w[j-1];
			j -= 1;
		}
		S[j] = i;
	}

	for(int k = 1; k <= j; ++k){
		for(int l = S[k-1] + 1; l <= S[k]; ++l){
			output[l-1] = yprime[k];
		}
	}

	free(yprime);
	free(w);
	free(S);
}


/* Test isoReg()
int main(){
	int n = 1e6;
	double scale = 0.5;

	double *input = (double *) malloc(n * sizeof(double));
	double *output = (double *) malloc(n * sizeof(double));

	//char *input_filename = "isoReg_input.csv";
	//isotonic_regression(n, input_filename, output);

	//Simulate U(0, 1)
	for(int i = 0; i < n; ++i){
		input[i] = i*1.0/n + scale * (-1 + ((double) rand () / RAND_MAX) * 2);
	}
	printf("input:");
	for(int i = 0; i < n; ++i) printf("%f ", input[i]);

	clock_t tic = clock();

	iso_reg(n, input, output);
	printf("\noutput:");
	for(int i = 0; i < n; ++i) printf("%f ", output[i]);

	clock_t toc = clock();
	printf("\nElapsed: %f seconds", (double)(toc - tic) / CLOCKS_PER_SEC);

	free(input);
	free(output);

	return 0;
}
*/


/***********************************************************/
/* The functions below are deprecated */
/***********************************************************/

/**
* @brief         Solving: min{||x-y||^2_2: x_1<=...<=x_n}
*                (Implementation of pool adjacent violators algorithm on linear graph)
* @param n       number of observations
* @param y       response vector
* @param output  allocated space for the output
* @return void
*/

/*
void iso_reg_deprecated(int n, double *y, double *output)
{
    int i;
    int j;
    int block;
    int *groups;
    int max_grp;
    double grp_avg;
    double thr = 1e-16;
    bool stop;
    bool stop_inner;

    block = 0;
    groups = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; ++i) groups[i] = i;
    stop = false;

    while(!stop){
        if(any(n, groups, block + 1)){
            if(mean_group(n, y, groups, block) <= mean_group(n, y, groups, block + 1) + thr){
                block += 1;
            } else {
                for(i = 0; i < n; ++i){
                    if(groups[i] > block) groups[i] -= 1;
                }
                stop_inner = false;
                while(!stop_inner){
                    if(any(n, groups, block - 1)){
                        if(mean_group(n, y, groups, block-1) > mean_group(n, y, groups, block) + thr){
                            for(i = 0; i < n; ++i){
                                if(groups[i] >= block) groups[i] -= 1;
                            }
                            block -= 1;
                        } else {
                            stop_inner = true;
                        }
                    } else {
                        stop_inner = true;
                    }
                }
            }
        } else {
            stop = true;
        }
    }

    max_grp = max_array(n, groups);

    for(i = 0; i <= max_grp; ++i){
        grp_avg = mean_group(n, y, groups, i);
        for(j = 0; j < n; ++j){
            if(groups[j] == i){
                output[j] = grp_avg;
            }
        }
    }

    free(groups);
}
*/


/**
* @brief            check if an int element exist in an int array
*
* @param n          number of elements in an array
* @param element
* @param array
* @return bool
*/

/*
bool any(int n, int *array, int element){
    for(int i = 0; i < n; ++i){
        if(array[i] == element) return true;
    }
    return false;
}
*/


/**
* @brief            find maximum elements in an integer array
*
* @param n          number of elements in array
* @param array
* @return int
*/

/*
int max_array(int n, int *array){
    int max = INT_MIN;
    for(int i = 0; i < n; ++i){
        if(array[i] > max) max = array[i];
    }
    return max;
}
*/


/**
* @brief            take the mean of a particular group
*
* @param n          number of elements
* @param y          response vector
* @param groups     group vector
* @param block      block number
* @return double
*/

/*
double mean_group(int n, double *y, int *groups, int block)
{
    double grp_avg;
    int grp_cnt = 0;
    double grp_sum = 0;

    for(int i = 0; i < n; ++i){
        if(groups[i] == block){
            grp_sum += y[i];
            grp_cnt += 1;
        }
    }
    grp_avg = grp_sum / grp_cnt;
    return grp_avg;
}
*/


/*
void isotonic_regression(int n, char* input_filename, double* output){
	double *y = (double *) malloc(n * sizeof(double));
	double *x = (double *) malloc(n * sizeof(double));
	double *w = (double *) malloc(n * sizeof(double));
	read_csv(n, y, x, w, input_filename);
	iso_reg(n, y, output);
	free(y);
	free(x);
	free(w);
}
*/



