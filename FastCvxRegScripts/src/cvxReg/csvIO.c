/**
 * @file csvIO.c
 * @author Ang Li
 * @date 2017-11-6
 * @brief functions for reading input from csv and writing output into csv file
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "csvIO.h"

/**
* @brief              Create a csv, and write a 1d array into it
* @param filename     csv filename
* @param arr          1d array
* @param n            length of the 1d array
*
* @return void
*/
void create_csv(char *filename, double arr[], int n){
	printf("\nCreating %s file",filename);
	FILE *fp;

	fp=fopen(filename,"w+");
	fprintf(fp,"Output\n");

	for(int i = 0;i < n; ++i){
		if(i == n-1) fprintf(fp,"%f", arr[i]);
		else{
			fprintf(fp,"%f\n", arr[i]);
		}
	}

	fclose(fp);
	printf("\n%s file created",filename);
	return;
}


/**
* @brief              Read 3-column data table (response, x, weight) from a csv file, and write into three 1-d arrays
* @param n            number of data (rows)
* @param y, x, w      response, x, weight
* @param filename     csv filename
* @return void
*/
void read_csv(int n, double y[], double x[], double w[], char *filename){
	FILE *fp = fopen(filename, "r");
	char *record;
	char *line;
	char buffer[150];
	int i = 0;

	if(fp == NULL){
		printf("\n %s file opening failed", filename);
	    return;
	}

	while((line = fgets(buffer, sizeof(buffer), fp)) != NULL && (i < n + 1)){
		if(i > 0){
			record = strtok(line, ",");
			y[i-1] = atof(record);

			record = strtok(NULL, ",");
			x[i-1] = atof(record);

			record = strtok(NULL, ",");
			w[i-1] = atof(record);
		}
		++i;
	}
}

