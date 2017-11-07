/**
 * @file csvIO.h
 * @author Ang Li
 * @date 2017-11-6
 * @brief functions for reading input from csv and writing output into csv file
 */

#ifndef CSVIO_H
#define CSVIO_H

/* Create a csv, and write a 1d array into it */
void create_csv(char *filename, double arr[], int n);

/* Read 3-column data table (response, x, weight) from a csv file, and write into three 1-d arrays */
void read_csv(int n, double y[], double x[], double w[], char *filename);

#endif
