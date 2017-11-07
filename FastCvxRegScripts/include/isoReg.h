/**
 * @file isoReg.h
 * @author Ang Li
 * @date 2017-11-6
 * @brief isotonic regression function, based on the pool adjacent violators algorithm (PAVA)
 */

#ifndef ISOREG_H
#define ISOREG_H

/* 1-d isotonic regression (linear order) */
void iso_reg(int n, double *y, double *output);

#endif
