
/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define N 64
#define M 1.0
#define W 1.0
#define Delta 1.5
#define N_SWEEP 2500
#define t_M_max 30
#define t_fisico_max 10
#define D_BIN 25
#define N_BIN 100

#if defined MAIN_PROGRAM
  #define EXTERN
#else
  #define EXTERN extern
#endif

EXTERN double xx[N];

#undef EXTERN

#endif

