#ifndef FUNCTIONS_GENERAL
#define FUNCTIONS_GENERAL

#include "phg.h"
#include "unit.h"
#include "options.h"
#include "coefficient.h"
#include "functions_aly.h"
#define INTERFACE BDRY_USER0

extern void Read_pqr(char *fn_pqr);
extern void Get_h(GRID *g, FLOAT *h_max, FLOAT *h_min);
extern void Cal_energy(DOF *u1_h);
extern void cal_error(XFEM_INFO *xi, DOF *u1_h, DOF *u2_h, DOF *u1_old, DOF *u2_old,
                    FLOAT *L2_norm, FLOAT *H1_norm, FLOAT *L2_error, FLOAT *H1_error);
extern char* Creat_Domin(FLOAT *protein_len, FLOAT protein_rmax);

#endif