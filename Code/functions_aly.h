#ifndef FUNCTIONS_ALY
#define FUNCTIONS_ALY

#include "phg.h"
#include "unit.h"
#include "options.h"
#include "coefficient.h"
#include "IPFEMPB.h"

extern void ls_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
extern void ls_grad_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);

extern void func_jD(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
extern void func_jN(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);
extern void func_g1D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
extern void func_g2D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
extern void func_f1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
extern void func_f2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);

extern void func_u1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
extern void func_u2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);


#endif