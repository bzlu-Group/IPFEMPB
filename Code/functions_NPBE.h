#ifndef FUNCTIONS_NPBE
#define FUNCTIONS_NPBE

#include "phg.h"
#include "unit.h"
#include "coefficient.h"
#include "options.h"
#include "functions_LPBE.h"
#include "functions_aly.h"

void Solve_NPBE(XFEM_INFO *xi, DOF *u1_h, DOF *u2_h, DOF *u1h_old, DOF *u2h_old, INT printtype);
#endif