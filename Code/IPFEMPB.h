#ifndef DEMO_H
#define DEMO_H

#include "phg.h"
#include "unit.h"
#include "options.h"
#include "coefficient.h"
#include "functions_aly.h"
#include "functions_general.h"
#include "functions_LPBE.h"
#include "functions_NPBE.h"

extern GRID *g;
extern DOF *u1_h;
extern DOF *u2_h;
extern DOF *ls;
extern DOF *ls_grad;
extern XFEM_INFO *xi;

#endif