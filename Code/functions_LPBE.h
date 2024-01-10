#ifndef FUNCTIONS_LPBE
#define FUNCTIONS_LPBE

#include "phg.h"
#include "unit.h"
#include "options.h"
#include "coefficient.h"
#include "functions_aly.h"
#include "functions_general.h"
#include "functions_LPBE.h"

extern void DoFace_LPBE(SOLVER *solver, XFEM_INFO *xi,
	                    QCACHE *qc,  ELEMENT *e,  int face,  int dof_id,  FLOAT coe,
	                    QCACHE *qc1, ELEMENT *e1, int face1, int dof_id1, FLOAT coe1);
extern void Solve_LPBE(XFEM_INFO *xi, DOF *u1_h, DOF *u2_h, INT PrintType);


#endif