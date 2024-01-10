#include <math.h>
#include <string.h>

#include "functions_NPBE.h"
static DOF *u1n = NULL, *u2n = NULL, *u1n_grad = NULL, *u2n_grad = NULL;
static INT Face_tmp = 0;
static int Q_f = 0, Q_gD = 0, Q_jD = 0, Q_jN = 0;
static int Q_gcosh = 0, Q_gsinh = 0, Q_gradun = 0, Q_un = 0;
static int Q_un_jump = 0, Q_un_grad_jump = 0, Q_un_avg = 0, Q_un_grad_avg = 0;
static DOF *cosh_un = NULL, *sinh_un = NULL, *u1n_jump = NULL, *u2n_jump = NULL, *u1n_grad_jump = NULL, *u2n_grad_jump = NULL;
static DOF *u1n_avg = NULL, *u2n_avg = NULL, *u1n_grad_avg = NULL, *u2n_grad_avg = NULL;
/*************************************************jump & average functions****************************************************************************************/
static void lambdafunc_sinh_du(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
	
	assert(u2n != NULL);
	FLOAT tmp_u2h;
	phgDofEval(u2n, e, lambda, &tmp_u2h);
	if(fabs(tmp_u2h) >= 1000){
		FLOAT x, y, z;
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgError(-1, "element: %d, loc = (%lf, %lf, %lf), tmp_u2h = %e\n",e->index, (double)x, (double)y, (double)z,(double)tmp_u2h);
		tmp_u2h = 100;
	}
	*values = 0.5 * (Exp(tmp_u2h) - Exp(-tmp_u2h));
    return;
}
static void lambdafunc_cosh_du(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
	
	assert(u2n != NULL);
	FLOAT tmp_u2h;
	phgDofEval(u2n, e, lambda, &tmp_u2h);
	//phgPrintf("%lf\n",(double)tmp_u2h);
	if(fabs(tmp_u2h) >= 1000){
		phgError(-1, "tmp_u2h = %e\n", (double)tmp_u2h);
		tmp_u2h = 100;
	}
	*values = 0.5 * (Exp(tmp_u2h) + Exp(-tmp_u2h));
    return;
}
static void lambdafunc_u1n_jump(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u1n, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u1n)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u1n, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u1n); i++){
			values[i] -= values1[i];
		}
		if (e->generation < e1->generation){
		    for(i = 0; i < DofDim(u1n); i++){
			values[i] = -values[i];
		    }
		}
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index)){
		    for(i = 0; i < DofDim(u1n); i++){
			values[i] = -values[i];
		    }
		}	
        return;
	}
}
static void lambdafunc_u2n_jump(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u2n, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u2n)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u2n, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u2n); i++){
			values[i] -= values1[i];
		}
		if (e->generation < e1->generation){
		    for(i = 0; i < DofDim(u2n); i++){
			values[i] = -values[i];
		    }
		}
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index)){
		    for(i = 0; i < DofDim(u2n); i++){
			values[i] = -values[i];
		    }
		}	
	}
}
static void lambdafunc_u1n_grad_jump(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u1n_grad, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u1n_grad)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u1n_grad, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u1n_grad); i++){
			values[i] -= values1[i];
		}
		if (e->generation < e1->generation){
		    for(i = 0; i < DofDim(u1n_grad); i++){
			values[i] = -values[i];
		    }
		}
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index)){
		    for(i = 0; i < DofDim(u1n_grad); i++){
			values[i] = -values[i];
		    }
		}	
	}
}
static void lambdafunc_u2n_grad_jump(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u2n_grad, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u2n_grad)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u2n_grad, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u2n_grad); i++){
			values[i] -= values1[i];
			//values[i] *= Epsilon_s;
		}
		if (e->generation < e1->generation){
		    for(i = 0; i < DofDim(u2n_grad); i++){
			values[i] = -values[i];
			//values[i] *= Epsilon_s;
		    }
		}
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index)){
		    for(i = 0; i < DofDim(u2n_grad); i++){
			values[i] = -values[i];
			//values[i] *= Epsilon_s;
		    }
		}	 
	}
}
static void lambdafunc_u1n_avg(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u1n, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u1n)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u1n, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u1n); i++){
			values[i] += values1[i];
			values[i] /= 2.0;
		}
	}
}
static void lambdafunc_u2n_avg(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u2n, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u2n)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u2n, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u2n); i++){
			values[i] += values1[i];
			values[i] /= 2.0;
		}
	}
}
static void lambdafunc_u1n_grad_avg(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u1n_grad, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u1n_grad)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u1n_grad, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u1n_grad); i++){
			values[i] += values1[i];
			values[i] /= 2.0;
		}
	}
}
static void lambdafunc_u2n_grad_avg(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    ELEMENT *e1 = phgGetNeighbour(g, e, Face_tmp);
	phgDofEval(u2n_grad, e, lambda, values);
	if(e1 != NULL){
		FLOAT values1[DofDim(u2n_grad)], x, y, z, lambda1[Dim + 1];
		phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		phgGeomXYZ2Lambda(g, e1, x, y, z, lambda1);
		phgDofEval(u2n_grad, e1, lambda1, values1);
		int i;
		for(i = 0; i < DofDim(u2n_grad); i++){
			values[i] += values1[i];
			values[i] /= 2.0;
		}
	}
}
static void Creat_lambdaFEM(){
    cosh_un = phgDofNew(g, DOF_ANALYTIC, 1, "cosh_du", DofLambdaFunction);
	phgDofSetLambdaFunction(cosh_un, lambdafunc_cosh_du);
	sinh_un = phgDofNew(g, DOF_ANALYTIC, 1, "sinh_du", DofLambdaFunction);
	phgDofSetLambdaFunction(sinh_un, lambdafunc_sinh_du);

	u1n_jump = phgDofNew(g, DOF_ANALYTIC, DofDim(u1n), "u1n_jump", DofLambdaFunction);
	phgDofSetLambdaFunction(u1n_jump, lambdafunc_u1n_jump);
	u2n_jump = phgDofNew(g, DOF_ANALYTIC, DofDim(u2n), "u2n_jump", DofLambdaFunction);
	phgDofSetLambdaFunction(u2n_jump, lambdafunc_u2n_jump);
	u1n_grad_jump = phgDofNew(g, DOF_ANALYTIC, DofDim(u1n_grad), "u1n_grad_jump", DofLambdaFunction);
	phgDofSetLambdaFunction(u1n_grad_jump, lambdafunc_u1n_grad_jump);
	u2n_grad_jump = phgDofNew(g, DOF_ANALYTIC, DofDim(u2n_grad), "u2n_grad_jump", DofLambdaFunction);
	phgDofSetLambdaFunction(u2n_grad_jump, lambdafunc_u2n_grad_jump);

	u1n_avg = phgDofNew(g, DOF_ANALYTIC, DofDim(u1n), "u1n_avg", DofLambdaFunction);
	phgDofSetLambdaFunction(u1n_avg, lambdafunc_u1n_avg);
	u2n_avg = phgDofNew(g, DOF_ANALYTIC, DofDim(u2n), "u2n_avg", DofLambdaFunction);
	phgDofSetLambdaFunction(u2n_avg, lambdafunc_u2n_avg);
	u1n_grad_avg = phgDofNew(g, DOF_ANALYTIC, DofDim(u1n_grad), "u1n_grad_avg", DofLambdaFunction);
	phgDofSetLambdaFunction(u1n_grad_avg, lambdafunc_u1n_grad_avg);
	u2n_grad_avg = phgDofNew(g, DOF_ANALYTIC, DofDim(u2n_grad), "u2n_grad_avg", DofLambdaFunction);
	phgDofSetLambdaFunction(u2n_grad_avg, lambdafunc_u2n_grad_avg);
    return;
}
static void Free_lambdaFEM(){
    phgDofFree(&cosh_un);
	phgDofFree(&sinh_un);

	phgDofFree(&u1n_jump);
	phgDofFree(&u2n_jump);
	phgDofFree(&u1n_grad_jump);
	phgDofFree(&u2n_grad_jump);

	phgDofFree(&u1n_avg);
	phgDofFree(&u2n_avg);
	phgDofFree(&u1n_grad_avg);
	phgDofFree(&u2n_grad_avg);
    return;
}
/**********************************************************NPBE functions**************************************************************************************/
static void 
DoFace_NPBE(SOLVER *solver, XFEM_INFO *xi,
	        QCACHE *qc,  ELEMENT *e,  int face,  int dof_id,  FLOAT coe,
	        QCACHE *qc1, ELEMENT *e1, int face1, int dof_id1, FLOAT coe1, int Case)
/* This function computes integrations on a surface area, which can be:
 *
 *   1) an element face, in this case:
 * 	e != e1 && dof_id == dof_id1, e and e1 are both non interface elements
 *
 *   2) part of an element face, in this case:
 *   	e != e1 && dof_id == dof_id1, e and e1 are both interface elements
 *
 *   3) \Gamma \cap e, in this case:
 *	e == e1, dof_id = 0, dof_id1 = 1, and e is an interface element
 *
 *   4) an element face which is entirely contained in \Gamma, in this case:
 *	e1 != NULL, e1 != e, mark(e) < 0, mark(e1) > 0, dof_id = 0, dof_id1 = 1
 *
 * For cases 1) and 2), the normal vector stored in qc is the outward normal
 * vector of e, and e1 can be NULL (boundary face).
 *
 * For the cases 3) and 4) the normal vectors stored in qc points from
 * \Omega^- to \Omega^+.
 */
  /*Case: 2: Newton; 3: Residual*/
{   
    assert(Case == 2);
    DOF *u_h = qc->fe, *u1_h = qc1->fe;
    GRID *g = u_h->g;
    int n = 0, n1 = 0, p = 0;   /* p = polynomial order */
    INT I = 0, J = 0;
    int i = 0, j = 0;
    FLOAT val = 0.0;
    FLOAT G0 = 0.0, G1 = 0.0, h = 0.0, a = 0.0;     /* G0 = gamma0*p^2/h, G1 = gamma1*h/p^2 */
    BTYPE bdry_flag = 0;

    if (phgQCGetNP(qc) == 0)
	return;

    n = qc->qd->n_bas(qc, e->index);
    p = DofTypeOrder(u_h, e);
    h = e == e1 ? phgGeomGetDiameter(g, e) :
		  phgGeomGetFaceDiameter(g, e, face);
    if (e1 == NULL) {
	/* boundary face */
    bdry_flag = DIRICHLET;
	n1 = 0;
	G0 = coe * Gamma0_XFEM * p * p / h;
	G1 = coe * Gamma1_XFEM * h / (p * (FLOAT)p);
	/* RHS */
	for (i = 0; i < n; i++) {
	    I = phgSolverMapE2G(solver, dof_id, e, i);
        val = 0.0;
	    if (bdry_flag == DIRICHLET) {	/* Dirichlet boundary */
		/* -\beta\int_\Gamma_D g_D (A\grad v).n */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD,   PROJ_NONE, 0,
				qc, e->index, face, Q_GRAD, PROJ_DOT,  i)
				* coe;
		val = -Theta_XFEM * a;
		/* G0 \int_\Gamma_D g_D v */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD, PROJ_NONE, 0,
				qc, e->index, face, Q_BAS, PROJ_NONE, i);
		val += G0 * a;
	    }
	    phgSolverAddGlobalRHSEntry(solver, I, val);
	}
    }
    else {
	bdry_flag = INTERIOR;
    if(dof_id1 != dof_id){
		bdry_flag = (bdry_flag | INTERFACE);
	}
	n1 = qc1->qd->n_bas(qc1, e1->index);
	i = DofTypeOrder(u1_h, e1);
	if (p < i)
	    p = i;
	G0 = ((coe + coe1) / 2.0) * Gamma0_XFEM * p * (FLOAT)p / h;
	G1 = ((coe + coe1) / 2.0) * Gamma1_XFEM * h / (p * (FLOAT)p);
    }
    /* The following macros depend on the assertion */

#define Sel(i, o, o1)	(i < n ? o : o1)
#define Qc(i)	Sel(i, qc, qc1)
#define Ele(i)	Sel(i, e, e1)
#define Fac(i)	Sel(i, face, face1)
#define Bas(i)	Sel(i, i, i - n)
#define Dof(i)	Sel(i, dof_id, dof_id1)
#define Coe(i)	Sel(i, coe, coe1)
#define Eid(i)	Ele(i)->index
#define Quad(fid1, proj1, i1, fid2, proj2, i2) phgQCIntegrateFace( \
		Qc(i1), Eid(i1), Fac(i1), fid1, proj1, Bas(i1), \
		Qc(i2), Eid(i2), Fac(i2), fid2, proj2, Bas(i2))

    /* loop on {basis funcs in e} \cup {basis funcs in e1} */

for (i = 0; i < n + n1; i++) {
	I = phgSolverMapE2G(solver, Dof(i), Ele(i), Bas(i));
	/* loop on {basis funcs in e} \cup {basis funcs in e1} */
	for (j = 0; j < n + n1; j++) {
		
	    J = phgSolverMapE2G(solver, Dof(j), Ele(j), Bas(j));
	    /* skip jumps for interior face and continuous element */
	    if (DofFESpace(u_h) == FE_H1 && e != e1 && e1 != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark >= 0) {
		continue;
	    }

	    /*-----------------------------------------------------------------
	     * Note if the normal vector is reversed, the sign of [.] changes,
	     * while the sign of {.} is unaffected.
	     *----------------------------------------------------------------*/

	    val = 0.0;

	    if (bdry_flag != NEUMANN) {

		/* -\int {A\grad u}.n [v] (i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_GRAD,  PROJ_DOT, j, Q_BAS, PROJ_NONE, i) * Coe(j);
		if (bdry_flag & INTERIOR)
		    a *= (i < n ? 0.5 : -0.5);
		val = -a;

		/* -\int \beta [u]{A\grad v}.n (note: i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_BAS, PROJ_NONE, j, Q_GRAD,  PROJ_DOT, i) * Coe(i);
		if (bdry_flag & INTERIOR)
		    a *= (j < n ? 0.5 : -0.5);
		val += -Theta_XFEM * a;

		/* \int G0 [u][v] (i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_BAS, PROJ_NONE, j, Q_BAS, PROJ_NONE, i);
		if ((bdry_flag  & INTERIOR) && (i < n) != (j < n))
		    a = -a;
		val += G0 * a;
	    }

	    if (bdry_flag != DIRICHLET) {
		/* \int G1 [A\grad u].n [A\grad v].n */
		a = Quad(Q_GRAD,  PROJ_DOT, j, Q_GRAD,  PROJ_DOT, i)
			* Coe(j) * Coe(i);
		if ((bdry_flag  & INTERIOR) && (i < n) != (j < n))
		    a = -a;
		val += G1 * a;
	    }
	    phgSolverAddGlobalMatrixEntry(solver, I, J, val);
		
	}

    if(Case == 2){
	/*a(un,v)*/
	val = 0.0;
    if(bdry_flag & INTERFACE){
		FLOAT a1, a2;
		a1 = phgQCIntegrateFace(qc, e->index, face, Q_gradun,  PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_gradun,  PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a = (Epsilon_p * a1 + Epsilon_s * a2) * (i < n ? 0.5 : -0.5);
		val = -a;

		a1 = phgQCIntegrateFace(qc, e->index, face, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a = a1 - a2;
		val += -Theta_XFEM * 0.5 * Coe(i) * a;

		a1 = phgQCIntegrateFace(qc, e->index, face, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a = a1- a2;
	    val += G0 * (i < n ? a : -a);

		a1 = phgQCIntegrateFace(qc, e->index, face, Q_gradun, PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_gradun, PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a = Epsilon_p * a1 - Epsilon_s * a2;
	    val += Coe(i) * G1 * (i < n ? a : -a);
	}
	else{
		if (DofFESpace(u_h) == FE_H1 && e != e1 && e1 != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark >= 0) {
		assert(bdry_flag == INTERIOR && dof_id == dof_id1);
		continue;
	    }
		assert(dof_id1 == dof_id);
		FLOAT coe_tmp;
		coe_tmp = dof_id == 0 ? Epsilon_p : Epsilon_s;
		Face_tmp = face;
		assert(Face_tmp != -1);
	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_grad_avg,  PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i)) * coe_tmp;
	    val = i < n ? -a : a;

	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_jump,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		if(bdry_flag & INTERIOR){
			a = a * 0.5;
		}
	    val += -Theta_XFEM * Coe(i) * a;

	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_jump,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	    val += G0 * (i < n ? a : -a);

		if(bdry_flag & INTERIOR){
	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_grad_jump, PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i)) * coe_tmp;
	    val += Coe(i) * G1 * (i < n ? a : -a);
		}
	}
	phgSolverAddGlobalRHSEntry(solver, I, -val);
	
	}

	if (e1 == NULL || dof_id == dof_id1)
	    continue;		/* non interface face */

	/* this face is part of the interface, apply jump conditions */
	assert(dof_id == 0 && dof_id1 == 1);
	assert((e == e1 && xi->info[e->index].mark == 0) ||
	       (xi->info[e->index].mark < 0 && xi->info[e1->index].mark > 0));
    val = 0.0;

	FLOAT tmp_val1, tmp_val2, tmp_val3, tmp_val4;
	/* \int [A*gN.n] {v} */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,  PROJ_DOT,  0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	tmp_val1 = a * 0.5;
	val = tmp_val1;

	/* G1 \int [A*gN.n] [(A\grad v).n] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,   PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	tmp_val4 = Coe(i) * G1 * (i < n ? a : -a);
	val += tmp_val4;

	/* -beta \int [gD] {(A\grad v).n} (func[1] := gD) */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,   PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	tmp_val2 = -Theta_XFEM * 0.5 * Coe(i) * a;
	val += tmp_val2;

	/* G0 \int [gD] [v] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	tmp_val3 = G0 * (i < n ? a : -a);
	val += tmp_val3;

    phgSolverAddGlobalRHSEntry(solver, I, val);
}

#undef Sel
#undef Qc
#undef Ele
#undef Fac
#undef Bas
#undef Dof
#undef Coe
#undef Quad
}

static void Build_linear_system_NPBE(XFEM_INFO *xi, SOLVER *solver, DOF *u1_h, DOF *u2_h, int Case){
/*Case: 1: Picard; 2: Newton; 3: Residual*/
    assert(Case >= 1 && Case <= 2);
    GRID *g = xi->ls->g;
    QCACHE *qm, *qp;
    ELEMENT *e;
	Creat_lambdaFEM();

    qm = phgXFEMNewQC(xi, u1_h, 0);
    Q_f = phgQCAddXYZFunction(qm, func_f1, 1);			/* f */
	Q_gD = phgQCAddXYZFunction(qm, func_g1D, 1);	
    Q_jD = phgQCAddXYZFunction(qm, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qm, func_jN, 3);		/* jN */
	Q_un = phgQCAddFEFunction(qm, u1n);
	Q_gradun = phgQCAddFEFunction(qm, u1n_grad);
	Q_un_jump = phgQCAddFEFunction(qm, u1n_jump);
	Q_un_grad_jump = phgQCAddFEFunction(qm, u1n_grad_jump);
    Q_un_avg = phgQCAddFEFunction(qm, u1n_avg);
	Q_un_grad_avg = phgQCAddFEFunction(qm, u1n_grad_avg);

    qp = phgXFEMNewQC(xi, u2_h, 1);
    /* The following must be in exactly the same order as above */
    Q_f = phgQCAddXYZFunction(qp, func_f2, 1);			/* f  */
    Q_gD = phgQCAddXYZFunction(qp, func_g2D, 1);	
    Q_jD = phgQCAddXYZFunction(qp, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qp, func_jN, 3);		/* jN */
	Q_un = phgQCAddFEFunction(qp, u2n);
	Q_gradun = phgQCAddFEFunction(qp, u2n_grad);
	Q_un_jump = phgQCAddFEFunction(qp, u2n_jump);
	Q_un_grad_jump = phgQCAddFEFunction(qp, u2n_grad_jump);
	Q_un_avg = phgQCAddFEFunction(qp, u2n_avg);
	Q_un_grad_avg = phgQCAddFEFunction(qp, u2n_grad_avg);
	Q_gcosh = phgQCAddFECoefficient(qp, cosh_un, Q_BAS);
	Q_gsinh = phgQCAddFEFunction(qp, sinh_un);
	
#if USE_OMP
#pragma omp parallel for private(e, Face_tmp)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {
	FLOAT *rule_m = NULL, *rule_0 = NULL, *rule_p = NULL;
	int N = qm->qd->n_bas(qm, e->index);
	INT I = 0, J = 0, eno = e->index;
	int i = 0, j = 0, k = 0, face = 0, order = 0;
	FLOAT val = 0.0;

	/* construct quadrature rules for the subdomains and the interface */
	i = DofTypeOrder(u1_h, e);
	j = DofTypeOrder(u2_h, e);
	order = 2 * (i > j ? i : j);
	if (xi->info[e->index].mark == 0)
	    phgQuadInterface(xi->ls, xi->ls_grad, e, order + Extra_orders,
				&rule_m, &rule_0, &rule_p);
	else if (xi->info[e->index].mark < 0)
	    phgQuadInterface(NULL, xi->ls_grad, e, order /*- 2*/,
			     &rule_m, NULL, NULL);
	else if (xi->info[e->index].mark > 0)
	    phgQuadInterface(NULL, xi->ls_grad, e, order/* - 2*/,
			     NULL, NULL, &rule_p);
    
	/*==================== volume integrals ====================*/
	phgQCSetRule(qm, rule_m, -1.);
	phgQCSetRule(qp, rule_p, -1.);
    
	for (k = 0; k < 2; k++) {
	    QCACHE *q;
	    if ((k == 0 && xi->info[e->index].mark > 0) || (k == 1 && xi->info[e->index].mark < 0)){
            continue;
		} 
	    q = k == 0 ? qm : qp;
	    for (i = 0; i < N; i++) {
		I = phgSolverMapE2G(solver, k, e, i);
		for (j = 0; j < N; j++) {
		    /* \int_T A\grad u1 . \grad v */
		    J = phgSolverMapE2G(solver, k, e, j);
		    val = phgQCIntegrate(q, eno, Q_GRAD, j, q, eno, Q_GRAD, i) * (k == 0 ? Epsilon_p : Epsilon_s);
            if(Case == 2){
                if(k == 1){
                val += phgQCIntegrate(q, eno, Q_gcosh, j, q, eno, Q_BAS, i) * a_noline;
			    }
            }
		    phgSolverAddGlobalMatrixEntry(solver, I, J, val);
		}
		
		/* \int_T f1 v */
		val = phgQCIntegrate(q, eno, Q_BAS, i, q, eno, Q_f, 0);
        if(Case == 1){
            if(k == 1){
			    val -= phgQCIntegrate(q, eno, Q_gsinh, 0, q, eno, Q_BAS, i) * a_noline;
		    }
        }
        else if(Case == 2){
            val -= phgQCIntegrate(q, eno, Q_gradun, 0, q, eno, Q_GRAD, i) * (k == 0 ? Epsilon_p : Epsilon_s);
		    if(k == 1){
			    val -= phgQCIntegrate(q, eno, Q_gsinh, 0, q, eno, Q_BAS, i) * a_noline;
		    }
        }
		phgSolverAddGlobalRHSEntry(solver, I, val);
	    }	
	}
	phgFree(rule_m);
	phgFree(rule_p);
	rule_m = rule_p = NULL;
    
	/*==================== face integrals ====================*/
	if (xi->info[e->index].mark == 0) {	
        /* the interface */
	    phgQCSetRule(qm, rule_0, -1.);
	    phgQCSetRule(qp, rule_0, -1.);
        if(Case == 1){
            DoFace_LPBE(solver, xi, qm, e, -1, /*dof_id*/0,  Epsilon_p, qp, e, -1, /*dof_id1*/1, Epsilon_s);
        }
        else{
	        DoFace_NPBE(solver, xi, qm, e, -1, /*dof_id*/0,  Epsilon_p, qp, e, -1, /*dof_id1*/1, Epsilon_s, Case);
        }
	}
	phgFree(rule_0);
	rule_0 = NULL;  
	/* the faces of the element */
	for (face = 0; face < NFace; face++){
	    FLOAT nv[Dim];
	    ELEMENT *e1 = phgGetNeighbour(g, e, face);
	    int face1 = -1;
	    if (e1 != NULL) {
		/* a face is only processed by the smaller of the two
		 * neighbouring elements, and by the one with smaller
		 * global index if the two elements are of the same size,
		 * to avoid double counting or redundant computation */
		if (e->generation < e1->generation)
		    continue;
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index))
		    continue; /* process each interior face just once */
		face1 = phgOppositeFace(g, e, face, e1);
	    }

	    phgGeomGetFaceOutNormal(g, e, face, nv);

	    if (xi->info[e->index].mark == 0)
		/* the face might be cut into 2 parts by the interface */
		phgQuadInterfaceFace(xi->ls, xi->ls_grad, e, face, order + Extra_orders, &rule_m, NULL, &rule_p);
	    else if (xi->info[e->index].mark < 0)
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face, order, &rule_m, NULL, NULL);
	    else
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face, order, NULL, NULL, &rule_p);
	    /* Note: e1 is either NULL or valid (since gen(e) >= gen(e1)) */
	    if (e1 != NULL && xi->info[e->index].mark * xi->info[e1->index].mark < 0){
		/* The face is flat and is part of the interface
		 * (e and e1 in different subdomains) */
		if (xi->info[e->index].mark < 0){
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    phgQCSetRule(qp, rule_m, -1.);
		    phgQCSetConstantNormal(qp, nv);
            if(Case == 1){
                DoFace_LPBE(solver, xi, qm, e,  face, 0, Epsilon_p, qp, e1, face1, 1, Epsilon_s);
            }
            else{
	            DoFace_NPBE(solver, xi, qm, e,  face, 0, Epsilon_p,qp, e1, face1, 1, Epsilon_s, Case);
            }
		}
		else{
		    nv[0] = -nv[0]; nv[1] = -nv[1]; nv[2] = -nv[2];
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    phgQCSetRule(qm, rule_p, -1.);
		    phgQCSetConstantNormal(qm, nv);
            if(Case == 1){
                DoFace_LPBE(solver, xi, qm, e1, face1, 0, Epsilon_p, qp, e,  face, 1, Epsilon_s);
            }
            else{
	            DoFace_NPBE(solver, xi, qm, e1, face1, 0, Epsilon_p, qp, e,  face, 1, Epsilon_s, Case);
            }
		}
	    }
	    else {
		/* e and e1 in a same subdomain 
         * element face */
		if (xi->info[e->index].mark <= 0) {
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
            if(Case == 1){
                DoFace_LPBE(solver, xi, qm, e,  face,  0, Epsilon_p, qm, e1, face1, 0, Epsilon_p);
            }
            else{
	            DoFace_NPBE(solver, xi, qm, e,  face,  0, Epsilon_p, qm, e1, face1, 0, Epsilon_p, Case);
            }
		}
		if (xi->info[e->index].mark >= 0) {
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
            if(Case == 1){
                DoFace_LPBE(solver, xi, qp, e,  face,  1, Epsilon_s, qp, e1, face1, 1, Epsilon_s);
            }
            else{
	            DoFace_NPBE(solver, xi, qp, e,  face,  1, Epsilon_s, qp, e1, face1, 1, Epsilon_s, Case);
            }
		}
	    }
	    phgFree(rule_m);
	    phgFree(rule_p);
	    rule_m = rule_p = NULL;
	}
    }ForAllElementsEnd

	Free_lambdaFEM();
    phgQCFree(&qm);
    phgQCFree(&qp);
    /* set the diagonal entry of empty rows to 1 */
    phgXFEMProcessEmptyRows(solver);
    if (dump_solver){
        phgSolverDumpMATLAB_(solver, "A", "b", NULL, TRUE);
	}
}

/************************************************Residual functions******************************************************************************************/
static void
Cal_face(SOLVER *solver, XFEM_INFO *xi, VEC *F,
	QCACHE *qc,  ELEMENT *e,  int face,  int dof_id,  FLOAT coeff1,
	QCACHE *qc1, ELEMENT *e1, int face1, int dof_id1, FLOAT coeff2)
{
    GRID *g;
    DOF *u_h, *u1_h;
    int n = 0, n1 = 0, p = 0;   /* p = polynomial order */
    INT I;
    int i;
    FLOAT val = 0.0;
    FLOAT G0, G1, h, a;     /* G0 = gamma0*p^2/h, G1 = gamma1*h/p^2 */
    BTYPE bdry_flag;

    if (phgQCGetNP(qc) == 0)
	return;

    u_h = qc->fe;
    u1_h = qc1->fe;
    g = u_h->g;

    n = qc->qd->n_bas(qc, e->index);
    p = DofTypeOrder(u_h, e);

    h = e == e1 ? phgGeomGetDiameter(g, e) :
		  phgGeomGetFaceDiameter(g, e, face);
    if(e1 != NULL){
		n1 = qc1->qd->n_bas(qc1, e1->index);
	}
	else{
		n1 = 0;
	}
	int N = n+n1;
	FLOAT values[N];
	INT indices[N];
	bzero(values, N*sizeof(FLOAT));
    if (e1 == NULL) {
	/* boundary face */
    bdry_flag = DIRICHLET;
	n1 = 0;
	G0 = coeff1 * Gamma0_XFEM * p * p / h;
	G1 = coeff1 * Gamma1_XFEM * h / (p * (FLOAT)p);
	/* RHS */
	for (i = 0; i < n; i++) {
	    I = phgSolverMapE2G(solver, dof_id, e, i);
        indices[i] = I;
		val = 0.0;
	    if (bdry_flag == DIRICHLET) {	/* Dirichlet boundary */
		/* -\beta\int_\Gamma_D g_D (A\grad v).n */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD,   PROJ_NONE, 0,
				qc, e->index, face, Q_GRAD, PROJ_DOT,  i) * coeff1;
		val = -(-Theta_XFEM * a);
		/* G0 \int_\Gamma_D g_D v */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD, PROJ_NONE, 0,
				qc, e->index, face, Q_BAS, PROJ_NONE, i);
		val -= G0 * a;
	    }
		values[i] = val;
	}

	phgVecAddGlobalEntries(F, 0, N, indices, values);

    }
    else {
	/*Interior Face*/
	bdry_flag = INTERIOR;
	n1 = qc1->qd->n_bas(qc1, e1->index);
	i = DofTypeOrder(u1_h, e1);
	if (p < i)
	    p = i;
	G0 = ((coeff1 + coeff2) / 2) * Gamma0_XFEM * p * (FLOAT)p / h;
	G1 = ((coeff1 + coeff2) / 2) * Gamma1_XFEM * h / (p * (FLOAT)p);
    }

    /* The following macros depend on the assertion */

#define Sel(i, o, o1)	(i < n ? o : o1)
#define Qc(i)	Sel(i, qc, qc1)
#define Ele(i)	Sel(i, e, e1)
#define Fac(i)	Sel(i, face, face1)
#define Bas(i)	Sel(i, i, i - n)
#define Dof(i)	Sel(i, dof_id, dof_id1)
#define Coe(i)	Sel(i, coeff1, coeff2)
#define Eid(i)	Ele(i)->index
#define Quad(fid1, proj1, i1, fid2, proj2, i2) phgQCIntegrateFace( \
		Qc(i1), Eid(i1), Fac(i1), fid1, proj1, Bas(i1), \
		Qc(i2), Eid(i2), Fac(i2), fid2, proj2, Bas(i2))

    /* loop on {basis funcs in e} \cup {basis funcs in e1} */
int Case_face = 0;
if(e == e1){
	/*interface in element*/
	Case_face = 1;
}
else if(e1 != NULL && xi->info[e->index].mark * xi->info[e1->index].mark < 0){
	/*flat interface*/
    Case_face = 2;
}
else if(xi->info[e->index].mark < 0){
	/*\Omega- element face*/
	Case_face = 3;
}
else if(xi->info[e->index].mark > 0){
	/*\Omega+ element face*/
	Case_face = 4;
}
else if(xi->info[e->index].mark == 0){
	/*\Omega+- element face*/
	Case_face = 5;
}
else{
	phgPrintf("no that case in Do_face\n");
	exit(0);
}

    for (i = 0; i < n + n1; i++) {
	I = phgSolverMapE2G(solver, Dof(i), Ele(i), Bas(i));
	indices[i] = I;
	/* loop on {basis funcs in e} \cup {basis funcs in e1} */
    if(Case_face == 1 || Case_face == 2){
		/*interface*/
		FLOAT a1 = 0.0, a2 = 0.0;
		val = 0.0;
		a1 = phgQCIntegrateFace(qc, e->index, face, Q_gradun,  PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_gradun,  PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a = (Epsilon_p * a1 + Epsilon_s * a2) * (i < n ? 0.5 : -0.5);
		val = -a;

		a1 = phgQCIntegrateFace(qc, e->index, face, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a = a1 - a2;
		val += -Theta_XFEM * 0.5 * Coe(i) * a;

		a1 = phgQCIntegrateFace(qc, e->index, face, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_un,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
		a = a1- a2;
	    val += G0 * (i < n ? a : -a);

		a1 = phgQCIntegrateFace(qc, e->index, face, Q_gradun, PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a2 = phgQCIntegrateFace(qc1, e1->index, face1, Q_gradun, PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		a = Epsilon_p * a1 - Epsilon_s * a2;
	    val += Coe(i) * G1 * (i < n ? a : -a);
	}
	else{
		if (DofFESpace(u_h) == FE_H1 && e != e1 && e1 != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark >= 0) {
		assert(bdry_flag == INTERIOR && dof_id == dof_id1);
		continue;
	    }
		//assert(Case_face == 5);
		assert(dof_id1 == dof_id);
		FLOAT coe_tmp;
		coe_tmp = dof_id == 0 ? Epsilon_p : Epsilon_s;
		Face_tmp = face;
		assert(Face_tmp != -1);
		val = 0.0;
	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_grad_avg,  PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i)) * coe_tmp;
	    val = i < n ? -a : a;

	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_jump,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
		if(bdry_flag == INTERIOR){
			a = a * 0.5;
		}
	    val += -Theta_XFEM * Coe(i) * a;

	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_jump,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	    val += G0 * (i < n ? a : -a);

		if(bdry_flag == INTERIOR){
	    a = phgQCIntegrateFace(qc, e->index, face, Q_un_grad_jump, PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i)) * coe_tmp;
	    val += Coe(i) * G1 * (i < n ? a : -a);
		}
	}

    values[i] = val;

	if (e1 == NULL || dof_id == dof_id1)
	    continue;		/* non interface face */

	/* this face is part of the interface, apply jump conditions */
	assert(dof_id == 0 && dof_id1 == 1);
	assert((e == e1 && xi->info[e->index].mark == 0) ||
	       (xi->info[e->index].mark < 0 && xi->info[e1->index].mark > 0));
    val = 0.0;
	/* \int [A*gN.n] {v} */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,  PROJ_DOT,  0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	val = a * 0.5;

	/* G1 \int [A*gN.n] [(A\grad v).n] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,   PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	val += Coe(i) * G1 * (i < n ? a : -a);

	/* -beta \int [gD] {(A\grad v).n} (func[1] := gD) */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,   PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	val += -Theta_XFEM * 0.5 * Coe(i) * a;

	/* G0 \int [gD] [v] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	val += G0 * (i < n ? a : -a);

    //phgSolverAddGlobalRHSEntry(solver, I, val);
	values[i] += -val;
    }
	phgVecAddGlobalEntries(F, 0, N, indices, values);

#undef Sel
#undef Qc
#undef Ele
#undef Fac
#undef Bas
#undef Dof
#undef Coe
#undef Quad
}


static void 
Assemble(XFEM_INFO *xi, SOLVER *solver, DOF *u1, DOF *u2, VEC *F){

    GRID *g = xi->ls->g;
    QCACHE *qm, *qp;
    ELEMENT *e;
	DOF *f1_h, *f2_h;
	f1_h = phgDofNew(u1->g, DOF_ANALYTIC, 1, "f1_h", func_f1);
	f2_h = phgDofNew(u2->g, DOF_ANALYTIC, 1, "f2_h", func_f2);
	Creat_lambdaFEM();

    qm = phgXFEMNewQC(xi, u1_h, 0);
    Q_f = phgQCAddFEFunction(qm, f1_h);			/* f */
	Q_gD = phgQCAddXYZFunction(qm, func_g1D, 1);	
    Q_jD = phgQCAddXYZFunction(qm, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qm, func_jN, 3);		/* jN */
	Q_un = phgQCAddFEFunction(qm, u1n);
	Q_gradun = phgQCAddFEFunction(qm, u1n_grad);
	Q_un_jump = phgQCAddFEFunction(qm, u1n_jump);
	Q_un_grad_jump = phgQCAddFEFunction(qm, u1n_grad_jump);
    Q_un_avg = phgQCAddFEFunction(qm, u1n_avg);
	Q_un_grad_avg = phgQCAddFEFunction(qm, u1n_grad_avg);
        

    qp = phgXFEMNewQC(xi, u2_h, 1);
    /* The following must be in exactly the same order as above */
    Q_f = phgQCAddFEFunction(qp, f2_h);			/* f  */
    Q_gD = phgQCAddXYZFunction(qp, func_g2D, 1);	
    Q_jD = phgQCAddXYZFunction(qp, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qp, func_jN, 3);		/* jN */
	Q_un = phgQCAddFEFunction(qp, u2n);
	Q_gradun = phgQCAddFEFunction(qp, u2n_grad);
	Q_un_jump = phgQCAddFEFunction(qp, u2n_jump);
	Q_un_grad_jump = phgQCAddFEFunction(qp, u2n_grad_jump);
	Q_un_avg = phgQCAddFEFunction(qp, u2n_avg);
	Q_un_grad_avg = phgQCAddFEFunction(qp, u2n_grad_avg);
	Q_gsinh = phgQCAddFEFunction(qp, sinh_un);
	
#if USE_OMP
#pragma omp parallel for private(e, Face_tmp)
#endif	/* USE_OMP */
    
    ForAllElementsBegin(g, e) {
	phgVecDisassemble(F);
	FLOAT *rule_m = NULL, *rule_0 = NULL, *rule_p = NULL;
	int N = qm->qd->n_bas(qm, e->index);
	INT I, eno = e->index;
	int i, k, face, order;
	FLOAT val;
	FLOAT values[N];
	INT indices[N];
	bzero(values, N*sizeof(FLOAT));

	/* construct quadrature rules for the subdomains and the interface */
	i = DofTypeOrder(u1_h, e);
	k = DofTypeOrder(u2_h, e);
	order = 2 * (i > k ? i : k);
	if (xi->info[e->index].mark == 0)
	    phgQuadInterface(xi->ls, xi->ls_grad, e, order + Extra_orders,
				&rule_m, &rule_0, &rule_p);
	else if (xi->info[e->index].mark < 0)
	    phgQuadInterface(NULL, xi->ls_grad, e, order /*- 2*/,
			     &rule_m, NULL, NULL);
	else if (xi->info[e->index].mark > 0)
	    phgQuadInterface(NULL, xi->ls_grad, e, order/* - 2*/,
			     NULL, NULL, &rule_p);  
	/*==================== volume integrals ====================*/
	phgQCSetRule(qm, rule_m, -1.);
	phgQCSetRule(qp, rule_p, -1.);
	for (k = 0; k < 2; k++) {

	    QCACHE *q;
	    if ((k == 0 && xi->info[e->index].mark > 0) || (k == 1 && xi->info[e->index].mark < 0)){
            continue;
		}
		
	    q = k == 0 ? qm : qp;
	    for (i = 0; i < N; i++) {
		val = 0.0;
		I = phgSolverMapE2G(solver, k, e, i);
		indices[i] = I;
		/* - \int_T f1 v */
		val = -phgQCIntegrate(q, eno, Q_BAS, i, q, eno, Q_f, 0);
		/*\epsilon \int_T \grad un \grad v*/
        val += phgQCIntegrate(q, eno, Q_gradun, 0, q, eno, Q_GRAD, i) * (k == 0 ? Epsilon_p : Epsilon_s);
		if(k == 1){
			/*\Kappa^2 \int_{T \cup \Omega+} sinh(un) v*/
			val += phgQCIntegrate(q, eno, Q_gsinh, 0, q, eno, Q_BAS, i) * a_noline;
		}
		values[i] = val;
	    }
		phgVecAddGlobalEntries(F, 0, N, indices, values);
	}
	phgFree(rule_m);
	phgFree(rule_p);
	rule_m = rule_p = NULL;
	/*==================== face integrals ====================*/

	if (xi->info[e->index].mark == 0) {	/* the interface */
	    phgQCSetRule(qm, rule_0, -1.);
	    phgQCSetRule(qp, rule_0, -1.);
	    Cal_face(solver, xi, F, qm, e, -1, /*dof_id*/0,  Epsilon_p,
								qp, e, -1, /*dof_id1*/1, Epsilon_s);
	}

	phgFree(rule_0);
	rule_0 = NULL;
	/* the faces of the element */
	for (face = 0; face < NFace; face++) {
	    FLOAT nv[Dim];
	    ELEMENT *e1 = phgGetNeighbour(g, e, face);
	    int face1 = -1;
	    if (e1 != NULL) {
		/* a face is only processed by the smaller of the two
		 * neighbouring elements, and by the one with smaller
		 * global index if the two elements are of the same size,
		 * to avoid double counting or redundant computation */
		if (e->generation < e1->generation)
		    continue;
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index))
		    continue; /* process each interior face just once */
		face1 = phgOppositeFace(g, e, face, e1);
	    }
	    phgGeomGetFaceOutNormal(g, e, face, nv);

	    if (xi->info[e->index].mark == 0)
		/* the face might be cut into 2 parts by the interface */
		phgQuadInterfaceFace(xi->ls, xi->ls_grad, e, face,
				     order + Extra_orders,
				     &rule_m, NULL, &rule_p);
	    else if (xi->info[e->index].mark < 0)
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face, order,
				     &rule_m, NULL, NULL);
	    else
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face, order,
				     NULL, NULL, &rule_p);
	    /* Note: e1 is either NULL or valid (since gen(e) >= gen(e1)) */
	    if (e1 != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark < 0) { 
		/* The face is flat and is part of the interface
		 * (e and e1 in different subdomains) */
		if (xi->info[e->index].mark < 0) {
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    phgQCSetRule(qp, rule_m, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    Cal_face(solver, xi, F, qm, e,  face,  0, Epsilon_p,
									qp, e1, face1, 1, Epsilon_s);
		}
		else {
		    nv[0] = -nv[0]; nv[1] = -nv[1]; nv[2] = -nv[2];
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    phgQCSetRule(qm, rule_p, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    Cal_face(solver, xi, F, qm, e1, face1, 0, Epsilon_p,
									qp, e,  face,  1, Epsilon_s);
		}
	    }
	    else {  
		/* e and e1 in a same subdomain */
		if (xi->info[e->index].mark <= 0) {
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    Cal_face(solver, xi, F, qm, e,  face,  0, Epsilon_p,
										   qm, e1, face1, 0, Epsilon_p);
		}
		if (xi->info[e->index].mark >= 0) {
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    Cal_face(solver, xi, F, qp, e,  face,  1, Epsilon_s,
										   qp, e1, face1, 1, Epsilon_s);
		}
	    }
	    phgFree(rule_m);
	    phgFree(rule_p);
	    rule_m = rule_p = NULL;
	}
	//phgVecAssemble(F);
    }ForAllElementsEnd
    
    Free_lambdaFEM();
    phgQCFree(&qm);
    phgQCFree(&qp);
    
    phgDofFree(&f1_h);
	phgDofFree(&f2_h);
}

static void 
Evaluation_Residual(XFEM_INFO *xi, DOF *u1, DOF *u2, FLOAT *value){
	assert(u1 != NULL && u2 != NULL);
	VEC *F;
	MAP *map;
	DOF *u1_grad, *u2_grad;
    SOLVER *solver;

	u1_grad = phgDofGradient(u1, NULL, NULL, "u1_grad");
	u2_grad = phgDofGradient(u2, NULL, NULL, "u2_grad");

	solver = phgSolverCreate(SOLVER_DEFAULT, u1, u2, NULL);
	solver->mat->handle_bdry_eqns = FALSE;
    map = solver->rhs->map;
    
	F = phgMapCreateVec(map, 1);

	DOF *u1n_old, *u2n_old, *u1n_grad_old, *u2n_grad_old;

	u1n_old = u1n;
	u2n_old = u2n;
	u1n_grad_old = u1n_grad;
	u2n_grad_old = u2n_grad;

	u1n = u1;
	u2n = u2;
    u1n_grad = u1_grad;
	u2n_grad = u2_grad;
	
	Assemble(xi, solver, u1, u2, F);
	
	u1n = u1n_old;
	u2n = u2n_old;
	u1n_grad = u1n_grad_old;
	u2n_grad = u2n_grad_old;

	*value = phgVecNorm2(F, 0, NULL);
	*value = 0.5 * (*value) * (*value);

	phgVecDestroy(&F);
	phgSolverDestroy(&solver);
	phgDofFree(&u1_grad);
	phgDofFree(&u2_grad);
	return;
}

/*****************************************************************Solve_NPBE****************************************************************************************/
void Solve_NPBE(XFEM_INFO *xi, DOF *u1_h, DOF *u2_h, DOF *u1h_old, DOF *u2h_old, INT printtype){
    FLOAT L2err = 1.0, H1err = 1.0;
	FLOAT val_res = 0.0, val_old = 0.0;
	int iterate_tol = 1, iterate_Newton = 1;
	SOLVER *solver = NULL;
	double t_tol = phgGetTime(NULL), t_ite;
	DOF  *u1_old = NULL, *u2_old = NULL;

	phgPrintf("\nSolve the linear PB equation as initial value\n");
	/*Solve the linear PB equation as initial value*/
	
    Solve_LPBE(xi, u1_h, u2_h, printtype - 1);

	L2err = 1.0;
	H1err = 1.0;
	
	u1_old = phgDofCopy(u1_h, NULL, NULL, "u1_old");
	u2_old = phgDofCopy(u2_h, NULL, NULL, "u2_old");
	u1n = u1_old;
	u2n = u2_old;
	if(u1n_grad != NULL){
		phgDofFree(&u1n_grad);
	}
	u1n_grad = phgDofGradient(u1n, NULL, NULL, "u1n_grad");
	if(u2n_grad != NULL){
		phgDofFree(&u2n_grad);
	}
	u2n_grad = phgDofGradient(u2n, NULL, NULL, "u2n_grad");

	Evaluation_Residual(xi, u1_old, u2_old, &val_old);
	phgPrintf("res_init = %e\n", (double)val_old);
	if(val_old < tolerance_Newton){
		return;
	}

/************************Newton*************************/
	DOF *du1, *du2;
	FLOAT h = 1e-5, rho = 0.618, gamma_init = 1.0, gamma_min = 0.05, gamma = 1.0;
	du1 = phgDofNew(g, u1_h->type, 1, "du1", DofInterpolation);
	du2 = phgDofNew(g, u2_h->type, 1, "du2", DofInterpolation);
    
	while(iterate_tol < iterate_Max){
		t_ite = phgGetTime(NULL);
		phgPrintf("\n==================%d Newton, %d tol iterate=====================\n", iterate_Newton, iterate_tol);


		solver = phgSolverCreate(SOLVER_DEFAULT, du1, du2, NULL);
	    solver->mat->handle_bdry_eqns = FALSE;

        phgPrintf("Building linear equations ...\n");

		u1n = u1_old;
	    u2n = u2_old;
	    if(u1n_grad != NULL){
			phgDofFree(&u1n_grad);
	    }
		u1n_grad = phgDofGradient(u1n, NULL, NULL, "u1n_grad");
		if(u2n_grad != NULL){
			phgDofFree(&u2n_grad);
		}
		u2n_grad = phgDofGradient(u2n, NULL, NULL, "u2n_grad");
		
		Build_linear_system_NPBE(xi, solver, du1, du2, 2);

        phgPrintf("Solving linear equations ...\n");
		phgXFEMSetInitialSolution(xi, du1, du2);
	    phgSolverSolve(solver, TRUE, du1, du2, NULL);
	    phgXFEMUpdateSolution(xi, du1, du2);

	    phgSolverDestroy(&solver);
        // DOF *tmp_u1, *tmp_u2;
		// tmp_u1 = phgDofCopy(u1_h, NULL, NULL, "tmp_u1");
		// tmp_u2 = phgDofCopy(u2_h, NULL, NULL, "tmp_u2");
		// phgDofAXPY(gamma, du1, &tmp_u1);
		// phgDofAXPY(gamma, du2, &tmp_u2);
		// //phgPrintf("success!\n");
		// Evaluation_Residual(xi, tmp_u1, tmp_u2, &val_res);
		// //phgPrintf("gamma update = %lf\n", (double)gamma);
		// //Evaluation_Residual(xi, u1_old, u2_old, &val_old);
		// phgPrintf("res_update = %e, res_uk = %e, res_update - (1 - 2 * h * gamma * res_uk) = %e\n",(double)val_res, (double)val_old, (double)(val_res - (1 - 2 * h * gamma) * val_old));
		// if((val_res > (1 - 2 * h * gamma) * val_old) && (gamma > gamma_min)){
        //     gamma = rho * gamma;
		// 	phgPrintf("gamma update = %lf\n", (double)gamma);
		// }
		// val_old = val_res;

		gamma = gamma_init;
		DOF *tmp_u1, *tmp_u2;
		tmp_u1 = phgDofCopy(u1_h, NULL, NULL, "tmp_u1");
		tmp_u2 = phgDofCopy(u2_h, NULL, NULL, "tmp_u2");
		phgDofAXPY(gamma, du1, &tmp_u1);
		phgDofAXPY(gamma, du2, &tmp_u2);
		Evaluation_Residual(xi, tmp_u1, tmp_u2, &val_res);
		phgPrintf("gamma = %lf, res(u^k + gamma du) = %e\n", (double)gamma, (double)val_res);
		while((val_res > (1 - 2 * h * gamma) * val_old) && (gamma > gamma_min)){
			gamma = gamma * rho;
			phgDofFree(&tmp_u1);
			phgDofFree(&tmp_u2);
			tmp_u1 = phgDofCopy(u1_h, NULL, NULL, "tmp_u1");
			tmp_u2 = phgDofCopy(u2_h, NULL, NULL, "tmp_u2");
			phgDofAXPY(gamma, du1, &tmp_u1);
			phgDofAXPY(gamma, du2, &tmp_u2);
			Evaluation_Residual(xi, tmp_u1, tmp_u2, &val_res);
			phgPrintf("gamma = %lf, res(u^k + gamma du) = %e\n", (double)gamma, (double)val_res);
		}
		phgDofAXPY(gamma, du1, &u1_h);
		phgDofAXPY(gamma, du2, &u2_h);
		phgDofFree(&tmp_u1);
		phgDofFree(&tmp_u2);

		if((val_res < tolerance_Newton)){
        	cal_error(xi, u1_h, u2_h, u1_old, u2_old, NULL, NULL, &L2err, &H1err);
        	phgPrintf("L2err: %0.10le; H1err: %0.10le\n", (double)L2err, (double)H1err);
			if(L2err < tolerance_Newton){
				phgPrintf("After %d step, Newton iterate converge!\n", iterate_Newton);
				break;
			}
		}
		
		phgDofFree(&u1_old);
		phgDofFree(&u2_old);
		u1_old = phgDofCopy(u1_h, NULL, NULL, "u1_old");
		u2_old = phgDofCopy(u2_h, NULL, NULL, "u2_old");
		iterate_Newton++;
		iterate_tol++;
		phgPrintf("use time:  %0.2lgs\n", phgGetTime(NULL) - t_ite);
	}

	if(iterate_tol >= iterate_Max && (L2err >= tolerance_Newton || val_res >= tolerance_Newton)){
		phgPrintf("\nNewton not converge!!\n");
		exit(0);
	}

	else{
		phgPrintf("\nAfter %d Newton itegrate, Nolinear PB Equation converge! use %0.2lgs\n", iterate_Newton, (phgGetTime(NULL) - t_tol));
	}
	phgDofFree(&u1_old);
	phgDofFree(&u2_old);
	phgDofFree(&du1);
	phgDofFree(&du2);
	phgDofFree(&u1n_grad);
	phgDofFree(&u2n_grad);
	return;
}