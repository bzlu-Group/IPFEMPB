#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#include"phg.h"
#include"functions_LPBE.h"
static int Q_f, Q_gD, Q_jD, Q_jN;

void DoFace_LPBE(SOLVER *solver, XFEM_INFO *xi,
				QCACHE *qc,  ELEMENT *e,  int face,  int dof_id,  FLOAT coe,
				QCACHE *qc1, ELEMENT *e1, int face1, int dof_id1, FLOAT coe1){
	
	GRID *g;
    DOF *u_h, *u1_h;
    int n, n1, p;   /* p = polynomial order */
    INT I, J;
    int i, j;
    FLOAT val;
    FLOAT G0, G1, h, a; /* G0 = coeff*gamma0*p^2/h, G1 = coeff*Gamma1_XFEM*h/p^2 */
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

    if (e1 == NULL) {
	/* boundary face */
	bdry_flag = DIRICHLET;
	n1 = 0;	
	G0 = coe * Gamma0_XFEM * p * p / h;
	G1 = coe * Gamma1_XFEM * h / (p * (FLOAT)p);
	/* RHS */
	for (i = 0; i < n; i++) {
	    I = phgSolverMapE2G(solver, dof_id, e, i);
	    if (bdry_flag == DIRICHLET){	/* Dirichlet boundary */
		/* -\Theta_XFEM\int_\Gamma_D g_D (A\grad v).n */
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
	n1 = qc->qd->n_bas(qc1, e1->index);
	i = DofTypeOrder(u1_h, e1);
	if (p < i)
	    p = i;
	G0 = 0.5 * (coe + coe1) * Gamma0_XFEM * p * (FLOAT)p / h;
	G1 = 0.5 * (coe + coe1) * Gamma1_XFEM * h / (p * (FLOAT)p);
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
		xi->info[e->index].mark * xi->info[e1->index].mark >= 0){
			assert(bdry_flag == INTERIOR && dof_id == dof_id1);
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

		/* -\int \Theta_XFEM [u]{A\grad v}.n (note: i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_BAS, PROJ_NONE, j, Q_GRAD,  PROJ_DOT, i) * Coe(i);
		if (bdry_flag & INTERIOR)
		    a *= (j < n ? 0.5 : -0.5);
		val += -Theta_XFEM * a;

		/* \int G0 [u][v] (i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_BAS, PROJ_NONE, j, Q_BAS, PROJ_NONE, i);
		if ((bdry_flag & INTERIOR) && (i < n) != (j < n))
		    a = -a;
		val += G0 * a;
	    }

	    if (bdry_flag != DIRICHLET) {
		/* \int G1 [A\grad u].n [A\grad v].n */
		a = Quad(Q_GRAD,  PROJ_DOT, j, Q_GRAD,  PROJ_DOT, i)
			* Coe(j) * Coe(i);
		if (bdry_flag & INTERIOR && (i < n) != (j < n))
		    a = -a;
		val += G1 * a;
	    }
	    phgSolverAddGlobalMatrixEntry(solver, I, J, val);
	}

	if (e1 == NULL || dof_id == dof_id1)
	    continue;		/* non interface face */

	/* this face is part of the interface, apply jump conditions */
	assert(dof_id == 0 && dof_id1 == 1);
	assert((e == e1 && xi->info[e->index].mark == 0) ||
	       (xi->info[e->index].mark < 0 && xi->info[e1->index].mark > 0));

	/* \int [A*gN.n] {v} */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,  PROJ_DOT,  0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	val = a * 0.5;

	/* G1 \int [A*gN.n] [(A\grad v).n] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,   PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	val += Coe(i) * G1 * (i < n ? a : -a);

	/* -Theta_XFEM \int [gD] {(A\grad v).n} (func[1] := gD) */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,   PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	val += -Theta_XFEM * 0.5 * Coe(i) * a;

	/* G0 \int [gD] [v] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	val += G0 * (i < n ? a : -a);

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

static void Build_linear_system_LPBE(XFEM_INFO *xi, SOLVER *solver, DOF *u1_h, DOF *u2_h){
	GRID *g = xi->ls->g;
    QCACHE *qm, *qp;
    ELEMENT *e;
    assert(u1_h->dim == 1 && u2_h->dim == 1);

    qm = phgXFEMNewQC(xi, u1_h, 0);
    Q_f = phgQCAddXYZFunction(qm, func_f1, 1);/* f */
    Q_gD = phgQCAddXYZFunction(qm, func_g1D, 1);		/* gD */
    Q_jD = phgQCAddXYZFunction(qm, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qm, func_jN, 3);		/* jN */

    qp = phgXFEMNewQC(xi, u2_h, 1);
    /* The following must be in exactly the same order as above */
    Q_f = phgQCAddXYZFunction(qp, func_f2, 1);	/* f  */
    Q_gD = phgQCAddXYZFunction(qp, func_g2D, 1);		/* gD */
    Q_jD = phgQCAddXYZFunction(qp, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qp, func_jN, 3);		/* jN */
	
#if USE_OMP
#pragma omp parallel for private(e)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {
	FLOAT *rule_m = NULL, *rule_0 = NULL, *rule_p = NULL;
	int N = qm->qd->n_bas(qm, e->index);
	INT I, J, eno = e->index;
	int i, j, k, face, order;
	FLOAT val;
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
	    if ((k == 0 && xi->info[e->index].mark > 0) ||
		(k == 1 && xi->info[e->index].mark < 0))
		continue;
	    q = k == 0 ? qm : qp;
	    for (i = 0; i < N; i++) {
		I = phgSolverMapE2G(solver, k, e, i);
		for (j = 0; j <= i; j++) {
		    /* \int_T A\grad u1 . \grad v */
		    J = phgSolverMapE2G(solver, k, e, j);
		    val = phgQCIntegrate(q, eno, Q_GRAD, j, q, eno, Q_GRAD, i) * (k == 0 ? Epsilon_p : Epsilon_s);
			if(k == 1){
				val += phgQCIntegrate(q, eno, Q_BAS, j, q, eno, Q_BAS, i) * a_noline;
			}
		    phgSolverAddGlobalMatrixEntry(solver, I, J, val); 
		    if (i != j)
				phgSolverAddGlobalMatrixEntry(solver, J, I, val);

		}
		/* \int_T f1 v */
		val = phgQCIntegrate(q, eno, Q_BAS, i, q, eno, Q_f, 0);
		phgSolverAddGlobalRHSEntry(solver, I, val);
	    }
	}
	phgFree(rule_m);
	phgFree(rule_p);
	rule_m = rule_p = NULL;
	/*==================== face integrals ====================*/

	if (xi->info[e->index].mark == 0) {	/* the interface */
	    phgQCSetRule(qm, rule_0, -1.);
	    phgQCSetRule(qp, rule_0, -1.);
	    DoFace_LPBE(solver, xi, qm, e, -1, /*dof_id*/0,  Epsilon_p,
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
		    DoFace_LPBE(solver, xi, qm, e,  face,  0, Epsilon_p,
					qp, e1, face1, 1, Epsilon_s);
		}
		else {
		    nv[0] = -nv[0]; nv[1] = -nv[1]; nv[2] = -nv[2];
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    phgQCSetRule(qm, rule_p, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    DoFace_LPBE(solver, xi, qm, e1, face1, 0, Epsilon_p,
					qp, e,  face,  1, Epsilon_s);
		}
	    }
	    else {
		/* e and e1 in a same subdomain */
		if (xi->info[e->index].mark <= 0) {
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    DoFace_LPBE(solver, xi, qm, e,  face,  0, Epsilon_p,
					qm, e1, face1, 0, Epsilon_p);
		}
		if (xi->info[e->index].mark >= 0) {
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    DoFace_LPBE(solver, xi, qp, e,  face,  1, Epsilon_s,
					qp, e1, face1, 1, Epsilon_s);
		}
	    }
	    phgFree(rule_m);
	    phgFree(rule_p);
	    rule_m = rule_p = NULL;
	}
    } ForAllElementsEnd

    phgQCFree(&qm);
    phgQCFree(&qp);

    /* set the diagonal entry of empty rows to 1 */
    phgXFEMProcessEmptyRows(solver);
}

void Solve_LPBE(XFEM_INFO *xi, DOF *u1_h, DOF *u2_h, INT printtype){
	/*Solve LPBE*/
    SOLVER *solver;
    
	//phgOptionsSetOptions("-solver_symmetry=spd");

	solver = phgSolverCreate(SOLVER_DEFAULT, u1_h, u2_h, NULL);
	solver->mat->handle_bdry_eqns = FALSE;
    phgPrintf("Building linear equations ...\n");
	Build_linear_system_LPBE(xi, solver, u1_h, u2_h);
  
    phgPrintf("Solving linear equations ...\n");

    phgXFEMSetInitialSolution(xi, u1_h, u2_h);
	phgSolverSolve(solver, TRUE, u1_h, u2_h, NULL);
	phgXFEMUpdateSolution(xi, u1_h, u2_h);
    
	phgSolverDestroy(&solver);
}