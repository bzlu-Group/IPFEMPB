#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#include"phg.h"
#include"functions_general.h"


char* 
Creat_Domin(FLOAT *protein_len, FLOAT protein_rmax){
	/*Create a Calculation Area*/
	if(fn_mesh != NULL){
		return fn_mesh;
	}
	if(!use_aly){
	int i = 0;
	BOOLEAN user = TRUE;
	FLOAT length = 0.0;
	for(i = 0; i < 3; i++){
		if(Domain_length < (FLOAT)((int)(length + protein_rmax) + 6)){
			user = FALSE;
		}
		length = Max(length, protein_len[i]);
	}
	if(!user){
		if(Domain_length > 0.0){
			phgPrintf("<<<<<<<<<<<< warning >>>>>>>>>>>>>\n");
			phgPrintf("The current computational region cannot cover the entire protein or boundary too close to protein.\nThe program will use the default computational region = [%lf, %lf]^3\n",
				(double)-Domain_length, (double)Domain_length);
		}
		Domain_length = (FLOAT)((int)(length + protein_rmax) + 6);
	}
		phgPrintf("\n");
	}
	else{
		Domain_length = 2.0;
	}
	static char mesh_name[100];
	sprintf(mesh_name, "cube_%d.mesh", (int)(Domain_length));
	fn_mesh = mesh_name;

	FILE *f_out = fopen(mesh_name, "w");

	fprintf(f_out, "MeshVersionFormatted 1\n");
	fprintf(f_out, "Dimension\n3\n\n");
	fprintf(f_out, "Vertices\n8\n\n");
	fprintf(f_out, "%lf %lf %lf 2\n", (double)-Domain_length, (double)-Domain_length, (double)-Domain_length);
	fprintf(f_out, "%lf %lf %lf 0\n", (double)Domain_length, (double)-Domain_length, (double)Domain_length);
	fprintf(f_out, "%lf %lf %lf 2\n", (double)Domain_length, (double)Domain_length, (double)-Domain_length);
	fprintf(f_out, "%lf %lf %lf 0\n", (double)-Domain_length, (double)Domain_length, (double)Domain_length);
	fprintf(f_out, "%lf %lf %lf 0\n", (double)-Domain_length, (double)-Domain_length, (double)Domain_length);
	fprintf(f_out, "%lf %lf %lf 2\n", (double)Domain_length, (double)-Domain_length, (double)-Domain_length);
	fprintf(f_out, "%lf %lf %lf 0\n", (double)Domain_length, (double)Domain_length, (double)Domain_length);
	fprintf(f_out, "%lf %lf %lf 2\n", (double)-Domain_length, (double)Domain_length, (double)-Domain_length);
	fprintf(f_out, "\nHexahedra\n1\n1 6 3 8 5 2 7 4 0\n\n");
	fprintf(f_out, "# quadrilaterals are used to specify boundary types (default: 0)\n");
	fprintf(f_out, "\nQuadrilaterals\n1\n5 2 7 4 1\n\n");
	fprintf(f_out, "END\n");
	fclose(f_out);
	return mesh_name;
}

void Read_pqr(char *fn_pqr){
	FLOAT t = phgGetTime(NULL);
    FILE *f_in;
    FLOAT charge=0.0;
	FLOAT sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    f_in = fopen(fn_pqr, "r");
	int i,j;
    if (f_in == NULL){
        phgError(1, "Can not open file %s.\n", fn_pqr);
    }
    
    N_m = 0;
	assert(atoms == NULL);
    if(!(atoms = (ATOM *)phgAlloc(Max_N_m * sizeof(*atoms)))){
        phgError(-1, "Error in memory allocation for fixed_chargess!\n");
		return;
	}
	double x = 0.0, y = 0.0, z = 0.0, zi = 0.0, r = 0.0;
	char s[5][100];
	//int number = 0;
	phgPrintf("****************************Begin read %s *********************************\n", fn_pqr);
    while(fscanf(f_in, "%s", s[0]) != EOF) {
		if(strcmp(s[0], "ATOM") != 0) continue;
		fscanf(f_in, "%s%s%s%s%lf%lf%lf%lf%lf",s[1],s[2],s[3],s[4],&x,&y,&z,&zi,&r);
		if(r < 1e-3){
			phgPrintf("warning %d atom: r = %lf\n", N_m, r);
		}
		if(N_m >= (Max_N_m - 1)){
			phgError(-1, "Too many atoms! %s%d\n", __FILE__, __LINE__);
		}
		atoms[N_m].number = N_m + 1;
		atoms[N_m].loc[0] = (FLOAT)x;
		atoms[N_m].loc[1] = (FLOAT)y;
		atoms[N_m].loc[2] = (FLOAT)z;
		atoms[N_m].zi = (FLOAT)zi;
		atoms[N_m].r = (FLOAT)r;
		charge += zi;
		sum_x += x;
		sum_y += y;
		sum_z += z;
		N_m ++;
    }
	sum_x /= N_m;
	sum_y /= N_m;
	sum_z /= N_m;

    fclose(f_in);
	for(i=0; i<N_m; i++){
		atoms[i].loc[0] -= sum_x;
		atoms[i].loc[1] -= sum_y;
		atoms[i].loc[2] -= sum_z;
		for(j = 0; j < 3; j++){
			Protein_MaxXYZ[j] = Max(Protein_MaxXYZ[j], Fabs((double)atoms[i].loc[j]));
		}
		Protein_MaxR = Max(Protein_MaxR, Fabs((double)atoms[i].r));
	}
	phgPrintf("N_m = %d, total_charge = %lf, max_txyz = (%lf, %lf, %lf), Protein_MaxR = %lf\n",N_m, (double)charge,
	          (double)Protein_MaxXYZ[0], (double)Protein_MaxXYZ[1], (double)Protein_MaxXYZ[2], (double)Protein_MaxR);
	phgPrintf("****************************End read %s  use %0.4lfs*********************************\n", fn_pqr, phgGetTime(NULL)-t);
}

void 
Get_h(GRID *g, FLOAT *h_max, FLOAT *h_min){
    /* get h_max and h_min */
    FLOAT d;
    ELEMENT *e;
    int i = 0;

	ForAllElements(g, e) {
	    FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
	    for (i = 0; i < Dim; i++) {
		d = corners[1][i] - corners[0][i];
		if (*h_max < d)
		    *h_max = d;
		if (*h_min > d)
		    *h_min = d;
	    }
	}

#if USE_MPI
	if (g->nprocs > 1) {
	    double tmp0[2] = {-*h_min, *h_max}, tmp[2];
	    MPI_Reduce(tmp0, tmp, 2, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	    *h_min = -tmp[0];
	    *h_max = tmp[1];
	}
	
#endif	/* USE_MPI */
return;
}

void 
Cal_energy(DOF *u1_h){
	/*calculate solvation energy*/
    int i;
    COORD pos[N_m];
    FLOAT u_r[N_m],  sum = 0.0;
    DOF *u1_r;
	FLOAT kBT = K_B * T * Na / 4.184e3; 
		
	for(i = 0; i < N_m; i++){
        pos[i][0] = atoms[i].loc[0];
        pos[i][1] = atoms[i].loc[1];
        pos[i][2] = atoms[i].loc[2];
    }

    u1_r = phgDofCopy(u1_h, NULL, NULL, NULL);
    phgInterGridDofEval(u1_r, N_m, pos, u_r, -1);
    sum = 0.0;
    for(i = 0; i < N_m; i++){
        sum += u_r[i] * atoms[i].zi;
    }
    sum *= 0.5;
    Energy_sol = sum * kBT;
    phgPrintf("\nSolvation energy for numetical solution is = %f kcal/mol\n", (double)(sum * kBT));
    phgDofFree(&u1_r);
}

void 
cal_error(XFEM_INFO *xi, DOF *u1_h, DOF *u2_h, DOF *u1_old, DOF *u2_old, FLOAT *L2_norm, FLOAT *H1_norm, FLOAT *L2_error, FLOAT *H1_error){

    DOF *error1, *error2, *gerror1, *gerror2, *gu_h;
    FLOAT L2norm, H1norm, L2err, H1err;
    
    /* error := u, projected to the FE space */
    error1 = phgDofCopy(u1_old, NULL, NULL, "error1");
	error2 = phgDofCopy(u2_old, NULL, NULL, "error2");
	L2norm = Sqrt(phgXFEMDot(xi, error1, error2, error1, error2));
	/* gerror := Grad(u), projected to the FE space */
    gerror1 = phgDofGradient(error1, NULL, NULL, "gerror1");
	gerror2 = phgDofGradient(error2, NULL, NULL, "gerror2");
	H1norm = Sqrt(phgXFEMDot(xi, gerror1, gerror2, gerror1, gerror2));
	H1norm += L2norm;

	int realornot = 0;
	FLOAT old_l2 = 0.0, old_h1 = 0.0;
	FLOAT d;
	d = Sqrt(phgXFEMDot(xi, u1_h, u2_h, u1_h, u2_h));
	if (L2norm < d){
		old_l2 = L2norm;
        L2norm = d;
		realornot = 1;
	}
	DOF *u1h_grad, *u2h_grad;
	u1h_grad = phgDofGradient(u1_h, NULL, NULL, "u1h_grad");
	u2h_grad = phgDofGradient(u2_h, NULL, NULL, "u2h_grad");
	d = L2norm + Sqrt(phgXFEMDot(xi, u1h_grad, u2h_grad, u1h_grad, u2h_grad));
	if (H1norm < d){
		old_h1 = H1norm;
	    H1norm = d;
		realornot = 1;
	}
	phgDofFree(&u1h_grad);
	phgDofFree(&u2h_grad);

	/* error := u_h - error */
	phgDofAXPY(-1.0, u1_h, &error1);
	phgDofAXPY(-1.0, u2_h, &error2);
	/* L2err := |error| */
	L2err = Sqrt(phgXFEMDot(xi, error1, error2, error1, error2));
    
	/* gerror := Grad(u_h) - gerror */
	gu_h = phgDofGradient(u1_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, gu_h, &gerror1);
	phgDofFree(&gu_h);
	gu_h = phgDofGradient(u2_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, gu_h, &gerror2);
	phgDofFree(&gu_h);
	/* H1err := L2err + |gerror| */
	H1err = Sqrt(phgXFEMDot(xi, gerror1, gerror2, gerror1, gerror2));
	H1err += L2err;
	if(L2_norm != NULL){
		*L2_norm = L2norm;
	}
	if(H1_norm != NULL){
		*H1_norm = H1norm;
	}
	if(L2_error != NULL){
		*L2_error = L2err / (L2norm == 0. ? 1.0 : (double)L2norm);
	}
	if(H1_error != NULL){
		*H1_error = H1err / (H1norm == 0. ? 1.0 : (double)H1norm);
	}
	phgDofFree(&error1);
	phgDofFree(&gerror1);
	phgDofFree(&error2);
	phgDofFree(&gerror2);
	return;
}