#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "phg.h"
#include "IPFEMPB.h"
#define OUT_LOG 1


GRID *g = NULL;
DOF *u1_h = NULL;
DOF *u2_h = NULL;
DOF *ls = NULL;
DOF *ls_grad = NULL;
XFEM_INFO *xi = NULL;

static DOF *u1_old = NULL;
static DOF *u2_old = NULL;

int main(int argc,char *argv[]){

    int level = 0, total_level = 0;
    double t = 0.0, t_tol = phgGetTime(NULL), t_level = 0.0;
	BOOLEAN remove_IP = FALSE;

    /*Load phg options*/
    OptionsRegister();
    
	/*=============================PHG=============================================*/
    phgInit(&argc, &argv);

    /*Initialize and Print*/
    Initialize();
	char name[50];
	if(!use_aly){
	char tmp;
	int i = 0;
	for(i = 0; i < 50; i++){
		tmp = fn_pqr[i];
		if(tmp != '.'){
			name[i] = tmp;
		}
		else{
			name[i] = '\0';
			break;
		}
	}
	}
	else{
		sprintf(name, "%s%s", "test_", PB_type);
	}
#if OUT_LOG
    char log_name[50];
	sprintf(log_name, "%s%s", name, ".log");
	//remove(log_name);
	FILE *fp_out = freopen(log_name,"w",stdout);
    if(NULL == fp_out)
    {
        printf("reopen failed\n");
        return -1;
    }
#endif
	PrintConstant();
	PrintOptions();

    g = phgNewGrid(-1);

	if(fn_mesh == NULL){
		fn_mesh = Creat_Domin(Protein_MaxXYZ, Protein_MaxR);
		remove_IP = TRUE;
	}
	
    if (!phgImport(g, fn_mesh, TRUE)){
        phgError(1, "can't read file \"%s\".\n", fn_mesh);
	}
	if(box_size > 0.0){
		refine0 = Max(refine0, 3 * ((int)(log(2*Domain_length/box_size) / log(2.0)) + 3));
	}
    while(refine0 > 0){
		//phgPrintf("%d\n", __LINE__);
		level = refine0 > refine_step ? refine_step : refine0;
		/* uniform refinement */ 
		phgRefineAllElements(g, level);
    	phgBalanceGrid(g, 1.1, 1, NULL, 1.0);
		refine0 -= level;
		total_level += level;
	}
	
	FLOAT h_max = 0.0, h_min = 0.0;
	Get_h(g, &h_max, &h_min);
	phgPrintf("After pre refine: %"dFMT" element%s. h_max = %lf\n", g->nelem_global, g->nelem_global > 1 ? "s" : "", (double)h_max);

    /* project the levelset function to a FE space */
    assert(ls_order >= 1);
	ls = phgDofNew(g, DOF_Pn[ls_order], 1, "ls", ls_func);
    ls_grad = phgDofNew(g, DOF_Pn[ls_order - 1], 3, "ls_grad", ls_grad_func);
	phgDofSetDataByFunction(ls, ls_func);
	phgDofSetDataByFunction(ls_grad, ls_grad_func);

	u1_h = phgDofNew(g, DOF_DEFAULT, 1, "u1_h", DofInterpolation);
    u2_h = phgDofNew(g, DOF_DEFAULT, 1, "u2_h", DofInterpolation);
    phgPrintf("Initialization time %0.4lg\n", phgGetTime(NULL) - t_tol);

	t_level = phgGetTime(NULL);
    t = phgGetTime(NULL);
	phgPrintf("Resolving/merging interface elements: \n");

	xi = phgXFEMInit(ls, ls_grad, ls_order, 2 * DOF_DEFAULT->order + Extra_orders);

	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	phgSetupHalo(g, HALO_FACE);
	phgPrintf("%d proc%s, %"dFMT" element%s, LIF %0.2lf, refine time: %0.4lg\n",g->nprocs, g->nprocs > 1 ? "s" : "", 
	          g->nelem_global, g->nelem_global > 1 ? "s" : "", (double)g->lif, phgGetTime(NULL) - t);
	phgPrintf("    Elements by region: -:%"dFMT", 0:%"dFMT", +:%"dFMT"\n", xi->ecnts[0], xi->ecnts[1], xi->ecnts[2]);

	phgDofSetDataByFunction(ls, ls_func);
	phgDofSetDataByFunction(ls_grad, ls_grad_func);

	if(remove_IP){
		remove(fn_mesh);
	}

	if(PB_Type == 2){
		Solve_NPBE(xi, u1_h, u2_h, NULL, NULL, PrintType);
	}
	else if(PB_Type == 1){
		Solve_LPBE(xi, u1_h, u2_h, PrintType);
	}
	else{
		phgError(-1, "Wrong PB_type!\n");
	}

	if(!use_aly){
    	Cal_energy(u1_h);
	}
	else{
		FLOAT L2_norm = 0.0, H1_norm = 0.0, L2_error = 0.0, H1_error = 0.0;
		u1_old = phgDofNew(g, DOF_Pn[DOF_DEFAULT->order + 2], 1, "u1_aly", func_u1);
		u2_old = phgDofNew(g, DOF_Pn[DOF_DEFAULT->order + 2], 1, "u2_aly", func_u2);
		phgPrintf("\nerror of ||u_h - u||\n");
		cal_error(xi, u1_h, u2_h, u1_old, u2_old, &L2_norm, &H1_norm, &L2_error, &H1_error);
		phgPrintf("L2err: %0.10le; H1err: %0.10le; L2Norm: %0.10le, H1Norm: %0.10le\n",
		    (double)L2_error, (double)H1_error, (double)L2_norm, (double)H1_norm);
		phgDofFree(&u1_old);
		phgDofFree(&u2_old);
	}
	size_t mem_peak = 0.0;
	phgMemoryUsage(g, &mem_peak);
	phgPrintf("mem: %0.2lfGB; total level time: %0.2lgs\n", (double)mem_peak / (1024. * 1024. * 1024.), phgGetTime(NULL) - t_level);

    /*export vtk*/
	if(!use_aly && vtk){
		ELEMENT *e;
		char name_vtk[50];
		sprintf(name_vtk, "%s%s", name, ".vtk");
		ForAllElements(g, e)
			e->region_mark = xi->info[e->index].mark;
		//phgDofAXPBY(0.0, u2_h, (1/Beta), &u2_h);
		u2_h->name = "phi";
		phgExportVTK(g, name_vtk, u2_h, ls, NULL);
		phgPrintf("\"%s\" created.\n", name_vtk);
	}
#if OUT_LOG
	fclose(fp_out);
#endif
	phgXFEMFree(&xi);
    
    phgDofFree(&u1_h);
    phgDofFree(&u2_h);
    phgDofFree(&ls);
    phgDofFree(&ls_grad);
    phgFreeGrid(&g);
    phgFinalize(); 
    return 0;
}