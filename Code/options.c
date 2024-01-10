#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "phg.h"

#include "unit.h"
#include "options.h"

/*options for User*/
FLOAT Epsilon_s = 80.0;/*Relative dielectric constant of solent*/
FLOAT Epsilon_p = 2.0;/*Relative dielectric constant of solute*/
FLOAT C_b = 0.10;/*ionic concentration mol/L*/
FLOAT T = 298.15;/*Temperature： k*/
char *fn_pqr = NULL;/*pqr file*/
char *PB_type = NULL;
INT PB_Type = 1;/*1: LPBE; 2: NPBE. */
FLOAT Domain_length = -1.0;/*Domain: [-Domain_lenth, Domain_lenth]^3*/
FLOAT box_size = -1.0;
BOOLEAN vtk = FALSE;/*If or not output vtk-file*/

/*other parameters*/
INT Max_N_m = 10000;/*max Number of atoms*/
FLOAT Am = 1e-10;/*unit length*/
FLOAT c_Gauss = 1.0, d_Gauss = 1.0;/*parameters of Gaussian surface*/
INT ls_order = 2;/*projection order of level set function*/
FLOAT tolerance_Newton = 1e-6;
INT iterate_Max = 1000;
FLOAT Relax_coefficient = 1.0;
INT Theta_XFEM = 1;/*Coefficient of symmetry term - interface*/
FLOAT Gamma0_XFEM = 10.0;/*Jump coefficient of order 0 - interface*/
FLOAT Gamma1_XFEM = 0.1;/*Jump coefficient of order 1 - interface*/
INT Extra_orders = 3;/*extra integral order of interface*/
BOOLEAN dump_solver = FALSE;
INT refine0 = 12, refine_step = 3;
INT PrintType = 3;
char *fn_mesh = NULL;/*mesh file*/
BOOLEAN use_aly = FALSE;/*Only used in the box-sphere analytic solution algorithm.*/

void
OptionsRegister(){
    phgOptionsRegisterFloat("-epsilon_solute", "the relative permittivity of solute.", &Epsilon_p);
    phgOptionsRegisterFloat("-epsilon_solvent", "the relative permittivity of solvent.", &Epsilon_s);
    phgOptionsRegisterFloat("-concentration", "the ionic concentration. unit: mol/L", &C_b);
    phgOptionsRegisterFloat("-temp", "the temperature. unit: K", &T);
    phgOptionsRegisterFilename("-pqr", "the PQR file", &fn_pqr);
    phgOptionsRegisterFilename("-PB_type", "the type of PBE (LPBE or NPBE) to be solved", &PB_type);
    phgOptionsRegisterFloat("-box_length", "denotes the size of the calculation domain $[ -box_length, box_length]^3$", &Domain_length);
    phgOptionsRegisterFloat("-box_size", "the size of the initial mesh", &box_size);
    phgOptionsRegisterInt("-ls_order", "the order of polynomial interpolation for Gaussian molecular surface ", &ls_order);
    phgOptionsRegisterNoArg("-vtk", "If or not output vtk-file", &vtk);

    phgOptionsRegisterNoArg("-use_aly", "only used in the box-sphere analytic solution algorithm", &use_aly);
	phgOptionsRegisterFilename("-fn_mesh", "mesh file", &fn_mesh);
    phgOptionsRegisterInt("-refine0", "The number of times the initial grid is uniformly refined, \
                        each time for dichotomous refinement, it is recommended to choose a multiple of 3", &refine0);
    phgOptionsPreset("-dof_type P1 -solver gmres -gmres_type flgmres_r -gmres_pc_type solver -gmres_pc_opts {-solver asm -asm_overlap 5} -lgmres_augk 7");
    phgOptionsPreset("-xfem_ctol 1.0");
    phgOptionsPreset("-xfem_dbg_level=0");
}

void PrintOptions(){
    phgPrintf("\n****************************OPTIONS*****************************\n");
    phgPrintf("epsilon_solute = %lf, epsilon_solvent = %lf\n",  (double)Epsilon_p, (double)Epsilon_s);
    phgPrintf("concentration = %lf mol/L\n",(double)C_b);
    phgPrintf("temperature = %lf K\n",(double)T);
    phgPrintf("pqr = %s\n", fn_pqr);
    phgPrintf("solve %s\n", PB_type);
    phgPrintf("\n*****************************************************************\n");
}