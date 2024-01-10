#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdbool.h>
#include "phg.h"
#include "unit.h"
#include "coefficient.h"

/*options for User*/
extern FLOAT C_b;/*ionic concentration mol/L*/
extern char *fn_pqr;/*pqr file*/
extern FLOAT Domain_length;/*Domain: [-Domain_lenth, Domain_lenth]^3*/
extern INT PB_Type;/*1: LPBE; 2: NPBE. */
extern BOOLEAN use_Cube;/*if use acceleration skill*/
extern BOOLEAN vtk;/*If or not output vtk-file*/
extern char *fn_mesh;/*mesh file*/
extern BOOLEAN use_aly;/*Only used in the box-sphere analytic solution algorithm.*/
extern FLOAT T;
extern FLOAT box_size;
extern char *PB_type;

/*other parameters*/
extern INT Max_N_m;/*max Number of atoms*/
extern FLOAT Epsilon_s;/*Relative dielectric constant of solent*/
extern FLOAT Epsilon_p;/*Relative dielectric constant of solute*/
extern FLOAT Am;/*unit length*/
extern FLOAT c_Gauss, d_Gauss;/*parameters of Gaussian surface*/
extern INT ls_order;/*projection order of level set function*/
extern FLOAT tolerance_Newton;
extern INT iterate_Max;
extern FLOAT Relax_coefficient;
extern INT Theta_XFEM;/*Coefficient of symmetry term - interface*/
extern FLOAT Gamma0_XFEM;/*Jump coefficient of order 0 - interface*/
extern FLOAT Gamma1_XFEM;/*Jump coefficient of order 1 - interface*/
extern INT Extra_orders;/*extra integral order of interface*/
extern BOOLEAN dump_solver;
extern INT refine0, refine_step;
extern INT PrintType;


extern void OptionsRegister();
extern void PrintOptions();

#endif
