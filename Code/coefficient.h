#ifndef COEFFICIENT_H 
#define COEFFICIENT_H

#include "phg.h"
#include "unit.h"
#include "options.h"
#include "functions_general.h"

typedef struct{

    INT number;
    FLOAT zi;/*valence*/
    FLOAT loc[3];/*location*/
	FLOAT r;/*radius*/

}ATOM;

#define PAI _F(3.14159265358979323846264338327950288) /*pai*/
#define E _F(2.71828182845904523536) /*e*/
#define Epsilon_0  (8.8541878e-12)/*Vacuum Permittivity： kg-1 m-3 s4 A2*/
#define E_c (1.6021765e-19)/*Unit charge： s A*/

#define K_B (1.3806488e-23)/*Boltzmann constant： kg m2 s-2 k-1*/
#define Na (6.02214129e23) /*Avogadro constant*/

extern FLOAT Beta;/*Dimensionless parameter Beta = E_c/(K_B*T)*/
extern FLOAT Alpha;/*Dimensionless parameter Alpha = E_c*E_c/(Epsilon_0*K_B*T)*/
extern FLOAT Kappa;/*Debye-Huckel parameter Kappa = (2*C_b*E_c*E_c/(K_B*T*Epsilon_s*Epsilon_0))^(1/2)*/
extern FLOAT a_noline;/*Dimensionless parameter a_noline = 2 * C_b * Alpha  * Na * 1000.0 * Am * Am * Am*/
extern FLOAT I_s;/*ionic strength I_s = 0.5 * \sum (c_i * z_i * z_i)*/
extern FLOAT Energy_sol;/*Solvation energy*/
extern INT N_m;/*Number of atoms*/
extern FLOAT Protein_MaxXYZ[3];
extern FLOAT Protein_MaxR;

extern ATOM *atoms;

void Initialize();

void PrintConstant();

#endif
