#include <math.h>
#include <string.h>


#include "coefficient.h"


FLOAT Beta = 0.0;
FLOAT Alpha = 0.0;
FLOAT Kappa = 0.0;
FLOAT a_noline = 0.0;
FLOAT I_s = 0.0;
FLOAT Energy_sol = 0.0;
ATOM *atoms = NULL;
INT N_m = 0;/*Number of atoms*/
FLOAT Protein_MaxXYZ[3] = {0.0, 0.0, 0.0};
FLOAT Protein_MaxR = 0.0;

void
Initialize(){
    I_s = C_b;/*unit: mol/L*/
    Beta = E_c / (K_B * T);/*e_c / (K_B * T) unit: v-1 = kg-1 m-2 s3 A 38*/
    Alpha = Beta * E_c / (Epsilon_0 * Am);/*e_c * e_c / (epsilon_0 * K_B * T * Am) unit: Am  7039*/
    Kappa = Sqrt(2.0 * I_s * Alpha * Am * Na * 1000 / Epsilon_s);/*(2*I_s*E_c*E_c/(K_B*T*Epsilon_s*Epsilon_0))^(1/2) unit: Am-1 e-1*/
    a_noline = 2 * C_b * Alpha  * Na * 1000.0 * Am * Am * Am;
    if(!use_aly){
        Read_pqr(fn_pqr);
    }
    if(PB_type == NULL || !strcmp(PB_type, "LPBE")){
        PB_Type = 1;
        PB_type = "LPBE";
        
    }
    else if(!strcmp(PB_type, "NPBE")){
        PB_Type = 2;
    }
    else{
        phgError(-1, "\n Wrong for options PB_type, please type LPBE or NPEB! \n\n ");
    }
}

void 
PrintConstant(){
    phgPrintf("\n****************************CONSTANT*******************************************\n");
    phgPrintf("Number of ion species              %d\n",2);
    phgPrintf("ion[%d].q  %+4.2f e    ion[%d].c  %lf M\n",0, 1.0, 0, (double)C_b);
    phgPrintf("ion[%d].q  %+4.2f e    ion[%d].c  %lf M\n",1, -1.0, 1, (double)C_b);
    phgPrintf("Ionic strength Is                  %lf M\n\n",(double)I_s); 
    phgPrintf("Avogadro constant                  %.4e M-1\n",(double)Na);
	phgPrintf("Elementary charge                  %.4e C\n",(double)E_c);
	phgPrintf("Boltzmann constant                 %.4e J/K\n",(double)K_B);
	phgPrintf("vacuum permittivity                %.4e C^2/(J*K)\n",(double)Epsilon_0);
	phgPrintf("tolerance                          %.4e\n",(double)tolerance_Newton);
	// phgPrintf("Temperature                        %f K\n",(double)T);
	// phgPrintf("biomolecule relative permittivity  %f\n",(double)Epsilon_p);
	// phgPrintf("solvent relative permittivity      %f\n",(double)Epsilon_s);
	phgPrintf("Debye-Huckel                       %.4e m-1\n", (double)Kappa); 
    phgPrintf("\n");
}