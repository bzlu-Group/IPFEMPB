#include<stdio.h>
#include<string.h>
#include<math.h>

#include "phg.h"
#include "functions_aly.h"

FLOAT xc = 0.0, yc = 0.0, zc = 0.0, rc = 1.0;

void
ls_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if(!use_aly){
        *value = 0.0;
        int i = 0;
	    FLOAT tmp;
        for(i = 0; i < N_m; i++){
		    tmp =   (x - atoms[i].loc[0]) * (x - atoms[i].loc[0]) + 
			        (y - atoms[i].loc[1]) * (y - atoms[i].loc[1]) + 
				    (z - atoms[i].loc[2]) * (z - atoms[i].loc[2]);
		    *value += Exp(-d_Gauss * (tmp - (atoms[i].r * atoms[i].r)));
	    }
	    *value = c_Gauss - *value;
        return;
    }
    else{
        *value = (x - xc) * (x - xc) + (y -yc) * (y-yc) + (z - zc) * (z-zc) - rc * rc;
        return;
    }
}

void
ls_grad_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    if(!use_aly){
    values[0] = 0;
	values[1] = 0;
	values[2] = 0;
	
    int i = 0;
	FLOAT tmp1, tmp2;
    for(i=0; i<N_m; i++){
		tmp1 = (x-atoms[i].loc[0]) * (x-atoms[i].loc[0]) + 
			    (y-atoms[i].loc[1]) * (y-atoms[i].loc[1]) + 
				(z-atoms[i].loc[2]) * (z-atoms[i].loc[2]);
		tmp2 = Exp(-1.0 * d_Gauss * (tmp1 - (atoms[i].r * atoms[i].r)));
		values[0] += 2 * d_Gauss * (x - atoms[i].loc[0]) * tmp2;
		values[1] += 2 * d_Gauss * (y - atoms[i].loc[1]) * tmp2;
		values[2] += 2 * d_Gauss * (z - atoms[i].loc[2]) * tmp2;
	}
	return;
    }
    else{
        x -= xc; y -= yc; z -= zc;
        values[0] = x + x;
        values[1] = y + y;
        values[2] = z + z;
    }
}

void
func_u1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    assert(use_aly);
    if(PB_Type == 1){
        /*LPBE*/
        FLOAT R = rc * Am;
        FLOAT q = E_c;
        FLOAT phi = q * (1 / (Epsilon_s * (1 + Kappa * R)) - 1 / Epsilon_p) / (4 * PAI *Epsilon_0 * R);
        *value  = phi * Beta;
        return; 
    }
    else{
        /*NPBE*/
        *value = x * x + x * y + z;
        return;
    }
}

void
func_u2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    assert(use_aly);
    if(PB_Type == 1){
        /*LPBE*/
         x -= xc; y -= yc; z -= zc;
        FLOAT r = Sqrt(x * x + y * y + z * z) * Am;
        if(r < Am * Am){
            *value = 0.0;
            return;
        }
        FLOAT R = rc * Am;
        FLOAT q = E_c;
        FLOAT phi = q *  Pow(E, Kappa * R) * Pow(E, -Kappa * r) / (4 * PAI * Epsilon_0 * Epsilon_s * r * (1 + Kappa * R));
        *value  = phi * Beta;
        return; 
    }
    else{
        /*NPBE*/
        *value = Sin(x + 2*y + 3*z);
        return;
    }
}

static void
/*func_G*/
func_us(FLOAT x, FLOAT y, FLOAT z, FLOAT *value){

    assert(N_m > 0);

    INT i = 0;
    FLOAT sum = 0.0;
	FLOAT r;
    for(i = 0; i < N_m; i++){
		r =   (x-atoms[i].loc[0]) * (x-atoms[i].loc[0]) + 
			  (y-atoms[i].loc[1]) * (y-atoms[i].loc[1]) + 
			  (z-atoms[i].loc[2]) * (z-atoms[i].loc[2]);
        if(r < 1e-6){
            r = 1e-6;
        }
        sum = sum + atoms[i].zi / Sqrt(r);
    }
    *value = Alpha * sum / (4 * PAI * Epsilon_p);
    return;
}


void
func_jD(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{   
    if(!use_aly){
    func_us(x,y,z,value);
    *value = - *value;
    return;
    }
    else{
        FLOAT u2;
        func_u1(x, y, z, value);
        func_u2(x, y, z, &u2);
        *value -= u2;
        return;
    }
}

static void
func_grad_u1(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{   
    assert(use_aly);
    if(PB_Type == 2){
        /*NPBE*/
        values[0] = 2 * x + y;
        values[1] = x;
        values[2] = 1;
        return;
    }
    if(PB_Type == 1){
        /*LPBE*/
        values[0] = 0.0;
        values[1] = 0.0;
        values[2] = 0.0;
        return;
    }
}

static void
func_grad_u2(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{   
    assert(use_aly);
    if(PB_Type == 2){
        /*NPBE*/
        values[0] = Cos(x + 2*y + 3*z);
        values[1] = 2 * Cos(x + 2*y + 3*z);
        values[2] = 3 * Cos(x + 2*y + 3*z);
        return;
    }
    if(PB_Type == 1){
        /*LPBE*/
        x -= xc; y -= yc; z -= zc;
        FLOAT r = Sqrt(x * x + y * y + z * z) * Am;
        if(fabs(r) < Am * Am){
            values[0] = 0.0;
            values[1] = 0.0;
            values[2] = 0.0;
            //phgWarning("grad_u2: r= %e\n", (double)r);
            return;
        }
        FLOAT R = rc * Am;
        FLOAT q = E_c;
        FLOAT coe = q * Exp(Kappa * R) * (Kappa * r + 1) * (Exp(-Kappa * r))/ (Pow(r,3) * 4 * PAI * Epsilon_0 * Epsilon_s * (1 + Kappa * R));
        values[0] = coe * x * Beta * Am * Am;
        values[1] = coe * y * Beta * Am * Am;
        values[2] = coe * z * Beta * Am * Am;
        return;
    }
}

static void
/*func_G_grad*/
func_grad_us(FLOAT x, FLOAT y, FLOAT z, FLOAT *values){
    assert(N_m > 0);
    INT i = 1;
    FLOAT r;
    FLOAT sums[3] = {0.0, 0.0 ,0.0};
    for(i = 0; i < N_m; i++){
        r = Sqrt((x-atoms[i].loc[0]) * (x-atoms[i].loc[0]) + \
			    (y-atoms[i].loc[1]) * (y-atoms[i].loc[1]) + \
				(z-atoms[i].loc[2]) * (z-atoms[i].loc[2]));
        if(r <= 1e-6){
            r = 1e-6;
        }
        sums[0] += -atoms[i].zi *(x-atoms[i].loc[0]) / (r * r * r);
        sums[1] += -atoms[i].zi *(y-atoms[i].loc[1]) / (r * r * r);
        sums[2] += -atoms[i].zi *(z-atoms[i].loc[2]) / (r * r * r);
    }

    values[0] = Alpha * sums[0] / (4 * PAI * Epsilon_p);
    values[1] = Alpha * sums[1] / (4 * PAI * Epsilon_p);
    values[2] = Alpha * sums[2] / (4 * PAI * Epsilon_p);
    return;
}

void
func_jN(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    if(!use_aly){
    FLOAT tmp_values[3];
    func_grad_us(x,y,z,tmp_values);
    values[0] = -Epsilon_p * tmp_values[0];
    values[1] = -Epsilon_p * tmp_values[1];
    values[2] = -Epsilon_p * tmp_values[2];
    }
    else{
        FLOAT un2[Dim] = {0.0, 0.0, 0.0};
        func_grad_u1(x, y, z, values);
        func_grad_u2(x, y, z, un2);
        values[0] = Epsilon_p * values[0] - Epsilon_s * un2[0];
        values[1] = Epsilon_p * values[1] - Epsilon_s * un2[1];
        values[2] = Epsilon_p * values[2] - Epsilon_s * un2[2];
        return;
    }
}

void
func_g1D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{   
    *value = 0.0;
    return;
}

void
func_g2D(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{   
    if(!use_aly){
    INT i = 1;
    FLOAT r;
    FLOAT sum = 0.0;
    for(i = 0; i < N_m; i++){
        r = Sqrt((x-atoms[i].loc[0]) * (x-atoms[i].loc[0]) + 
			    (y-atoms[i].loc[1]) * (y-atoms[i].loc[1]) + 
				(z-atoms[i].loc[2]) * (z-atoms[i].loc[2]));
        sum += atoms[i].zi * Exp(-Kappa * r * Am) / r;
    }
    *value = Alpha * sum / (4 * PAI * Epsilon_s);
    return;
    }
    else{
        func_u2(x,y,z, value);
        return;
    }
}

void
func_f1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if(!use_aly){
        *value = 0.0;
        return;
    }
    else{
        if(PB_Type == 1){
            *value = 0.0;
            return;
        }
        else if(PB_Type == 2){
            *value = - 2.0 * Epsilon_p;
            return;
        }
    }
}

void
func_f2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if(!use_aly){
        *value = 0.0;
        return;
    }
    else{
        if(PB_Type == 1){
            *value = 0.0;
            return;
        }
        else if(PB_Type == 2){
            FLOAT tmp_value;
            func_u2(x,y,z,&tmp_value);
            *value = 14.0 * Sin(x + 2 * y + 3 * z) * Epsilon_s + a_noline * sinh(tmp_value);
            return;
        }
    }
}
