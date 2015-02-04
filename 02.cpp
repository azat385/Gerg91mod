#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//const double counter_global = 4294967295;
const double  CONST_R = 8.31451, Tc = 293.15,Pc = 0.101325; // молярная масса сухого воздуха стандартного состава
// Переменные

double
KK,
ROwork,
k,
U,
mu,
Qmax,
Qmin,
m_P,
m_T,
kappa,
P,
T,
ROc,
Xa,
Xy;

/*double stepppp(double x , double n) 
{
       int n_int;
       if (fabs(2.0 - n) < 1E-9) 
          return  x * x;
       if (fabs(3.0 - n) < 1E-9) 
          return x * x * x;
       if (x > 0) 
          return exp(log(x) * n);
       n_int = (int)round(n);
       if (fabs(n_int - n) > 1E-9) 
	   { 	
		printf("#ERR# pow()");
        return 0;
	    }
        if (!(n_int && 1)) 
           return exp(log(fabs(x)) * n);
        else 
         return (-1) * exp(log(fabs(x)) * n); 
}
*/


double gerg91m(double ROc , double Xa , double Xy , double P , double T)
{
double 
Bm = 0,
Cm= 0,
Xe= 0,
B1= 0,
B2= 0,
B23= 0,
B3= 0,
C1= 0,
C2= 0,
C3= 0,
C223= 0,
C233= 0,
Bb= 0,
Cc= 0,
H= 0,
Me= 0,
A0= 0,
A1= 0,
A2= 0,
B0= 0,
C0= 0,
b= 0,
Zc= 0,
Z= 0,
T2= 0,
Xe2= 0.0,
Xe3= 0,
Xa2= 0,
Xa3= 0,
Xy2= 0,
Xy3= 0,
p3= 0,
tmp;
	
	Xe = 1 - Xa - Xy;
	Zc = 1 - pow((0.0741*ROc-0.006-0.063*Xa-0.0575*Xy), 2);
	Me = ((double)24.05525*Zc*ROc-28.0135*Xa-44.01*Xy)/Xe;

	H = 128.64 + 47.479*Me;

	T2 = pow(T, 2);
	Xe2 = pow(Xe, 2);
	Xe3 = pow(Xe, 3);
	Xa2 = pow(Xa, 2);
	Xa3 = pow(Xa, 3);
	Xy2 = pow(Xy, 2);
	Xy3 = pow(Xy, 3);

	B1 = -0.425468+2.865E-3*T-4.62073E-6*T2+(8.77118E-4-5.56281E-6*T+8.8151E-9*T2)*H+
		(-8.24747E-7+4.31436E-9*T-6.08319E-12*T2)*pow(H,2);

	B2 = -0.1446+7.4091E-4*T-9.1195E-7*T2;
	B23 = -0.339693+1.61176E-3*T-2.04429E-6*T2;
	B3 = -0.86834+4.0376E-3*T-5.1657E-6*T2;

	Bb = 0.72 + 1.875E-5 * pow((320 - T), 2);

	Bm = Xe2*B1+Xe*Xa*Bb*(B1+B2)-1.73*Xe*Xy*pow(B1*B3,0.5)+Xa2*B2+2*Xa*Xy*B23+Xy2*B3;

	C1 = -0.302488+1.95861E-3*T-3.16302E-6*T2+(6.46422E-4-4.22876E-6*T+6.88157E-9*T2)*H+
		(-3.32805E-7+2.2316E-9*T-3.67713E-12*T2)*pow(H,2);
	C2 = 7.8498E-3-3.9895E-5*T+6.1187E-8*T2;
	C3 = 2.0513E-3+3.4888E-5*T-8.3703E-8*T2;
	C223 = 5.52066E-3-1.68609E-5*T+1.57169E-8*T2;
	C233 = 3.58783E-3+8.06674E-6*T-3.25798E-8*T2;
	Cc = 0.92 + 0.0013 * (T - 270);

	p3 = 1/3;
	Cm = Xe3*C1+3*Xe2*Xa*Cc*pow(C1*C1*C2,p3)+2.76*Xe2*Xy*pow(C1*C1*C3,p3) +
		3*Xe*Xa2*Cc*pow(C1*C2*C2,p3)+6.6*Xe*Xa*Xy*pow(C1*C2*C3,p3)+
		2.76*Xe*Xy2*pow(C1*C3*C3,p3) + Xa3*C2+3*Xa2*Xy*C223+3*Xa*Xy2*C233+Xy3*C3;

	b = 1E3*P/(2.7715*T);
	B0 = b*Bm;
	C0 = b*b*Cm;
	A0 = 1+1.5*(B0+C0);
	A1 = 1+B0;
	A2 = pow((A0-pow((A0*A0-A1*A1*A1),0.5)),p3);

	Z = (1+A2+A1/A2)/3;
	printf("FUNCTION:input: ROc=%f\t Xy=%f\t Xa=%f\t T=%f\t P=%f\t \n",ROc, Xy,Xa,T,P);
	return Z/Zc;
}


double calc_ROwork(void)
{
	return (ROc*P*Tc)/(Pc*T*KK);
}


double calc_kappa(void)
{ 
	return 1.556*(1+0.074*Xa)-3.9E-4*T*(1-0.68*Xa)-0.208*ROc+pow(P/T,1.43)*
		(384*(1-Xa)*pow(P/T,0.8)+26.4*Xa);
}


double calc_U(void)
{ 
	return 18.591*pow(T*k*KK/ROc,0.5);
}


double calc_mu(void)
{
double
Mt,
Cm,
Pp,
Tp,
Ppk,
Tpk;
	Mt = 3.24*(pow(T,0.5)+1.37-9.09*pow(ROc,0.125))/(pow(ROc,0.5)+2.08-1.5*(Xa+Xy));

	if(P < 0.5) return Mt;
	
	Ppk = 2.9585*(1.608-0.05994*ROc+Xy-0.392*Xa);
	Tpk = 88.25*(0.9915+1.759*ROc-Xy-1.681*Xa);
	
	Pp = P/Ppk;
	Tp = T/Tpk;
	
	Cm = 1+(Pp*Pp)/(30*(Tp-1));
	
	return Mt*Cm;
}


double calc_Qmax(void)
{
	return 92.819*(0.51447*ROc+0.05603-0.65689*Xa-Xy);
}


double calc_Qmin(void)
{
	return  85.453*(0.52190*ROc+0.04242-0.65197*Xa-Xy);
}



int main(void) 
{

	int t1c = 0,
		t2c = 0,
		t3c = 0,
		t1p = 0,
		t1p2 = 0,
		t1p3 = 0,
		i = 0,
		j = 0,
		k = 0,
		count1 = 0,
		count2 = 0,
		count3 = 0,
		x = 0;
	double p = 0,
				temp = 0, 
				flow_c_z = 0,
				flow_c_1zp = 0,
				flow_c_z1 = 0,
				flow_c_z2= 0,
				temp1 = 0,
				temp2 = 0,
				temp3 = 0,
				tempr = 0,
				tempr1 = 0,
				tempr2 = 0,
				tempr3 = 0,
				ps = 0,
				ps2 = 0,
				ps3 = 0,
				c = 0, 
				g = 0,
				h1sc = 0,
				h2sc = 0,
				h1s2c = 0,
				h2s2c = 0,
				h1s3c = 0,
				h2s3c = 0,		 
				h1s = 0,
				h2s = 0,
				h1s2 = 0,
				h2s2 = 0,
				h1s3 = 0,
				h2s3 = 0,
				dens = 0,
				dens1 = 0,
				dens2 = 0,
				dens3 = 0,
				dens1o = 0,
				dens2o = 0,
				dens3o = 0,
				flow_u_1= 0,
				flow_u_1s = 0,
				flow_c_1 = 0,
				flow_c_1s = 0,
				flow_u_2 = 0,
				flow_u_2s =0,
				flow_u_1s2 = 0,
				flow_c_2 = 0,
				flow_c_2s = 0,
				K = 0,
				dens_std= 0,
				dens_work = 0,
				K_Factor_1=3600,
				K_Factor_2=3600,
				data[9],
	RO,
	Xy = 0.0097,
	Xa = 0.0013;
	
	

	/*while (TRUE)  {
		request_resource(IO_SYSTEM);
		key = dbase(MODBUS, 1002);
		p = read_val(0, 40011);
		temp = read_val(0, 40013);
		key2 = dbase(MODBUS, 1003);
		release_resource(IO_SYSTEM);
		if (key2 == TRUE)  {
			K_Factor_1 =  read_val(0,  40017);
			K_Factor_2 =  read_val(0,  40019);
			request_resource(IO_SYSTEM);
			databaseWrite(MODBUS, 1003, FALSE);
			release_resource(IO_SYSTEM);
		}


	
	    P = 5;
		T = 5 + 273.15;
		dens_std = 0.78;
*/        
	    P = 5.65;
		T = 7.0 + 273.15;
		dens_std = 0.7781;	
		
		ROc = dens_std;
		KK = gerg91m(ROc, Xa, Xy, P, T);
	    K = KK;
	    //write_val(0, 40049, K);
		dens_work = calc_ROwork();
		
		printf("K=%f  dens_work=%f\n ",K,dens_work);
		//write_val(0, 40051, dens_work);
/*		flow_c_z = read_val(0, 30011);
		if (flow_c_z > flow_c_1zp)  flow_c_1 =  flow_c_z2 + (flow_c_z)/K_Factor_1; 
			else  {
				flow_c_z2 = flow_c_z2 + (counter_global + flow_c_z)/K_Factor_1;
				flow_c_1 =  flow_c_z2;
			}
		flow_c_1zp = flow_c_z;
		flow_u_1 = read_val(0, 49011);
		flow_u_1s = (flow_u_1*p*Tc)/(Pc*(temp + 273.15)*K);
		flow_c_1s += flow_u_1s;
		write_val(0, 40115, flow_u_1*3600);
		write_val(0, 40117, flow_u_1s*3600);
		write_val(0, 40123, flow_c_1);
		write_val(0, 40125, flow_c_1s);
*/        
	/* k = calc_kappa;

	U = calc_U;
	mu = calc_mu;
	Qmax = calc_Qmax;
	Qmin = calc_Qmin;
 */
	//kappa = k;
/*	if (t3c != t1p3 ) 
				{
				data[0] = ps3/k;
				data[1] = temp3/k;
				data[2] = h1s3c;
				data[3] = h1s3;
				data[4] = h2s3c;
				data[5] = h2s3;
				data[6] = dens3/k;
				data[7] = dens3o/k;
				data[8] = tempr3/j;
				if(count3 > 12) {
						for(x = 0; x < 11; x++) {
							c =  read_val(0, 40131 + 21*x);
							write_val(0, 40131 + 21*x + 21, c);
						}
						count3 = 12;
				}
				request_resource(IO_SYSTEM);
				databaseWrite(MODBUS, count3*20 + 40131, now.month);
				databaseWrite(MODBUS, count3*20 + 1 + 40131, now.year);
				release_resource(IO_SYSTEM);
				for(x = 0 ; x < 9; x++) {
					write_val(0, 40131 + count3*20 + 2 + x, data[x]);
				}
				count3++;
				ps3 = 0;
				temp3 = 0;
				h1s3c = 0;
				h2s3c = 0;
				h1s3 = 0;
				h2s3 = 0;
				ps = 0;
				temp1 = 0;
				h1sc = 0;
				h2sc = 0;
				h1s = 0;
				h2s = 0;
				ps2 = 0;
				temp2 = 0;
				h1s2c = 0;
				h2s2c = 0;
				h1s2 = 0;
				h2s2 = 0;
				dens1 = 0;
				dens1o = 0;
				dens2 = 0;
				dens2o = 0;
				dens3 = 0;
				dens3o = 0;
				tempr1 = 0;
				tempr2 = 0;
				tempr3 = 0;
				i = 0;
				j = 0;
				t1p3 = t3c;
				t1p2 = t2c;
				t1p = t1c;
				k = 0;
				request_resource(IO_SYSTEM);
				databaseWrite(MODBUS, 1006, TRUE);
				release_resource(IO_SYSTEM);
				}
			ps3 += p;
			temp3 += temp;
			h1s3c += flow_u_1s;
			h2s3c += flow_u_2s;
			h1s3 += flow_u_1;
			h2s3 += flow_u_2;
			dens3 += dens_work;
			dens3o += dens_std;
			tempr3 += tempr;
			k++;
			g = temp3/k;
			write_val(0, 40087, g);
			g = ps3/k;
			write_val(0, 40089, g);
			g = h1s3c;
			write_val(0, 40091, g);
			g = h2s3c;
			write_val(0, 40093, g);
			g = h1s3;
			write_val(0, 40095, g);
			g = h2s3;
			write_val(0, 40097, g);
			g = dens3/k;
			write_val(0, 40099, g);
			g = dens3o/k;
			write_val(0, 40101, g);
			if (t2c != t1p2) 
				{
				data[0] = ps2/j;
				data[1] = temp2/j;
				data[2] = h1s2c;
				data[3] = h1s2;
				data[4] = h2s2c;
				data[5] = h2s2;
				data[6] = dens2/j;
				data[7] = dens2o/j;
				data[8] = tempr2/j;
				if(count2 > 30) {
						for(x = 0; x < 29; x++) {
							c =  read_val(0, 40381 + 21*x);
							write_val(0, 40381 + 21*x + 21, c);
						}
						count2 = 30;
				}
				request_resource(IO_SYSTEM);
				databaseWrite(MODBUS, count2*21  +  40381, now.day);
				databaseWrite(MODBUS, count2*21 + 1 + 40381, now.month);
				databaseWrite(MODBUS, count2*21 + 2 + 40381, now.year);
				release_resource(IO_SYSTEM);
				for(x = 0 ; x < 9; x++) {
					write_val(0, 40381 + count2*21 + 3 + x, data[x]);
				}
				count2++;
				ps2 = 0;
				temp2= 0;
				h1s2c = 0;
				h2s2c = 0;
				h1s2 = 0;
				h2s2 = 0;
				ps = 0;
				temp1 = 0;
				h1s = 0;
				h2s = 0;
				h1sc = 0;
				h2sc = 0;
				dens1 = 0;
				dens1o = 0;
				dens2 = 0;
				dens2o = 0;
				tempr1 = 0;
				tempr2 = 0;
				i = 0;
				t1p2 = t2c;
				t1p = t1c;
				j = 0;
				request_resource(IO_SYSTEM);
				databaseWrite(MODBUS, 1005, TRUE);
				release_resource(IO_SYSTEM);
				}
			ps2 += p;
			temp2 += temp;
			h1s2c += flow_u_1s;
			h2s2c += flow_u_2s;
			h1s2 += flow_u_1;
			h2s2 += flow_u_2;
			dens2 += dens_work;
			dens2o += dens_std;
			tempr2 += tempr;
			j++;
			g = ps2/j;
			write_val(0, 40071, g);
			g = temp2/j;
			write_val(0, 40073, g);
			g = h1s2c;
			write_val(0, 40075, g);
			g = h2s2c;
			write_val(0, 40077, g);
			g = h1s2;
			write_val(0, 40079, g);
			g = h2s2;
			write_val(0, 40081, g);
			g = dens2/j;
			write_val(0, 40083, g);
			g = dens2o/j;
			write_val(0, 40085, g);
			if (t1c != t1p) 
				{
				data[0] = ps/i;
				data[1] = temp1/i;
				data[2] = h1s;
				data[3] = h1sc;
				data[4] = h2s;
				data[5] = h2sc;
				data[6] = dens1/i;
				data[7] = dens1o/i;
				data[8] = tempr1/j;
				if(count1 > 24) {
						for(x = 0; x < 23; x++) {
							c =  read_val(0, 41031 + 22*x);
							write_val(0, 41031 + 22*x + 22, c);
						}
						count1 = 24;
				}
				request_resource(IO_SYSTEM);
				databaseWrite(MODBUS, count1*22 + 41031, now.hour);
				databaseWrite(MODBUS, count1*22 + 1 + 41031, now.day);
				databaseWrite(MODBUS, count1*22 + 2 + 41031, now.month);
				databaseWrite(MODBUS, count1*22+ 3 + 41031, now.year);
				release_resource(IO_SYSTEM);
				for(x = 0 ; x < 9; x++) {
					write_val(0, 41031 + count1*22 + 4 + x, data[x]);
				}
				count1++;
				ps = 0;
				temp1 = 0;
				h1s = 0;
				h2s = 0;
				h1sc = 0;
				h2sc = 0;
				dens1 = 0;
				dens1o = 0;
				t1p = t1c;
				tempr1 = 0;
				i = 0;
				request_resource(IO_SYSTEM);
				databaseWrite(MODBUS, 1004, TRUE);
				release_resource(IO_SYSTEM);
				}
			ps += p;
			temp1 += temp;
			h1sc += flow_u_1s;
			h2sc += flow_u_2s;
			h1s += flow_u_1;
			h2s += flow_u_2;
			dens1 += dens_work;
			dens1o += dens_std;
			tempr1 += tempr;
			i++;	
			g = ps/i;
			write_val(0, 40055, g);
			g = temp1/i;
			write_val(0, 40057, g);
			g = h1sc;
			write_val(0, 40059, g);
			g = h2sc;
			write_val(0, 40061, g);
			g = h1s;
			write_val(0, 40063, g);
			g = h2s;
			write_val(0, 40065, g);
			g = dens1/i;
			write_val(0, 40067, g);
			g = dens1o/i;
			write_val(0, 40069, g);			 
			
            */
			
			
            return 0;			
}
