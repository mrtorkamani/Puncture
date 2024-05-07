/* TP_CoordTransf.c */

#include "TwoPunctures.h"


void AB_To_RT(int nvar, double A, double B, double *R, double *T,
	       derivs *U)
				 /* R: radial coordinate                    T: \theta coordinate             */
				 /*On Entrance: U.d0[]=U[]; U.d1[] =U[]_A;  U.d2[] =U[]_B;  U.d3[] =U[]_3;   */
				 /*                          U.d11[]=U[]_AA; U.d12[]=U[]_AB; U.d13[]=U[]_A3; */
				 /*                          U.d22[]=U[]_BB; U.d23[]=U[]_B3; U.d33[]=U[]_33; */
				 /* At Exit:     U.d0[]=U[]; U.d1[] =U[]_R;  U.d2[] =U[]_T;  U.d3[] =U[]_3;  */
				 /*                          U.d11[]=U[]_RR; U.d12[]=U[]_RT; U.d13[]=U[]_r3; */
				 /*                          U.d22[]=U[]_TT; U.d23[]=U[]_T3; U.d33[]=U[]_33; */
{
						double A_R, A_RR, B_T, B_TT, tworm, tworm2, tworm4, par_m;
						int ivar;

						par_m = params_get_real("par_m");

					 *R = par_m/2 *1/(1-1/A);
					 *T = Pih * (B-1);

					 tworm =(2*(*R)-par_m);
					 tworm2 = tworm*tworm;
					 tworm4 = tworm2*tworm2;
					 A_R = - 2 * par_m/tworm2;
					 A_RR = 8*tworm *par_m/tworm4;
					 B_T = 1/Pih;
					 for (ivar = 0; ivar < nvar; ivar++) {
				     U->d11[ivar] = A_R * A_R * U->d11[ivar] + A_RR * U->d1[ivar];
				     U->d12[ivar] = A_R * B_T * U->d12[ivar];
				     U->d13[ivar] = A_R * U->d13[ivar];
				     U->d22[ivar] = B_T * B_T * U->d22[ivar];
				     U->d23[ivar] = B_T * U->d23[ivar];
				     U->d1[ivar] = A_R * U->d1[ivar];
				     U->d2[ivar] = B_T * U->d2[ivar];
				   }
}

void RT3_To_xyz(int nvar, double R, double T, double phi,
		 double *x ,double *y, double *z, derivs *U)
{
			 int jvar;
		   double
		     sin_phi = sin (phi),
		     cos_phi = cos (phi),
		     sin2_phi = sin_phi * sin_phi,
		     cos2_phi = cos_phi * cos_phi,
		     sin_2phi = 2 * sin_phi * cos_phi,
		     cos_2phi = cos2_phi - sin2_phi, r_inv = 1 / R, r_inv2 = r_inv * r_inv,
				 sin_T = sin (T),
				 cos_T = cos (T),
				 sin2_T= sin_T * sin_T,
				 cos2_T= sin_T * sin_T,
				 sin_2T = 2 * sin_T * cos_T,
				 cos_2T=cos2_T - sin2_T;

				 *x = R * cos_T;
				 *y = R * sin_T * cos_phi;
				 *z = R * sin_T * sin_phi;
				 for (jvar = 0; jvar < nvar; jvar++) {

			     double U_R = U->d1[jvar], U_T = U->d2[jvar], U_3 = U->d3[jvar],
			       U_RR = U->d11[jvar], U_RT = U->d12[jvar], U_R3 = U->d13[jvar],
			       U_TT = U->d22[jvar], U_T3 = U->d23[jvar], U_33 = U->d33[jvar];

			     U->d1[jvar] = U_R * cos_T - U_T * r_inv * sin_T;		/* U_x*/
			     U->d2[jvar] = U_R * sin_T * cos_phi + U_T* r_inv * cos_T * cos_phi - U_3 * r_inv * sin_phi / sin_T;	/* U_y*/
			     U->d3[jvar] = U_R * sin_T * sin_phi + U_T * r_inv * cos_T * sin_phi + U_3 * r_inv * cos_phi / sin_T;	/* U_z*/
			     U->d11[jvar] = U_RR *cos2_T + U_R * 2 * sin2_T * r_inv - U_RT * sin_2T * r_inv + U_T * sin_2T * r_inv2 + U_TT * sin2_T * r_inv2;		/* U_xx*/
			     U->d12[jvar] = U_RR * sin_T * cos_T * cos_phi + U_RT * r_inv * cos_2T * cos_phi - U_R * sin_T * cos_T * cos_phi * r_inv    /* U_xy*/
					 	 - U_TT * cos_T * sin_T * cos_phi * r_inv2 + U_T3 * sin_phi * r_inv - U_R3 * cos_T * sin_phi * r_inv/sin_T
						 - U_T * cos_2T * cos_phi * r_inv2;
			     U->d13[jvar] = U_RR * sin_T * cos_T * sin_phi + U_RT * r_inv * cos_2T * sin_phi - U_R * sin_T * cos_T * sin_phi * r_inv    /* U_xz*/
					 	 - U_TT * cos_T * sin_T * sin_phi * r_inv2 - U_T3 * cos_phi * r_inv2 + U_R3 * cos_T * cos_phi * r_inv/sin_T
						 - U_T * cos_2T * sin_phi * r_inv2;
			     U->d22[jvar] = U_RR * sin_2T * cos2_phi + U_RT * sin_2T * cos2_phi *r_inv + U_R * (cos2_phi * cos2_T + sin2_phi) * r_inv /* U_yy */
					   + U_TT * cos2_phi * cos2_T *r_inv2 - U_T3 * sin_2phi * cos_T * r_inv2/sin_T + U_33 * sin2_phi * r_inv2/sin2_T
						 + U_3 * sin_2phi * r_inv2/sin2_T - U_R3 * sin_2phi * r_inv - U_T * ((cos_T * cos2_phi) *(2 * sin2_T +1)-cos_T)*r_inv2/sin_T;
			     U->d23[jvar] = U_RR * 0.5 * sin_2phi * cos2_T + U_RT * 0.5 * sin_2T * sin_2phi * r_inv + U_R * 0.5 * sin_2phi * sin2_T * r_inv
					 	 + U_TT * 0.5 * sin_2phi * cos2_T * r_inv2 + U_T3 * cos_T * cos_2phi * r_inv2/sin_T - U_3 * cos_2phi * r_inv2/sin2_T
						 + U_R3 * cos_2phi * r_inv - U_33 * 0.5 * sin_2phi * r_inv2/sin2_T - U_T * 0.5 * sin_2phi * cos_T * (2 * sin2_T +1) * r_inv2/sin_T;
			     U->d33[jvar] = U_RR * sin2_T * sin2_phi + U_RT * sin2_phi * sin_2T * r_inv + U_R * sin2_phi * cos2_T * r_inv  /* U_zz */
					 	 + U_TT * sin2_phi * cos2_T * r_inv2 + U_T3 * sin_2phi * cos_T * r_inv2/sin_T + U_33 * cos2_phi * r_inv2/sin2_T
						 - U_3 * sin_2phi * r_inv2/sin2_T + U_R3 * sin_2phi * r_inv + U_T * (cos2_phi * cos_T - 2 * sin2_phi * sin2_T * cos_T) * r_inv2/sin_T;
			   }
}
