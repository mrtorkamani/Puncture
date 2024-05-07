#include <stdio.h>
#include <math.h>

#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164	/* Pi/2*/
#define Piq 0.78539816339744830961566084582	/* Pi/4*/

static inline double min (double const x, double const y)
{
  return x<y ? x : y;
}

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
						double A_R, A_RR, B_T, B_TT, tworm, tworm2, tworm4, par_m = 1;
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


int main() {
  double par_b = 2;
  double x = 1000;
  double y = 120;
  double z = 120;
  double xs, ys, zs, rs2, phi, X, R, A, B, aux1, aux2, result, Ui;

  xs = x / par_b;
  ys = y / par_b;
  zs = z / par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;

  A = 2 * tanh (0.5 * X) - 1;
  B = tan (0.5 * R - Piq);

  printf("A= %f \t B= %f \n",A,B );

  AB_To_RT (A,B);

  double x1, x2;
  x1 = par_b * (A*A + 1) / (A*A -1) * 2*B/(1+B*B);
  x2 = par_b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * cos(phi);

  printf("X=%f \t R= %f \n", X,R);

  return 0;
}
