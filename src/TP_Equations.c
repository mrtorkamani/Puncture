/* TP_Equations.c */

#include "TwoPunctures.h"

double BY_KKofxyz (double x, double y, double z)
{

  double par_P[3];
  double par_S[3];

  par_P[0] = params_get_real("par_P1");
  par_P[1] = params_get_real("par_P2");
  par_P[2] = params_get_real("par_P3");

  par_S[0] = params_get_real("par_S1");
  par_S[1] = params_get_real("par_S2");
  par_S[2] = params_get_real("par_S3");
  int i, j;
  double r, r2, r3, n_P,
    Aij, AijAij, n[3], n_S[3];

  r2 = x * x + y * y + z * z;
  r = sqrt (r2);
  r3 = r * r2;

  n[0] = x / r;
  n[1] = y / r;
  n[2] = z / r;


  /* dot product: n_P = n.P; */
  n_P = 0;
  for (i = 0; i < 3; i++) {
    n_P += n[i] * par_P[i];
  }

  /* cross product: n_S[i] = [n x S]_i;*/
  n_S[0] = n[1] * par_S[2] - n[2] * par_S[1];
  n_S[1] = n[2] * par_S[0] - n[0] * par_S[2];
  n_S[2] = n[0] * par_S[1] - n[1] * par_S[0];

  AijAij = 0;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Bowen-York-Curvature :*/
      Aij =
	+ 1.5 * (par_P[i] * n[j] + par_P[j] * n[i]
                 + n_P * n[i] * n[j]) / r2
	- 3.0 * (n_S[i] * n[j] + n_S[j] * n[i]) / r3;
      if (i == j)
	Aij -= +1.5 * (n_P / r2);
      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}

void BY_Aijofxyz (double x, double y, double z, double Aij[3][3])
{


  double TP_epsilon = params_get_real("TP_epsilon");
  double TP_Tiny = params_get_real("TP_Tiny");

  double par_P[3];
  double par_S[3];

  par_P[0] = params_get_real("par_P1");
  par_P[1] = params_get_real("par_P2");
  par_P[2] = params_get_real("par_P3");

  par_S[0] = params_get_real("par_S1");
  par_S[1] = params_get_real("par_S2");
  par_S[2] = params_get_real("par_S3");


  int i, j;
  double r, r2, r3, n_P,
    n[3], n_S[3];

  r2 = x * x + y * y + z * z;
  r2 = sqrt (pow (r2, 2) + pow (TP_epsilon, 4));
  if (r2 < pow(TP_Tiny,2))
    r2 = pow(TP_Tiny,2);
  r = sqrt (r2);
  r3 = r * r2;

  n[0] = x / r;
  n[1] = y / r;
  n[2] = z / r;

   /* dot product: n_P = n.P; */
  n_P = 0;
  for (i = 0; i < 3; i++) {
    n_P += n[i] * par_P[i];
  }

  /* cross product: n_S[i] = [n x S]_i;*/
  n_S[0] = n[1] * par_S[2] - n[2] * par_S[1];
  n_S[1] = n[2] * par_S[0] - n[0] * par_S[2];
  n_S[2] = n[0] * par_S[1] - n[1] * par_S[0];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Bowen-York-Curvature :*/
      Aij[i][j] =
        + 1.5 * (par_P[i] * n[j] + par_P[j] * n[i]
		 + n_P * n[i] * n[j]) / r2
	- 3.0 * (n_S[i] * n[j] + n_S[j] * n[i]) / r3;
      if (i == j)
	Aij[i][j] -= +1.5 * (n_P / r2);
    }
  }
}

void NonLinEquations (double rho_adm,
       double A, double B, double T, double R,
       double x, double phi,
       double y, double z, derivs *U, double *values)
{
  double par_m = params_get_real("par_m");

  double r1, psi, psi2, psi4, psi7;
  double mu;

  r1 = sqrt (x * x + y * y + z * z);

  psi = 1. + 0.5 * par_m / r1 + U->d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi * psi2 * psi4;

  values[0] = U->d11[0] + U->d22[0] + U->d33[0] + 0.125 * BY_KKofxyz (x, y, z) / psi7
    + 2.0 * Pi / psi2/psi * rho_adm;
}

void LinEquations (double A, double B, double T, double R,
		   double x, double phi,
		   double y, double z, derivs *dU, derivs *U, double *values)
{
  double par_m = params_get_real("par_m");

  double r1, psi, psi2, psi4, psi8;

  r1 = sqrt (x  * x + y * y + z * z);

  psi = 1. + 0.5 * par_m / r1 + U->d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;

  values[0] = dU->d11[0] + dU->d22[0] + dU->d33[0] - 0.875 * BY_KKofxyz (x, y, z) / psi8 * dU->d0[0];
}
