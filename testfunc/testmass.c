#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164	/* Pi/2*/
#define Piq 0.78539816339744830961566084582	/* Pi/4*/

static inline double min (double const x, double const y)
{
  return x<y ? x : y;
}


int main() {
  /*
  double par_b = .00001;
  double x = 10.0;
  double y = 12.0;
  double z = 12.0;
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

  double x1, x2, x3;
  x1 = par_b * (A*A + 1) / (A*A -1) * 2*B/(1+B*B);
  x2 = par_b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * cos(phi);
  x3 = par_b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * sin(phi);

  printf("x=%f \t y= %f \t z= %f \n", x1,x2,x3);
//////////////////////////////////////////////////////////////////////////////////////////////////////////
  double X, R, A=.5, B=.5, phi=1;
  double At = 0.5 * (A + 1);

  X = 2 * atanh (At);
  R = Pih + 2 * atan (B);

  double C_c2, U_cb, U_CB, x, r,y,z, par_b = 1.0;
  gsl_complex C,C1, C_c, C_cc, c, c_C, c_CC, U_c, U_cc, U_C, U_CC;
  int ivar;

  C = gsl_complex_rect (X,R);

  c = gsl_complex_mul_real (gsl_complex_cosh (C), par_b);	/ c=b*cosh(C)/
  C1 = gsl_complex_inverse (C);

  //c = gsl_complex_mul_real(gsl_complex_add(C , C1),par_b/2);

  x = GSL_REAL(c);
  r = GSL_IMAG(c);

  y = r * cos (phi);
  z = r * sin (phi);

  printf("x =%f y =%f z =%f \n",x,y,z);

  double x1, x2, x3;
  x1 = par_b * (A*A + 1) / (A*A -1) * 2*B/(1+B*B);
  x2 = par_b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * cos(phi);
  x3 = par_b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * sin(phi);

  printf("x=%f \t y= %f \t z= %f \n", x1,x2,x3);
  */

  double eta, eps,r ,x,y,z;
  double A=0,B=1,phi=1,b=1;
  gsl_complex kesi,c;
  eta = 2 * atanh (A);
  eps = Pih + 2* atan (B);

  kesi = gsl_complex_rect (eps,eta);
  c = gsl_complex_mul_real(gsl_complex_cosh (kesi),b);

  //x = GSL_REAL(c);
  //r = GSL_IMAG(c);
  x = b* cosh (eps) * cos (eta);
  r = b* sinh (eps) * sin (eta);
  y = r * cos (phi);
  z = r * sin (phi);
  printf("x =%f y =%f z =%f \n",x,y,z);

  double x1, x2, x3;
  x1 = b * (A*A + 1) / (A*A -1) * 2*B/(1+B*B);
  x2 = b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * cos(phi);
  x3 = b * 2*A/(1-A*A) * (1-B*B)/(1+B*B) * sin(phi);

  printf("x=%f \t y= %f \t z= %f \n", x1,x2,x3);
  return 0;
}
