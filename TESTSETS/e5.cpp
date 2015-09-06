#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"

int main ()
{
	ifstream in("params.stiff");
	ifstream in2("params.nonstiff");

	StiffIntegratorT stiffT(in);
	
	cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
	
	stiffT.Integrate();
	
	NonStiffIntegratorT nonstiffT(in2);

	cout << "\n\n*******Problem integrated with DOPRI5*******\n\n";
	
	nonstiffT.Integrate();

	return (0);
}

void fEval(double x, double *y, double *f)
{

	double prod1, prod2, prod3, prod4;

	prod1 = 7.89e-10*y[0];
	prod2 = 1.1e7*y[0]*y[2];
	prod3 = 1.13e9*y[1]*y[2];
	prod4 = 1.13e3*y[3];
	f[0] = -prod1 - prod2;
	f[1] = prod1 - prod3;
	f[3] = prod2 - prod4;
	f[2] = f[1] - f[3];

	return;

} // fEval

void Jacobian(double x, double *y, double **J)
{
 	double a = 7.89e-10, b = 1.1e7, cm = 1.13e9, c = 1.13e3;

	J[0][0] = -a - b*y[2];
	J[0][1] = 0.0;
	J[0][2] = -b*y[0];
	J[0][3] = 0.0;
	J[1][0] = a;
	J[1][1] = -cm*y[2];
	J[1][2] = -cm*y[1];
	J[1][3] = 0.0;
	J[2][0] = a - b*y[2];
	J[2][1] = -cm*y[2];
	J[2][2] = -b*y[0] - cm*y[1];
	J[2][3] = c;
	J[3][0] = b*y[2];
	J[3][1] = 0.0;
	J[3][2] = b*y[0];
	J[3][3] = -c;

	return;

} // Jacobian

void Mass(double **M)
{


} // Mass
