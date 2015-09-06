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
	f[0] = -0.04*y[0] + 1.0e4*y[1]*y[2];
	f[2] = 3.0e7*y[1]*y[1];
	f[1] = -f[0] - f[2];

	return;
} // fEval

void Jacobian(double x, double *y, double **J)
{
	double prod1, prod2, prod3;

	prod1 = 1.0e4*y[1];
	prod2 = 1.0e4*y[2];
	prod3 = 6.0e7*y[1];

	J[0][0] = -0.04;
	J[0][1] = prod2;
	J[0][2] = prod1;

	J[1][0] = 0.04;
	J[1][1] = -prod2 - prod3;
	J[1][2] = -prod1;

	J[2][0] = 0.0;
	J[2][1] = prod3;
	J[2][2] = 0.0;

	return;

} // Jacobian

void Mass(double **M)
{

} // Mass
