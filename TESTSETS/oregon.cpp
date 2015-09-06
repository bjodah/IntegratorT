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
	f[0] = 77.27*(y[1] + y[0]*(1.0 - 8.375e-6*y[0] - y[1]));
	f[1] = (y[2] - (1.0 + y[0])*y[1])/77.27;
	f[2] = 0.161*(y[0] - y[2]);

} // fEval

void Jacobian(double x, double *y, double **J)
{
	J[0][0] = 77.27*(1.0 - 2.0*8.375e-6*y[0] - y[1]);
	J[0][1] = 77.27*(1.0 - y[0]);
	J[0][2] = 0.0;

	J[1][0] = -y[1]/77.27;
	J[1][1] = -(1.0 + y[0])/77.27;
	J[1][2] = 1.0/77.27;

	J[2][0] = 0.161;
	J[2][1] = 0.0;
	J[2][2] = -0.161;

	return;

} // Jacobian

void Mass(double **M)
{


} // Mass
