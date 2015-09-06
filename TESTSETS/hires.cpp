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
	f[0] = -1.71*y[0] + 0.43*y[1] + 8.32*y[2] + 0.0007;
	f[1] = 1.71*y[0] - 8.75*y[1];
	f[2] = -10.03*y[2] + 0.43*y[3] + 0.035*y[4];
	f[3] = 8.32*y[1] + 1.71*y[2] - 1.12*y[3];
	f[4] = -1.745*y[4] + 0.43*y[5] + 0.43*y[6];
	f[5] = -280.0*y[5]*y[7] + 0.69*y[3] + 1.71*y[4] - 0.43*y[5] + 0.69*y[6];
	f[6] = 280.0*y[5]*y[7] - 1.81*y[6];
	f[7] = -f[6];

} // fEval

void Jacobian(double x, double *y, double **J)
{
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			J[i][j] = 0.0;

	J[0][0] = -1.71;
	J[0][1] = 0.43;
	J[0][2] = 8.32;

	J[1][0] = 1.71;
	J[1][1] = -8.75;

	J[2][2] = -10.03;
	J[2][3] = 0.43;
	J[2][4] = 0.035;

	J[3][1] = 8.32;
	J[3][2] = 1.71;
	J[3][3] = -1.12;

	J[4][4] = -1.745;
	J[4][5] = 0.43;
	J[4][6] = 0.43;

	J[5][3] = 0.69;
	J[5][4] = 1.71;
	J[5][5] = -0.43 - 280.0*y[7];
	J[5][6] = 0.69;
	J[5][7] = -280.0*y[5];

	J[6][5] = 280.0*y[7];
	J[6][6] = -1.81;
	J[6][7] = 280.0*y[5];

	J[7][5] = -280.0*y[7];
	J[7][6] = 1.81;
	J[7][7] = -280.0*y[5];

	return;

} // Jacobian

void Mass(double **M)
{


} // Mass
