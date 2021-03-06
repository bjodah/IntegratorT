// main.cpp for Amplifier problem

#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"

int main ()
{
	// dimension of problem
	const int n(8);
	// initial values for y
	double y[n] = {0.0, 6.0, 3.0, 3.0, 6.0, 3.0, 3.0, 0.0};
	// initial value for x
	double xbeg(0.0);
	// final value for x
	const double xend(0.03);
	// interval of x for printing output
	double dx(0.0025);
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	double *rtoler = new double(1.0e-5);
	// absolute tolerance
	double *atoler = new double(1.0e-11);
	// use SolutionOutput routine
	const int iout(1);
	// analytical Jacobian function provided
	const int ijac(1);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(1);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(2);
	// Mass matrix routine is provided
	const int imas(1);
	// number of non-zero rows below main diagonal of mass matrix
	int mlmas(1);
	// number of non-zero rows above main diagonal of mass matrix
	int mumas(1);
	
// Use default values (see header files) for these parameters:
	double hinit(0.0);
	double hmax(0.0);
	int nmax(0);
	double uround(0.0), safe(0.0), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);

	StiffIntegratorT stiffT(n, y, xbeg, xend, dx, itoler, rtoler, atoler,
		iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
		mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
		m1, m2, hess, fnewt, quot1, quot2, thet);

	cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
	
	stiffT.Integrate();
	
	// print statistics
	cout << "fcn = " << stiffT.NumFunction() <<
		    " jac = " << stiffT.NumJacobian() <<
			" step = " << stiffT.NumStep() <<
			" accpt = " << stiffT.NumAccept() <<
			" rejct = " << stiffT.NumReject() <<
			" dec = " << stiffT.NumDecomp() <<
			" sol = " << stiffT.NumSol() << endl;

	cout << "\n\n*******Problem integrated with DOPRI5*******\n\n";

	double y2[n] = {0.0, 6.0, 3.0, 3.0, 6.0, 3.0, 3.0, 0.0};
	const int iout2 = 2;
	
	// Use default values (see header files) for these parameters:
	hinit = 0.0;
	double beta = 0.0;
	int nstiff = 0;
	int nrdens = n;
	unsigned *icont = NULL;

	NonStiffIntegratorT nonstiffT(n, y2, xbeg, xend, dx, nrdens, itoler, rtoler,
		atoler, iout2, hinit, hmax, nmax, uround, safe, facl, facr, beta,
		nstiff, icont);

	nonstiffT.Integrate();

	// print statistics
	cout << "fcn = " << nonstiffT.NumFunction() <<
		    " step = " << nonstiffT.NumStep() <<
			" accpt = " << nonstiffT.NumAccept() <<
			" rejct = " << nonstiffT.NumReject() << endl;

	delete rtoler;
	delete atoler;

	return (0);
}

void Function(double x, double *y, double *f)
{
	double ue, ub, uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
	double w, uet, fac1, fac2;

	ue = 0.1;
	ub = 6.0;
	uf = 0.026;
	alpha = 0.99;
	beta = 1.0e-6;
	r0 = 1000.0;
	r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = 9000.0;
	w = 2.0*3.141592654*100.0;
	uet = ue*sin(w*x);
	fac1 = beta*(exp((y[3] - y[2])/uf) - 1.0);
	fac2 = beta*(exp((y[6] - y[5])/uf) - 1.0);

	f[0] = y[0]/r9;
	f[1] = (y[1] - ub)/r8 + alpha*fac1;
	f[2] = y[2]/r7 - fac1;
	f[3] = y[3]/r5 + (y[3] - ub)/r6 + (1.0 - alpha)*fac1;
	f[4] = (y[4] - ub)/r4 + alpha*fac2;
	f[5] = y[5]/r3 - fac2;
	f[6] = y[6]/r1 + (y[6] - ub)/r2 + (1.0 - alpha)*fac2;
	f[7] = (y[7] - uet)/r0;

	return;

} // Function

void Jacobian(double x, double *y, double **J)
{

	double uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
	double fac14, fac27;

	uf = 0.026;
	alpha = 0.99;
	beta = 1.0e-6;
	r0 = 1000.0;
	r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = 9000.0;
	fac14 = beta*exp((y[3] - y[2])/uf)/uf;
	fac27 = beta*exp((y[6] - y[5])/uf)/uf;

	for (int j = 0; j < 8; j++) {
		J[0][j] = 0.0;
		J[1][j] = 0.0;
		J[3][j] = 0.0;
	}
	J[2][0] = 1.0/r9;
	J[2][1] = 1.0/r8;
	J[1][2] = -alpha*fac14;
	J[0][3] = alpha*fac14;
	J[2][2] = 1.0/r7 + fac14;
	J[1][3] = -fac14;
	J[2][3] = 1.0/r5 + 1.0/r6 + (1.0 - alpha)*fac14;
	J[3][2] = -(1.0 - alpha)*fac14;
	J[2][4] = 1.0/r4;
	J[1][5] = -alpha*fac27;
	J[0][6] = alpha*fac27;
	J[2][5] = 1.0/r3 + fac27;
	J[1][6] = -fac27;
	J[2][6] = 1.0/r1 + 1.0/r2 + (1.0 - alpha)*fac27;
	J[3][5] = -(1.0 - alpha)*fac27;
	J[2][7] = 1.0/r0;

	return;

} // Jacobian

void Mass(double **M)
{
 	const double c1 = 1.0e-6;
	const double c2 = 2.0e-6;
	const double c3 = 3.0e-6;
	const double c4 = 4.0e-6;
	const double c5 = 5.0e-6;

	for (int i = 0; i < 8; i++) {
		M[0][i] = 0.0;
		M[2][i] = 0.0;
	}

	M[1][0] = -c5;
	M[0][1] = c5;
	M[2][0] = c5;
	M[1][1] = -c5;
	M[1][2] = -c4;
	M[1][3] = -c3;
	M[0][4] = c3;
	M[2][3] = c3;
	M[1][4] = -c3;
	M[1][5] = -c2;
	M[1][6] = -c1;
	M[0][7] = c1;
	M[2][6] = c1;
	M[1][7] = -c1;

	return;

} // Mass
