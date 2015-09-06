#include <radau5.h>

void fbruss(const int n, double x, double *y, double *fy, double *rpar,
	int *ipar);

void jbruss(const int n, double x, double *y, double **J, double *rpar,
	int *ipar);

void mbruss(const int n, double **am, double *rpar, int *ipar);

int solout(int nr, double xold, double x, double *y, double *cont, const int n,
			double *rpar, int *ipar, double hsol);

int main ()
{

	const int nd = 1000;

	double *y = new double[nd];

	double work[9];
	int iwork[17];

	const double pi = 3.14159265358979324;
	int n = nd/2;

	double rpar[1];
	int ipar[1];

	rpar[0] = 0.0;
	ipar[0] = 0;

	//initial value of x
	double x = 0.0;
	//endpoint of integration
	const double xend = 10.0;

	for (int i = 0; i < n; i++) {
		y[2*i + 1] = 3.0;
		y[2*i] = 1.0 + 0.5*sin(2.0*pi*double(i+1)/double(n+1));
	}

	// compute the jacobian analytically
	const int ijac = 1;
	// jacobian is a banded matrix
	int mljac = 2, mujac = 2;

	// differential equation is in implicit form
	int imas = 0, mlmas = 0, mumas = 0;

	// output routine is not used during integration
	const int iout = 0;

	// required tolerance
	int itol = 0;
	double rtol = 1.0e-10;
	double atol = rtol;

	// initial step size
	double h = 1.0e-6;

	// set default values
	for(int i = 0; i < 9; i++) work[i] = 0.0;
	for(int i = 0; i < 17; i++) iwork[i] = 0;

	int idid;

	idid = radau5(nd, fbruss, x, y, xend, &h, &rtol, &atol, itol, jbruss,
			ijac, mljac, mujac, mbruss, imas, mlmas, mumas, solout,
			iout, work, iwork, rpar, ipar);

	if (idid < 0) {
		cout << " Computation failed " << endl;
		return (0);
	}

	// print final solution
	cout << "\tt = " << xend;
	for (int i = 0; i < 2; i++) cout << "\ty[" << i << "] = " << y[i] << "\t";
	cout << endl;

	// print statistics
	cout << "rtol = " << rtol << endl;
	cout << "fcn = " << iwork[10] << " jac = " << iwork[11] << " step = " <<
			iwork[12] << " accpt = " << iwork[13] << " rejct = " << iwork[14] <<
			" dec = " << iwork[15] << " sol = " << iwork[16] << endl;

	delete [] y;

	return (0);
}

void fbruss(const int nd, double x, double *y, double *f, double *rpar,
	int *ipar)
{
	int n, i = 0, iu, iv;
	double ui, vi, uim, vim, uip, vip, prod;
	double gamma, usdelq;

	n = nd/2;
	usdelq = (double(n+1))*(double(n+1));
	gamma = 0.02*usdelq;

	iu = 2*i;
	iv = 2*i + 1;
	ui = y[iu];
	vi = y[iv];
	uim = 1.0;
	vim = 3.0;
	uip = y[iu+2];
	vip = y[iv+2];
	prod = ui*ui*vi;
	f[iu] = 1.0 + prod - 4.0*ui + gamma*(uim - 2.0*ui + uip);
	f[iv] = 3.0*ui - prod + gamma*(vim - 2.0*vi + vip);

	for (i = 1; i < n-1; i++) {
		iu = 2*i;
		iv = 2*i + 1;
		ui = y[iu];
		vi = y[iv];
		uim = y[iu-2];
		vim = y[iv-2];
		uip = y[iu+2];
		vip = y[iv+2];
		prod = ui*ui*vi;
		f[iu] = 1.0 + prod - 4.0*ui + gamma*(uim - 2.0*ui + uip);
		f[iv] = 3.0*ui - prod + gamma*(vim - 2.0*vi + vip);
	}

	i = n - 1;
	iu = 2*i;
	iv = 2*i + 1;
	ui = y[iu];
	vi = y[iv];
	uim = y[iu-2];
	vim = y[iv-2];
	uip = 1.0;
	vip = 3.0;
	prod = ui*ui*vi;
	f[iu] = 1.0 + prod - 4.0*ui + gamma*(uim - 2.0*ui + uip);
	f[iv] = 3.0*ui - prod + gamma*(vim - 2.0*vi + vip);

	return;
}


/* writes the jacobian (NOT transposed!!)*/
void jbruss(const int nd, double x, double *y, double **J, double *rpar,
	int *ipar)
{

	int n, iu, iv;
	double ui, vi, uivi, ui2;
	double gamma, gamma2, usdelq;

	n = nd/2;
	usdelq = (double(n+1))*(double(n+1));
	gamma = 0.02*usdelq;
	gamma2 = 2.0*gamma;

	for (int i = 0; i < n; i++) {
		iu = 2*i;
		iv = 2*i + 1;
		ui = y[iu];
		vi = y[iv];
		uivi = ui*vi;
		ui2 = ui*ui;
		J[2][iu] = 2.0*uivi - 4.0 - gamma2;
		J[1][iv] = ui2;
		J[3][iu] = 3.0 - 2.0*uivi;
		J[2][iv] = -ui2 - gamma2;
		J[1][iu] = 0.0;
		J[3][iv] = 0.0;
	}
	for (int i = 0; i < 2*n - 2; i++) {
		J[0][i+2] = gamma;
		J[4][i] = gamma;
	}

	return;

}

/* matrix M for the amplifier problem */
void mbruss(const int n, double **M, double *rpar, int *ipar)
{

}


/* this is the solout function for radau5.
   As we use the continuous output function, we need
   a special version calling contr5 */
int solout(int nr, double xold, double x, double *y, double *cont,
			const int n, double *rpar, int *ipar, double hsol)
{

 /*       SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(2),Y(3),NR-1
 99     FORMAT(1X,'X =',F11.4,'    Y =',3E18.10,'    NSTEP =',I4)
        RETURN
        END*/

	double dt = 0.0025;

//#define EVERY_STEP

#ifdef EVERY_STEP
	cout << "Step " << nr << ":\tt = " << d;
	for (int i = 0; i < 2; i++) cout << "\ty[" << i << "] = " << yd[i] << "  ";
	cout << endl;
#else
	static double d;
	double *yd = new double[n];

	if (nr == 1) d = xold;

	while (d < x) {
		if ((xold <= d) && (x >= d)) {
			cout << "Step " << nr - 1 << ": t = " << d;
			for (int i = 0; i < n; i++) {
				yd[i] = contr5(i, d, cont, n, x, hsol);
				cout << "\ty[" << i << "] = " << yd[i];
			}
			cout << endl;
			d += dt;
		}
	}
	delete [] yd;
#endif
	return 0;

}
