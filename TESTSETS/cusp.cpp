#include <radau5.h>

void fcusp(const int n, double x, double *y, double *fy, double *rpar,
	int *ipar);

void jcusp(const int n, double x, double *y, double **J, double *rpar,
	int *ipar);

void mcusp(const int n, double **am, double *rpar, int *ipar);

int solout(int nr, double xold, double x, double *y, double *cont, const int n,
			double *rpar, int *ipar, double hsol);

int main ()
{

	const int nd = 96;

	double *y = new double[nd];

	double work[9];
	int iwork[17];

	int nnerv = nd/3;

	// compute the jacobian numerically
	const int ijac = 0;
	// jacobian is a banded matrix
	int mljac = 3, mujac = 3;

	// differential equation is in explicit form
	int imas = 0, mlmas = 0, mumas = 0;

	// output routine is not used during integration
	const int iout = 0;

	double rpar[1];
	int ipar[1];
	
	rpar[0] = 0.0;
	ipar[0] = 0;

	//initial value of x
	double x = 0.0;
	//endpoint of integration
	const double xend = 1.1;

	double del = 2.0*3.14159265358979324/nnerv;

	for (int i = 0; i < nnerv; i++) {
		y[3*i] = 0.0;
		y[3*i+1] = -2.0*cos(double(i+1)*del);
		y[3*i+2] = 2.0*sin(double(i+1)*del);
	}

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

	idid = radau5(nd, fcusp, x, y, xend, &h, &rtol, &atol, itol, jcusp,
			ijac, mljac, mujac, mcusp, imas, mlmas, mumas, solout,
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

void fcusp(const int nd, double t, double *y, double *f, double *rpar,
	int *ipar)
{

	int nnerv = nd/3;
	double diffus = double(nnerv*nnerv)/144.0;

	double x, a, b, xright, aright, bright, xleft, aleft, bleft;
	double xdot, adot, bdot, u, v;

	for (int i = 0; i < nnerv; i++) {
		x = y[3*i];
		a = y[3*i+1];
		b = y[3*i+2];
		if (i == 0) {
			xright = y[3*nnerv-3];
			aright = y[3*nnerv-2];
			bright = y[3*nnerv-1];
		}
		else {
			xright = y[3*i-3];
			aright = y[3*i-2];
			bright = y[3*i-1];
		}
		if (i == nnerv-1) {
			xleft = y[0];
			aleft = y[1];
			bleft = y[2];
		}
		else {
			xleft = y[3*i+3];
			aleft = y[3*i+4];
			bleft = y[3*i+5];
		}
		xdot = -10000.0*(b + x*(a + x*x));
		u = (x - 0.7)*(x - 1.3);
		v = u/(u + 0.1);
		adot = b + 0.07*v;
		bdot = ((1.0 - a*a)*b - a) - 0.4*x + 0.035*v;
		f[3*i] = xdot + diffus*(xleft - 2.0*x + xright);
		f[3*i+1] = adot + diffus*(aleft - 2.0*a + aright);
		f[3*i+2] = bdot + diffus*(bleft - 2.0*b + bright);
	}

	return;

}


/* dummy routine */
void jcusp(const int nd, double x, double *y, double **J, double *rpar,
	int *ipar)
{

}

/* dummy routine */
void mcusp(const int n, double **M, double *rpar, int *ipar)
{

}


/* this is the solout function for radau5.
   As we use the continuous output function, we need
   a special version calling contr5 */
int solout(int nr, double xold, double x, double *y, double *cont,
			const int n, double *rpar, int *ipar, double hsol)
{

/*        SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(2),NR-1
 99     FORMAT(1X,'X =',F6.2,'    Y =',2F18.7,'    NSTEP =',I4)
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
