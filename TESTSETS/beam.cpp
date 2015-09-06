#include <radau5.h>

void fbeam(const int n, double x, double *y, double *fy, double *rpar,
	int *ipar);

void jbeam(const int n, double x, double *y, double **J, double *rpar,
	int *ipar);

void mbeam(const int n, double **am, double *rpar, int *ipar);

int solout(int nr, double xold, double x, double *y, double *cont, const int n,
			double *rpar, int *ipar, double hsol);

int main ()
{
	const int nd = 80;

	// differential equation is in explicit form
	int imas = 0;
	int mlmas = 0, mumas = 0;
	// output routine is not used during integration
	int iout = 0;

	// required tolerance
	int itol = 0;
	double rtol = 1.0e-10;
	double atol = rtol;

	double work[9];
	int iwork[17];

	//initial value of x
	double x = 0.0;

	//endpoint of integration
	const double xend = 5.0;

	int idid;
	double rpar[1];
	int ipar[1];
	
	rpar[0] = 0.0;
	ipar[0] = 0;

	double *y = new double[nd];

	// set initial values for y
	for (int i = 0; i < nd; i++) y[i] = 0.0;

	// initial step size
	double h = 1.0e-3;

	// set default values
	for(int i = 0; i < 9; i++) work[i] = 0.0;
	for(int i = 0; i < 17; i++) iwork[i] = 0;

	work[2] = 0.1;
	work[3] = 0.3;
	work[4] = 0.99;
	work[5] = 2.0;

	// compute the jacobian numerically
	const int ijac = 0;
	int mujac = 0;

	// second order option
	iwork[8] = nd/2;
	int mljac = nd/2;

	idid = radau5(nd, fbeam, x, y, xend, &h, &rtol, &atol, itol, jbeam,
			ijac, mljac, mujac, mbeam, imas, mlmas, mumas, solout,
			iout, work, iwork, rpar, ipar);

	if (idid < 0) {
		cout << " Computation failed " << endl;
		return (0);
	}

	// print final solution
	cout << "\tt = " << xend;
	for (int i = 0; i < 2; i++) cout << "\ty[" << i << "] = " << y[i] << "  ";
	cout << endl;

	// print statistics
	cout << "rtol = " << rtol << endl;
	cout << "fcn = " << iwork[10] << " jac = " << iwork[11] << " step = " <<
			iwork[12] << " accpt = " << iwork[13] << " rejct = " << iwork[14] <<
			" dec = " << iwork[15] << " sol = " << iwork[16] << endl;

	delete [] y;

	return (0);
}

void fbeam(const int nd, double t, double *th, double *f, double *rpar,
	int *ipar)
{

	int n = nd/2;
	int nsq = n*n;
	int nquatr = nsq*nsq;

//	double u[n], v[n], w[n], alpha[n], beta[n], sth[n], cth[n];
	double u[40], v[40], w[40], alpha[40], beta[40], sth[40], cth[40];
	double term1, term2, q, fx, fy, fabsol;

// calcul th[i] and sin and cos

	for (int i = 1; i < n; i++) {
		sth[i] = sin(th[i] - th[i-1]);
		cth[i] = cos(th[i] - th[i-1]);
	}

// calcul du cote droit du systeme lineaire

	if (t > 3.14159265358979324) {
	// case t > pi
		// i = 0
		term1 = (-3.0*th[0] + th[1])*nquatr;
		v[0] = term1;
		// i = 1,...,n-2
		for (int i = 1; i < n-1; i++) {
			term1 = (th[i-1] - 2.0*th[i] + th[i+1])*nquatr;
			v[i] = term1;
		}
		// i = n-1
		term1 = (th[n-2] - th[n-1])*nquatr;
		v[n-1] = term1;
	}
	else {
	// case t <= pi
		fabsol = 1.5*sin(t)*sin(t);
		fx = -fabsol;
		fy = fabsol;
		// i = 0
		term1 = (-3.0*th[0] + th[1])*nquatr;
		term2 = nsq*(fy*cos(th[0]) - fx*sin(th[0]));
		v[0] = term1 + term2;
		// i = 1,...,n-2
		for (int i = 1; i < n-1; i++) {
			term1 = (th[i-1] - 2.0*th[i] + th[i+1])*nquatr;
			term2 = nsq*(fy*cos(th[i]) - fx*sin(th[i]));
			v[i] = term1 + term2;
		}
		// i = n-1
		term1 = (th[n-2] - th[n-1])*nquatr;
		term2 = nsq*(fy*cos(th[n-1]) - fx*sin(th[n-1]));
		v[n-1] = term1 + term2;
	}

// compute product dv = w
	w[0] = sth[1]*v[1];
	for (int i = 1; i < n-1; i++)
		w[i] = -sth[i]*v[i-1] + sth[i+1]*v[i+1];
	w[n-1] = -sth[n-1]*v[n-2];

// term 3
	for (int i = 0; i < n; i++)
		w[i] += th[n+i]*th[n+i];

// solve system cw = w
	alpha[0] = 1.0;
	for (int i = 1; i < n; i++) {
		alpha[i] = 2.0;
		beta[i-1] = -cth[i];
	}
	alpha[n-1] = 3.0;
	for (int i = n-2; i >= 0; i--) {
		q = beta[i]/alpha[i+1];
		w[i] -= w[i+1]*q;
		alpha[i] -= beta[i]*q;
	}
	w[0] = w[0]/alpha[0];
	for (int i = 1; i < n; i++)
		w[i] = (w[i] - beta[i-1]*w[i-1])/alpha[i];

// compute u = cv + dw
	u[0] = v[0] - cth[1]*v[1] + sth[1]*w[1];
	for (int i = 1; i < n-1; i++)
		u[i] = 2.0*v[i] - cth[i]*v[i-1] - cth[i+1]*v[i+1]
				-sth[i]*w[i-1] + sth[i+1]*w[i+1];
	u[n-1] = 3.0*v[n-1] - cth[n-1]*v[n-2] - sth[n-1]*w[n-2];

// put derivatives in right place
	for (int i = 0; i < n; i++) {
		f[i] = th[n+i];
		f[n+i] = u[i];
	}

	return;

}

// dummy function, calculates the Jacobian numerically
void jbeam(const int n, double x, double *y, double **J, double *rpar,
	int *ipar)
{
}

// dummy function, the equation is explicit
void mbeam(const int n, double **M, double *rpar, int *ipar)
{
}


/* this is the solout function for radau5.
   As we use the continuous output function, we need
   a special version calling contr5 */
int solout(int nr, double xold, double x, double *y, double *cont,
			const int n, double *rpar, int *ipar, double hsol)
{

  /*      SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N,RPAR,IPAR,IRTRN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CONT(LRC)
           WRITE (6,99) X,Y(1),Y(2),NR-1
 99     FORMAT(1X,'X =',F6.2,'    Y =',2F18.7,'    NSTEP =',I4)
        RETURN
        END */

	double dt = 0.1;

//#define EVERY_STEP

#ifdef EVERY_STEP
	cout << "Step " << nr << ":\tt = " << d;
	for (int i = 0; i < n; i++) cout << "\ty[" << i << "] = " << yd[i] << "  ";
	cout << endl;
#else
	static double d;
	double *yd = new double[n];

	if (nr == 1) d = xold;

	while (d < x) {
		if ((xold <= d) && (x >= d)) {
			cout << "Step " << nr << ": t = " << d;
			for (int i = 0; i < 2; i++) {
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

