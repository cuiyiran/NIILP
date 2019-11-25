/* CGLS.c */
#include "mex.h"
#include "blas.h"

// compressed column storage data structure
typedef struct sparseCCS
{
    double *AC;  /* numerical values, size nzmax */
    mwIndex *ia; /* row indices */
    mwIndex *jp; /* column pointers (n+1) */
    mwSize m;    /* number of rows */
    mwSize n;    /* number of columns */
} ccs;

double eps = 1.0e-6, ieps = 1.0e-1, one = 1.0, zero = 0.0;
mwSize maxnin = 1;


// How to use
void usage()
{
    mexPrintf("CGNE: CGNE method\n");
    mexPrintf("This Matlab-MEX function is for the minimum-norm solution of.\n");
    mexPrintf("linear systems Ax=b.\n");

    mexPrintf("  [x, relres, iter] = CGNE(A', b, tol, maxit);\n\n");
    mexPrintf("  valuable | size | remark \n");
    mexPrintf("  A'        m-by-n   coefficient matrix. must be sparse array.\n");
    mexPrintf("  ** REMARK:   matrix A must be TRANSPOSED ** \n");
    mexPrintf("  b         n-by-1   right-hand side vector\n");
    mexPrintf("  tol       scalar   tolerance for stopping criterion.\n");
    mexPrintf("  maxit     scalar   maximum number of iterations.\n");
    mexPrintf("  x         m-by-1   resulting approximate solution.\n");
    mexPrintf("  relres   iter-by-1 relative residual history.\n");
    mexPrintf("  iter      scalar   number of iterations required for convergence.\n");
}

// NE-SSOR inner iterations
void NESSOR(const ccs *A, double *rhs, double *Aei, double *u, double *x)
{
	double d, *AC;
	mwIndex i, *ia, inc1 = 1, j, *jp, k, k1, k2, l;
	mwSize m, n;

	AC = A->AC;
	ia = A->ia;
	jp = A->jp;
	m  = A->m;
	n  = A->n;

	for (j=0; j<n; ++j) u[j] = zero;
	for (i=0; i<m; ++i) x[i] = zero;

	for (k=0; k<maxnin; ++k) {

		for (j=0; j<n; ++j) {
			k1 = jp[j];
			k2 = jp[j+1];
			for (d=zero, l=k1; l<k2; ++l) d += AC[l]*x[ia[l]];
			d = (rhs[j] - d) * Aei[j];
			for (l=k1; l<k2; ++l) x[ia[l]] += d*AC[l];
			u[j] = u[j] + d;
		}

		j = n;
		while (j--) {
			k1 = jp[j];
			k2 = jp[j+1];
			for (d=zero, l=k1; l<k2; ++l) d += AC[l]*x[ia[l]];
			d = (rhs[j] - d) * Aei[j];
			for (l=k1; l<k2; ++l) x[ia[l]] += d*AC[l];
			u[j] = u[j] + d;
		}

	}
}


// CGNE method
void CGNE(const ccs *A, double *b, mwIndex maxit, double *iter, double *relres, double *x){

	double *p, *r, *w, *z, *Aei, *AC;
	double alpha, beta, gmm, gmm_prv, inprod, invnrmb, nrmb, nrmr, tmp, Tol;
	mwIndex i, *ia, inc1 = 1, j, *jp, k, k1, k2, l;
	mwSize m, n;

	AC = A->AC;
	ia = A->ia;
	jp = A->jp;
	m  = A->m;
	n  = A->n;

	// Allocate p
	if ((p = (double *)mxMalloc(sizeof(double) * m)) == NULL) {
		mexErrMsgTxt("Failed to allocate p");
	}

	// Allocate r
	if ((r = (double *)mxMalloc(sizeof(double) * n)) == NULL) {
		mexErrMsgTxt("Failed to allocate r");
	}

	// Allocate s
	if ((w = (double *)mxMalloc(sizeof(double) * m)) == NULL) {
		mexErrMsgTxt("Failed to allocate s");
	}

	// Allocate Ap
	if ((z = (double *)mxMalloc(sizeof(double) * n)) == NULL) {
		mexErrMsgTxt("Failed to allocate Ap");
	}

	// Allocate Aei
	if ((Aei = (double *)mxMalloc(sizeof(double) * n)) == NULL) {
		mexErrMsgTxt("Failed to allocate Aei");
	}

	for (j=0; j<n; ++j) {
		k1 = jp[j];
		k2 = jp[j+1];
		for (inprod=zero, l=k1; l<k2; ++l) inprod += AC[l]*AC[l];
		if (inprod > zero) {
			Aei[j] = one / sqrt(inprod);
			for (l=k1; l<k2; ++l) AC[l] = AC[l]*Aei[j];
		} else {
			mexPrintf("%.15e\n", AC[j]);
			mexErrMsgTxt("'warning: ||aj|| = 0");
		}
		b[j] = Aei[j] * b[j];
        Aei[j] = one; //cui: assign ||a(i,:)|| = 1
	}

  	for (i=0; i<m; ++i) x[i] = zero;

  	for (j=0; j<n; ++j) r[j] = b[j];

  	nrmb = dnrm2(&n, b, &inc1);

  	invnrmb = one / nrmb;

	Tol = eps * nrmb;

   	NESSOR(A, r, Aei, z, p);

   	for (gmm=zero, j=0; j<n; ++j) gmm += r[j] * z[j];
	gmm_prv = gmm;

	for (k=0; k<maxit; ++k) {

		tmp = dnrm2(&m, p, &inc1);
		tmp = tmp * tmp;
		alpha = gmm / tmp;

		daxpy(&m, &alpha, p, &inc1, x, &inc1);

		// w = A x
		for (j=0; j<n; ++j) {
			k1 = jp[j];
			k2 = jp[j+1];
			for (tmp=zero, l=k1; l<k2; ++l) tmp += AC[l]*p[ia[l]];
			z[j] = tmp;
		}

		tmp = -alpha;
		daxpy(&n, &tmp, z, &inc1, r, &inc1);

		nrmr = dnrm2(&n, r, &inc1);
		relres[k] = nrmr * invnrmb;

		if (nrmr < Tol) {
			iter[0] = (double)(k + 1);
			mxFree(p);
			mxFree(r);
			mxFree(w);
			mxFree(z);
			mxFree(Aei);
			return;
		}

		NESSOR(A, r, Aei, z, w);

		for (gmm=zero, j=0; j<n; ++j) gmm += r[j] * z[j];

		beta = gmm / gmm_prv;
		gmm_prv = gmm;

		for (i=0; i<m; ++i) p[i] = w[i] + beta*p[i];

	}

	mexPrintf("Failed to converge.\n");
	iter[0] = (double)(k);

	mxFree(p);
	mxFree(r);
	mxFree(w);
	mxFree(z);
	mxFree(Aei);

	return;

}

/* form sparse matrix data structure */
ccs *form_ccs(ccs *A, const mxArray *Amat)
{
    A->jp = (mwIndex *)mxGetJc(Amat);
    A->ia = (mwIndex *)mxGetIr(Amat);
    A->m = mxGetM(Amat);
    A->n = mxGetN(Amat);
    A->AC = mxGetPr(Amat);
    return (A) ;
}


// Main
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ccs *A, Amat;
    double *b, *iter, *relres, *x;
    mwIndex maxit;
    mwSize m, n;

	// Check the number of input arguments
    if(nrhs > 4){
    	usage();
        mexWarnMsgTxt("Too many inputs. Ignored extras.\n");
    }

    // Check the number of output arguments
    if(nlhs > 4){
        mexWarnMsgTxt("Too many outputs. Ignored extras.");
    }

    // Check the number of input arguments
    if (nrhs < 1) {
        usage();
        mexErrMsgTxt("Input A.");
    } else if (nrhs < 2) {
        usage();
        mexErrMsgTxt("Input b.");
    }

	// Check the 1st argument
    if (!mxIsSparse(prhs[0]))  {
        usage();
        mexErrMsgTxt("1st input argument must be a sparse array.");
    } else if (mxIsComplex(prhs[0])) {
        mexErrMsgTxt("1st input argument must be a real array.");
    }

    A = form_ccs(&Amat, prhs[0]);
    m = A->m;
    n = A->n;

	// Check the 2nd argument
    if (mxGetM(prhs[1]) != n) {
    	usage();
    	mexErrMsgTxt("The length of b is not the numer of columns of A'.");
    }

    b = mxGetPr(prhs[1]);

	// Check the 3rd argument
    // Set eps
    if (nrhs < 3) {
        mexPrintf("Default: stopping criterion is set to 1e-6.\n");
    } else {
    	if (mxIsComplex(prhs[2]) || mxGetM(prhs[2])*mxGetN(prhs[2])!=1) {
    		usage();
    		mexErrMsgTxt("3nd argument must be a scalar");
    	} else {
    		eps = *(double *)mxGetPr(prhs[2]);
    		if (eps<zero || eps>=one) {
    			usage();
    			mexErrMsgTxt("3nd argument should be positive and less than or equal to 1.");
    		}
    	}
    }

	// Check the 4th argument
	// Set maxit
  	if (nrhs < 4) {
    	maxit = n;
		// mexPrintf("Default: max number of iterations is set to the number of columns.\n");
	} else {
      if (mxIsComplex(prhs[3]) || mxGetM(prhs[3])*mxGetN(prhs[3])!=1) {
   		usage();
    	mexErrMsgTxt("4th argument must be a scalar");
    } else {
   		maxit = (mwIndex)*mxGetPr(prhs[3]);
   		if (maxit < 1) {
   			usage();
   			mexErrMsgTxt("4th argument must be a positive scalar");
         }
      }
   }

	plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(maxit, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

   x = mxGetPr(plhs[0]);
   relres = mxGetPr(plhs[1]);
   iter = mxGetPr(plhs[2]);

    CGNE(A, b, maxit, iter, relres, x);

    // Reshape relres
    mxSetPr(plhs[1], relres);
    mxSetM(plhs[1], (mwSize)(iter[0]));
    mxSetN(plhs[1], 1);

}
