#include <R.h>
#include <stdio.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>
#include <R_ext/Lapack.h>

// ----------------------------------------------------------
// Update betas until convergence, i.e., the difference
// between two consecutive solutions is smaller than a threshold
// for a set of decreasing lambda: l0>l1>...>lk
// Using the whole matrix XtX as input stored in object XL passed as vector
//
//      n:         Number of beta parameters
//      XtX:        Crossprod matrix X'X in vector form, stacked by columns
//      Xty:       Crossprod vector X'y in vector form
//      q:         Number of lambdas
//      lambda:    Vector of lambdas for which the betas will be calculated
//      a:         alpha value in the Elastic-net problem
//      maxtole:   Maximum value between two consecutive solutions for beta to be accepted for convergence
//      maxsteps:  Number of iterations to run before the updating stops
// ----------------------------------------------------------
SEXP updatebeta_lambda(SEXP n, SEXP XtX, SEXP Xty, SEXP q, SEXP lambda, SEXP a, SEXP maxtole, SEXP maxsteps, SEXP echo)
{
    double *pXtX, *pB, *pXty, *plambda, *b, *currfit;
    double bOLS, bNew;
    double L1, L2, alpha, maxTol,maxdiff;
    int j, k, np, maxIter,iter,nlambda, verbose;
    int inc=1;
    double delta;
    SEXP list;

    np=INTEGER_VALUE(n);
    nlambda=INTEGER_VALUE(q);
    maxIter=INTEGER_VALUE(maxsteps);
    verbose=asLogical(echo);
    alpha=NUMERIC_VALUE(a);
    maxTol=NUMERIC_VALUE(maxtole);

    PROTECT(XtX=AS_NUMERIC(XtX));
    pXtX=NUMERIC_POINTER(XtX);

    PROTECT(Xty=AS_NUMERIC(Xty));
    pXty=NUMERIC_POINTER(Xty);

    PROTECT(lambda=AS_NUMERIC(lambda));
    plambda=NUMERIC_POINTER(lambda);

    SEXP B = PROTECT(allocMatrix(REALSXP, nlambda, np));
    pB=NUMERIC_POINTER(B);

    b=(double *) R_alloc(np, sizeof(double));
    currfit=(double *) R_alloc(np, sizeof(double));

    memset(currfit,0, sizeof(double)*np);  // Initialize all currentfit to zero
    memset(b,0, sizeof(double)*np);  // Initialize all coefficients to zero

    for(k=0; k<nlambda; k++)
    {
        L1=alpha*plambda[k];
        L2=(1-alpha)*plambda[k];
        iter=0;
        maxdiff=INFINITY;    // Set to INF to enter to the WHILE
        while(iter<maxIter && maxdiff>maxTol)
        {
            iter++;
            maxdiff=0;
            for(j=0; j<np; j++)
            {
                //bOLS=(pXty[j] - currfit)/pXtX[np*j + j];
                bOLS=pXty[j] - currfit[j];
                if(fabs(bOLS) > L1){
                    bNew=sign(bOLS)*(fabs(bOLS)-L1)/(1+L2); // (pXtX[np*j + j]+L2)
                }else{
                    bNew=0;
                }

                delta = bNew-b[j];
                if(fabs(delta)>0){
                    // update the current fit for all variables: cf=cf+XtX[j]*(bNew-b[j])
                    F77_NAME(daxpy)(&np, &delta, pXtX + j*np, &inc, currfit, &inc);
                    currfit[j] -= delta; //delta*pXtX[np*j + j]

                    if(fabs(delta)>maxdiff){
                        maxdiff=fabs(delta);
                    }
                }
                b[j]=bNew;
            }
        }
        if(verbose){
            printf(" lambda[%d]=\t%11.9f.  nIters=%7d.  maxError=%11.9f\n",k+1,plambda[k],iter,maxdiff);
            if(maxdiff>maxTol){
              printf("    Warning: The process did not converge after %d iterations for lambda[%d]=%11.9f\n",maxIter,k+1,plambda[k]);
            }
        }
        // memcpy(pB + np*k, b, sizeof(double)*np);
        F77_NAME(dcopy)(&np, b, &inc, pB+k, &nlambda);
    }

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 1));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, B);

    UNPROTECT(5);

    return(list);
}

// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP gaussian_kernel(SEXP n, SEXP XtX, SEXP h)
{
    double *pXtX, *pd;
    double bw;
    int np;
    int i,j;

    np=INTEGER_VALUE(n);
    bw=NUMERIC_VALUE(h);

    PROTECT(XtX=AS_NUMERIC(XtX));
    pXtX=NUMERIC_POINTER(XtX);

    pd=(double *) R_alloc(np, sizeof(double));

    // Diagonal values
    for(j=0; j<np; j++)
    {
        pd[j]=pXtX[np*j + j];
        pXtX[np*j + j]=1;
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          pXtX[np*j + i]=exp(-bw*(pd[i] + pd[j] -2*pXtX[np*j + i]));
          pXtX[np*i + j]=exp(-bw*(pd[i] + pd[j] -2*pXtX[np*i + j]));
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}
// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP laplacian_kernel(SEXP n, SEXP XtX, SEXP h)
{
    double *pXtX, *pd;
    double bw;
    int np;
    int i,j;

    np=INTEGER_VALUE(n);
    bw=NUMERIC_VALUE(h);

    PROTECT(XtX=AS_NUMERIC(XtX));
    pXtX=NUMERIC_POINTER(XtX);

    pd=(double *) R_alloc(np, sizeof(double));

    // Diagonal values
    for(j=0; j<np; j++)
    {
        pd[j]=pXtX[np*j + j];
        pXtX[np*j + j]=1;
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          pXtX[np*j + i]=exp(-bw*sqrt(pd[i] + pd[j] -2*pXtX[np*j + i]));
          pXtX[np*i + j]=exp(-bw*sqrt(pd[i] + pd[j] -2*pXtX[np*i + j]));
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP polynomial_kernel(SEXP n, SEXP XtX,SEXP a, SEXP b)
{
    double *pXtX;
    double aa, bb;
    int np;
    int i,j;

    np=INTEGER_VALUE(n);
    aa=NUMERIC_VALUE(a);
    bb=NUMERIC_VALUE(b);

    PROTECT(XtX=AS_NUMERIC(XtX));
    pXtX=NUMERIC_POINTER(XtX);

    // Diagonal values
    for(j=0; j<np; j++)
    {
        pXtX[np*j + j]=pow(aa*pXtX[np*j + j] + 1, bb);
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          pXtX[np*j + i]=pow(aa*pXtX[np*j + i] + 1, bb);
          pXtX[np*i + j]=pow(aa*pXtX[np*i + j] + 1, bb);
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP crossprod2distance(SEXP n, SEXP XtX)
{
    double *pXtX, *pd;
    int np;
    int i,j;

    np=INTEGER_VALUE(n);

    PROTECT(XtX=AS_NUMERIC(XtX));
    pXtX=NUMERIC_POINTER(XtX);

    pd=(double *) R_alloc(np, sizeof(double));

    // Diagonal values
    for(j=0; j<np; j++)
    {
        pd[j]=pXtX[np*j + j];
        pXtX[np*j + j]=0;
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          pXtX[np*j + i]=pd[i] + pd[j] -2*pXtX[np*j + i];
          pXtX[np*i + j]=pd[i] + pd[j] -2*pXtX[np*i + j];
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// Scale XtX matrix to have all diagonal elements equal to one.
// Scaling is perform by dividing each entry x_ij by
// x_ij = x_ij/(sqrt(x_ii)*sqrt(x_jj))
// where sqrt(x_ii) is the SD of the variable i. These values are returned
//
//      n:         Number of variable (columns in XtX)
//      XtX:       Crossprod matrix X'X in vector form, stacked by columns
// ----------------------------------------------------------
SEXP scaleXtX(SEXP n, SEXP XtX)
{
    double *pXtX, *psdx;
    int np;
    int i,j;

    np=INTEGER_VALUE(n);

    PROTECT(XtX=AS_NUMERIC(XtX));
    pXtX=NUMERIC_POINTER(XtX);

    psdx=(double *) R_alloc(np, sizeof(double));

    // Get standard deviations
    for(i=0; i<np; i++)
    {
        psdx[i]=sqrt(pXtX[np*i + i]);
        pXtX[np*i + i]=1;
    }

    for(j=0; j<np-1; j++)
    {
        for(i=j+1; i<np; i++)
        {
          pXtX[np*j + i]=pXtX[np*j + i]/(psdx[j]*psdx[i]);
          pXtX[np*i + j]=pXtX[np*i + j]/(psdx[j]*psdx[i]);
        }
    }

    // Creating a NULL variable with 1 elements:
    SEXP out = PROTECT(allocVector(NILSXP, 1));

    UNPROTECT(2);

    return(out);
}

// ----------------------------------------------------------
// Return the index of variables that are very correlated
// output a 0-1 vector indicating whether the variable has
// a correlation of at most a threshold
//
//      COV:         Covariance matrix, stacked by columns
//      maxVal:      Maximum correlation allowed
// ----------------------------------------------------------
SEXP getCorrelated(SEXP n, SEXP COV, SEXP maxVal)
{
    double *pCOV, *psdx;
    double maxCor, corre;
    int *pout;
    int np,isok, cont;
    int i,j;
    SEXP list;

    np=INTEGER_VALUE(n);
    maxCor=NUMERIC_VALUE(maxVal);

    PROTECT(COV=AS_NUMERIC(COV));
    pCOV=NUMERIC_POINTER(COV);

    SEXP out = PROTECT(allocVector(INTSXP, np));
    pout=INTEGER_POINTER(out);

    psdx=(double *) R_alloc(np, sizeof(double));

    // Get standard deviations
    for(i=0; i<np; i++)
    {
        psdx[i]=sqrt(pCOV[np*i + i]);
    }

    cont=0;
    for(j=0; j<np-1; j++)
    {
        if(psdx[j]>0)
        {
            isok=1;
            for(i=j+1; i<np; i++)
            {
                if(psdx[i]>0){
                    corre=pCOV[np*j + i]/(psdx[j]*psdx[i]);
                    if(corre>maxCor){
                        isok=0;
                        break;
                    }
                }
            }

            if(isok==1){
              pout[cont]=j+1;
              cont++;
            }
        }
    }

    if(psdx[np-1]>0){
        pout[cont]=np;
        cont++;
    }

    SEXP le = PROTECT(ScalarInteger(cont));

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, out);
    SET_VECTOR_ELT(list, 1, le);

    UNPROTECT(4);

    return(list);
}

// ----------------------------------------------------------
// the p x p-1 matrix R has been formed from a
// p x p upper-triangular matrix by deleting column k
// integer p,k,n,nz; double precision R[p,p-1], z[p,nz]
// ----------------------------------------------------------
SEXP delete_col(SEXP R, SEXP p0, SEXP k0, SEXP z, SEXP nz0)
{
    double *pR, *pz;
    int p, k, nz;
    int i,j;
    double a, b, c, s, tau;
    SEXP list;

    p=INTEGER_VALUE(p0);
    k=INTEGER_VALUE(k0);
    nz=INTEGER_VALUE(nz0);

    PROTECT(R=AS_NUMERIC(R));
    pR=NUMERIC_POINTER(R);

    PROTECT(z=AS_NUMERIC(z));
    pz=NUMERIC_POINTER(z);

    for(i=k-1; i<p-1; i++) // loop thats goes j=i+1,...,p-1
    {
        a=pR[p*i + i];
        b=pR[p*i + i + 1];
        if(b!=0)
        {
          // Compute the rotation
	         if(fabs(b)>fabs(a)){
             tau = -a/b;
             s = 1/sqrt(1 + tau*tau);
             c = s * tau;
           }else{
             tau = -b/a;
             c = 1/sqrt(1 + tau*tau);
             s = c * tau;
           }

           // update r and z
           pR[p*i + i]=c*a - s*b;
           pR[p*i + i + 1]=s*a + c*b;

           for(j=i+1; j<p-1; j++)  // loop thats goes j=i+1,...,p-1
           {
             a=pR[p*j + i];
             b=pR[p*j + i + 1];
             pR[p*j + i]=c*a - s*b;
             pR[p*j + i + 1]=s*a + c*b;
           }
           for(j=0; j<nz; j++)  // loop thats goes j=1,...,nz
           {
             a=pz[p*j + i];
             b=pz[p*j + i + 1];
             pz[p*j + i]=c*a - s*b;
             pz[p*j + i + 1]=s*a + c*b;
           }
        }
    }

    // Creating a list with 1 vector elements:
    PROTECT(list = allocVector(VECSXP, 2));

    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, R);
    SET_VECTOR_ELT(list, 1, z);

    UNPROTECT(3);

    return(list);
}

// ----------------------------------------------------------
// ----------------------------------------------------------

SEXP writeBinFile(SEXP filename, SEXP n, SEXP p, SEXP size, SEXP X, SEXP echo)
{
    FILE *f=NULL;
    int i, j;
    int nrows, ncols, sizevar ;
    int inc=1;
    double *pX;
    double *linedouble;
    float valuefloat;
    SEXP list;

    nrows=INTEGER_VALUE(n);
    ncols=INTEGER_VALUE(p);
    sizevar=INTEGER_VALUE(size);

    PROTECT(X=AS_NUMERIC(X));
    pX=NUMERIC_POINTER(X);

    linedouble=(double *) R_alloc(ncols, sizeof(double));

    f=fopen(CHAR(STRING_ELT(filename,0)),"wb");
    fwrite(&nrows,4, 1 , f);
    fwrite(&ncols,4, 1 , f);
    fwrite(&sizevar,4, 1 , f);

    // Write lines
    for(i=0; i<nrows; i++)
    {
        if(sizevar==4)
        {
            for(j=0; j<ncols; j++)
            {
                valuefloat = pX[nrows*j + i];
                fwrite(&valuefloat,sizevar, 1 , f);
            }
        }else{
            F77_NAME(dcopy)(&ncols, pX+i, &nrows,linedouble, &inc);
            fwrite(linedouble,sizevar, ncols , f);
        }

    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 3));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, ScalarInteger(nrows));
    SET_VECTOR_ELT(list, 1, ScalarInteger(ncols));
    SET_VECTOR_ELT(list, 2, ScalarInteger(sizevar));

    UNPROTECT(2);

    return(list);
}

// ----------------------------------------------------------
// ----------------------------------------------------------
SEXP readBinFile(SEXP filename, SEXP nsetRow, SEXP nsetCol, SEXP setRow, SEXP setCol)
{
    FILE *f=NULL;
    int i, j;
    int *psetRow, *psetCol;
    int nrows, ncols, sizevar, lsrow, lscol;
    int n, p;
    double *pX;
    double *linedouble;
    float *linefloat;
    SEXP list;

    lsrow=INTEGER_VALUE(nsetRow);
    lscol=INTEGER_VALUE(nsetCol);

    PROTECT(setRow=AS_INTEGER(setRow));
    psetRow=INTEGER_POINTER(setRow);

    PROTECT(setCol=AS_INTEGER(setCol));
    psetCol=INTEGER_POINTER(setCol);

    f=fopen(CHAR(STRING_ELT(filename,0)),"rb");
    fread(&nrows, 4, 1, f);
    fread(&ncols, 4, 1, f);
    fread(&sizevar, 4, 1, f);

    linedouble=(double *) R_alloc(ncols, sizeof(double));
    linefloat=(float *) R_alloc(ncols, sizeof(float));

    n=lsrow > 0 ? lsrow : nrows;
    p=lscol > 0 ? lscol : ncols;

    SEXP X=PROTECT(allocMatrix(REALSXP, n, p));
    pX=NUMERIC_POINTER(X);

    // Read lines
    for(i=0; i<n; i++)
    {
        if(lsrow > 0){
            fseek(f, 12 + ncols*sizevar*(psetRow[i]-1), SEEK_SET);
        }

        if(sizevar==4){
            fread(linefloat,sizevar,ncols,f);
        }else{
            fread(linedouble,sizevar,ncols,f);
        }

        for(j=0; j<p; j++)
        {
            if(lscol>0)
            {
                if(sizevar==4){
                    linedouble[psetCol[j]-1]=linefloat[psetCol[j]-1];
                }

                pX[n*j + i]=linedouble[psetCol[j]-1];
            }else{
                if(sizevar==4){
                    linedouble[j]=linefloat[j];
                }
                pX[n*j + i]=linedouble[j];
            }
        }

    }

    fclose(f);

    PROTECT(list = allocVector(VECSXP, 4));
    // Attaching outputs to list:
    SET_VECTOR_ELT(list, 0, ScalarInteger(n));
    SET_VECTOR_ELT(list, 1, ScalarInteger(p));
    SET_VECTOR_ELT(list, 2, ScalarInteger(sizevar));
    SET_VECTOR_ELT(list, 3, X);

    UNPROTECT(4);
    return(list);
}
