//
// Created by giacomo on 01/06/19.
//

#include <malloc.h>
#include "dgsol_seq_lib.h"

extern void dginitial(int *pInt, double *pDouble, int *pInt1, int *ind, int *jInd, double *lbd, double *ubd, int *seed, int *sp);
extern void dgopt(int *pInt, double *pDouble, double *opt, int *pInt1, int *ind, int *jInd, double *lbd, double *ubd, int *sp, double *wkSp);
extern void dgerr(int *pInt, double *pDouble, double *opt, int *pInt1, int *ind, int *jInd, double *lbd, double *ubd, double *err);

double* dgsol_j(int numDs, int numAtoms, int *iInd, int *jInd, double *lbd, double *ubd, int* seed, struct dgsol_summary* summary) {
    int n, k, lintWkSp, ldblWkSp, *intWkSp;
    double *xOpt, *fOpt, *dErr, *dblWkSp;

    for (k = 0; k < numDs; k ++) {
        lbd [k] = lbd [k] * lbd [k];
        ubd [k] = ubd [k] * ubd [k];
    }

    /*
    Allocate memory.
    */
    n = 3*numAtoms;
    lintWkSp = 2*numAtoms > numDs ? 2*numAtoms : numDs;

    /*
    In general ldblWkSp is 2*(m+1)*n + 2*m + 1, where m is the number
    of vectors used by VMLM.
    */
    ldblWkSp = 22*n + 21;

    fOpt = (double *) malloc (sizeof(double));
    dErr = (double *) malloc (3*sizeof(double));
    xOpt = (double *) malloc (n*sizeof(double));
    intWkSp = (int *) malloc (lintWkSp*sizeof(int));
    dblWkSp = (double *) malloc (ldblWkSp*sizeof(double));

    dginitial(&numAtoms,xOpt,&numDs,iInd,jInd,lbd,ubd,seed,intWkSp);
    dgopt(&numAtoms,xOpt,fOpt,&numDs,iInd,jInd,lbd,ubd,intWkSp,dblWkSp);
    dgerr(&numAtoms,xOpt,fOpt,&numDs,iInd,jInd,lbd,ubd,dErr);

    summary->f_err = *fOpt;
    summary->derr_min = dErr[0];
    summary->derr_avg = dErr[1];
    summary->derr_max = dErr[2];
    free(dErr);
    free(fOpt);
    free(intWkSp);
    free(dblWkSp);

    return fOpt;
}
