#ifndef __LINALGEXT_H__
#define __LINALGEXT_H__

extern void dgesdd_(const char *jobz, int *m, int *n, double *a, int *lda, double *s,
		    double *u, int *ldu, double *vt, int *ldvt,
		    double *work, int *lwork, int *iwork, int *info);

void linalg_dgesdd(double **, int, int, double *, double *, double **);
#endif
