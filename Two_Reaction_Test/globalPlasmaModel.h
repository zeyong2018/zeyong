#pragma once
void targetFunc(int n, double *y, double *F);
double Ew_cacl(int n, double *y);
double k1_cacl(int n, double *y);
void Jacobian_F(int n, double *y, double **jacobian);
void gaussElimination(int n, double **A, double *b, double *x);
double ** dmatrix(int m, int n);
void Jacobian_G(int n, double *y, double dh, double **jacobian);
void targetFunc_ODE(int n, double *y, double *y_old, double dh, double *F);