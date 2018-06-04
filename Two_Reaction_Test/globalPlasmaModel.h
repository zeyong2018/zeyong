#pragma once
void targetFunc(int n, double *y, double *F);
double Ew_cacl(int n, double *y);
double k1_cacl(int n, double *y);
void Jacobian_F(int n, double *y, double **jacobian);
void gaussElimination(int n, double **A, double *b, double *x);
void dmatrix(double **A, int n);