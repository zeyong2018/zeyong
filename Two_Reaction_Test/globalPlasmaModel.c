#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "System.h"
#include <math.h>
#include "globalPlasmaModel.h"

void targetFunc(int n, double *y, double *F)
{
	double k1, n_ar, n_e, n_arplus,V, Ew;	
	n_ar = *(y + 0);
	n_arplus = *(y + 1);
	n_e = n_arplus;

	V = PI * pow(R, 2.)*L;
	k1 = k1_cacl(n, y);
	Ew = Ew_cacl(n, y);
	*(F + 0) = qin / V - k1 * n_ar*n_e - k_pump * n_ar + k2 * n_arplus;
	*(F + 1) = -k2 * n_arplus + k1 * n_ar*n_e;
	*(F + 2) = p_abs / V - Ec1 * k1*n_ar*n_e - Ew * k2*n_arplus;
}

double Ew_cacl(int n, double *y)
{
	double Te;
	Te = *(y + 2);

	double Ew;
	Ew = 2.*Te + Te / 2.*(1 + log(m_arplus / 2 / PI / me));
	return Ew;
}

double k1_cacl(int n, double *y)
{
	double Te;
	Te = *(y + 2);

	double k1;
	k1 = (7.93e-14)*exp(-18.9 / (Te+TINY));
	return k1;
}

void Jacobian_F(int n, double *y, double **jacobian)
{
	double n_ar, n_arplus, n_e, Te, k1, Ew;
	n_ar = *(y + 0);
	n_arplus = *(y + 1);
	n_e = n_arplus;
	Te = *(y + 2);
	k1 = k1_cacl(n, y);
	Ew = Ew_cacl(n, y);

	double dk1dte, dewdte;
	dk1dte = (7.93e-14)*exp(-18.9 / (Te+TINY))*18.9 / pow(Te+TINY, 2.);
	dewdte = 2. + 1. / 2 * (1 + log(m_arplus / 2 / PI / me));

	jacobian[0][0] = -k1 * n_e - k_pump;
	jacobian[0][1] = -k1 * n_ar + k2;
	jacobian[0][2] = -n_ar * n_e * dk1dte;
	jacobian[1][0] = k1 * n_e;
	jacobian[1][1] = -k2 + k1 * n_ar;
	jacobian[1][2] = n_ar * n_e * dk1dte;
	jacobian[2][0] = -Ec1 * k1*n_e;
	jacobian[2][1] = -Ec1 * k1*n_ar - Ew * k2;
	jacobian[2][2] = -Ec1 * n_ar*n_e * dk1dte - k2 * n_arplus* dewdte;
}

/*
	gaussElimination() solves Ax=b
*/
void gaussElimination(int n, double **A, double *b, double *x)
{
	
	// Extended Matrix
	double **A_extend;
	A_extend = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		A_extend[i] = (double *)malloc((n + 1) * sizeof(double));
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A_extend[i][j] = A[i][j];
		}
		A_extend[i][n] = b[i];
	}

	// Gauss Elimination
	int p;
	double temp;
	for (int i = 0; i < n-1; i++)
	{
		p = -1;
		for (int j = i; j < n; j++)
		{
			if (fabs(A_extend[j][i]) > 0.)
			{
				p = j;
				break;
			}
		}
		if (p==-1)
		{
			printf("No unique solution exist!\n");
			exit(EXIT_FAILURE);
		}
		else
		{
			if (p!=i)
			{
				for (int k = 0; k < n + 1; k++)
				{
					temp = A_extend[i][k];
					A_extend[i][k] = A_extend[p][k];
					A_extend[p][k] = temp;

				}
			}
		}

		// Elimination
		for (int j = i+1; j < n; j++)
		{
			temp = A_extend[j][i] / A_extend[i][i];
			for (int  k= i; k < n+1; k++)
			{
				A_extend[j][k] -= A_extend[i][k] * temp;
			}
		}

	}
	if (fabs(A_extend[n-1][n-1])<TINY)
	{
		printf("No unique solution exist!\n");
		exit(EXIT_FAILURE);
	}

	// Back Substitution
	x[n-1] = A_extend[n - 1][n] / A_extend[n - 1][n - 1];

	for (int i = n-2; i >=0; i--)
	{
		temp = 0.;
		for (int j = i+1; j < n; j++)
		{
			temp += A_extend[i][j] * x[j];
		}
		x[i] = (A_extend[i][n] - temp) / A_extend[i][i];
	}

}

void dmatrix(double **A, int n) 
{
	A = (double **)malloc(3 * sizeof(double *));
	for (int i = 0; i < 3; i++)
	{
		A[i] = (double *)malloc(3 * sizeof(double));
	}
}
