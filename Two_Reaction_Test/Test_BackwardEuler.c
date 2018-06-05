#define  _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "globalPlasmaModel.h"
#include "stdlib.h"
#include "System.h"

void Test_BackwardEuler(void)
{
	qin = 7.1728e18;
	L = 0.14;
	R = 0.15;
	k_pump = 3.4088;
	k2 = 2.6582e4;
	ee = 1.6022E-19;
	p_abs = 700. / ee;
	Ec1 = 57.7141;
	me = 9.1095E-31;
	m_arplus = 6.6335209e-26 - me;

	double y_prev[3] = { 7.1728e18, 7.1728e15, 2.};
	double y_temp[3] = { 0,0,0 };
	double y_iter_old[3] = { 0.,0,0 };
	double **Jacobian;
	int n = 3;
	Jacobian = dmatrix(3, 3);
	double G[3] = { 0,0,0 };

	double dh;
	int N_STEP = 2000000;
	dh = 2. / N_STEP;
	FILE *f_out;
	f_out = fopen("Output.csv", "w");
	fprintf(f_out, "%12s,%12s,%12s,%12s\n", "t", "n_ar", "n_ar+", "Te");
	double y[3] = { 7.1728e18, 7.1728e15, 2. }; // initial guess for first step
	int N_MAX = 1000;
	double t,sum;
	fprintf(f_out, "%12.4e,%12.4e,%12.4e,%12.4e\n", 0., y_prev[0], y_prev[1], y_prev[2]);
	for (int i = 0; i < N_STEP; i++)
	{
		t = (i + 1)*dh;
		for (int j = 0; j < N_MAX; j++)
		{
			targetFunc_ODE(n, y, y_prev, dh, G);
			Jacobian_G(n, y, dh, Jacobian);
			gaussElimination(n, Jacobian, G, y_temp);
			sum = 0.;
			for (int k = 0; k < n; k++)
			{
				y_iter_old[k] = y[k];
				y[k] -= y_temp[k];
				sum += fabs((y[k] - y_iter_old[k]) / (y[k] + TINY));
			}
			if (sum<=3.e-3)
			{
				break;
			}


		}
		printf("%12.4e\t", t);
		for (int j = 0; j < n; j++)
		{
			y_prev[j] = y[j];
			printf("%12.4e\t", y[j]);
		}
		printf("\n");
		if (i%10000==0)
		{
			fprintf(f_out, "%12.4e,%12.4e,%12.4e,%12.4e\n", t, y_prev[0], y_prev[1], y_prev[2]);
		}
	}

	// Check Residual
	targetFunc_ODE(n, y, y_prev, dh, G);
	for (int i = 0; i < n; i++)
	{
		printf("%12.4e\t", G[i]);
	}
	printf("\n");

	for (int i = 0; i < 3; i++)
	{
		free(Jacobian[i]);
	}
	free(Jacobian);
	fclose(f_out);
}