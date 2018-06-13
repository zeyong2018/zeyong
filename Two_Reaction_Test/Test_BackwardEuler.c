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

	double n1, n2, t0;
	n1 = 7.1728e18;
	n2 = 7.1728e10;
	t0 = 2.;
	double y_prev[3] = { n1,n2,t0};
	double y_temp[3] = { 0,0,0 };
	double y_iter_old[3] = { 0.,0,0 };
	double **Jacobian;
	int n = 3;
	Jacobian = dmatrix(3, 3);
	double G[3] = { 0,0,0 };

	double dh=1.e-6;
	int N_STEP = ceil(5. / dh);
	FILE *f_out;
	f_out = fopen("Output.csv", "w");
	fprintf(f_out, "%12s,%12s,%12s,%12s\n", "t", "n_ar", "n_ar+", "Te");
	double y[3] = { n1,n2,t0 }; // initial guess for first step
	int N_MAX = 1000;
	double t,sum;
	fprintf(f_out, "%12.4e,%12.4e,%12.4e,%12.4e\n", 0., y_prev[0], y_prev[1], y_prev[2]);
	double Te;
	for (int i = 0; i < N_STEP; i++)
	{
		t = (i + 1)*dh;
		for (int j = 0; j < N_MAX; j++)
		{
			printf("=====================================\n");
			//Jacobian_F(n, y, Jacobian);
			targetFunc_ODE(n, y, y_prev, dh, G);
			Jacobian_G(n, y, dh, Jacobian);
			Te = y[2];
			printf("K for i=%d, j=%d\n", i, j);
			printf("%12.4e\t%12.4e\n", (7.93e-14)*exp(-18.9 / (Te + TINY)), 2.6582E4);
			printf("\n");
			printf("DKDTE for i=%d, j=%d\n", i, j);
			printf("%12.4e\t%12.4e\n", (7.93e-14)*exp(-18.9 / (Te + TINY))*18.9 / pow(Te + TINY, 2.), 0.);
			printf("\n");
			printf("KE for i=%d, j=%d\n", i, j);
			printf("%12.4e\t%12.4e\n", (7.93e-14)*exp(-18.9 / (Te + TINY))*57.7141, 2.6582*(2. + 1. / 2 * (1 + log(m_arplus / 2 / PI / me)))*Te);
			printf("\n");
			printf("DKEDTE for i=%d, j=%d\n", i, j);
			printf("%12.4e\t%12.4e\n", (7.93e-14)*exp(-18.9 / (Te + TINY))*18.9 / pow(Te + TINY, 2.)*57.7141, 2.6582*(2. + 1. / 2 * (1 + log(m_arplus / 2 / PI / me))));
			printf("\n");
			

			printf("Jacobian for i=%d, j=%d\n", i, j);
			for (int k = 0; k < 3; k++)
			{
				for (int l = 0; l < 3; l++)
				{
					printf("%12.4e\t", Jacobian[k][l]);
				}
				printf("\n");
			}
			printf("G for i=%d, j=%d\n", i, j);
			for (int k = 0; k < 3; k++)
			{
				printf("%12.4e\t", G[k]);
			}
			printf("\n");
			gaussElimination(n, Jacobian, G, y_temp);
			
			printf("y_temp for i=%d, j=%d\n", i, j);
			for (int k = 0; k < 3; k++)
			{
				printf("%12.4e\t", y_temp[k]);
			}
			printf("\n");


			sum = 0.;
			for (int k = 0; k < n; k++)
			{
				y_iter_old[k] = y[k];
				y[k] -= y_temp[k];
				printf("y_new=%12.4e\t,y_old=%12.4e\t,diff=%12.4e\n", y[k], y_iter_old[k], (y[k] - y_iter_old[k]) / (y[k] + TINY));
				sum += fabs((y[k] - y_iter_old[k]) / (y[k] + TINY));
				printf("sum=%12.4e\n", sum);
			}
			printf("y for i=%d, j=%d\n", i, j);
			for (int k = 0; k < 3; k++)
			{
				printf("%12.4e\t", y[k]);
			}
			printf("sum=%12.4e\n", sum);
			printf("\n");
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