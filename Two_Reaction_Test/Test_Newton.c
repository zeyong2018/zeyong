
#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "globalPlasmaModel.h"
#include "stdlib.h"
#include "System.h"

void Test_Newton(void)
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

	//double y[3] = { 2.1239e20, 2.1092e17, 2.9294 };
	double y[3] = { 2.1239e20*1.3, 2.1092e17*1.3, 2.9294*0.8 };
	double y_temp[3] = { 0,0,0 };
	double **Jacobian;
	int n = 3;
	Jacobian = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++)
	{
		Jacobian[i] = (double *)malloc(n * sizeof(double));
	}
	double F[3] = { 0,0,0 };

	int N_MAX = 1000;
	for (int i = 0; i < N_MAX; i++)
	{
		targetFunc(n, y, F);
		Jacobian_F(n, y, Jacobian);
		gaussElimination(n, Jacobian, F, y_temp);
		for (int j = 0; j < n; j++)
		{
			y[j] -= y_temp[j];
			printf("%12.4e\t", y[j]);
		}
		printf("\n");
	}
	targetFunc(n, y, F);
	for (int i = 0; i < n; i++)
	{
		printf("%12.4e\t", F[i]);
	}
	printf("\n");

	free(Jacobian);
}