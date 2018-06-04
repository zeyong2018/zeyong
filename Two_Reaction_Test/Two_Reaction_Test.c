// Two_Reaction_Test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "globalPlasmaModel.h"
#include "stdlib.h"
#include "System.h"
#define _CRT_SECURE_NO_WARNINGS 

double qin, L, R, k_pump, k2, p_abs, Ec1, m_arplus, me, ee;

int main()
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

	double y_prev[3] = { 2.1239e20, 2.1092e17, 2.9294 };
	double y_temp[3] = { 0,0,0 };
	double **Jacobian;
	int n = 3;
	Jacobian = dmatrix(3, 3);
	double G[3] = { 0,0,0 };

	double dh = 1.e-6;

	double y[3] = { 2.1239e20, 2.1092e17, 2.9294 }; // initial guess for first step
	int N_MAX = 1000;
	double t;
	for (int i = 0; i < 200000; i++)
	{
		t = (i + 1)*dh;
		for (int j = 0; j < N_MAX; j++)
		{
			targetFunc_ODE(n, y, y_prev, dh, G);
			Jacobian_G(n, y, dh, Jacobian);
			gaussElimination(n, Jacobian, G, y_temp);
			for (int k = 0; k < n; k++)
			{
				y[k] -= y_temp[k];
			}
			
		}
		printf("%12.4e\t", t);
		for (int j = 0; j < n; j++)
		{
			y_prev[j] = y[j];
			printf("%12.4e\t", y[j]);
		}
		printf("\n");

	}

	// Check Residual
	targetFunc_ODE(n, y, y_prev, dh, G);
	for (int i = 0; i < n; i++)
	{
		printf("%12.4e\t", G[i]);
	}
	printf("\n");

	free(Jacobian);

    return 0;
}

