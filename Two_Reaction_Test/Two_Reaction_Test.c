// Two_Reaction_Test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "globalPlasmaModel.h"
#include "stdlib.h"
#include "System.h"

double qin, L, R, k_pump, k2, p_abs, Ec1, m_arplus, me, ee;

int main()
{
	Test_BackwardEuler();
    return 0;
}

