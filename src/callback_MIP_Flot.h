#pragma once

//======================================================================================
// fichier callback.h : ce fichier contient uniquement les callbacks pour le MIP MIP_flot.h
//
// ATTENTION : comme les fonctions callbacks sont des macros on doit mettre la déclaration
// ET la définition ensemble dans ce fichier 
//======================================================================================

#include "util.h"

#include <ilcplex/ilocplex.h>
#include <string>
ILOSTLBEGIN



// N = nb jobs (1...N) [0 = source, inutile ici]
//pop(i) = population du job i = 1..N (attention il faut un vecteur qui commence a 1 et donc decaler de 1 le vecteur _pop de l'instance _ins qui n'a pas de source
ILOLAZYCONSTRAINTCALLBACK7(TangenteLazyCallback, int, N, IloNumVarArray, start, IloNumVarArray, end, IloNumVarArray, rate, IloArray<IloNumVarArray>, x,
	IloNumArray, pop, long &, nbCut)
{

	IloEnv env = getEnv();

	//cout << " ----- TangenteLazyCallback ----- " << endl;


	//on stocke la solution courante
	IloNum startCour;
	IloNum endCour;
	IloNum rateCour;


	for (int i = 1; i <= N; ++i)
	{
		startCour = getValue(start[i]);
		endCour = getValue(end[i]);
		rateCour = getValue(rate[i]);

		if (endCour - startCour < pop[i] / rateCour - EPSILON_CUT_CPLEX)
		{
			double w = rateCour;
			double w2 = w * w;

			add(end[i] - start[i] + pop[i] * rate[i] / w2 >= 2 * pop[i] / w - EPSILON_CUT_CPLEX).end();
			nbCut++;
		}
	}

	//je pense que ce callback n'est appele que sur des sol entieres, on verifie par curiosite :
	//for (int i = 1; i <= N; ++i)
	//{
	//	for (int j = 1; j <= N; ++j)
	//	{
	//		if (i != j)
	//		{
	//			if (getValue(x[i][j]) > EPSILON && getValue(x[i][j]) < 1 - EPSILON)//IloNumIsInteger donne faux pour des valeurs < 10^-9
	//			{
	//				cout << "lazy : sol non entiere" << endl;
	//				cout << "x[" << i << "][" << j << "]=" << getValue(x[i][j]) << endl;
	//			}
	//		}
	//	}
	//}


	//cout << " ----- Fin TangenteLazyCallback ----- " << endl;

	return;

} 





ILOUSERCUTCALLBACK7(TangenteUserCallback, int, N, IloNumVarArray, start, IloNumVarArray, end, IloNumVarArray, rate, IloArray<IloNumVarArray>, x,
	IloNumArray, pop, long &, nbCut)
{

	IloEnv env = getEnv();


	//on stocke la solution courante
	IloNum startCour;
	IloNum endCour;
	IloNum rateCour;


	for (int i = 1; i <= N; ++i)
	{
		startCour = getValue(start[i]);
		endCour = getValue(end[i]);
		rateCour = getValue(rate[i]);

		if (endCour - startCour < pop[i] / rateCour - EPSILON_CUT_CPLEX)
		{
			double w = rateCour;
			double w2 = w * w;

			add(end[i] - start[i] + pop[i] * rate[i] / w2 >= 2 * pop[i] / w).end();
			nbCut++;
		}
	}

	//je pense que ce callback n'est appele que sur des sol fractionnaires, on verifie par curiosite :
	//bool entier = true;
	//for (int i = 1; i <= N; ++i)
	//{
	//	for (int j = 1; j <= N; ++j)
	//	{
	//		if (i != j)
	//		{
	//			if (getValue(x[i][j]) > EPSILON && getValue(x[i][j]) < 1 - EPSILON)//IloNumIsInteger donne faux pour des valeurs < 10^-9
	//			{
	//				entier = false;
	//			}
	//		}
	//	}
	//}
	//if (entier)
	//	cout << "user callback : sol entiere" << endl;



	return;


}
