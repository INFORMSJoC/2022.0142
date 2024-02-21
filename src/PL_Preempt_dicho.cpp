#include "PL_Preempt_dicho.h"


#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


double PL_Preempt_dicho::run(double precision)
{

	double bestMargeCible = -1;


	double margeCibleMin = 0;
	double margeCibleMax = calculMargeCibleMax();


	vector<double> date;

	//boucle principale dicho
	while (abs(margeCibleMax - margeCibleMin) > precision)
	{
		double margeCible = (margeCibleMin + margeCibleMax) / 2.0;

		//tri les dates ES _ LFi-mergecible
		date = trierDate(margeCible);

		//lance un PL; retourne true si margeCible est realisable
		bool ok = isMargeRealisable(margeCible, date);

		if (ok)
		{
			margeCibleMin = margeCible;
			if (margeCible > bestMargeCible)
				bestMargeCible = margeCible;

		}
		else
			margeCibleMax = margeCible;
	}

	return bestMargeCible;
}



//calcul la meilleure marge possible compte tenu des population
double PL_Preempt_dicho::calculMargeCibleMax()
{
	double margeMax = TIME_INFINITY;

	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		double m = _ins->_LF[i] - (_ins->_ES[i] + static_cast<double>(_ins->_pop[i]) / _ins->_debitMax[i]);
		if (margeMax > m)
			margeMax = m;
	}

	return margeMax;
}


//on ordonne l'ensemble des dates ESi, LFi-margeCible par ordre croissant
//renvoie le vecteur des valeurs triees
vector<double> PL_Preempt_dicho::trierDate(double margeCible)
{
	vector<double> date;
	date.reserve(_ins->_nbJob * 2);

	//ajoute le vecteur _ES dans date
	copy(_ins->_ES.begin(), _ins->_ES.end(), back_inserter(date));
	transform(_ins->_LF.begin(), _ins->_LF.end(), back_inserter(date), [&](double lf) {return lf - margeCible; });

	sort(date.begin(), date.end());
	auto it = unique(date.begin(), date.end());

	date.resize(distance(date.begin(), it));

	return date;
}


//on lance un PL pour savoir si toutes les pop peuvent evacuer avec une marge >= margeCible  
bool PL_Preempt_dicho::isMargeRealisable(double margeCible, vector<double> & date)
{
	bool res = true;
	int nbDate = static_cast<int> (date.size());

	//==========================================================
	//env

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	//==========================================================
	//variables

	IloArray<IloNumVarArray> debit(env, _ins->_nbJob); //debit[i][k] = debit du job i entre les dates date[k] et date[k+1]

	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		debit[i] = IloNumVarArray(env, nbDate - 1);

		for (int k = 0; k < nbDate - 1; ++k)
		{
			stringstream ss;
			ss << "d_" << i << "_" << k;
			debit[i][k] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, ss.str().c_str());
		}
	}//fin for i 

	//==========================================================
	//contraintes

	//1. le debit ne depasse pas le debit max
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		for (int k = 0; k < nbDate - 1; ++k)
		{
			model.add(debit[i][k] <= _ins->_debitMax[i]);
		}

	}//fin for i 


	//2. prise en compte des fenetres de temps [ES,LF-marge]
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		for (int k = 0; k < nbDate - 1; ++k)
		{
			if (_ins->_ES[i] >= date[k + 1] - EPSILON || _ins->_LF[i] - margeCible <= date[k] + EPSILON)
				model.add(debit[i][k] == 0);
		}

	}//fin for i 

	//3. on evacue tout le monde
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		IloExpr expr(env);
		for (int k = 0; k < nbDate - 1; ++k)
		{
			expr += debit[i][k] * (date[k + 1] - date[k]);
		}

		
		model.add(expr == _ins->_pop[i]);
		expr.end();
	}//fin for i 



	for (int a = 0; a < _ins->_nbArc; ++a)
	{

		int nbJobSurA = static_cast<int> (_ins->_arcToJob[a].size());

		for (int k = 0; k < nbDate - 1; ++k)
		{

			IloExpr expr(env);

			for (int j = 0; j < nbJobSurA; ++j)
			{
				int i = _ins->_arcToJob[a][j];
				expr += debit[i][k];
			}
			model.add(expr <= _ins->_cap[a]);
			expr.end();
		}


	}//fin for a

	//===========================================================================
	// resolution


	cplex.setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);

	cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 0.00001);
	//cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 0.0001);
	cplex.setOut(env.getNullStream());

	//cplex.exportModel("PL_dicho.lp");
	cplex.solve();

	//cout << "marge = " << margeCible << " "<< cplex.getStatus() << endl;

	if (cplex.getStatus() == IloAlgorithm::Infeasible)
		res = false;

	cplex.end();
	model.end();
	env.end();

	return res;
}