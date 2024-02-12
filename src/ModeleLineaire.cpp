#include "ModeleLineaire.h"


//allocation des variables	
void ModeleLineaire::allocVar()
{
	int N = static_cast<int> (_expNet->_arc.size());

	_flot = IloNumVarArray(*_env, _expNet->_arc.size());

	for (int i = 0; i < N; ++i)
	{
		stringstream ss;
		ss << "F_" << _expNet->getArcName(i);
		_flot[i] = IloNumVar(*_env, 0, IloInfinity, IloNumVar::Float, ss.str().c_str());

	}//fin for i 

}


// si on veut mettre une fonction objectif on met obj a true,
// si obj = false c juste un pb de faisabilite

void ModeleLineaire::creationCtr_Obj(bool obj)
{

	//===========================================================================
	//===========================================================================

	int Nnodes = _ins->_nbNode;
	int T = _ins->_TMAX;

	IloExpr exprObj(*_env);


	//================================
	//1. objective function


	if (obj)
	{
		//vieil obj...
		//for (int i = _expNet->_firstRouteExpArc; i <= _expNet->_lastRouteExpArc; ++i)
		//{

		//	exprObj += _expNet->_arc[i]._cost * _flot[i];
		//}

		for (int i = _expNet->_firstRouteExpArc; i <= _expNet->_lastRouteExpArc; ++i)
		{
			int nodeIdOrig = _expNet->_arc[i]._orig._nodeIndex;
			int nodeIdDest = _expNet->_arc[i]._dest._nodeIndex;

			if (nodeIdDest >= _ins->_firstSafeNode && nodeIdDest <= _ins->_lastSafeNode)
			{
				int t = _expNet->_arc[i]._dest._time;
				exprObj += t * _flot[i];
			}
			
		}

		_mod->add(IloMinimize(*_env, exprObj));
	}



	//================================
	//2. constraint


	//================================
	//2.1. population

	for (int i = _expNet->_firstInputExpArc; i <= _expNet->_lastInputExpArc; ++i)
	{
		int nodeId = _expNet->_arc[i]._dest._nodeIndex;
		int population = _ins->_pop[nodeId];
		stringstream ss;
		ss << "pop_" << i;

		IloRange range(*_env, population, _flot[i], population, ss.str().c_str());
		_mod->add(range);
	}



	//================================
	//2.2. capacity safe node

	for (int i = _expNet->_firstOutputExpArc; i <= _expNet->_lastOutputExpArc; ++i)
	{
		int nodeId = _expNet->_arc[i]._orig._nodeIndex;
		int cap = _ins->_capNode[nodeId];

		stringstream ss;
		ss << "capSafe_" << i;

		IloRange range(*_env, _flot[i], cap, ss.str().c_str());

		_mod->add(range);

	}


	//================================
	//2.3. capacity transit node
	for (int i = _expNet->_firstStationExpArc; i <= _expNet->_lastStationExpArc; ++i)
	{
		int nodeId = _expNet->_arc[i]._orig._nodeIndex;
		int cap = _ins->_capNode[nodeId];


		stringstream ss;
		ss << "capTransit_" << i;

		IloRange range(*_env, _flot[i], cap, ss.str().c_str());

		_mod->add(range);

	}

	//================================
	//2.4. capacity route arc
	for (int i = _expNet->_firstRouteExpArc; i <= _expNet->_lastRouteExpArc; ++i)
	{
		int nodeIdOrig = _expNet->_arc[i]._orig._nodeIndex;
		int nodeIdDest = _expNet->_arc[i]._dest._nodeIndex;

		int cap = _ins->_capArc[nodeIdOrig][nodeIdDest];

		stringstream ss;
		ss << "capRoute_" << i;

		IloRange range(*_env, _flot[i], cap, ss.str().c_str());

		_mod->add(range);
	}

	//================================
	//2.5. conservation du flot

	for (size_t i = 1; i <= Nnodes; ++i)
	{

		for (size_t t = 0; t <= T; ++t)
		{
			size_t sizeIn = _expNet->_arcIn[i][t].size();
			size_t sizeOut = _expNet->_arcOut[i][t].size();

			IloExpr exprIn(*_env);
			IloExpr exprOut(*_env);


			if (sizeIn > 0 || sizeOut > 0)
			{
				if (sizeIn > 0)
				{
					for (size_t k = 0; k < sizeIn; ++k)
					{

						exprIn += _flot[_expNet->_arcIn[i][t][k]];

					}
				}

				if (sizeOut > 0)
				{
					for (size_t k = 0; k < sizeOut; ++k)
					{

						exprOut += _flot[_expNet->_arcOut[i][t][k]];
					}
				}

				stringstream ss;
				ss << "Kirch_" << i << "_" << t;

				IloRange range(*_env, 0, exprIn - exprOut, 0, ss.str().c_str());

				_mod->add(range);
			}
		}
	}

}


void ModeleLineaire::creationModele(bool obj)
{

	allocVar();




	creationCtr_Obj(obj);
}


double ModeleLineaire::resolutionIterative()
{
	double res = -1;

	//=================================================================================
	// FIXER LES PARAMETRES

	//_cplex->exportModel("modele.lp");

	_cplex->setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	_cplex->setParam(IloCplex::Param::ClockType, CLOCK_TYPE_CPLEX);
	_cplex->setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_CPLEX);
	_cplex->setParam(IloCplex::Param::MIP::Limits::TreeMemory, RAM_LIMIT_CPLEX);

	_cplex->setOut(_env->getNullStream());

	

	//=================================================================================
	// RESOLUTION

	vector<ExpArcValue>  solution;//derniere solution trouvee
	vector<ExpArcValue>  bestSol;//derniere solution trouvee
	int bestValMin = -1;


	// 1. resol du probleme init
	_cplex->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Network);
	_cplex->solve();

	// 2. tant qu on peut tuer des arcs, on tue et on reoptimize => on utilise le dual simplex

	_cplex->setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);

	while (_cplex->getStatus() == IloAlgorithm::Optimal)
	{
		solution.clear();

		//2.1 recherche de l arc a qui nous minimise dead - t

		int valMin = INT_INFINITY;
		int indArc = -1;
		
		//les deadline sont sur les arcs route
		for (int i = _expNet->_firstRouteExpArc; i <= _expNet->_lastRouteExpArc; ++i)
		{
			if (_cplex->isExtracted(_flot[i]) && _cplex->getValue(_flot[i]) > EPSILON)
			{
				ExpArcValue A = _expNet->getArcInfo(i);

				//on ne s interesse qu aux acrs transverses
				if ( (A._orig != A._dest) && (A._orig != 0) && (A._dest != _ins->_nbNode+1) )
				{
					int val = _ins->_dead[A._orig][A._dest] - A._timeOrig;

					if (val < valMin)
					{
						valMin = val;
						indArc = i;
					}

					//stocke sol car on s arretera sur une infaisabilite et la denirer sol faisable est la sol optimale

					A._value = _cplex->getValue(_flot[i]);
					//A.afficheSiRoute();
					solution.push_back(A);
				}
			
			}
		}//for i

		//on sauvegarde la solution si elle est meilleure (ça arrive que la solution se degrade : on interdit des arcs et Cplex peut en trouver d'autres plus mauvais)

		if (valMin > bestValMin)
			bestSol = solution;

		//2.2 ajout de la ctr


		//cout << "plus petite marge cplex = " << valMin << endl;

		//_expNet->traceSol(solution);
		//_expNet->afficheMaxMin(solution);

		ExpArcValue A = _expNet->getArcInfo(indArc);

		//on interdit tous les arcs transverses de la forme (i,t+k) -> (j,t+k+dist_ij), k >= 0
		interdireArcs(A._orig, A._dest, A._timeOrig);


		_cplex->solve();

		//cout << _cplex->getStatus() << endl;

	}//fin while

	//=================================================================================
	// AFFICHAGE SOLUTION ET COUT


	if (!bestSol.empty())
	{
		_expNet->traceSol(bestSol);

		res = _expNet->calculMaxMin(bestSol);
	}


	return res;
}

//on interdit tous les arcs transverses de la forme (i,t+k) -> (j,t+k+dist_ij), k >= 0
void ModeleLineaire::interdireArcs(int i, int j, int t)
{
	//il faut parcourir tous les arcs (c est couteux => faire un stockage plus malin !!)

	for (int k = _expNet->_firstRouteExpArc; k <= _expNet->_lastRouteExpArc; ++k)
	{
		const ExpandedNetwork::ExpArc & A = _expNet->_arc[k];
			
		if ( A._orig._nodeIndex == i && A._dest._nodeIndex == j && A._orig._time >= t)
			_mod->add(_flot[k] == 0);
	}
}

double ModeleLineaire::resolutionSimple()
{
	double res = -1;
	//=================================================================================
	// FIXER LES PARAMETRES

	_cplex->exportModel("modele.lp");

	_cplex->setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	_cplex->setParam(IloCplex::Param::ClockType, CLOCK_TYPE_CPLEX);
	_cplex->setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_CPLEX);
	_cplex->setParam(IloCplex::Param::MIP::Limits::TreeMemory, RAM_LIMIT_CPLEX);

	_cplex->setOut(_env->getNullStream());



	//=================================================================================
	// RESOLUTION

	vector<ExpArcValue>  solution;//derniere solution trouvee

	// 1. resol du probleme init
	_cplex->solve();

	//solution
	if (_cplex->getStatus() == IloAlgorithm::Optimal)
	{
		res = _cplex->getObjValue();

		for (int i = 0; i < _expNet->_arc.size(); ++i)
		{
			if (_cplex->isExtracted(_flot[i]) && _cplex->getValue(_flot[i]) > EPSILON)
			{
				ExpArcValue A = _expNet->getArcInfo(i);
				A._value = _cplex->getValue(_flot[i]);
				solution.push_back(A);

			}
		}//for i

	}//fin if

	//=================================================================================
	// AFFICHAGE SOLUTION ET COUT


	if (!solution.empty())
	{
		_expNet->traceSol(solution);

	}


	return res;
}

double ModeleLineaire::creeEtResout(bool obj)
{
	creationModele(obj);
	
	return resolutionIterative();
	//return resolutionSimple();
}