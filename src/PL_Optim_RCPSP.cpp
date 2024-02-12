#include "PL_Optim_RCPSP.h"


//allocation des variables	

void PL_Optim_RCPSP::allocVar(const vector<pair<int, int> > & arcActif)
{
	
	//========================================================
	// 1. alloc : (_nbSommet + 2) * (_nbSommet + 2) * _nbArc
	_deltaFlot = NumVarArray3(*_env, _nbSommet + 2);//+2 = source et puits

	for (int i = 0; i <= _nbSommet + 1; ++i)
	{

		_deltaFlot[i] = NumVarArray2(*_env, _nbSommet + 2);//+2 = source et puits

		for (int j = 0; j <= _nbSommet + 1; ++j)
		{
			_deltaFlot[i][j] = IloNumVarArray(*_env, _nbArc);

			//on ne cree que les deltaFlot[i][j][e] si e est une ressource de i et j (i.e. e appartient au chemin de i et j)
			for (int k = 0; k < _graph->_listeArcPartage[i][j].size(); ++k)
			{
				int o = _graph->_listeArcPartage[i][j][k]._orig;
				int d = _graph->_listeArcPartage[i][j][k]._dest;

				int e = _graph->_arc2Ind[o][d];

				stringstream ss;
				ss << "F_" << i << "_" << j << "_" << e;
				_deltaFlot[i][j][e] = IloNumVar(*_env, -IloInfinity, IloInfinity, IloNumVar::Float, ss.str().c_str());

			}
		}

	}//fin for i 



    //========================================================
	
	_deltaDebit = IloNumVarArray(*_env, _nbSommet + 1);

	for (int i = 1; i <= _nbSommet; ++i)
	{
		stringstream ss;
		ss << "D_" << i;
		_deltaDebit[i] = IloNumVar(*_env, -IloInfinity, IloInfinity, IloNumVar::Float, ss.str().c_str());
	}

}


//on cree le modele sachant que les seuls arcs actifs sont ceux stockes dans arcActif : 
//le flot sera nul pour les autres arcs (les variables sont crees par commodite mais non incluses dans les ctr)
void PL_Optim_RCPSP::creationModele(const vector< double> & debitMin, const vector<pair<int, int> > & arcActif, const vector<double> & lambda)
{
	
	//vecteur pour savoir si un couple est dans arcActif
	vector<vector<bool> > isActif(_nbSommet + 2, vector<bool>(_nbSommet + 2, false));

	//1. la variation de flot ne doit pas rendre le flot négatif

	for (int k = 0; k < arcActif.size(); ++k)
	{
		int i = arcActif[k].first;
		int j = arcActif[k].second;
		isActif[i][j] = true; //init isActif ici 


		for (int k = 0; k < _graph->_listeArcPartage[i][j].size(); ++k)
		{
			int o = _graph->_listeArcPartage[i][j][k]._orig;
			int d = _graph->_listeArcPartage[i][j][k]._dest;

			int e = _graph->_arc2Ind[o][d];

			stringstream ss;
			ss << "flot_" << i << "_" << j << "_" << _graph->_arcs[e]._orig << "_" << _graph->_arcs[e]._dest;

			
			IloRange range(*_env, - _rcpspAlgo->_flot[i][j][e], _deltaFlot[i][j][e] , +IloInfinity, ss.str().c_str());
			_mod->add(range);
			
		}
	}

	//2. debit = somme des flots pour chaque arc e

	for (int i = 1; i <= _nbSommet; ++i)//kirchhoff pour tout sommet sauf source et puits
	{
		//pour chaque e dans le chemin de i

		for (int e = 0; e < _nbArc; ++e)
		{
			int o = _graph->_arcs[e]._orig;
			int d = _graph->_arcs[e]._dest;


			if (_ins->appartientChemin(i, o, d) >= 0)
			{
				IloExpr exprIn(*_env);
				IloExpr exprOut(*_env);

				for (int j = 0; j <= _nbSommet + 1; ++j)
				{
					if (_graph->isDansArcPartage(i, j, o, d))
					{
						if (isActif[j][i])
							exprIn += _deltaFlot[j][i][e];

						if (isActif[i][j])
							exprOut += _deltaFlot[i][j][e];
					}
				}

				stringstream ss1;
				ss1 << "Kirchhoff_" << i << "_" << e;
				IloRange range1(*_env, 0, exprIn - exprOut, 0, ss1.str().c_str());
				_mod->add(range1);

				stringstream ss2;
				ss2 << "flotEtDebit" << i << "_" << e;
				IloRange range2(*_env, 0, exprIn - _deltaDebit[i], 0, ss2.str().c_str());
				_mod->add(range2);

			}


		}
	}


	//3. pour chaque arc e : pas de creation de flot


	for (int e = 0; e < _nbArc; ++e)
	{

		IloExpr exprOut(*_env);

		//on parcourt tous les jobs j qui ont e comme ressource (i.e. e appartient au chemin de j)
		for (int j = 1; j <= _nbSommet; ++j)
		{
			int o = _graph->_arcs[e]._orig;
			int d = _graph->_arcs[e]._dest;

		
			if ( _ins->appartientChemin(j, o, d) >= 0 && isActif[0][j] )
				exprOut += _deltaFlot[0][j][e];
		}

		//on ajoute le puits
		exprOut += _deltaFlot[0][_nbSommet+1][e];

		stringstream sse;
		sse << "source" << "_" << e;

		IloRange range1(*_env, 0, exprOut, 0, sse.str().c_str());
		_mod->add(range1);


		//IloRange range2(*_env, 0, exprIn, 0, "puits");
		//_mod->add(range2);
	}


	//4. debit borne
	for (int i = 1; i <= _nbSommet; ++i)
	{
		stringstream ss;
		ss << "debit_" << i;
		IloRange range(*_env, debitMin[i] - _rcpspAlgo->_debit[i], _deltaDebit[i] , _ins->getDebitMax(i) - _rcpspAlgo->_debit[i], ss.str().c_str());
		_mod->add(range);
	}

	//5. on encadre la quantite "gradient projete" par le pas
	// 7/05/2019 => on max au lieu de d'imposer un pas
	IloExpr expr(*_env);
	
	for (int i = 1; i <= _nbSommet; ++i)
	{
		double debitCarre = _rcpspAlgo->_debit[i] * _rcpspAlgo->_debit[i];
		expr += _deltaDebit[i] * lambda[i] * _ins->_pop[i]  / debitCarre; 
	}

	//IloRange range(*_env, pas, expr, 2* pas, "pas" );
	//IloRange range(*_env, pas, expr, +IloInfinity, "pas");
	//_mod->add(range);

	_mod->add(IloMaximize(*_env,expr));
}


//cree le modele et resout
bool PL_Optim_RCPSP::creeEtResout(const vector<double> & debitMin, const vector<pair<int, int> > & arcActif,const vector<double> & lambda)
{
	
	bool ok = false;

	//=================================================================================
	// CREATION MODELE	

	allocVar(arcActif);
	creationModele(debitMin, arcActif,  lambda);



	//=================================================================================
	// FIXER LES PARAMETRES

	//_cplex->exportModel("PL_rcpsp.lp");

	_cplex->setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	_cplex->setParam(IloCplex::Param::ClockType, CLOCK_TYPE_CPLEX);
	_cplex->setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_CPLEX);
	_cplex->setParam(IloCplex::Param::MIP::Limits::TreeMemory, RAM_LIMIT_CPLEX);

	_cplex->setOut(_env->getNullStream());


	//=================================================================================
	// RESOLUTION

	// 1. resol du probleme init
	
	_cplex->solve();


	if (_cplex->getCplexStatus() == IloCplex::CplexStatus::Optimal || _cplex->getCplexStatus() == IloCplex::CplexStatus::Feasible)
		ok = true;


	return ok;
} 


//affiche les var non nulles
void PL_Optim_RCPSP::afficheSol( )
{
	double val;

	cout << "PL_Optim_RCPSP : " << endl;

	//on n'affiche pas avec source et puits sinon c'est illisible
	for (int i = 0; i <= _nbSommet ; ++i)
	{
		for (int j = 1; j <= _nbSommet+1 ; ++j)
		{
			for (int k = 0; k < _graph->_listeArcPartage[i][j].size(); ++k)
			{
				int o = _graph->_listeArcPartage[i][j][k]._orig;
				int d = _graph->_listeArcPartage[i][j][k]._dest;
				int e = _graph->_arc2Ind[o][d];

				val = getDeltaFlot(i, j, e);

				/*
				if ( val > EPSILON || val < -EPSILON ) 
					cout << "deltaFlot(" << i << "," << j << "," << o << " -> " << d << ") = " << val << endl;*/
			}
		}
	}//fin for i 


	//========================================================

	int cpt = 0;
	for (int i = 1; i <= _nbSommet; ++i)
	{
		val = getDeltaDebit(i);
		if (val > EPSILON || val < -EPSILON)
		{
			cout << "deltaDebit(" << i << ") = " << val << endl;
			cpt++;
		}
	}
	if (cpt == 0)
		cout << "pas de modif de debit !! " << endl;

}