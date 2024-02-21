#include "PL_RCPSP_fixeDead.h"
#include <queue>


//allocation des variables	
void PL_RCPSP_fixeDead::allocVar()
{
	//========================================================
	// 1. alloc : (_nbSommet + 2) * (_nbSommet + 2) * _nbArc
	_flot = NumVarArray3(*_env, _nbSommet + 2);//+2 = source et puits

	for (int i = 0; i <= _nbSommet + 1; ++i)
	{

		_flot[i] = NumVarArray2(*_env, _nbSommet + 2);//+2 = source et puits

		for (int j = 0; j <= _nbSommet + 1; ++j)
		{
			_flot[i][j] = IloNumVarArray(*_env, _nbArc);

			for (int e = 0; e < _nbArc; ++e)
			{
				stringstream ss;
				ss << "F_" << i << "_" << j << "_" << e;
				_flot[i][j][e] = IloNumVar(*_env, 0, IloInfinity, IloNumVar::Float, ss.str().c_str());
			}
		}

	}//fin for i 

	//========================================================

	_debit = IloNumVarArray(*_env, _nbSommet + 1);

	for (int i = 1; i <= _nbSommet; ++i)
	{
		stringstream ss;
		ss << "D_" << i;
		_debit[i] = IloNumVar(*_env, 0, IloInfinity, IloNumVar::Float, ss.str().c_str());
	}

}


//flot autorise de i vers j uniquement si i << j dans ordre
//ordre ne doit pas contenir source et puits
void PL_RCPSP_fixeDead::creationModele( vector<int> & ordre )
{

	int puits = _nbSommet + 1;


	//1. kirchhoff

	for (int i = 1; i <= _graph->_nbSommet; ++i)//pour chaque vraie action
	{
		
		//pour chaque arc du chemin de i
		int prec = _ins->_chemins[i][0]._som;
		for (int k = 1; k < _ins->_chemins[i].size(); ++k)
		{

			IloExpr exprIn(*_env);
			IloExpr exprOut(*_env);


			int cour = _ins->_chemins[i][k]._som;
			int e = _graph->_arc2Ind[prec][cour];

			//pour chaque j qui contient e dans son chemin
			for (int j = 1; j <= _graph->_nbSommet ; ++j)//pour chaque action+source+puits
			{
				// source et puits sont forcement dans le chemin
				if (i!=j &&  _ins->appartientChemin(j, prec, cour) >= 0 )//retourn -1 si prec,cour pas dans chemin j
				{
					exprIn += _flot[j][i][e];
					exprOut += _flot[i][j][e];
				}
			}

			exprIn += _flot[0][i][e]; //source entre
			exprOut += _flot[i][_graph->_nbSommet + 1][e]; // puits


			prec = cour;

			stringstream ss, ss2;
			ss << "kirch.IN_" << i << "_" << e ;
			ss2 << "kirch.OUT_" << i << "_" << e;

			IloRange range(*_env, 0, exprIn - _debit[i], 0, ss.str().c_str());
			_mod->add(range);

			IloRange range2(*_env, 0, exprOut - _debit[i], 0, ss2.str().c_str());
			_mod->add(range2);

		}//fin for k 

	}

	//2. capacite

	for (int e = 0; e < _nbArc; ++e)
	{
		IloExpr exprIn(*_env);
		IloExpr exprOut(*_env);

		stringstream ss;
		ss << "CAP_" << e;

		int orig = _graph->_arcs[e]._orig;
		int dest = _graph->_arcs[e]._dest;

		for (int i = 1; i <= _graph->_nbSommet; ++i)
		{
			if (_ins->appartientChemin(i, orig, dest) >= 0) //retourne -1 si (prec,cour) pas dans chemin j
			{
				exprIn += _flot[0][i][e];
				exprOut += _flot[i][puits][e];
			}
		}
		exprIn += _flot[0][puits][e];
		exprOut += _flot[0][puits][e];


		IloRange range(*_env, _ins->_capArc[orig][dest], exprIn, _ins->_capArc[orig][dest], ss.str().c_str());
		_mod->add(range);//ce qui sort de la source = capa arc

		_mod->add(exprIn - exprOut == 0);//facultatif ? : ce qui entre au puits = ce qui sort de la source

	}


	//3. debit  max 
	//attention il ne faut pas imposer de débit min car on le debit min est calcule avec les deadlines actuelles
	//qu on va changer

	
	for (int i = 1; i <= _nbSommet; ++i)
	{
		stringstream ss;
		ss << "debitMax" << i ;

		double debitMax = _ins->getDebitMax(i);

		IloRange range(*_env, 0,  _debit[i], debitMax, ss.str().c_str());
		_mod->add(range);
	}


	//4. ordre

	for (int e = 0; e < _nbArc; ++e)
	{
		for (int i = 0; i < _nbSommet; ++i)
		{
			int som1 = ordre[i];
			_mod->add(_flot[som1][0][e] == 0);//source avant tout le monde
			_mod->add(_flot[puits][som1][e] == 0);//puits ap tout le monde
			
			for (int j = i+1; j < _nbSommet; ++j)
			{
				int som2 = ordre[j];
				_mod->add(_flot[som2][som1][e] == 0);//flot dans le sens som1 -> som2 autorise uniquement

			}

		}
	}


	//objectif
	IloExpr expr(*_env);

	for (int i = 1; i <= _nbSommet; ++i)
	{
		expr += _debit[i] * _ins->_pop[i];
	}
	

	
	_mod->add(IloMaximize(*_env, expr));
	expr.end();
}



//cree le modele et resout
bool PL_RCPSP_fixeDead::creeEtResout()
{
	bool ok = false;

	//===============================================
	// construction ordre alea (sans source ni puits)
	
	vector<int> ordre;

	/*for (int i = 1; i <= _nbSommet; ++i)
		ordre[i - 1] = i;

	random_shuffle(ordre.begin(), ordre.end());
	*/

	//on met sigma dans l'ordre des plus petits chemins
	vector < pair<double, int> > ordre2;
	for (int i = 1; i <= _nbSommet; ++i)
	{
		double l = _ins->getLongueurCh(i);
		ordre2.push_back({l,i});
	}
	sort(ordre2.begin(), ordre2.end());

	for (int i = 0; i < _nbSommet; ++i)
	{
		
		ordre.push_back(ordre2[i].second);
	}

	//=================================================================================
	// CREATION MODELE	

	allocVar();
	creationModele(ordre);

	//=================================================================================
	// FIXER LES PARAMETRES

	_cplex->exportModel("PL_fixeDead.lp");

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

void PL_RCPSP_fixeDead::afficheSol()
{
	cout << "obj PL FIX = " << _cplex->getObjValue() << endl;

	for (int i = 0; i <= _nbSommet + 1; ++i)
	{
		for (int j = 0; j <= _nbSommet + 1; ++j)
		{
			for (int e = 0; e < _nbArc; ++e)
			{
				if (_cplex->isExtracted(_flot[i][j][e]) && _cplex->getValue(_flot[i][j][e]))
				{
					cout << "flot[" << i <<"," << j << "," << e << "]" << " = " << _cplex->getValue(_flot[i][j][e]) << endl;
				}
			}
		}

	}//fin for i 


	

	for ( int i = 1; i <= _nbSommet; ++i )
	{
		if ( _cplex->isExtracted(_debit[i]) && _cplex->getValue(_debit[i]) )
		{
			cout << "debit[" << i << "]" << " = " << _cplex->getValue(_debit[i]) << endl;
		}
	}

}

void PL_RCPSP_fixeDead::calculeDebutEtFin( )
{
	
	int S = static_cast<int>(_graph->_arcs.size());

	queue<int> file;
	file.push(0);//source

	//==========================================================
	// dates de debut

	//1. reinit les debuts a 0
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		_debut[i] = -1;
		_fin[i] = -1;
	}


	//2. boucle tant qu on augmente des dates 

	while (!file.empty())
	{
		int cour = file.front();
		file.pop();

		//on cherche les suivants : ceux a qui cour transmet du flot
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			if ( existeFlot(cour,j))//cour -> j 
			{

				double deb = 0;
				if (cour != 0)
					deb = _debut[cour] + _graph->getTLcond(cour, j, _cplex->getValue(_debit[cour]));

				if (deb > _debut[j] + EPSILON)
				{
					_debut[j] = deb;
					file.push(j);
				}
			}
		}

	}//fin while


	//==========================================================
	// dates de fin

	for (int i = 1; i <= _nbSommet; ++i)
	{
		_fin[i] = _debut[i] + _ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / _cplex->getValue(_debit[i]) - 1;
	}


}

//renvoie vrai si flot de i vers j pour au moins une ressource
bool PL_RCPSP_fixeDead::existeFlot(int i, int j)
{
	bool exist = false;

	for (int e = 0; e < _nbArc; ++e)
	{
		if (_cplex->isExtracted(_flot[i][j][e]) &&  _cplex->getValue(_flot[i][j][e]) > EPSILON)
			exist = true;
	}

	return exist;
}

//dessine le graphe uniquement avec les arcs dont le flot pour la ressource arc e est non nulle
void PL_RCPSP_fixeDead::traceFlotArc(int e)
{
	string name = "flot";

	ofstream ficOut("flot.txt", ios::out);// ouverture fichier de sortie

	//premiere ligne du fichier de sortie
	ficOut << "digraph G {" << endl;

	//graphe horizontal :
	ficOut << "rankdir=LR" << endl;

	for (int i = 0; i <= _nbSommet + 1; ++i)
	{
		for (int j = 0; j <= _nbSommet + 1; ++j)
		{

			if (_cplex->isExtracted(_flot[i][j][e]) && _cplex->getValue(_flot[i][j][e]) > EPSILON)
			{
				ficOut << i << " -> " << j << "[label=" << _cplex->getValue(_flot[i][j][e]) << "]" << endl;
			}
		}

	}//fin for i 



	ficOut << "}" << endl;
	ficOut.close();

	stringstream ligCom;


#if defined(_WIN32) || defined(_WIN64)
	ligCom << "dot.exe -Tpdf -o" << name << ".pdf " << name << ".txt";
#endif

	system(ligCom.str().c_str());
}