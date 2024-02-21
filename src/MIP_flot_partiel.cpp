#include "MIP_flot_partiel.h"


//retourne borne inf, borne sup
double MIP_flot_partiel::creeEtResout(const vector<pair<int, int>> & arcSupport)
{
	double res = -1;
	
	vector<double> startCour(_ins->_nbJob + 1);
	vector<double> endCour(_ins->_nbJob + 1);
	vector<double> rateCour(_ins->_nbJob + 1);

	allocVar();
	creationCtrLineaire(arcSupport);
	creationFctObj();

	creationCtrTangente(true, 5);

	
	_cplex->setOut(_env->getNullStream());
	//_cplex->setParam(IloCplex::Param::Simplex::Display, 2);
	//_cplex->exportModel("mod.lp");

	_cplex->solve();

	

	bool fin = false;
	_cplex->setParam(IloCplex::RootAlg, IloCplex::Dual);

	while (!fin && _cplex->getStatus() == IloAlgorithm::Status::Optimal)
	{
		fin = true; //fin passe a faux si on ajoute des ctr violees
		
		// 1. on recupere les valeurs (on ne peut plus recuperer la valeur d'une variable apres l'appel a _mod->add
		for (int i = 1; i <= _ins->_nbJob; ++i)
		{
			if (_existe[i - 1])
			{
				startCour[i] = _cplex->getValue(_start[i]);
				endCour[i] = _cplex->getValue(_end[i]);
				rateCour[i] = _cplex->getValue(_rate[i]);
			}
		}

		// 2. RECHERCHE DES CTR TANGENTES VIOLEES
		for (int i = 1; i <= _ins->_nbJob; ++i)
		{
			if (_existe[i - 1])
			{

				if (endCour[i] - startCour[i] < (double)(_ins->_pop[i-1]) / rateCour[i] - EPSILON)
				{
					double w = rateCour[i];
					double w2 = w * w;

					_mod->add(_end[i] - _start[i] + _ins->_pop[i - 1] * _rate[i] / w2 >= 2 * (double)(_ins->_pop[i - 1]) / w);
					fin = false;
				}
			}
		}
		if (!fin)
		{

			_cplex->solve();
			//_cplex->writeSolution("mipPartiel.sol");
			
		}
	}

	if (_cplex->getStatus() == IloAlgorithm::Status::Optimal)
	{
		res = _cplex->getObjValue();

		if (!verifSol())
			stopProg("MIP_flot_partiel : sol partiel non realisable");
	}
	//afficheSol();

	//_cplex->writeSolution("mipPartiel.sol");
	return res;
}


//allocation des variables
void MIP_flot_partiel::allocVar()
{
	int N = _ins->_nbJob + 2;//+2 => on ajoute source et puits
	int puits = N - 1;
	int K = _ins->_nbArc;

	//1. allocation des variables flots
	_flot = IloArray< IloArray< IloNumVarArray > >(*_env, N);

	for (int i = 0; i < N; ++i)
	{
		_flot[i] = IloArray<IloNumVarArray>(*_env, N);

		for (int j = 0; j < N; ++j)
		{
			_flot[i][j] = IloNumVarArray(*_env, K);

			for (int k = 0; k < K; ++k)
			{
				stringstream ss;
				ss << "f_" << i << "_" << j << "_" << k;
				_flot[i][j][k] = IloNumVar(*_env, 0, IloInfinity, IloNumVar::Float, ss.str().c_str());
			}
		}

	}//fin for i 


	//2. allocation des tableaux 1D : debut / fin / vitesse
	_start = IloNumVarArray(*_env, N);
	for (int i = 1; i < puits; ++i)
	{
		stringstream ss;
		ss << "S_" << i;
		_start[i] = IloNumVar(*_env, _ins->_ES[i - 1], IloInfinity, IloNumVar::Float, ss.str().c_str());
	}

	_end = IloNumVarArray(*_env, N);
	for (int i = 1; i < puits; ++i)
	{
		stringstream ss;
		ss << "C_" << i;
		_end[i] = IloNumVar(*_env, 0, _ins->_LF[i - 1], IloNumVar::Float, ss.str().c_str());
	}

	_rate = IloNumVarArray(*_env, N);
	for (int i = 1; i < puits; ++i)
	{
		stringstream ss;
		ss << "r_" << i;
		_rate[i] = IloNumVar(*_env, _ins->_debitMin[i - 1], _ins->_debitMax[i - 1], IloNumVar::Float, ss.str().c_str());

	}

	//3. allocation de la var pour la fonction objectif
	_obj = IloNumVar(*_env, -IloInfinity, IloInfinity, IloNumVar::Float, "OBJ");
}


//on cree les contraintes lineaires : toutes sauf "fin(j) = deb(j) + pop(j)/rate(j)"
void MIP_flot_partiel::creationCtrLineaire(const vector<pair<int, int>> & arcSupport)
{
	int N = _ins->_nbJob + 2;//+2 => on ajoute source (0) et puits (N-1)
	int puits = N - 1;
	int K = _ins->_nbArc;

	vector<vector<int>> arcSucc(N);//arc[i] contient la liste des j tel que i -> j est un arcSupport
	vector<vector<int>> arcPrec(N);//arc[i] contient la liste des j tel que j -> i est un arcSupport

	//attenion dans l'instances les "vrais jobs" commencent à 0 !!

	//1. obj = min des marges
	for (int i = 1; i < puits; ++i)
	{
		if (_existe[i - 1])
			_mod->add(_obj - _ins->_LF[i - 1] + _end[i] <= 0);
	}



	//2. bornes sur les dates de fin avec le debit max
	// les autres bornes sont fixees quand on declare les var
	for (int i = 1; i < puits; ++i)
	{
		if (_existe[i - 1])
			_mod->add(_end[i] - _start[i] >= static_cast<double>(_ins->_pop[i - 1]) / _ins->_debitMax[i - 1]);
	}

	//3. contraintes de precedence
	for (auto p : arcSupport)
	{
		int i = p.first, j = p.second;
		_mod->add(_start[j + 1] - _end[i + 1] >= 0);
		arcSucc[i].push_back(j);
		arcPrec[j].push_back(i);
	}




	//4. contraintes de flot 

	//4.2 ctr de ressources

	for (int k = 0; k < K; ++k)
	{
		IloExpr expr1(*_env), expr2(*_env);

		for (int i = 1; i <= puits; ++i)
		{
			if (i == puits || (_existe[i - 1] && _ins->_utilise[i - 1][k]))
				expr1 += _flot[0][i][k];
		}

		for (int i = 0; i < puits; ++i)
		{
			if (i == 0 || (_existe[i - 1] && _ins->_utilise[i - 1][k]))
				expr2 += _flot[i][puits][k];
		}

		_mod->add(expr1 - expr2 == 0); //ce qui sort de la source = ce qui entre au puits
		_mod->add(expr1 == _ins->_cap[k]);

		expr1.end();
		expr2.end();
	}

	//4.3 loi de kirchhoff

	for (int k = 0; k < K; ++k)
	{

		for (int i = 1; i < puits; ++i)
		{
			//si i emprunte l'arc k
			if (_existe[i - 1] && _ins->_utilise[i - 1][k])
			{
				IloExpr exprIn(*_env), exprOut(*_env);

				exprIn += _flot[0][i][k];
				for (int j : arcPrec[i - 1]) //j = pred  de i dans les arcs support
				{
					if (_ins->_utilise[j][k])
					{
						j++; // decalage de un a cause de la source
						if (i != j)
							exprIn += _flot[j][i][k];
					}
				}

				exprOut += _flot[i][puits][k];
				for (int j : arcSucc[i - 1]) //j = succ de i dans les arcs support
				{
					if (_ins->_utilise[j][k])
					{
						j++; // decalage de un a cause de la source
						if (i != j)
							exprOut += _flot[i][j][k];
					}
				}

				_mod->add(exprIn - exprOut == 0);
				_mod->add(exprIn - _rate[i] == 0);
				exprIn.end();
				exprOut.end();
			}
		}
	}

}



//on cree la fonction objectif
void MIP_flot_partiel::creationFctObj()
{
	_mod->add(IloMaximize(*_env, _obj));
}




//affiche la sol. (suppose que _cplex->solve a ete appele !!)
void MIP_flot_partiel::afficheSol()
{
	cout << "=============================================" << endl;
	cout << "statut = " << _cplex->getStatus() << endl;
	cout << "=============================================" << endl;

	for (int i = 1; i <= _ins->_nbJob; ++i)
	{
		if (_existe[i-1])
			cout << "r" << i << " = " << _cplex->getValue(_rate[i]) << endl;
	}

	for (int i = 1; i <= _ins->_nbJob; ++i)
	{
		if (_existe[i - 1]) 
		{
			cout << "Start_" << i << " = " << _cplex->getValue(_start[i]) << " --- ";
			cout << "End_" << i << " = " << _cplex->getValue(_end[i]) << endl;
		}
	}


	for (int k = 0; k < _ins->_nbArc; ++k)
	{
		double sum = 0;
		for (int i = 1; i <= _ins->_nbJob; ++i)
		{
			if (_existe[i - 1])
			{
				if (_cplex->isExtracted(_flot[0][i][k]))
					sum += _cplex->getValue(_flot[0][i][k]);
			}
		}

		cout << "sum " << k << " = " << sum << endl;
	}

}


//on cree nbPoints tangentes qui linearisent (partiellement)"fin(j) = deb(j) + pop(j)/rate(j)"
//si uniforme = true alors on en crée nbPoint espacés régulièrement, sinon on en crée nbPoint non uniformément espacées
//les autres tangentes seront générées dynamiquement (dans le lazy ctr callback)
void MIP_flot_partiel::creationCtrTangente(bool uniforme, int nbPoint)
{
	//ATTENTION : decalage d'un entre les vrais jobs cplex (1..N) et les jobs dans l'instance (0..N-1)
	//a cause de l'ajout de la source / du puits pour modeliser le flot dans cplex

	int N = _ins->_nbJob;
	double S = 1; // S = 1 + 1/2 + 1/3 +...+ 1/(nbPoint-1) utile si pas non constant

	if (!uniforme)
	{
		for (int k = 2; k < nbPoint; ++k)
			S += 1.0 / k;
	}

	for (int i = 1; i <= N; ++i)
	{
		if (_existe[i - 1])
		{
			//on calcul y tel que debitMin = DebitMax - y * S
			// => si pas non constant alors pas = y/(nbPoint - 1 - k)
			double y = (_ins->_debitMax[i - 1] - _ins->_debitMin[i - 1]) / S;

			double pas = (_ins->_debitMax[i - 1] - _ins->_debitMin[i - 1]) / (nbPoint - 1);
			if (!uniforme)
				pas = y / (nbPoint - 1);

			double x = _ins->_debitMin[i - 1]; //abscisse auquel on calcul la tangente

			for (int k = 0; k < nbPoint; ++k)
			{

				double a = _ins->_pop[i - 1] / x / x; //pop / x^2
				_mod->add(_end[i] - _start[i] + a * _rate[i] >= 2 * _ins->_pop[i - 1] / x);

				//abscise pour la prochaine tangente
				x += pas;

				if (!uniforme && k != nbPoint - 2)
					pas = y / (nbPoint - k - 2);
			}

		}
	}
}


//verif solution donnee par cplex
bool MIP_flot_partiel::verifSol()
{
	int N = _ins->_nbJob;
	bool ok = true;
	double obj = INT_INFINITY;
	vector<double> start(N), end(N), rate(N);

	//prerequis : on a une solution !
	if (_cplex->getStatus() != IloAlgorithm::Status::Optimal && _cplex->getStatus() != IloAlgorithm::Status::Feasible)
		return true;



	for (int i = 0; i < N; ++i)
	{
		if (_existe[i])
		{
			start[i] = _cplex->getValue(_start[i + 1]);
			end[i] = _cplex->getValue(_end[i + 1]);
			rate[i] = _cplex->getValue(_rate[i + 1]);

			//1. on commence apres la release date
			if (start[i] < _ins->_ES[i] - EPSILON_DEBIT)
			{
				cout << "start time < release : " << start[i] << " < " << _ins->_ES[i] << endl;
				ok = false;
			}

			//2. on finit avant la due date
			if (end[i] > _ins->_LF[i] + EPSILON_DEBIT)
			{
				cout << "end time > LF : " << end[i] << " > " << _ins->_LF[i] << endl;
				ok = false;
			}

			//3. on a end = start + pop / debit
			if (end[i] < start[i] + _ins->_pop[i] / rate[i] - EPSILON_DEBIT)
			{
				cout << setprecision(6) << fixed;
				cout << "end time <  start + pop/rate : " << end[i] << " < " << start[i] + _ins->_pop[i] / rate[i] << endl;
				ok = false;
			}

			//4. la fonction obj est ok
			obj = min(obj, _ins->_LF[i] - end[i]);
		}
	}

	//verif fct obj
	if (abs(obj - _cplex->getObjValue()) > EPSILON_DEBIT)
	{
		cout << "fonction obj = " << _cplex->getObjValue() << " -- recalculee = " << obj << endl;
		ok = false;
	}

	//verif capa des arcs
	ok = ok && verifDebit(start, end, rate);

	return ok;
}


bool MIP_flot_partiel::verifDebit(const vector<double> & dateDeb, const vector<double> & dateFin, const vector<double> & debit) const
{
	bool ok = true;
	int N = _ins->_nbJob;

	typedef struct triplet
	{
		double date;
		int job;
		int type; // 0 = date de fin, 1 = date de debut
	}triplet;

	vector<triplet> dates;

	for (int i = 0; i < N; ++i)
	{
		if (_existe[i])
		{
			dates.push_back( { dateDeb[i],i,1 });
			dates.push_back( { dateFin[i],i,0 });

		}
	}

	//ici on met EPSILON_DEBIT pour les dates sinon on rsique de faux positifs avec EPSILON_DATE
	sort(dates.begin(), dates.end(), [](const triplet & a, const triplet & b)
	{return a.date + EPSILON_DATE <= b.date || (abs(a.date - b.date) < EPSILON_DATE && a.type < b.type); });


#ifdef _VERIF_
	for (int i = 0; i < dates.size() - 1; ++i)
		if (dates[i].date > dates[i].date)
			stopProg("MIP_flot::verifDebit : les dates ne sont pas dans le bon ordre");
#endif

	vector <double> debitArc(_ins->_nbArc, 0);

	for (int i = 0; ok && i < dates.size(); ++i)
	{
		int type = dates[i].type;
		int job = dates[i].job;

		if (type == 1)
		{
			for (int a : _ins->_jobToArc[job])
			{
				debitArc[a] += debit[job];
				if (debitArc[a] > _ins->_cap[a] + EPSILON_DEBIT)
				{
					cout << "debit total = " << debitArc[a] << "capa = " << _ins->_cap[a] << endl;
					ok = false;
				}
			}
		}
		else //type = 0
		{
			for (int a : _ins->_jobToArc[job])
			{
				debitArc[a] -= debit[job];
#ifdef _VERIF_
				if (debitArc[a] < -EPSILON)
					stopProg("BranchBoundDicho::reconstruireSolCouranteAmeliore : debit devient < 0");
#endif
			}
		}

	}

	return ok;

}

