#include "MIP_flot.h"
#include "callback_MIP_Flot.h"



//retourne borne inf, borne sup
pair<double, double> MIP_flot::creeEtResout(const vector<double> & debut, const vector<double> & debit, const vector<vector<vector<double>>> & flot)
{
	pair<double, double>  res = { -1,-1 };	
	
	allocVar();

	creationCtrLineaire();
	creationFctObj();
	
	creationCtrTangente(true, 5);
	//creationCtrTangente_NonUniformeCarre(10);


	// ==== warm start ====
	if (debut.size() > 0)
	{
		warmStart1(debut, debit);
		warmStart2(debut, debit);
		warmStart3(debut, debit, flot);
		//warmStart4(flot);
		//_cplex->setParam(IloCplex::Param::Advance, 2);
	}


	_cplex->setParam(IloCplex::Param::ClockType, CLOCK_TYPE_CPLEX);
	_cplex->setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_CPLEX);
	//_cplex->setParam(IloCplex::Param::MIP::Limits::TreeMemory, RAM_LIMIT_CPLEX);

	_cplex->setParam(IloCplex::Param::MIP::Display, 0);
	//_cplex->setOut(_env->getNullStream());

	//  -----  on fixe quelques param pour le branch & cut   -------   
	_cplex->setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
	_cplex->setParam(IloCplex::Param::Threads, 1);
	_cplex->setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
	//--------------------------------------------------------------------
	
	//uniquement pour les instances avec des pb numeriques (non utilise ici)
	//_cplex->setParam(IloCplex::Param::Emphasis::Numerical, true);
	//_cplex->setParam(IloCplex::Param::MIP::Tolerances::Integrality, 0.001);//default 10-5
	//_cplex->setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.000001);//default 10-4
	_cplex->setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 0.0001);//default 10-6
	
	
	//  -----  integration des callbacks   -------   
	int N = _ins->_nbJob;
	IloNumArray pop (*_env, N + 1);
	for (int i = 1; i <= N; ++i)
		pop[i] = _ins->_pop[i - 1];

	_cplex->use(TangenteLazyCallback(*_env, N, _start, _end, _rate, _x, pop, _nbCutLazy));
	//_cplex->use(TangenteUserCallback(*_env, N, _start, _end, _rate, _x, pop, _nbCutUser));
	//--------------------------------------------------------------------

	//_cplex->exportModel("mipflot.lp");
	//_cplex->readMIPStarts("mipStart.sol");
	_cplex->solve();
	

	//ATTENTION : il faut enlever 1 a l'obj pour etre conforme avec l'article GEOSAFE et les autres algos
	if (_cplex->getStatus() == IloAlgorithm::Status::Optimal || _cplex->getStatus() == IloAlgorithm::Status::Feasible)
	{
		//_cplex->writeSolution("solCplex.txt");

#ifdef _VERIF_
		//if (!verifSol())
		//	stopProg("MIP_flot::creeEtResout : sol non faisable");
#endif
		res.first = _cplex->getObjValue() - 1; //obj du current best sol = BI
	}
	
	if (_cplex->getStatus() == IloAlgorithm::Status::Optimal)
		res.second = res.first;
	else
		res.second = _cplex->getBestObjValue() - 1; //borne sup (on maximise)

	cout << "MIP BI, BS = " << res.first << ", " << res.second << endl;

	//if (_cplex->getStatus() == IloAlgorithm::Status::Optimal || _cplex->getStatus() == IloAlgorithm::Status::Feasible)
	//	afficheSol();

	//_cplex->writeMIPStarts("mipStart.sol");

	return res;
}


//allocation des variables
void MIP_flot::allocVar()
{
	int N = _ins->_nbJob + 2;//+2 => on ajoute source et puits
	int puits = N - 1;
	int K = _ins->_nbArc;

	_x = IloArray< IloNumVarArray>(*_env, N); //on alloue mais ne concerne ni la source, ni le puits

	//1. allocation des variables binaires
	for (int i = 1; i < puits; ++i)
	{
		_x[i] = IloNumVarArray(*_env, N);

		for (int j = 1; j < puits; ++j)
		{
			stringstream ss;
			ss << "x_" << i <<  "_" << j;
			_x[i][j] = IloNumVar(*_env, 0, 1, IloNumVar::Bool, ss.str().c_str());
		}

	}//fin for i 

	//2. allocation des variables flots
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


	//3. allocation des tableaux 1D : debut / fin / vitesse
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

	//4. allocation de la var pour la fonction objectif
	_obj = IloNumVar(*_env, -IloInfinity, IloInfinity, IloNumVar::Float, "OBJ");
}


//on cree les contraintes lineaires : toutes sauf "fin(j) = deb(j) + pop(j)/rate(j)"
void MIP_flot::creationCtrLineaire()
{
	int N = _ins->_nbJob + 2;//+2 => on ajoute source (0) et puits (N-1)
	int puits = N - 1;
	int K = _ins->_nbArc;

	//attenion dans l'instances les "vrais jobs" commencent à 0 !!
	
	//1. obj = min des marges
	for (int i = 1; i < puits; ++i)
	{
		_mod->add(_obj - _ins->_LF[i - 1] + _end[i] <= 0);
	}



	//2. bornes sur les dates de fin avec le debit max
	// les autres bornes sont fixees quand on declare les var
	for (int i = 1; i < puits; ++i)
	{
		_mod->add(_end[i] - _start[i] >= static_cast<double>(_ins->_pop[i - 1]) / _ins->_debitMax[i - 1]);
	}

	//3. contraintes de precedence
	for (int i = 1; i < puits; ++i)
	{
		for (int j = i + 1; j < puits; ++j)
		{
			//4.1. i << j
			double M = _ins->_LF[j - 1] - _ins->_ES[i - 1];
			_mod->add(_start[i] - _end[j] - _x[j][i] * M >= -M);

			//4.2. j << i
			M = _ins->_LF[i - 1] - _ins->_ES[j - 1];
			_mod->add(_start[j] - _end[i] - _x[i][j] * M >= -M);

		} 
	}


	////4. contraintes de flot 

	////4.1 x_i_j = 0 => flot(i,j) = 0
	
	for (int i = 1; i < puits; ++i)
	{
		for (int j = 1; j < puits; ++j)
		{
			if (i != j)
			{
				for (int k = 0; k < K; ++k)
				{
					if (_ins->_utilise[i - 1][k] && _ins->_utilise[j - 1][k])
						_mod->add(_flot[i][j][k] - min(_ins->_debitMax[i - 1], _ins->_debitMax[j - 1]) * _x[i][j] <= 0);
				}
			}

		}
	}

	////4.2 ctr de ressources

	for (int k = 0; k < K; ++k)
	{
		IloExpr expr1(*_env), expr2(*_env);

		for (int i = 1; i <= puits; ++i)
		{
			if (i == puits || _ins->_utilise[i - 1][k])
				expr1 += _flot[0][i][k];
		}

		for (int i = 0; i < puits; ++i)
		{
			if (i == 0 || _ins->_utilise[i - 1][k])
				expr2 += _flot[i][puits][k];
		}

		_mod->add(expr1 - expr2 == 0);
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
			if (_ins->_utilise[i - 1][k])
			{
				IloExpr exprIn(*_env), exprOut(*_env);

				exprIn += _flot[0][i][k];
				for (int j : _ins->_arcToJob[k]) //j = job qui emprunte arc k
				{
					j++; // decalage de un a cause de la source
					if (i != j)
						exprIn += _flot[j][i][k];
				}

				exprOut += _flot[i][puits][k];
				for (int j : _ins->_arcToJob[k]) //j = job qui emprunte arc k
				{
					j++; // decalage de un a cause de la source
					if (i != j)
						exprOut += _flot[i][j][k];
				}

				_mod->add(exprIn - exprOut == 0);
				_mod->add(exprIn - _rate[i] == 0);
				exprIn.end();
				exprOut.end();
			}
		}
	}

}

// creation de la contrainte "fin(j) >= deb(j) + pop(j)/rate(j)"
// ne fonctionne pas : cplex n accepte que des variables binaires dans les ctr quad convexe
void MIP_flot::creationCtrConvexe()
{
	int N = _ins->_nbJob + 2; //+2 => on ajoute source (0) et puits (N-1)
	int puits = N - 1;
	
	
	for (int i = 1; i < puits; ++i)
	{
		_mod->add(_end[i]* _rate[i] - _start[i]* _rate[i]  >= _ins->_pop[i - 1]);
	}


}

//on cree la fonction objectif
void MIP_flot::creationFctObj()
{
	_mod->add(IloMaximize(*_env, _obj));
}



//
//
//int _nbJob; //nb jobs = nb de population a evacuer
//int _nbArc; //nb d'arcs du chemin
//
//vector<int> _job; //_job[i] = numero du job dans Istance
//
//vector<int> _pop;//_pop[j] = population (nombre de personnes a evacuer) pour j [0..._nbJob-1]
//vector<double> _ES;//_ES[j] = date au plus tot d'arrivee a la porte de j [0..._nbJob-1]
//vector<double> _LF;//_LF[j] = date au plus tard d'arrivee a la porte [0..._nbJob-1]
//double _ESmin; //ESmin = min des ES[j] = date au plus ou une population peut arriver a la porte
//vector<double> _debitMax; //debitMax[j] = debit max de j
//vector<double> _debitMin; //debitMin[j] = debit min de j (pour arriver avant sa deadline)
//
//vector<int> _cap; //_cap[i] = capacite de l arc i [0..._nbArc-1]
//vector<vector<int>> _arcToJob; //_arcToJob[i] = liste des jobs qui passent sur l'arc i
//vector<int> _jobToArcInit; //_jobToArcInit[j] = premier arc d'evacuation de j
//vector<vector<int>> _jobToArc; //_jobToArc[i] = liste des arcs empruntes par i
//
//vector<vector<bool> > _utilise; //utilise[i][a] = true si le job i utilise l'arc a


//affiche la sol. (suppose que _cplex->solve a ete appele !!)
void MIP_flot::afficheSol()
{
	cout << "=============================================" << endl;
	cout << "statut = " << _cplex->getStatus() << endl;
	cout << "=============================================" << endl;

	for (int i = 1; i <= _ins->_nbJob; ++i)
	{
		cout << "r" << i << " = " << _cplex->getValue(_rate[i]) << endl;
	}

	for (int i = 1; i <= _ins->_nbJob; ++i)
	{
		cout << "Start_" << i << " = " << _cplex->getValue(_start[i]) << " --- ";
		cout << "End_" << i << " = " << _cplex->getValue(_end[i]) << endl;
	}


	for (int k = 0; k < _ins->_nbArc; ++k)
	{
		double sum = 0;
		for (int i = 1; i <= _ins->_nbJob; ++i)
		{
			if (_cplex->isExtracted(_flot[0][i][k]))
				sum += _cplex->getValue(_flot[0][i][k]);
		}

		cout << "sum " << k << " = " << sum << endl;
	}

}

//retourne le nb de node explores par cplex
long MIP_flot::getExploredNode()
{

	return static_cast<long>(_cplex->getNnodes());
}



//on cree nbPoints tangentes qui linearisent (partiellement)"fin(j) = deb(j) + pop(j)/rate(j)"
//si uniforme = true alors on en crée nbPoint espacés régulièrement, sinon on en crée nbPoint non uniformément espacées
//les autres tangentes seront générées dynamiquement (dans le lazy ctr callback)
void MIP_flot::creationCtrTangente(bool uniforme, int nbPoint)
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
		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->_debitMax[i-1] - _ins->_debitMin[i-1]) / S;

		double pas = (_ins->_debitMax[i-1] - _ins->_debitMin[i-1]) / (nbPoint - 1);
		if (!uniforme)
			pas = y / (nbPoint - 1);

		double x = _ins->_debitMin[i-1]; //abscisse auquel on calcul la tangente

		for (int k = 0; k < nbPoint; ++k)
		{

			double a = _ins->_pop[i-1] / x / x; //pop / x^2
			_mod->add(_end[i] - _start[i] + a * _rate[i]  >= 2 * _ins->_pop[i - 1] / x - EPSILON_CUT_CPLEX);

			//abscise pour la prochaine tangente
			x += pas;

			if (!uniforme && k != nbPoint - 2)
				pas = y / (nbPoint - k - 2);
		}

	}
}


//idem creationCtrTangente mais on genere les points non uniformement 
// avec une formule differente qui depend du carre du nb de points
void MIP_flot::creationCtrTangente_NonUniformeCarre(int nbPoint)
{
	//ATTENTION : decalage d'un entre les vrais jobs cplex (1..N) et les jobs dans l'instance (0..N-1)
	//a cause de l'ajout de la source / du puits pour modeliser le flot dans cplex

	int N = _ins->_nbJob;
	int K2 = (nbPoint - 1)*(nbPoint - 1);

	for (int i = 1; i <= N; ++i)
	{
		double Dmin = _ins->_debitMin[i - 1];
		double Dmax = _ins->_debitMax[i - 1];

		double x = Dmin; //abscisse auquel on calcul la tangente
		
		for (int k = 1; k <= nbPoint; ++k)
		{

			double a = _ins->_pop[i - 1] / x / x; //pop / x^2
			_mod->add(_end[i] - _start[i] + a * _rate[i] >= 2 * _ins->_pop[i - 1] / x - EPSILON_CUT_CPLEX);

			//abscise pour la prochaine tangente
			int k2 = k * k;

			x = k2 * Dmax + (K2 - k2)*Dmin;
			x /= K2;
		}

	}
	

}



//verif solution donnee par cplex
bool MIP_flot::verifSol()
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
		start[i] = _cplex->getValue(_start[i+1]);
		end[i] = _cplex->getValue(_end[i+1]);
		rate[i] = _cplex->getValue(_rate[i+1]);

		//1. on commence apres la release date
		if (start[i] < _ins->_ES[i] - EPSILON_DEBIT)
		{
			cout << "start time < release : " << start[i]  << " < " << _ins->_ES[i] << endl;
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


bool MIP_flot::verifDebit(const vector<double> & dateDeb, const vector<double> & dateFin, const vector<double> & debit) const
{
	bool ok = true;
	int N = _ins->_nbJob;

	typedef struct triplet
	{
		double date;
		int job;
		int type; // 0 = date de fin, 1 = date de debut
	}triplet;

	vector<triplet> dates(2 * N);
	int cpt = 0;
	for (int i = 0; i < N; ++i)
	{
		dates[cpt] = { dateDeb[i],i,1 };
		dates[cpt + 1] = { dateFin[i],i,0 };
		cpt += 2;
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
	int N2 = 2 * N;
	for (int i = 0; ok && i < N2; ++i)
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


//fourni a cplex une solution initiale pour demarrer l'optimisation
//on donne les debuts, debits, fin et x_ij qui font 1
void MIP_flot::warmStart1(const vector<double> & debut, const vector<double> & debit)
{
	int N = _ins->_nbJob; //attention les jobs commencent a 0 dans l'instance et 1 dans le MIP (0 = source dans le MIP)


	IloNumVarArray startVar(*_env);
	IloNumArray startVal(*_env);

	for (int i = 0; i < N; ++i)
	{
		if (_cplex->isExtracted(_start[i + 1]))
		{
			startVar.add(_start[i + 1]);
			startVal.add(debut[i]);
		}
		else
			cout << "start" << i << endl;
		if (_cplex->isExtracted(_rate[i + 1]))
		{
			startVar.add(_rate[i + 1]);
			startVal.add(debit[i]);
		}
		else
			cout << "_rate" << i << endl;

		if (_cplex->isExtracted(_end[i + 1]))
		{
			startVar.add(_end[i + 1]);
			startVal.add(debut[i] + (double)(_ins->_pop[i]) / debit[i]);
		}
		else
			cout << "_end" << i << endl;

		for (int j = 0; j < N; ++j)
		{
			if (i != j)
			{
				if (_cplex->isExtracted(_x[i + 1][j + 1]) && debut[i] + (double)(_ins->_pop[i]) / debit[i] <= debut[j])
				{
					startVar.add(_x[i + 1][j + 1]);
					startVal.add(1);
				}
				if (!_cplex->isExtracted(_x[i + 1][j + 1]))
					cout << "_x" << i << " " << j << endl;
			}
		}
	}

	_cplex->addMIPStart(startVar, startVal);

	startVal.end();
	startVar.end();
}

//fourni a cplex une solution initiale partielle pour demarrer l'optimisation : 
//on ne fixe que les x_ij = 1 qui correspondent a fin(i) <= debut(j)
void MIP_flot::warmStart2(const vector<double> & debut, const vector<double> & debit)
{
	int N = _ins->_nbJob; //attention les jobs commencent a 0 dans l'instance et 1 dans le MIP (0 = source dans le MIP)


	IloNumVarArray startVar(*_env);
	IloNumArray startVal(*_env);

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			if (i != j)
			{
				//si i << j
				if (_cplex->isExtracted(_x[i + 1][j + 1]) && debut[i] + (double)(_ins->_pop[i]) / debit[i] <= debut[j])
				{
					startVar.add(_x[i + 1][j + 1]);
					startVal.add(1);
				}

			}
		}
	}

	_cplex->addMIPStart(startVar, startVal);

	startVal.end();
	startVar.end();
}

//fourni a cplex une solution initiale partielle pour demarrer l'optimisation : 
void MIP_flot::warmStart3(const vector<double> & debut, const vector<double> & debit, const vector<vector<vector<double>>> & flot)
{
	int N = _ins->_nbJob; //attention les jobs commencent a 0 dans l'instance et 1 dans le MIP (0 = source dans le MIP)


	IloNumVarArray startVar(*_env);
	IloNumArray startVal(*_env);

	for (int i = 0; i < N; ++i)
	{
		if (_cplex->isExtracted(_start[i + 1]))
		{
			startVar.add(_start[i + 1]);
			startVal.add(debut[i]);
		}

		if (_cplex->isExtracted(_rate[i + 1]))
		{
			startVar.add(_rate[i + 1]);
			startVal.add(debit[i]);
		}


		if (_cplex->isExtracted(_end[i + 1]))
		{
			startVar.add(_end[i + 1]);
			startVal.add(debut[i] + (double)(_ins->_pop[i]) / debit[i]);
		}

		
		for (int j = 0; j < N; ++j)
		{
			if (i != j)
			{
				double sumFlot = 0;
				for (int k = 0; k < _ins->_nbArc; ++k)
				{
					if (_cplex->isExtracted(_flot[i + 1][j + 1][k]))
					{
						startVar.add(_flot[i + 1][j + 1][k]);
						startVal.add(flot[i][j][k]);
					}
					sumFlot += flot[i][j][k];
				}

				if (sumFlot > 0)
				{
					startVar.add(_x[i + 1][j + 1]);
					startVal.add(1);
				}
			}
		}
	}

	_cplex->addMIPStart(startVar, startVal);

	startVal.end();
	startVar.end();
}

//fourni a cplex une solution initiale partielle pour demarrer l'optimisation : 
void MIP_flot::warmStart4(const vector<vector<vector<double>>> & flot)
{
	int N = _ins->_nbJob; //attention les jobs commencent a 0 dans l'instance et 1 dans le MIP (0 = source dans le MIP)


	IloNumVarArray startVar(*_env);
	IloNumArray startVal(*_env);

	for (int i = 0; i < N; ++i)
	{
		

		for (int j = 0; j < N; ++j)
		{
			if (i != j)
			{
				double sumFlot = 0;
				for (int k = 0; k < _ins->_nbArc; ++k)
				{
					if (_cplex->isExtracted(_flot[i + 1][j + 1][k]))
					{
						startVar.add(_flot[i + 1][j + 1][k]);
						startVal.add(flot[i][j][k]);
					}
					sumFlot += flot[i][j][k];
				}

				
			}
		}
	}

	_cplex->addMIPStart(startVar, startVal);

	startVal.end();
	startVar.end();
}