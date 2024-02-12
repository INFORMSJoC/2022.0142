#include "BranchBoundDicho.h"
#include "my_time.h"

#include <algorithm>


unsigned int BranchBoundDicho::Noeud::_nbNoeudTotal = 0;

//===========================================================================
// méthode principale : effectue le B&B

pair<double, double> BranchBoundDicho::run(double cpuMax)
{
	double margeEval = -1;
	my_time debut_cpu = give_time();
	

	//initialisation : dans le constructeur (utilisation du résultat de l'heuristique)

	_nonAmelioreb0 = 0;
	_nonAmelioreb1 = 0;
	_nonAmelioreb2 = 0;
	_nbBranchType1 = 0;
	_nbBranchType2ou3 = 0;
	_nbNoeudExplore = 0; 


	//=============================================================================
	//noeud racine

	Noeud racine;
	Noeud::_nbNoeudTotal++;

	racine._rel = _ins->_ES;
	racine._due = _ins->_LF;
	racine._pointMil_1.assign(_ins->_nbJob, -1);
	racine._pointMil_2.assign(_ins->_nbJob, -1);


	double margeRacine = TIME_INFINITY;
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		double val = -1;
		if (_ins->_jobToArcInit.size() != 0)
			val = _ins->_LF[i] - _ins->_ES[i] - static_cast<double>(_ins->_pop[i]) / _ins->_cap[_ins->_jobToArcInit[i]];
		else
		{
			double cap = INT_INFINITY;
			for (int e : _ins->_jobToArc[i])
				cap = min(cap, (double)(_ins->_cap[e]));
			val = _ins->_LF[i] - _ins->_ES[i] - static_cast<double>(_ins->_pop[i]) / cap;
		}

		margeRacine = min(margeRacine, val);
	}

	margeEval = evalNoeud(racine, margeRacine);

	//on s'arrete ici si on ne veut que le resultat des heuristiques a la racine
	//return{ -1.0,  margeEval };


	_coutRacine = margeEval;

#ifdef _VERIF_
	if (margeEval == -1)
		stopProg("noeud racine infaisable");
#endif

	

	//si margeEval > _bestMargeRealisable on continue => il est possible de trouver mieux que _bestMargeRealisable 
	if (abs(margeEval - _bestMargeRealisable) > EPSILON)
	{
		Noeud noeud = racine;

		//creer le noeud
		
		SignatureBranchement sign = calculParamBranchement();

		if (sign._job >= 0)
		{
			noeud._marge = margeEval;
			noeud._sign = sign;
			noeud._priorite = 0;//peu importe la priorite pour la racine

			_nbBranchType1 = (sign._type == 1) ? _nbBranchType1 + 1 : _nbBranchType1;
			_nbBranchType2ou3 = (sign._type != 1) ? _nbBranchType2ou3 + 1 : _nbBranchType2ou3;
			
#ifdef _VERIF_BRANCHEMENT_
			noeud._sol = _solCour;
			noeud._date = _dateCour;
#endif

			_listeNoeud.push(noeud);

		}
	}

	//======================================================================
	// boucle principale du B&B


	unsigned int seuil = 0;

	while (!_listeNoeud.empty() && give_time() - debut_cpu < cpuMax)
	{
		_nbNoeudExplore++;

		Noeud noeudCour = _listeNoeud.top(); //on copie car on va modifier la  file
		_listeNoeud.pop();

		//cout << noeudCour._priorite << endl;

		//si on depile un noeud plus mauvais que la solution courante, on passe au suivant (inutile de l'evaluer !)
		if (noeudCour._marge <= _bestMargeRealisable + EPSILON_DEBIT)
			continue;
		
		if (Noeud::_nbNoeudTotal > seuil)
		{
			//cout << "noeud " << noeudCour._id << "(" << noeudCour._marge << ") " << endl;
			//cout << "cpu = " << give_time() - debut_cpu << endl;
			seuil += 1000;
			//seuil += 10;
		}


		// on cree 3 fils (b=0...2)
		for (int b = 0; b < 3; ++b)
		{
			Noeud fils = noeudCour;

			//modifie le noeud fils en fonction de sa signature et de b (numero de branche)
			appliqueSignatureEtBranch(fils, b);


			//on evalue le noeud s il est capable de faire mieux que la sol. actuelle
			margeEval = -1;
			if (noeudCour._marge - EPSILON_COMP > _bestMargeRealisable )
				 margeEval = evalNoeud(fils, noeudCour._marge);


#ifdef _VERIF_BRANCHEMENT_
			if (margeEval > -EPSILON)
			{
				bool solIdem = false;
				if (fils._date.size() == _dateCour.size())
				{
					bool dateIdem = true;
					for (unsigned int k = 0; k < _dateCour.size(); ++k)
					{
						if (abs(_dateCour[k] - fils._date[k]) > EPSILON_DATE_VERIF)
							dateIdem = false;
					}

					if (dateIdem)
					{
						solIdem = true;
						for (int i = 0; i < _ins->_nbJob; ++i)
						{
							for (unsigned int k = 0; k < _dateCour.size(); ++k)
							{
								if (abs(_solCour[i][k] - fils._sol[i][k]) > EPSILON_DEBIT_VERIF)
									solIdem = false;
							}
						}
					}
				}
				if (solIdem)
					stopProg("BranchBoundDicho::run : le branchement ne coupe pas");
			}
#endif


			if (abs(margeEval - noeudCour._marge) < EPSILON)
			{
				_nonAmelioreb0 = (b == 0) ? _nonAmelioreb0 + 1 : _nonAmelioreb0;
				_nonAmelioreb1 = (b == 1) ? _nonAmelioreb1 + 1 : _nonAmelioreb1;
				_nonAmelioreb2 = (b == 2) ? _nonAmelioreb2 + 1 : _nonAmelioreb2;
			}

			//si l'évaluation du noeud fils donne une borne sup (sol. optimiste - on maximise) meilleure que la solution courante
			// on sauvegarde ce noeud
			if (margeEval > -EPSILON && margeEval > _bestMargeRealisable)
			{
				//calcule la signature de branchement pour la solution courante associee au vecteur date courant
				SignatureBranchement sign = calculParamBranchement();
				if (sign._job >= 0)
				{
					Noeud::_nbNoeudTotal++;

					fils._marge = margeEval;
					fils._sign = sign;
					fils._id = Noeud::_nbNoeudTotal;
					fils._idPere = noeudCour._id;
					fils._priorite = 0;//si on fait un parcours en prof, peu importe la priorite

#ifndef _PARCOURS_PROF_
					fils._priorite = margeEval;// +10000 * nbJobResolu();
					//fils._priorite = intervalMin();
#endif // !_PARCOURS_PROF_


					_nbBranchType1 = (sign._type == 1) ? _nbBranchType1 + 1 : _nbBranchType1;
					_nbBranchType2ou3 = (sign._type != 1) ? _nbBranchType2ou3 + 1 : _nbBranchType2ou3;

#ifdef _VERIF_BRANCHEMENT_
					fils._sol = _solCour;
					fils._date = _dateCour;
#endif

					//cout << "nouv = " << fils._id << "(" << fils._marge << ") de pere " << noeudCour._id << endl;


					_listeNoeud.push(fils); //si file priorite : insere par marge decroissante

				}

			}
		}

	}



	double BS = _bestMargeRealisable;
	while (!_listeNoeud.empty())
	{
		BS = max(BS, _listeNoeud.top()._marge);
		_listeNoeud.pop();
	}

	cout << "nombre total de noeuds = " << Noeud::_nbNoeudTotal << endl;


	
	return{ _bestMargeRealisable,  BS};
}


//==========================================================================
// sous fonction : allocation / initialisation des attributs de classe
void BranchBoundDicho::init()
{

	int N = _ins->_nbJob;
	_due_2.resize(N);
	_indRel.resize(N);
	_indDue_2.resize(N);
	_indMil_1.resize(N);
	_indMil_2.resize(N);

	_bestDebit.resize(N);
	_bestDebut.resize(N);
}


//verifie que la sol best (bestDebit ; bestDebut) est ok
bool BranchBoundDicho::verifieSolBest(const double marge) const
{
	bool ok = true;
	double margeRecalculee = TIME_INFINITY;

	vector<double> dateFin(_ins->_nbJob);

	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		dateFin[i] = _bestDebut[i] + static_cast<double>(_ins->_pop[i]) / _bestDebit[i];
		margeRecalculee = min(margeRecalculee, _ins->_LF[i] - dateFin[i]);
	}

	ok = verifDebit(_bestDebut, dateFin, _bestDebit);

	if (abs(margeRecalculee - marge) > EPSILON_DATE)
		ok = false;


	return ok;
}

//evalue un noeud avec une procedure dicho basee sur un PL
double BranchBoundDicho::evalNoeud(const Noeud & noeud, double margePere)
{
	double res = -1;

#ifdef _VERIF_
	if (margePere < _bestMargeRealisable - EPSILON_DEBIT)
	{
		cout << "margePere = " << margePere << endl;
		cout << "_bestMargeRealisable = " << _bestMargeRealisable << endl;
		stopProg(" BranchBoundDicho::evalNoeud : margePere < _bestMargeRealisable");
	}
#endif

	// resolutionPL()  stocke la solution  dans _solCour
	
	bool ok = resolutionPL(noeud, margePere);
	
	if (ok)
	{
		//la margepere donne une solution preemptive...
		res = margePere;

		//... on verifie si on peut la transformer en solution non preemptive
		if (reconstruireSolutionCourante(true, margePere))
		{
			_bestMargeRealisable = margePere;

			//cout << " ***** sterilisation : best = " << margePere << endl;
		}
		
	}		
	else
	{
		double margeMin = max(0.0, _bestMargeRealisable);
		
		ok = resolutionPL(noeud, margeMin);
		
		
		//les contraintes de branchement ne rendent pas irrealisable la best marge => on explore
		if (ok)
		{
			double margeSup = margePere;

			while (margeSup - margeMin > EPSILON_DEBIT) //mettre EPSILON_DEBIT (on change mommentanement pour répondre aux rapporteurs mars 2023
			{
				double margeCour = (margeMin + margeSup) / 2.0;
				ok = resolutionPL(noeud, margeCour);

				if (!ok)
				{
					margeSup = margeCour;
				}
				else
				{
					margeMin = margeCour;
					//tente de reconstruire une solution non preemptive a partir de la solution courante
					bool sterilisation = reconstruireSolutionCourante(true, margeCour);
					

					if (sterilisation)
					{
#ifdef _VERIF_
						if (_bestMargeRealisable > margeCour-EPSILON)
							stopProg("BranchBoundDicho::evalNoeud : le test de sterilisation doit etre lance uniquement si la solution cour est meilleure");
#endif
						_bestMargeRealisable = margeCour;
						//cout << " ***** sterilisation : best = " << _bestMargeRealisable << endl;
						//la solution reconstruite a ete sauvegardee dans bestSol
					}
				
				}

#ifdef _VERIF_
			if (margeSup < margeMin - EPSILON)
				stopProg(" BranchBoundDicho::evalNoeud : margeSup < margeMin");
#endif
			}

			res = margeMin;
		}

	}


	//==========================================================================================================
	//on essaie de reconstruire la meilleure solution trouvee a ce noeud (celle qui donne la marge preemptive stockee dans res)

	//on n'a pas sterilise le noeud  => on lance les heuristiques avec proba 1/3
	if (_bestMargeRealisable + EPSILON < res )//&& rand() % 3 == 0)
	{
		vector <int>  i_deb, i_fin;
		vector <double> duree, debit;
		vector <pair<int, int> > precedence;
		vector <vector<int> > cliques;


		calculInfoJob(i_deb, i_fin, duree, debit, false);
		calculPredClique(i_deb, i_fin, precedence, cliques);

		//_res_DebitFixe = reconstruireSolHeuristique(debit, duree);

		//les fonctions suivantes calcul ube solution heuristique et mettent a jour _bestMargeRealisable si besoin
		//_res_heurCut = reconstruireSolCourantePL_heur(i_deb, i_fin, duree, debit);
		//_res_tangente = reconstruireSolCourantePL_exact(precedence, cliques, TANGENTE, false);
		//_res_secante = reconstruireSolCourantePL_exact(precedence, cliques, SECANTE, false);
		
		//calcul une solution et la sauvegarde si elle est meilleure que la solution actuelle
		_res_cut = reconstruireSolCourantePL_exact(precedence, cliques, CUT, true);



		//cout << "res_heur = " << res_heur << " -- res_tangente = " << res_tangente << " --- res_secante = " << res_secante  << endl;
		//if (res_tangente > 0.1)
		//{
		//	cout << "res_heur = " << res_heur << " -- res_tangente = " << res_tangente << " --- res_secante = " << res_secante << " --- res_cut = " << res_cut << endl;
		//	cout << (res_tangente - res_secante) * 100 / res_secante << endl;
		//}
		//else
		//	cout << "echec" << endl;

//#ifdef _VERIF_
//			if (res_heur > res_cut + EPSILON_COMP)
//				stopProg("BranchBoundDicho::evalNoeud : res_heur > res_cut");
//
//			if (res_secante > res_cut + EPSILON_COMP)
//				stopProg("BranchBoundDicho::evalNoeud : res_secante > res_cut");
//
//			if (res_cut > res_tangente + EPSILON_COMP)
//				stopProg("BranchBoundDicho::evalNoeud : res_cut > res_tangente");
//#endif
	}



	return res;
}



//calcule les donnees extraites du noeud  pour construire le PL
void BranchBoundDicho::pretraitementPL(const Noeud & noeud, double marge)
{
	int N = _ins->_nbJob;

	_date.clear(); 
	_date.reserve(3 * N);


	//1. on remplit le vecteur date
	for (int i = 0; i < N; ++i)
	{
		_date.push_back(noeud._rel[i]);
		_due_2[i] = min(noeud._due[i], _ins->_LF[i] - marge);
		_date.push_back(_due_2[i]);
	}
	for (double val : noeud._LAdd)
		_date.push_back(val);

	//on trie on ordre croissant et on enleve les doublons
	sort(_date.begin(), _date.end());
	auto it = unique(_date.begin(), _date.end());
	_date.resize(distance(_date.begin(), it));


	//si on a 2 dates consecutives de différence < EPSILON_DATE alors on supprime la plus petite
	//attention ne pas supprimer la plus grande => sinon le calcul de _indRel et _indDue_2 devient faux

	for (int i = 0; i < _date.size() - 1; ++i)
	{
		if (_date[i + 1] - _date[i] < EPSILON_DATE)
		{
			_date.erase(_date.begin()+i);
			i--;
		}
	}


	//on arrondi les dates a 10-6

	//for_each(_date.begin(), _date.end(), [](double & x) {x /= EPSILON_ARRONDI;  x = ceil(x); x *= EPSILON_ARRONDI; });


#ifdef _VERIF_
	for (int i = 0; i < _date.size() - 1; ++i)
	{
		if (_date[i + 1] - _date[i] < EPSILON_DATE)
			stopProg("BranchBoundDicho::pretraitementPL : deux dates consecutives a moins de EPSILON_DATE");
	}
#endif

	//2. on remplit les indices
	double mil1 = 0, mil2 = 0;
	for (int i = 0; i < N; ++i)
	{
		// lower_bound => recherche dicho, retourne un pointeur sur la premiere occurence >= a la valeur cherchée (-v.begin() donne l'indice)
		// on est certain de trouver noeud._rel[i] et _due_2[i] dans date (inutile de verifier)
		_indRel[i] = static_cast<int>(lower_bound(_date.begin(), _date.end(), noeud._rel[i]) - _date.begin());
		_indDue_2[i] = static_cast<int>(lower_bound(_date.begin(), _date.end(), _due_2[i]) - _date.begin());

#ifdef _VERIF_
		//attention, on peut avoir supprime deux dates consecutives a moins de 10-4 et maintemant les 2 dates consecutives 
		//restantes sont a un peu plus que 10-4
		if (abs(_date[_indRel[i]] - noeud._rel[i]) > 10*EPSILON_DATE)
			stopProg("BranchBoundDicho::pretraitementPL : mauvais indice rel");
		if (abs(_date[_indDue_2[i]] - _due_2[i]) > 10*EPSILON_DATE)
			stopProg("BranchBoundDicho::pretraitementPL : mauvais indice due");
#endif

		//pour _pointMil_1[i] et _pointMil_2[i] on n'est pas sur => verifier
		//RQ : si _pointMil_1[i] / _pointMil_2[i] est défini alors normalement ils ont une date qui n'est ni la plus petite (car > rel) ni la plus grande (car < due)
		_indMil_1[i] = static_cast<int>(_date.size());
		mil1 = noeud._pointMil_1[i];



		if (mil1 != -1 && mil1 < _date[_date.size() - 1] )
		{
			if (mil1 < _date[0])
				_indMil_1[i] = 0;
			else
			{
				int k = static_cast<int>(lower_bound(_date.begin(), _date.end(), mil1) - _date.begin());

				if (mil1 < _date[k] - EPSILON)
					k--;
				_indMil_1[i] = k;

#ifdef _VERIF_
				if (k < 0 || k >= _date.size())
					stopProg("BranchBoundDicho::pretraitementPL : indice en dehors des bornes pour _pointMil_1");

				if (_date[k] > mil1 + EPSILON || mil1 >= _date[k + 1] - EPSILON)
					stopProg("BranchBoundDicho::pretraitementPL : mauvais indice pour _pointMil_1");
#endif
			}
		}

		_indMil_2[i] = 0;
		mil2 = noeud._pointMil_2[i];

		if (mil2 != -1 && mil2 > _date[0] + EPSILON)
		{
			if ( mil2 > _date[_date.size()-1] )
				_indMil_2[i] = static_cast<int>(_date.size() - 2);
			else
			{
				int k = static_cast<int>(lower_bound(_date.begin(), _date.end(), mil2 - EPSILON) - _date.begin());
				k--;//on veut mil2 \in ]date[k],date[k+1]]

				_indMil_2[i] = k;

#ifdef _VERIF_

				if (k < 0 || k >= _date.size())
					stopProg("BranchBoundDicho::pretraitementPL : indice en dehors des bornes pour _pointMil_2");

				if (_date[k] >= mil2 - EPSILON || mil2 > _date[k + 1] + EPSILON)
					stopProg("BranchBoundDicho::pretraitementPL : mauvais indice pour _pointMil_2");
#endif
			}

		}
	}





}


//on utilise un PL pour trouver une solution preemptive au noeud courant
//retourne true si PL faisable
bool BranchBoundDicho::resolutionPL(const Noeud & noeud, double marge)
{

	DEBUT:

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

#ifdef _VERIF_

	if (marge < -EPSILON)
		stopProg("BranchBoundDicho::resolutionPL : marge negative " );
#endif

	//1. pretraitement : construit les vecteurs dates et les vecteurs d'indices en fonction du noeud et de la marge en parametre
	pretraitementPL(noeud, marge);



	//2. construction du PL et resolution
	bool res = true;
	int nbDate = static_cast<int> (_date.size());
	int nbJob = _ins->_nbJob, nbArc = _ins->_nbArc;


	//==========================================================
	//variables

	IloArray<IloNumVarArray> debit(env, nbJob); //debit[i][k] = debit du job i entre les dates date[k] et date[k+1]

	for (int i = 0; i < nbJob; ++i)
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



	//1. prise en compte des fenetres de temps [ES,LF-marge] et contraintes de branchement
	for (int i = 0; i < nbJob; ++i)
	{
		
		for (int k = 0; k <= _indRel[i] - 1; ++k)
		{
			model.add(debit[i][k] == 0);
		}
		for (int k = _indDue_2[i]; k < nbDate - 1; ++k)
		{
			model.add(debit[i][k] == 0);
		}
		
		
		for (int k = _indMil_1[i]; k < nbDate - 2; ++k)
		{
			model.add(debit[i][k+1] <= debit[i][k]);
		}
		for (int k = 0; k <= _indMil_2[i]-1; ++k)
		{
			model.add(debit[i][k + 1] >= debit[i][k]);
		}

	}//fin for i 

	//2. on evacue tout le monde
	for (int i = 0; i < nbJob; ++i)
	{
		IloExpr expr(env);
		for (int k = 0; k < nbDate - 1; ++k)
		{
			expr += debit[i][k] * (_date[k + 1] - _date[k]);
		}

		model.add(expr == _ins->_pop[i]);
		expr.end();
	}//fin for i 


	//3. prise en compte des capacites des arcs
	
	for (int a = 0; a < nbArc; ++a)
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

	//normalement ce PL est resolu de maniere quasi immediate, si ce n'est pas le cas 
	// c'est qu'on a des pb numeriques et que cplex ne donnera pas de solution
	cplex.setParam(IloCplex::Param::TimeLimit, 10);



	//_cplex->setParam(IloCplex::Param::ClockType, CLOCK_TYPE_CPLEX);
	
	//_cplex->setParam(IloCplex::Param::MIP::Limits::TreeMemory, RAM_LIMIT_CPLEX);


	//cplex.setParam(IloCplex::Param::Emphasis::Numerical, true);
	//cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 0.0001);
	//cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 0.0001);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	//cplex.setError(env.getNullStream());

	//cplex.exportModel("PL_dicho.lp");



	//cout << "solve..." << endl;
	cplex.solve();
	 
	//cout << "solve OK" << endl;


	// CAS TRES ENNUYEUX : on n'arrive pas a conclure (pb de precision numerique) 
	// => on perturbe un peu la marge (d'une quantite < EPSILON_DEBIT) et on recommence
	if (cplex.getStatus() == IloAlgorithm::Unknown)
	{
		cplex.end();
		model.end();
		env.end();
		
		marge += EPSILON_DEBIT / 10;
		cout << "probleme de precision num => pas de reponse cplex, perturbation de la marge" << endl;
		goto DEBUT;
	}

	
	if (cplex.getStatus() == IloAlgorithm::Infeasible)
	{
		res = false;
	}
	else
	{
		//1. on recupere les dates courantes qui permettent de reconstruire la solution
		_dateCour = _date;

		//2. on recupere la solution qui devient la solution courante du B&B
		//solCour a volontairmeent une case de plus que debit pour toujours terminer sur un debit = 0 (evite les cas particuliers dans le calcul du branchement)
		_solCour.assign(nbJob, vector<double>(nbDate, 0));
		
		for (int i = 0; i < nbJob; ++i)
		{
			double prec = 0.0;
			for (int k = 0; k < nbDate-1; ++k)
			{
				
				if (cplex.isExtracted(debit[i][k]))
				{
					double val = cplex.getValue(debit[i][k]);

					if (val > EPSILON_DEBIT)
					{
						_solCour[i][k] = val;

						//pour eviter les problemes numerioques dans le branchement lies a 2 debits consecutifs "presque" egaux
						//si on a deux debits consecutifs "presque" egaux, on remplace le plus petit par le plus grand

						// = > finalement on garde la solution "telle quelle", 
//						double diff = abs(prec - val);
//						if (diff < EPSILON_DEBIT && prec != val)
//						{
//#ifdef _VERIF_
//							if (k == 0)
//								stopProg("BranchBoundDicho::resolutionPL : on ne devrait pas passer ici avec k = 0");
//#endif 
//
//							if (prec > val)
//								_solCour[i][k] = prec;
//							else
//								_solCour[i][k - 1] = val;
//						}
					}
					prec = val;
				}
			}
		}
	}

	cplex.end();
	model.end();

	env.end();
	return res;

}




//a partir de la solution courante _solCour preemptve (obtenue par le PL) on tente de reconstruire une solution non preemptive
//en faisant la moyenne des debits
//on passe la marge correspond a la sol. courante (utile pour faire des stats uniquement)
bool BranchBoundDicho::reconstruireSolutionCourante(bool stocke, const double marge)
{

	bool ok = true;
	int N = _ins->_nbJob;
	int nbDate = static_cast<int>(_dateCour.size());

	vector<int> deb(N, 0);
	vector<int> fin(N, nbDate);
	vector<double> duree(N);

	//1. on remplit les vecteurs deb, fin, duree
	for (int i = 0; i < N; ++i)
	{
		int k = 0;
		//1.1. calcul plus petit k tq _solCour[i][k] != 0 (existe forcement)
		while (_solCour[i][k] < EPSILON)
			k++;

		deb[i] = k;

		//1.2. calcul plus grand k tq _solCour[i][k] != 0 (existe forcement)
		k = nbDate - 2; //-2 car _solCour[i] index sur 0...nbDate-2
		while (_solCour[i][k] < EPSILON)
			k--;

		fin[i] = k;


		//1.3. duree du job i
		duree[i] = _dateCour[fin[i] + 1] - _dateCour[deb[i]];
	}

	//2. on construite la solution "moyenne"
	vector<vector<double>> solMoy(N, vector<double>(nbDate-1,0));

	for (int i = 0; i < N; ++i)
	{
		for (int k = 0; k < nbDate - 1; ++k)
		{
			if (deb[i] <= k && k <= fin[i])
				solMoy[i][k] = static_cast<double>(_ins->_pop[i]) / duree[i];
		}
	}

	//3. on teste si la solution moyenne est realisable
	for (int a = 0; ok && a < _ins->_nbArc; ++a)
	{
		int nbJobSurA = static_cast<int> (_ins->_arcToJob[a].size());

		for (int k = 0; k < nbDate - 1; ++k)
		{

			double sumDebit = 0;

			for (int j = 0; j < nbJobSurA; ++j)
			{
				int i = _ins->_arcToJob[a][j];
				sumDebit += solMoy[i][k];
			}

			if (sumDebit > _ins->_cap[a] + EPSILON_DEBIT)
				ok = false; //solution non realisable
				
			
		}
	}
	if (ok && _steril1 < marge)
		_steril1 = marge;

	if (ok )
	{
		if (stocke)
		{

			for (int i = 0; i < _ins->_nbJob; ++i)
			{
				_bestDebut[i] = _dateCour[deb[i]];
				_bestDebit[i] = solMoy[i][deb[i]];
				_hasSol = true;
				
			}
#ifdef _VERIF_
			if (!verifieSolBest(marge))
			{
				stopProg("BranchBoundDicho::reconstruireSolutionCourante : sol fausse");
			}
#endif 
		}
	}
	else //commenter si on veut les resultats pour steril2 a chauqe fois
	{
		//2eme tentative de reconstruction en faisant une moyenne "intelligente" 
		
		ok = reconstruireSolCouranteAmeliore(deb, fin, solMoy, marge, stocke);
		if (ok && _steril2 < marge)
			_steril2 = marge;

	}



	return ok;
}

bool BranchBoundDicho::reconstruireSolCouranteAmeliore(const vector<int> & i_deb, const vector<int> & i_fin, const vector<vector<double>> & solMoy, double marge, bool stocke)
{
	bool ok = true;
	bool modif = false;
	int N = _ins->_nbJob;
	int nbDate = static_cast<int>(_dateCour.size());


	//une solution va etre definie par son débit (constant), sa date de debut, sa date de fin
	vector<double> debit(N);
	vector<double> dateDeb(N);
	vector<double> dateFin(N);


	for (int i = 0; i < N; ++i)
	{
		//1.1. par defaut on garde les dates de debut, fin courant et le debit moyen
		dateDeb[i] = _dateCour[i_deb[i]];
		dateFin[i] = _dateCour[i_fin[i] + 1];
		debit[i] = solMoy[i][i_deb[i]];

		//1.2. on precalcule la date qui suit la date de debut, la date qui precede la date de fin et
		// les quantites de populations evacuees sur la premiere et dernier "marche"
		double dateApDeb = _dateCour[i_deb[i] + 1];
		double dateAvFin = _dateCour[i_fin[i]];

		double popDebut = (dateApDeb - dateDeb[i]) * _solCour[i][i_deb[i]];
		double popFin = (dateFin[i] - dateAvFin) * _solCour[i][i_fin[i]];

		//2. on regarde si on peut "ameliore"
		//2.1 la premiere et dernière marche ont un debit inf au debit moyen => on les raccourcis
		if (i_deb[i] + 2 <= i_fin[i] && _solCour[i][i_deb[i]] < debit[i] && _solCour[i][i_fin[i]] < debit[i])
		{
			debit[i] = static_cast<double>(_ins->_pop[i] - popDebut - popFin) / (dateAvFin - dateApDeb);
			dateDeb[i] = dateApDeb - popDebut / debit[i];
			dateFin[i] = dateAvFin + popFin / debit[i];
			modif = true;

		}
		else
		{
			if (i_deb[i] + 1 <= i_fin[i])
			{
				//2.2 seulement la premiere marche a un debit inf au debit moyen
				if (_solCour[i][i_deb[i]] < debit[i])
				{
					debit[i] = static_cast<double>(_ins->_pop[i] - popDebut) / (dateFin[i] - dateApDeb);
					dateDeb[i] = dateApDeb - popDebut / debit[i];
					modif = true;
				}
				else 
				{
					//2.3 seulement la derniere marche a un debit inferieur au debit moyen
					if (_solCour[i][i_fin[i]] < debit[i])
					{
						debit[i] = static_cast<double>(_ins->_pop[i] - popFin) / (dateAvFin - dateDeb[i]);
						dateFin[i] = dateAvFin + popFin / debit[i];
						modif = true;
					}
				}

			}

		}

#ifdef _VERIF_
		if (abs(_ins->_pop[i] - debit[i] * (dateFin[i]- dateDeb[i]) ) > EPSILON)
			stopProg("BranchBoundDicho::reconstruireSolCouranteAmeliore : probleme avec les debits");
#endif
	}

	//3. si la solution a ete modifiee alors on regarde si elle est faisable
	if (modif)
	{
		ok = verifDebit(dateDeb, dateFin, debit);

		if (ok && stocke)
		{
			for (int i = 0; i < _ins->_nbJob; ++i)
			{
				_bestDebut[i] = dateDeb[i];
				_bestDebit[i] = debit[i];
				_hasSol = true;

			}
#ifdef _VERIF_
			if (!verifieSolBest(marge))
			{
				stopProg("BranchBoundDicho::reconstruireSolutionCouranteAmeliore : sol fausse");
			}
			
#endif 
		}

//#ifdef _VERIF_
//		Solution sol(_insInit);
//		sol.init();
//
//		for (int i  = 0; i < N; ++i)
//		{
//			int j = _ins->_job[i];
//
//			sol._debit[j] = debit[i];
//			sol._debut[j] = dateDeb[i] - _ins->_ES[i];
//			sol._fin[j] = dateFin[i] - _ins->_ES[i];
//		}
//
//		if (sol.verifieDebit() != ok)
//		{
//			stopProg("BranchBoundDicho::reconstruireSolCouranteAmeliore : pb avec le test de realisabilité");
//		}
//#endif

	}//fin if modif

	return modif && ok;
}


bool BranchBoundDicho::verifDebit(const vector<double> & dateDeb, const vector<double> & dateFin, const vector<double> & debit) const
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
			stopProg("BranchBoundDicho::verifDebit : les dates ne sont pas dans le bon ordre");
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
					//cout << "debit total = " << debitArc[a] << "capa = " << _ins->_cap[a] << endl;
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


void BranchBoundDicho::calculInfoJob(vector<int> & i_deb, vector<int> & i_fin, vector<double> & duree, vector<double> & debit, bool reduireDebit )
{
	int K = static_cast<int> (_dateCour.size());
	int N = _ins->_nbJob;

	i_deb.assign(N, 0);
	i_fin.assign(N, K);
	duree.resize(N);
	debit.resize(N);


	for (int i = 0; i < N; ++i)
	{
		int k = 0;
		//1. calcul plus petit k tq _solCour[i][k] != 0 (existe forcement)
		while (_solCour[i][k] < EPSILON)
			k++;

		i_deb[i] = k;

		//2. calcul plus grand k tq _solCour[i][k] != 0 (existe forcement)
		k = K - 2; //-2 car _solCour[i] index sur 0...nbDate-2
		while (_solCour[i][k] < EPSILON)
			k--;

		i_fin[i] = k;


		//3. duree du job i
		duree[i] = _dateCour[i_fin[i] + 1] - _dateCour[i_deb[i]];

		//4 debit
		debit[i] = static_cast<double>(_ins->_pop[i]) / duree[i];




		//5. on reduit legerement le debit si besoin
		if (reduireDebit)
		{
			debit[i] -= debit[i] * 0.01;
			duree[i] = static_cast<double>(_ins->_pop[i]) / debit[i];



			//on verifie que le debit est ok par rapport a la capa de l'arc (probleme de precision numeriques avec cplex
			for (int a : _ins->_jobToArc[i])
			{
				if (debit[i] > _ins->_cap[a] + EPSILON_COMP)
				{


#ifdef _VERIF_
					if (debit[i] - _ins->_cap[a] > 0.001)//EPSILON_DEBIT)
					{
						cout << debit[i] - _ins->_cap[a];
						stopProg("  : BranchBoundDicho::calculInfoJob : debit > cap");
					}
#endif // _VERIF_


				}

			}
		}

	}


}



//calcul les jobs qui forment une clique (overlap) ou qui sont l'un avant l'autre en utilisant les debut et fin des jobs de la solution courante
void BranchBoundDicho::calculPredClique(const vector <int> & i_deb, const vector<int> & i_fin, vector <pair<int, int> > & precedence, vector <vector<int> > & clique)
{
	int K = static_cast<int> (_dateCour.size());
	int nbArc = static_cast<int> (_ins->_nbArc);

	vector<vector<int>> jobsDebutent(K); //jobDebutent[k] = ensemble des jobs qui debutent a la date k
	vector<vector<int>> jobsFinissent(K); //jobsFinissent[k] = ensemble des jobs qui finissent a la date k
	clique.resize(K);

	int N = _ins->_nbJob;



	//-----------------------------------------------------------------------
	//1. on repere les jobs qui debutent / finissent en k
	for (int i = 0; i < N; ++i)
	{
		jobsDebutent[i_deb[i]].push_back(i);
		jobsFinissent[i_fin[i] + 1].push_back(i);
	}

	//-------------------------------------------------------------------------
	//2. on calcule les cliques "overlap"

	clique[0] = jobsDebutent[0];

	for (int k = 1; k < K; ++k)
	{
		clique[k] = clique[k - 1];

		//2.1 on enleve les jobs qui se finissent en k
		for (int j : jobsFinissent[k])
		{
			auto it = find(clique[k].begin(), clique[k].end(), j);
			clique[k].erase(it);
		}

		//2.2 on ajoute ceux qui commencent en k
		for (int j : jobsDebutent[k])
			clique[k].push_back(j);
		
		//3. contraintes de precedence : les jobs qui finissent en k sont avant ceux qui commence en k, k+1 ...
		for (int i : jobsFinissent[k])
		{
			for (int l = k; l < K; ++l)
			{
				for (int j : jobsDebutent[l])
					precedence.push_back({ i,j });
			}
		}
	}


}

//idem en utilisant un PL basé sur l'ordre des jobs donne par la solution preemptive
double BranchBoundDicho::reconstruireSolCourantePL_heur(const vector<int> & i_deb, const vector<int> & i_fin, const vector<double> & duree, const vector<double> & debit)
{

	int K = static_cast<int> (_dateCour.size());
	int nbArc = static_cast<int> (_ins->_nbArc);

	vector<vector<int>> jobsDebutent(K); //jobDebutent[k] = ensemble des jobs qui debutent a la date k
	vector<vector<int>> jobsFinissent(K); //jobsFinissent[k] = ensemble des jobs qui finissent a la date k


	int N = _ins->_nbJob;



	//-----------------------------------------------------------------------
	//1. on repere les jobs qui debutent / finissent en k
	for (int i = 0; i < N; ++i)
	{
		jobsDebutent[i_deb[i]].push_back(i);
		jobsFinissent[i_fin[i] + 1].push_back(i);
	}

	//-------------------------------------------------------------------------
	//2. on calcule les cliques "overlap"
	vector<vector<int> > clique(K); //clique[k] = ensemble des jobs qui se croisent a la date k

	clique[0] = jobsDebutent[0];

	for (int k = 1; k < K; ++k)
	{
		clique[k] = clique[k - 1];

		//2.1 on enleve les jobs qui se finissent en k
		for (int j : jobsFinissent[k])
		{
			auto it = find(clique[k].begin(), clique[k].end(), j);
			clique[k].erase(it);
		}

		//2.2 on ajoute ceux qui commencent en k
		for (int j : jobsDebutent[k])
			clique[k].push_back(j);
		
	}

	//-----------------------------------------------------------------------
	//3. on utilise les cliques pour construire le PL
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);


	//chaque debit est multiplié par un coeff reducteur, on maximise le min de ce coeff
	IloNumVarArray coeff(env, N); //coeff[i] = coeff reducteur du job i 
	IloNumVar lambda(env, 0, 1, IloNumVar::Float, "lambda");

	//3.1 deifnition des variables
	for (int i = 0; i < N; ++i)
	{
		string name = "c_" + to_string( i );

		coeff[i] = IloNumVar(env, 0, 1, IloNumVar::Float, name.c_str());

		model.add(coeff[i] - lambda >= 0);
	}


	//3.2 capa des arcs
	for (int k = 0; k < K; ++k)
	{
		
		for (int a = 0; a < nbArc; ++a)
		{
			IloExpr expr(env);
			int nbJobSurA = static_cast<int> (_ins->_arcToJob[a].size());
			int cpt = 0;
			for (int j = 0; j < nbJobSurA; ++j)
			{
				int i = _ins->_arcToJob[a][j];

				//on ajoute la ctr pour les jobs dans clique[k] qui utilisent a
				if (find(clique[k].begin(), clique[k].end(), i) != clique[k].end())
				{
					expr += coeff[i] * debit[i];
					cpt++;
				}
			}
			if (cpt > 0)
			{
				//string name_ctr = "C" + to_string(k) + "_" + to_string(a);
				//IloRange ctr(env, 0, expr, _ins->_cap[a], name_ctr.c_str());
				//model.add(ctr);

				model.add(expr <= _ins->_cap[a]);
			}

			expr.end();
		}
	}//fin for a

	//3.3 objectif : avoir des coeff les plus grands possible => lambda le plus grand possible

	model.add(IloMaximize(env, lambda));


	//3.4 resolution
	cplex.setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	cplex.setOut(env.getNullStream());
	//cplex.exportModel("recons3.lp");

	//cout << "solve..." << endl;
	cplex.solve();

	//cout << "solve OK" << endl;

	if (cplex.getStatus() != IloAlgorithm::Optimal)
		stopProg("BranchBoundDicho::reconstruireSolCourantePL : devrait fournir la solution optimale");


	//----------------------------------------------------------------------------
	//4. on utilise le PL pour reconstruire une solution => calcul les nouveaux debut et fin avec les nouveau debit

	vector<double> nouvDebit(N);
	vector<double> longueur(N);
	vector<bool> marq(N, false);

	vector<double> debut(N,0);
	vector<double> fin(N);

	for (int i = 0; i < N; ++i)
	{
		debut[i] = _dateCour[i_deb[i]]; //debut au plus tot, maj avec Bellman dans la suite
		nouvDebit[i] = cplex.getValue(coeff[i]) * debit[i];
		longueur[i] = static_cast<double> (_ins->_pop[i]) / nouvDebit[i];
	}

	for (int k = 0; k < K; ++k)
	{
		for (int j : jobsFinissent[k])
			marq[j] = true;

		for (int j : jobsDebutent[k])
		{
			for (int i = 0; i < N; ++i)
				if (marq[i])
					debut[j] = max(debut[j], debut[i] + longueur[i]);

			fin[j] = debut[j] + longueur[j];
		}
	}


#ifdef _VERIF_

	for (int j = 0; j < N; ++j)
		if (abs(fin[j] - debut[j] - longueur[j]) > EPSILON_COMP)
			stopProg("BranchBoundDicho::reconstruireSolCourantePL : fin et debut non coherent");
#endif
	//-------------------------------------------------------------------------------------------------------
	//5. on calcule la marge ainsi obtenue => si elle est < 0 alors la solution construite n'est pas realisable

	double marge = TIME_INFINITY;


	for (int i = 0; i < N; ++i)
		marge = min(marge, _ins->_LF[i] - fin[i]);


	if (marge > 0 && _bestMargeRealisable < marge)
	{
		
		cout << "on gagne " << marge - _bestMargeRealisable << " avec le PL" << endl;
		_bestMargeRealisable = marge;
		cout << " ***** sterilisation PL = " << _bestMargeRealisable << endl;
		
	}

#ifdef _VERIF_
	if (_coutRacine != -1 && marge > _coutRacine + EPSILON)
		stopProg("BranchBoundDicho::reconstruireSolCourantePL : sol meilleure que BS");

	if (marge > 0)
	{
		if (!verifDebit(debut, fin, nouvDebit))
		{
			stopProg("BranchBoundDicho::reconstruireSolCourantePL : sol fausse");
		}
	}
#endif

	cplex.end();
	model.end();
	env.end();

	return marge;
}


//idem reconstruireSolCourantePL_heur : utilise  l'ordre des jobs donne par la solution preemptive pour constrire une solution
//mais on resout de maniere "quasi" exact en approximant la ctr non lineaire par ses tangentes (version optimiste) ou par ses secantes (version pessimiste)
// ou avec une genration de coupes (version exacte mais risque de degenerescence)
//precedence donne la liste des couples (i,j) avec i << j
//cliques donne l'ensemble des cliques
//pasConstant = true si o nutilise un pas constant, false sinon

double BranchBoundDicho::reconstruireSolCourantePL_exact(const vector<pair<int,int> > & precedence, const vector<vector<int> > & cliques,
	const TYPE_RELAX & relax, bool pasConstant)
{
	double res = -1.0;

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);


	int nbJob = _ins->_nbJob, nbArc = _ins->_nbArc;


	//==========================================================
	//variables

	IloNumVarArray debit(env, nbJob); //debit[i] = debit du job i 
	IloNumVarArray debut(env, nbJob); //debut[i] = debut du job i 
	IloNumVar margeMin = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, "M"); //marge minimum => quantite a maximiser

	for (int i = 0; i < nbJob; ++i)
	{
		string nomDebit = "D_" + to_string(i);
		string nomDebut = "T_" + to_string(i);

		debit[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, nomDebit.c_str());
		debut[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, nomDebut.c_str());

	}//fin for i 


	//==========================================================
	//contraintes temporelles

	
	int nbPoint = 5; // 5*_ins->_nbJob;//nombre de tangentes par job (pour ctr 2. et 3.)
	double S = 1; // S = 1 + 1/2 + 1/3 +...+ 1/(nbPoint-1) utile si pas non constant

	if (!pasConstant)
	{
		for (int k = 2; k < nbPoint; ++k)
			S += 1.0 / k;
	}


	for (int i = 0; i < nbJob; ++i)
	{

		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->_debitMax[i] - _ins->_debitMin[i]) / S;

		//1. debut au plus tot
		model.add(debut[i] >= _ins->_ES[i]);


		//2. la marge est plus petite ou egal que la fin des jobs : 
		//linearisation de la contrainte "dead(i) - debut(i) - pop(i)/debit(i) >= marge
		// => on approxime par la tangente
		
		double pas = (_ins->_debitMax[i] - _ins->_debitMin[i]) / (nbPoint - 1);
		if (!pasConstant)
			pas = y / (nbPoint - 1);

		double x = _ins->_debitMin[i]; //abscisse auquel on calcul la tangente
		double x2 = x + pas;

		for (int k = 0; k < nbPoint; ++k)
		{
			if (relax == SECANTE)
			{
				double a = _ins->_pop[i] * (1.0 / x - 1.0 / x2) / (x - x2);
				double b = _ins->_pop[i] * (1.0 / x - (1.0 - x / x2) / (x - x2));
				model.add(a * debit[i] + margeMin + debut[i] <= _ins->_LF[i] - b );
			}
			else //TANGENTE ou CUT (cut => on commence avec les tangentes, on ajoute les cut dynamiquement ensuite)
			{
				double a = _ins->_pop[i] / x / x; //pop / x^2
				model.add(a * debit[i] - margeMin - debut[i] >= 2 * _ins->_pop[i] / x - _ins->_LF[i]);
			}

			//abscise pour la prochaine tangente
			x += pas;
			x2 += pas;

			if (!pasConstant && k != nbPoint-2)
				pas = y / (nbPoint - k - 2);
		}

	}

	//3. constraintes de precedence

	for (pair<int, int> p : precedence)
	{
		//p.first << p.second

		int i = p.first, j = p.second;

		double pas = (_ins->_debitMax[i] - _ins->_debitMin[i]) / (nbPoint - 1);

		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->_debitMax[i] - _ins->_debitMin[i]) / S;

		if (!pasConstant)
			pas = y / (nbPoint - 1);

		double x = _ins->_debitMin[i]; //abscisse auquel on calcul la tangente
		double x2 = x + pas;

		for (int k = 0; k < nbPoint; ++k)
		{

			if (relax == SECANTE)
			{
				double a = _ins->_pop[i] * (1.0 / x - 1.0 / x2) / (x - x2);
				double b = _ins->_pop[i] * (1.0 / x - (1.0 - x / x2) / (x - x2));
				model.add(debut[j] - debut[i] - a * debit[i] >= b);
			}
			else //TANGENTE ou CUT (cut => on commence avec les tangentes, on ajoute les cut dynamiquement ensuite)
			{
				double a = _ins->_pop[i] / x / x; //pop / x^2
				model.add(debut[j] - debut[i] + a * debit[i] >= 2 * _ins->_pop[i] / x);
			}

			x += pas;
			x2 += pas;

			if (!pasConstant && k != nbPoint - 2)
				pas = y / (nbPoint - k - 2);
		}

	}

	//==========================================================
	//contraintes de capacité

	//4. debits min et max
	for (int i = 0; i < nbJob; ++i)
	{
		model.add(debit[i] >= _ins->_debitMin[i]);
		model.add(debit[i] <= _ins->_debitMax[i]);
	}



	//5. ctr de capacité pour chaque clique

	int nbClique = static_cast<int>(cliques.size());
	bool nonVide = false;


	for (int a = 0; a < _ins->_nbArc; ++a)
	{
		for (int c = 0; c < nbClique; ++c)
		{
			IloExpr expr(env);

			for (int i : cliques[c])
			{
				if (_ins->_utilise[i][a])
				{
					expr += debit[i];
					nonVide = true;
				}
			}

			if (nonVide)
				model.add(expr <= _ins->_cap[a]);
			nonVide = false;
		}
	}


	//6. fonction objectif

	model.add(IloObjective(env, margeMin, IloObjective::Maximize, "obj"));


	//cplex.exportModel("PL_tangente.lp"); 
	
	
	cplex.setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	cplex.setOut(env.getNullStream());
	//cplex.setParam();

	cplex.solve();

	if (cplex.getStatus() == IloAlgorithm::Optimal)
	{
		res = cplex.getObjValue();


		//================================================================================
		//7. si on est dans le cas ou on genere des coupes dynamiquement
		if (relax == CUT)
		{
			//cout << "cutting..." << endl;

			bool fin = false;
			int cpt = 0;
			_nbAppelCut++;


			//vecteurs pour stocker la solution
			IloNumArray sol_debut(env, nbJob);
			IloNumArray sol_debit(env, nbJob);

			while (!fin)
			{
				bool ajout_ctr = false;
				_nbCut++;
				cpt++;
				
				//au cas ou on boucle sans s'arreter => on renvoie qu'on ne trouve pas de solution (a priori n'arrive jamais)
				if (cpt > 1000)
					return -1.0;


				//on recupere la solution pour pouvoir ajouter plusieurs ctr d'un coup
				cplex.getValues(debut, sol_debut);
				cplex.getValues(debit, sol_debit);

				//7.1 test contrainte "marge <= fin"
				for (int i = 0; i < nbJob; ++i)
				{
					

					//on teste la contrainte initiale (non lineaire)...
					if (_ins->_LF[i] - sol_debut[i] - static_cast<double>(_ins->_pop[i]) / sol_debit[i] < res - EPSILON_CUT)
					{
						//... si elle n'est pas respectee alors on ajoute la tangente au point x = debit[i]

						double x = sol_debit[i];
						double a = _ins->_pop[i] / x / x; //pop / x^2
						model.add(a * debit[i] - margeMin - debut[i] >= 2 * _ins->_pop[i] / x - _ins->_LF[i]);
						ajout_ctr = true;
					}
				}

				//7.2 test contrainte de precedence
				for (pair<int, int> p : precedence)
				{
					//p.first << p.second

					int i = p.first, j = p.second;

					if (sol_debut[j] < sol_debut[i] + static_cast<double>(_ins->_pop[i]) / sol_debit[i] - EPSILON_CUT)
					{
						double x = sol_debit[i];
						double a = _ins->_pop[i] / x / x; //pop / x^2
						model.add(debut[j] - debut[i] + a * debit[i] >= 2 * _ins->_pop[i] / x);
						ajout_ctr = true;
					}
					

				}


				//7.3 on relance la resolution avec le simplex reseau
				fin = true;
				if (ajout_ctr)
				{
					cplex.solve();
					res = -1;
					
					if (cplex.getStatus() == IloAlgorithm::Optimal)
					{
						res = cplex.getObjValue();
						fin = false;


					}
				}

			}//fin whlie (!fin)

			 //cout << " ==== nb tours cut = " << cpt << endl;

			//REMARQUE : on pourrait aussi sauvegarder dans le cas des secantes, mais on sait que les cuts sont meilleures
			if (res > _bestMargeRealisable)
			{
				//cout << "PL cut / secante ameliore : " << "best cour = " << _bestMargeRealisable << " ; res = " << res << endl;
				_bestMargeRealisable = res;
				for (int k = 0; k < _ins->_nbJob; ++k)
				{
					_bestDebit[k] = sol_debit[k];
					_bestDebut[k] = sol_debut[k];
					_hasSol = true;
				}
#ifdef _VERIF_
				if (!verifieSolBest(res))
					stopProg("BranchBoundDicho::reconstruireSolCourantePL_exact : sol fausse !!");
#endif
			}

	
		}//fin if (relax == CUT)
		//=============================================================================
	}




	cplex.end();
	model.end();
	env.end();


	return res;
}


//utilise les debits moyens courants pour essayer de reconstruire une solution avec un algo
// glouton 
double BranchBoundDicho::reconstruireSolHeuristique(const vector<double> & debit, const vector<double> & longueur)
{
	double marge = -1;

	int N = _ins->_nbJob;
	int A = _ins->_nbArc;

	int nbJobRestant = N; //nb de jobs qu'il reste a planifier

	double tCour = 0; // date courante => date d'insertion du prochain job
	vector<bool> marq(N, false);
	vector<double> debut(N, -1);
	vector<double> fin(N);
	vector<bool> enCours(N, false); //enCours[i] = vrai si le job i est en cours en tCour
	vector<double> charge(A, 0); //charge[a] = quantité de ressource de type "arc a" consommée a la dtae tCour

	vector<int> jobCour;
	jobCour.reserve(N);

	bool fini = false;

	//1. init du tCour au min des ES
	tCour = _ins->_ESmin;


	while (!fini && nbJobRestant > 0)
	{
		jobCour.clear();

		//2. calcul des jobs qui peuvent commencer en tCour
		for (int i = 0; i < N; ++i)
		{
			double debit_i = debit[i];
			if (_ins->_ES[i] <= tCour + EPSILON && !marq[i])
			{
				//on verifie qu'il y a assez de ressource pour i : 
				bool ok = true;
				for (int a : _ins->_jobToArc[i] )
				{
					if (debit_i > _ins->_cap[a] - charge[a] + EPSILON_COMP)
					{
						ok = false;
						break;
					}
				}

				if(ok)
					jobCour.push_back(i);
			}
		}

		//3. si aucun job ne peut commencer en tCour

		if (jobCour.empty())
		{
			//3.1 calcul de la prochaine date interessante
			double next = TIME_INFINITY;
			for (int i = 0; i < N; ++i)
			{
				if (enCours[i])
					next = min(next, fin[i]);
				if (tCour < _ins->_ES[i] - EPSILON)
					next = min(next, _ins->_ES[i]);
			}
			tCour = next;

			//3.2 mise a jour des donnees
			for (int i = 0; i < N; ++i)
			{
				//jobs qui se finissent en tCour
				if (enCours[i] && fin[i] <= tCour + EPSILON)
				{
					enCours[i] = false;
					for (int a : _ins->_jobToArc[i])
						charge[a] -= debit[i];					
				}
					
			}
		}
		else //ensemble des jobs non vide => on en choisit un et on le place en tCour
		{
			//if (jobCour.size() > 1)
			//	cout << jobCour.size() << endl;

			//on choisit un job
			int j = choisir(tCour, jobCour, longueur, marq, debit);

			if (j == -1) //signal d'echec : au moins un job ne peut plus respecter sa deadline
				fini = true;
			else
			{
				marq[j] = true;
				debut[j] = tCour;
				fin[j] = tCour + longueur[j];
				enCours[j] = true;
				nbJobRestant--;
				for (int a : _ins->_jobToArc[j])
					charge[a] += debit[j];
			}
		}


	}//fin while

	if (nbJobRestant == 0) //on a place tous les jobs => calcul de la marge
	{
		marge = TIME_INFINITY;

		for (int i = 0; i < N; ++i)
		{
			marge = min(marge, _ins->_LF[i] - fin[i]);
		}


		if (marge > _bestMargeRealisable)
		{
			cout << " ***** amelioration avec heur = " << marge - _bestMargeRealisable << endl;

			_bestMargeRealisable = marge;
		}

#ifdef _VERIF_
		if (!verifDebit(debut, fin, debit))
			stopProg("BranchBoundDicho::reconstruireSolHeuristique : sol fausse");

		for (int i = 0; i < N; ++i)
			if (debut[i] < -EPSILON)
				stopProg("BranchBoundDicho::reconstruireSolHeuristique : des jobs non init");

#endif
	}

	return marge;
}


//on choisit un job parmi l'ensmeble jobs (utilise dans reconstruireSolHeuristique)
int BranchBoundDicho::choisir(double tCour, const vector<int> & jobs, const vector<double> & longueur, 
	const vector<bool> & marq, const vector<double> & debit)
{
	//cas partivculier : un seul job => on le retourne direct

	if (jobs.size() == 1)
		return jobs[0];


	//cas general : plusieurs jobs :

	double alpha = 1.0;

	//choix au hasard parmi les 2 meilleurs
	int iBest = -1;
	int iSecondBest = -1;
	double margeMin = TIME_INFINITY;
	double congestionMax = 0;

	for (int i : jobs)
	{
		double margeCour = _ins->_LF[i] - tCour - longueur[i];
		if (margeCour < -EPSILON)
			return -1; //echec


		//double congestion = 0;

		//for (int a : _ins->_jobToArc[i])
		//{
		//	double sum = 0;
		//	for (int k : _ins->_arcToJob[a])
		//		if (!marq[k])
		//			sum += debit[k];

		//	congestion = max(congestion, longueur[i] * sum/ _ins->_cap[a]);
		//}

		if (margeCour < margeMin)
		{
	
			margeMin = margeCour;
			iSecondBest = iBest;
			iBest = i;
		}

		//if (congestion > congestionMax)
		//{

		//	congestionMax = congestion;
		//	iSecondBest = iBest;
		//	iBest = i;
		//}

	}

	//on choisit le deuxieme meilleur avec proba 0.5
	//if (rand() % 2 == 1)
	//	iBest = iSecondBest;

	return iBest;
}

//modifie le noeud fils en fonction de sa signature et de b (numero de branche = 0...2)
void BranchBoundDicho::appliqueSignatureEtBranch(Noeud & fils, int b)
{
	const SignatureBranchement & S = fils._sign;


	switch (b)
	{
	case 0:
		fils._rel[S._job] = S._alpha;
		break;

	case 1:
		fils._due[S._job] = S._beta;
		break;

	case 2: //la troisieme branche depend du type
		if (fils._pointMil_1[S._job] != -1)
			fils._pointMil_1[S._job] = min(S._alpha, fils._pointMil_1[S._job]);
		else
			fils._pointMil_1[S._job] = S._alpha;

		if (fils._pointMil_2[S._job] != -1)
			fils._pointMil_2[S._job] = max(S._beta, fils._pointMil_2[S._job]);
		else
			fils._pointMil_2[S._job] = S._beta;

		switch (S._type)
		{
		case 1://rien a faire
			break;
		case 2:
			fils._LAdd.push_back(S._alpha);
			break;
		case 3:
			fils._LAdd.push_back(S._beta);
			break;
		default:
			stopProg("BranchBoundDicho::appliqueSignatureEtBranch :type impossible");

		}
		break;
	default:
		stopProg("BranchBoundDicho::appliqueSignatureEtBranch : valeur de b impossible");
	}


}


//calcule la signature de branchement pour la solution courante associee au vecteur date courant
BranchBoundDicho::SignatureBranchement BranchBoundDicho::calculParamBranchement()
{
	int N = _ins->_nbJob;
	SignatureBranchement bestSign;
	

	double bestCreux = -1; //
	double bestEscalier = -1;
	double bestType2ou3 = -1;


	for (int i = 0; i < N; ++i)
	{
		//1. on recherche les jobs en creux
		//pour les jobs en creux on calcule la profondeur du creux, on branche en priorité sur la plus grande profondeur
		SignatureBranchement signCour;
		
		//double valCour = calculCreux(i, signCour);
		
		double valCour = calculCreuxAmeliore(i, signCour);

		if (valCour > bestCreux)
		{
			bestSign = move(signCour);
			bestCreux = valCour;

		}
		//on calcule les jobs en escalier si on n'a pas de job en creux
		if (bestCreux < 0)
		{
			//double val = calculEscalierPositif(i, signCour);
			//if (val > bestEscalier)
			//{
			//	bestEscalier = val;
			//	bestSign = signCour;//pas de move ici car on reutilise signCour
			//}
			//val = calculEscalierNegatif(i, signCour);
			//if (val > bestEscalier)
			//{
			//	bestEscalier = val;
			//	bestSign = move(signCour);
			//}


			//if (bestEscalier > 0)
			//	stopProg("BranchBoundDicho::SignatureBranchement : escalier ne devrait pas trouve si creux ameliore ne trouve pas");

			//on ne trouve pas d'escalier non plus alors on cherche un branchement de type 2 ou 3
			if (bestEscalier < 0)
			{
				double val = 0.0;
				vector<double> solGroup, dateGroup;
				solGroup.reserve(_dateCour.size());
				dateGroup.reserve(_dateCour.size());
				//on regroupe les valeurs identiques consecutives dans _solCour[i] (evite de devoir traiter les cas = dans la recherche des param de branchement)
				regroupe(i, solGroup, dateGroup);

				val = calculBranchType2ou3(i, solGroup, dateGroup, signCour);
				if (val > bestType2ou3)
				{		
#ifdef _VERIF_
					if (signCour._alpha == 0 || signCour._beta == -1)
						stopProg("BranchBoundDicho::calculSignatureBranchCour : probleme type 2 ou 3");
#endif
					bestType2ou3 = val;
					bestSign = move(signCour);
				}
			}
		}

	}




	//cout << "branchement type " << bestSign._type << "(" << bestSign._alpha << "," << bestSign._beta << ")" << endl;


//#ifdef _VERIF_

	//REMARQUE : il est possible que la sol soit non faisable mais qu'on ne trouve pas de branchmement : 
	//ca arrive dans le cas où :
	//1. on a des escaliers mais les marches ont une différence d'hauteur < 2* EPSILON_DEBIT => cplex va surment faire la moyenne et le branchement n'aura pas d'effet a  EPSILON_DEBIT pres
	//2. on a un profil ++- ou +-- tel que si on calcul un nouvel alpha ou un nouveau beta on aura au moins un intervalle de temps < EPSILON_DATE => les dates seront fusionnes et au moins une branche n'aura pas d'effet

//	if (bestSign._job == -1 || bestSign._alpha == 0 || bestSign._beta == -1)
//	{
//		cout << "ERREUR ??? : PAS DE BRANCHEMENT TROUVE" << endl;
//		cout << "sol. courante realisable ? " << reconstruireSolutionCourante(false) << endl;
//
//		afficheSolCour();		
//	}
//#endif


	return bestSign;
}



//calcul le creux le plus profond pour le job i s'il y existe, renvoie sa profondeur (-1 si n'existe pas)
// et met a jour sign en consequence
double BranchBoundDicho::calculCreux(int i, BranchBoundDicho::SignatureBranchement & sign)
{

	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (la derniere date ets fictive pour avoir toujours solCour[K-1] = toujours 0)
	const vector<double> & sol = _solCour[i];


	int k = 0, kplus = -1, kmoins = -1;
	double bestVal = -1;

	double debitPrec = 0;
#ifdef _VERIF_
	if (sol[K - 1] > EPSILON)
		stopProg("BranchBoundDicho::calculCreux : sol doit finir par 0");
#endif

	//il peut y avoir plusieurs creux, on cherche le plus profond
	while (k < K)
	{
		//prochaine rupture negative
		while (k < K && sol[k] + EPSILON >= debitPrec )
		{
			debitPrec = sol[k];
			k++;
		}
		if (k < K)//on a trouve une rupture negative
		{
			kmoins = k;
			double hauteurCour = 0;

			//tant qu'on descend on ajoute la difference de debit a la profondeur
			while (k < K && sol[k] <= debitPrec + EPSILON)
			{
				hauteurCour += debitPrec - sol[k];
				debitPrec = sol[k];
				k++;
			}

			//on compte maintenant la remontee (si elle existe)
			int kSav = k;
			while (k < K && sol[k] + EPSILON >= debitPrec )
			{
				hauteurCour +=  sol[k] - debitPrec;
				debitPrec = sol[k];
				k++;
			}

			if (k > kSav) //on a trouvee une remontee => on a bien un creux
			{
				kplus = k - 1; //la derniere rupture positive etait sur le k precedent


				double surfaceCour = hauteurCour * (_dateCour[kplus] - _dateCour[kmoins]);
				if (surfaceCour > bestVal)
				//if (hauteurCour > bestVal)
				//if (_dateCour[kplus] - _dateCour[kmoins] > bestVal)
				{
					sign._alpha = _dateCour[kmoins];
					sign._beta = _dateCour[kplus];
					//bestVal = hauteurCour;
					//bestVal = _dateCour[kplus] - _dateCour[kmoins];
					bestVal = surfaceCour;
#ifdef _VERIF_
					verifieBranchementValide(i, kmoins, kplus);
					if (bestVal < 0)
						stopProg("BranchBoundDicho::calculCreux : best val < 0");
#endif
				}
			}

		}

	}//fin while principal

	sign._type = 1;
	sign._job = i;


#ifdef _VERIF_
	if (bestVal > -0.999 && sign._alpha == 0)
		stopProg("BranchBoundDicho::calculCreux : alpha = 0");
#endif

	return bestVal;

}


//calcul le branchement le plus large possible (on cherche le 1er et le dernier k | sol[k]!=0, si la solution n'est pas constante 
//entre k1+1 et k2 on utilise les alpha beta associes pour faire un branchmeent de type 1
// et met a jour sign en consequence
//contient le branhcmeent creux dans le sens ou si creuxAmeliore n'existe pas alors creux non plus
double BranchBoundDicho::calculCreuxAmeliore(int i, BranchBoundDicho::SignatureBranchement & sign)
{

	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (la derniere date ets fictive pour avoir toujours solCour[K-1] = toujours 0)
	const vector<double> & sol = _solCour[i];


	int k = 0, kplus = -1, kmoins = -1;
	double gap = -1, denivele = -1;

#ifdef _VERIF_
	if (sol[K - 1] > EPSILON)
		stopProg("BranchBoundDicho::calculCreux : sol doit finir par 0");
#endif



	while (sol[k] < EPSILON_DEBIT)
		k++;

	kmoins = k; //1er k | sol[k] != 0

	k = K - 1;
	while (sol[k] < EPSILON_DEBIT)
		k--;

	kplus = k; //dernier k | sol[k] != 0

	//il faut au moins un ecart de 2 pour que le branchement soit valide
	if (kmoins + 2 <= kplus)
	{
		//il faut aussi que la solution ne soit pas constante entre kmoins+1 et kplus ou qu'elle soit inferieure
	

		double deniveleInterne = 0;
		for (int j = kmoins + 1; j < kplus - 1; ++j)
		{
			deniveleInterne += abs(sol[j + 1] - sol[j]);
		}

		//on impose d'avoir au moins 2*EPSILON_DEBIT au milieu de alpha et beta, sinon la 3eme branche ne coupera pas suffisamment la solution
		double EPSILON_FORT = 2 * EPSILON_DEBIT;

		if ( deniveleInterne > EPSILON_FORT || sol[kmoins] > sol[kmoins + 1] + EPSILON_FORT || sol[kplus] > sol[kplus - 1] + EPSILON_FORT) //solution non constante => on peut brancher ou de debit inferieur entre kmoins et kplus
		{

			denivele = deniveleInterne;
		
			//on compte dans le premier denivele si c'est une descente (on veut favoriser les creux)
			//idem pour le dernier
			denivele += max(0.0, sol[kmoins] - sol[kmoins + 1]);
			denivele += max(0.0, sol[kplus] - sol[kplus - 1]);

			
			kmoins++;
			
			//on essaie d'avancer alpha et reculer beta (chacun leur tour pour ne pas favoriser un branche) pour que les branches 1 et 2 coupent +
			bool chg = true;
			while (chg && kmoins + 1 < kplus)
			{
				chg = false;
				if (sol[kmoins] > sol[kmoins + 1] + EPSILON_DEBIT)
				{
					kmoins++;
					chg = true;
				}
				if (kmoins + 1 < kplus && sol[kplus-1] > sol[kplus - 2] + EPSILON_DEBIT)
				{
					chg = true;
					kplus--;
				}
			}

			sign._alpha = _dateCour[kmoins];
			sign._beta = _dateCour[kplus];
			gap = sign._beta - sign._alpha;


#ifdef _VERIF_
			verifieBranchementValide(i, kmoins, kplus);
			if (gap < 0)
				stopProg("BranchBoundDicho::calculCreuxAmeliore : le gap doit etre positif");
#endif
		}
	}

	sign._type = 1;
	sign._job = i;




	return denivele;

}


//calcul l'escalier  le plus haut pour le job i s'il y existe, renvoie sa hauteur (-1 si n'existe pas)
// et met a jour sign en consequence
double BranchBoundDicho::calculEscalierPositif(int i, BranchBoundDicho::SignatureBranchement & sign)
{

	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (solCour[K-1] = toujours 0)
	const vector<double> & sol = _solCour[i];


	int k = 0, k1 = -1, k2 = -1; //k1 = debut de la deuxieme marche et k2 = debut de la derniere marche
	double bestHauteur = -1;

	double debitPrec = 0;
#ifdef _VERIF_
	if (sol[K - 1] > EPSILON_DEBIT)
		stopProg("BranchBoundDicho::calculCreux : sol doit finir par 0");
#endif

	//il peut y avoir plusieurs escaliers, on cherche le plus haut
	while (k < K)
	{

		//on avance k jusqu'a la prochaine rupture strictement positive
		while (k < K && debitPrec >= sol[k] - EPSILON_DEBIT)
		{
			debitPrec = sol[k];
			k++;
		}
		//k est au debut de la premiere marche => 2eme rupture positive
		k1 = k + 1;
		int cpt = 0;
		//on avance k tant qu'on monte
		while (k < K && debitPrec <= sol[k] + EPSILON_DEBIT)
		{
			if (debitPrec + EPSILON_DEBIT < sol[k])
				cpt++;
			debitPrec = sol[k];
			k++;
		}
		k2 = k - 1;

		if (k2 > k1 && cpt >= 3)//on a au moins 3 ruptures > 0
		{
			double h = sol[k2] - sol[k1 - 1];
			//rq : on ne compte pas la premiere marche sinon on risque d'avoir un plateau (+ =...= - )(et non un escalier)
			if (h > max(bestHauteur, EPSILON) )
			{
				bestHauteur = h;
				sign._alpha = _dateCour[k1];
				sign._beta = _dateCour[k2];

#ifdef _VERIF_
				verifieBranchementValide(i, k1, k2);
#endif
			}
		}


	}//fin while principal

	sign._type = 1;
	sign._job = i;

	return bestHauteur;
}


//calcul l'escalier  le plus haut pour le job i s'il y existe, renvoie sa hauteur (-1 si n'existe pas)
// et met a jour sign en consequence
//rq : fonction appelee uniquement si pas de creux => l'escalier calcule est faux si la solution contient des creux
double BranchBoundDicho::calculEscalierNegatif(int i, BranchBoundDicho::SignatureBranchement & sign)
{

	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (solCour[K-1] = toujours 0)
	const vector<double> & sol = _solCour[i];


	int k = 0, k1 = -1, k2 = -1; //k1 = debut de la deuxieme marche et k2 = debut de la derniere marche
	double bestHauteur = -1;

	double debitPrec = 0;
#ifdef _VERIF_
	if (sol[K - 1] > EPSILON_DEBIT)
		stopProg("BranchBoundDicho::calculCreux : sol doit finir par 0");
#endif

	//il peut y avoir plusieurs escaliers, on cherche le plus haut
	while (k < K)
	{

		//on avance k jusqu'a la prochaine rupture strictement negative
		//on passe aussi les paliers a 0
		while (k < K && (sol[k] < EPSILON_DEBIT || debitPrec <= sol[k] + EPSILON_DEBIT) )
		{
			debitPrec = sol[k];
			k++;
		}
		//k est au debut de la premiere marche
		k1 = k;

		//on avance k tant qu'on descend (ou egal et != 0)
		while (k < K && sol[k] > EPSILON_DEBIT && debitPrec >= sol[k] - EPSILON_DEBIT)
		{
			debitPrec = sol[k];
			k++;
		}
		k2 = k - 1;

		if (k2 > k1)//on a au moins 3 marches consecutives
		{
			double h = sol[k1 - 1] - sol[k2];
			if (h > max (bestHauteur,EPSILON) )
			{
				bestHauteur = h;
				sign._alpha = _dateCour[k1];
				sign._beta = _dateCour[k2];

#ifdef _VERIF_
				verifieBranchementValide(i, k1, k2);
#endif
			}
		}


	}//fin while principal

	sign._type = 1;
	sign._job = i;

	return bestHauteur;
}

//calcul les paramètres de branchmeent pour les jobs de type ++- ou +--
//on utilise les vecteurs sol et date en parametre dans lesquels les valeurs de debit identiques et consecutives ont ete regroupees
double BranchBoundDicho::calculBranchType2ou3(int i, const vector<double> & sol, const vector<double> & date, BranchBoundDicho::SignatureBranchement & sign)
{

	//==========================================
	// attention cette fonction utilise les valeurs "groupees" : sol et date passees en parametre
	// et non _solCour et _date
	//============================================

	//init de sorte que beta-alpha = -1 si pas de branchement trouve
	sign._beta = -1;
	sign._alpha = 0;

	int K = static_cast<int> (date.size()); //solGroup[k] = debit entre les dates dateGroup[k], dateGroup[k+1] (solGroup[K-1] = toujours 0)
	
	

	
	double debitPrec = 0, bestGap = -1;
	int k = 0;

	while (k < K - 2)
	{
		//on avance jusqu'a la premiere rupture > +
		while (k < K && debitPrec >= sol[k] - EPSILON)
		{
			debitPrec = sol[k];
			k++;
		}
		//k => premiere rupture +

		//type + + - ?
		double gap = -1;

		if (k + 2 < K && sol[k] + EPSILON_DEBIT < sol[k + 1] && sol[k + 1] > sol[k + 2] + EPSILON_DEBIT)
		{
			double newAlpha = (date[k + 1] * sol[k + 1] - sol[k] * (date[k + 1] - date[k])) / sol[k + 1];
			gap = min(newAlpha - date[k], date[k + 1] - newAlpha);

			//meme si le gap est < espilon_date, on cree cette nouvelle date, elle sera fusionnee avec la plus proche dans le pretraitement du PL
			if ( gap > bestGap )
			{
				bestGap = gap;

				sign._type = 2;
				sign._alpha = newAlpha;
				sign._beta = date[k + 1];

//#ifdef _VERIF_
//				if (abs(sign._alpha - date[k]) < EPSILON_DATE || abs(sign._beta - sign._alpha) < EPSILON_DATE)
//					stopProg(" BranchBoundDicho::calculBranchType2ou3 : type 2, 2 dates trop proche");
//#endif
			}
		}

		//type + - - ?
		if (k + 2 < K && sol[k + 2] + EPSILON_DEBIT < sol[k + 1] && sol[k + 1] + EPSILON_DEBIT < sol[k])
		{
			k++;
			double newBeta = ((date[k + 1] - date[k])*sol[k] + date[k] * sol[k - 1]) / sol[k - 1];
			gap = min(date[k + 1] - newBeta, newBeta - date[k]);

			//meme si le gap est < espilon_date, on cree cette nouvelle date, elle sera fusionnee avec la plus proche dans le pretraitement du PL
			if (gap > bestGap)
			{
				bestGap = gap;

				sign._type = 3;
				sign._alpha = date[k];
				sign._beta = newBeta;

//#ifdef _VERIF_
//				if (abs(date[k + 1] - sign._beta) < EPSILON || abs(sign._beta - sign._alpha) < EPSILON)
//					stopProg(" BranchBoundDicho::calculBranchType2ou3 : type 3, 2 dates trop proche");
//#endif
			}
		}
		debitPrec = sol[k];
		k++;
	}



	sign._job = i;

#ifdef _VERIF_
	if (sign._beta > 0 && sign._beta < sign._alpha)
		stopProg("BranchBoundDicho::calculBranchType2ou3 : alpha et beta pas dans le bon sens");
#endif


	//si le gap est plus petit que EPSILON_DATE alors les dates seront fusionnées => le branchement sera inutile
	if (bestGap <= EPSILON_DATE)
		bestGap = -1;

	return bestGap;
}



//on verifie que si on branche le job i suivant les indices k1 et k2 alors les 3 branches resultantes coupent la solution courante
void BranchBoundDicho::verifieBranchementValide(int i, int k1, int k2)
{
	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (solCour[K-1] = toujours 0)
	const vector<double> & sol = _solCour[i];

	//on verifie que les trois branches vont couper
	// branche 1 : le job commence après k1 (on coupe si actuellement il commence avant)
	int k = k1 - 1;
	while (k >= 0 && sol[k] < EPSILON_DEBIT_VERIF)
		k--;
	if (k == -1)
		stopProg("BranchBoundDicho::verifieBranchementValide : la branche 1 ne coupe pas");

	// branche 2 : le job finit avant k2 (on coupe si actuellement il finit après)
	k = k2;
	while (k < K && sol[k] < EPSILON_DEBIT_VERIF)
		k++;
	if (k == K)
		stopProg("BranchBoundDicho::verifieBranchementValide : la branche 2 ne coupe pas");


	// branche 3 : le job est max debit entre k1 et k2 
	//(on coupe si actuellement le debit est plus eleve que sol[k1] avant k1 ou plus eleve que sol[k2-1] apres k2 ou non constant entre k1 et k2)

	double val = sol[k1];
	k = k1 + 1;

	//on avance k tant que valeur sol(k) = val, la solution est constante sur [k1,k2] si on sort avec k >= k2 (sol[k2-1] va jusqua date[k2])
	while (k < K && abs(sol[k] - val) < EPSILON_DEBIT_VERIF)
		k++;

	val = 0;
	if (k1 > 0)
		val = sol[k1 - 1];
	if (val <= sol[k1] && sol[k2] <= sol[k2 - 1] && k >= k2)
		stopProg("BranchBoundDicho::verifieBranchementValide : la branche 3 ne coupe pas");

}


//on regroupe les valeurs identiques consecutives dans _solCour[i] (evite de devoir traiter les cas = dans la recherche des param de branchement)
void BranchBoundDicho::regroupe(int i, vector<double> & solGroup, vector<double> & dateGroup)
{
	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (solCour[K-1] = toujours 0)
	const vector<double> & sol = _solCour[i];
	solGroup.clear();
	dateGroup.clear();

	double prec = -1;
	for (int k = 0; k < K; ++k)
	{
		if (abs(sol[k] - prec) > EPSILON_DEBIT)
		{
			solGroup.push_back(sol[k]);
			dateGroup.push_back(_dateCour[k]);
		}
		prec = sol[k];
	}

}


//affiche la solution courante
void BranchBoundDicho::afficheSolCour()
{
	int N = _ins->_nbJob;
	int K = static_cast<int>( _dateCour.size() );

	for (int i = 0; i < N; ++i)
	{
		cout << "job " << i << " : " << endl;
		for (int k = 0; k < K; ++k)
		{
			if (_solCour[i][k] > EPSILON)
			{
				cout << "[" << _dateCour[k] << ";" << _dateCour[k + 1] << "] : " << _solCour[i][k] << endl;
			}
		}
	}
}


//retourne le nombre de jobs non preemtif dans la solution courantes
int BranchBoundDicho::nbJobResolu()
{

	int nbJob = 0;

	int K = static_cast<int> (_dateCour.size()); //solCour[k] = debit entre les dates _date[k], _date[k+1] (la derniere date ets fictive pour avoir toujours solCour[K-1] = toujours 0)
	
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		const vector<double> & sol = _solCour[i];

		int k = 0;
		while (k < K && sol[k] < EPSILON_DEBIT)
			++k;

		while (k < K-1 && abs(sol[k] - sol[k + 1]) < EPSILON_DEBIT)
			++k;

		//le job finit sur un plateau => ok
		if (k == K)
			nbJob++;
		else
		{
			k++;
			//le job a un plateau puis que des 0 => ok
			while (k < K && sol[k] < EPSILON_DEBIT)
				++k;

			if (k == K)
				nbJob++;
		}

	}
	return nbJob;
}


//retourne la taille du plus petit intervalle dans dateCour
double BranchBoundDicho::intervalMin()
{
	
	int K = static_cast<int>(_dateCour.size() - 1);

	double res = TIME_INFINITY;

	for (int k = 0; k < K; ++k)
		if (_dateCour[k + 1] - _dateCour[k] < res)
			res = _dateCour[k + 1] - _dateCour[k];

	return res;

}


void BranchBoundDicho::ecrireHeurNoeud(const string & nomInst)
{
	ofstream fic("resBBNode_"+nomInst+".txt", ios::out);

	//attention on ajoute "-1" pour etre coherent avec la formule dans l'article

	fic << setprecision(3) << fixed;
	fic << setw(30) << nomInst
		<< setw(10) << _steril1 -1
		<< setw(10) << _steril2 -1
		<< setw(10) << _res_DebitFixe -1 
		<< setw(10) << _res_heurCut -1
		<< setw(10) << _res_tangente -1 
		<< setw(10) << _res_secante -1
		<< setw(10) << _res_cut -1 << endl;


	fic.close();
}


//on calcule le taux de parallélisme sur le dernier arc pour la best solution
//renvoie le taux de parallelisme moyen sur le dernier arc et le taux max
pair<double,int> BranchBoundDicho::statParallelisme()
{

	int N = _ins->_nbJob;


	//1. on trie les dates par ordre croissant

	vector<double> date;
	date.reserve(2 * N);

	for (int i = 0; i < N; ++i)
	{
		date.push_back(_bestDebut[i]);
		date.push_back(_bestDebut[i] + static_cast<double>(_ins->_pop[i]) / _bestDebit[i]);

		//cout << _bestDebut[i] << " ====" << _bestDebut[i] + static_cast<double>(_ins->_pop[i]) / _bestDebit[i] << endl;
	}


	//on trie on ordre croissant et on enleve les doublons
	sort(date.begin(), date.end());
	auto it = unique(date.begin(), date.end());
	date.resize(distance(date.begin(), it));


	for (unsigned int k = 0; k < date.size()-1; ++k)
	{
		if (abs(date[k] - date[k + 1]) < EPSILON_DATE)
		{
			date.erase(date.begin() + k + 1);
			k--;
		}
	}

	//2. on construit les cliques

	int K = static_cast<int> (date.size());
	

	vector<vector<int>> jobsDebutent(K); //jobDebutent[k] = ensemble des jobs qui debutent a la date k
	vector<vector<int>> jobsFinissent(K); //jobsFinissent[k] = ensemble des jobs qui finissent a la date k



	//-----------------------------------------------------------------------
	//on repere les jobs qui debutent / finissent en k


	//attention : a cause des problemes numeriques on peut avoir 
	// un meme job qui repond positivment au test abs(_bestDebut[i] - date[k]) < 0.001
	//pour 2 dates !=

	vector<bool> is_deb_fait(N, false);
	vector<bool> is_fin_fait(N, false);

	for (int i = 0; i < N; ++i)
	{
		for (int k = 0; k < K; ++k)
		{
			if (!is_deb_fait[i] && abs(_bestDebut[i] - date[k]) < 5 * EPSILON_DATE) //attention : mettre un peu plus que EPSILON_DATE sinon on risque de rater des jobs
			{
				jobsDebutent[k].push_back(i);
				is_deb_fait[i] = true;
			}

			if (!is_fin_fait[i] && abs(_bestDebut[i] + static_cast<double>(_ins->_pop[i]) / _bestDebit[i] - date[k]) < 5 * EPSILON_DATE)
			{
				jobsFinissent[k].push_back(i);
				is_fin_fait[i] = true;
			}
		}
	}

	//-------------------------------------------------------------------------
	// on calcule les cliques "overlap"

	vector<vector<int> > clique(K); //clique[k] = ensemble des jobs qui se croisent a la date [k, k+1[

	clique[0] = jobsDebutent[0];

	for (int k = 1; k < K; ++k)
	{
		clique[k] = clique[k - 1];

		//2.1 on enleve les jobs qui se finissent en k
		for (int j : jobsFinissent[k])
		{
			auto it = find(clique[k].begin(), clique[k].end(), j);
			clique[k].erase(it);
		}

		//2.2 on ajoute ceux qui commencent en k
		for (int j : jobsDebutent[k])
			clique[k].push_back(j);

	}

	//=========================================
	//on fait les stats

	double paralMoyen = 0;
	int paralMax = 0;


	for (int k = 0; k < K-1; ++k)//RQ : clique[K-1].size = 0, aucun job ne peut commencer a la fin
	{
		paralMax = max( paralMax, static_cast<int>(clique[k].size()) );
		paralMoyen += clique[k].size() * (date[k + 1] - date[k]);

	}

	paralMoyen /= (date[K - 1] - date[0]);

	return { paralMoyen, paralMax };

}