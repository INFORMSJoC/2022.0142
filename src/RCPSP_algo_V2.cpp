#include "RCPSP_algo_V2.h"
#include "PL_Optim_RCPSP.h"
#include "BSmarge.h"
#include "util.h"
#include "MIP_flot_partiel.h"

#include <queue>
#include <algorithm>

#include <random>
#include <numeric>
#include <tuple>



//=====================================================================
// algo principal



RCPSP_Algo_V2::Resultat RCPSP_Algo_V2::GRASP(int maxRun, double precision)
{
	Resultat bestRes;
	int nbChangement = 0;
	int cpt = 0; //compte le nb de tours de boucle principal

	//  ============= INITIALISATION  =====================

	vector<vector<int>> coupleInterdit(_ins->_nbJob); //coupleInterdit[i] = liste des j tel que i ne doit pas etre avant j dans sigma

	// Tmax est la plus grande des deadlines
	double TMax = 0;
	for (double lf : _ins->_LF)
		TMax = max(TMax, lf);

	//marges minimales et maximales (en supposant que l'instance est realisable)
	double lambda_min = 0, lambda_max = TIME_INFINITY;
	for (int i = 0; i < _ins->_nbJob; ++i)
		lambda_max = min(lambda_max, _ins->_LF[i] - (_ins->_ES[i] + (double)(_ins->_pop[i]) / _ins->_debitMax[i]));

	//liste ordonnee des jobs suivant leur priorite
	vector<int> sigma;
	construireSigmaProba(sigma, coupleInterdit, 1);

	//  ============= BOUCLE PRINCIPALE : construction d'une solution realisable =====================

	bool first = true;
	double proba_sigma = 0.5;

	vector < pair<int, int> > precedence;
	vector < vector<int> > cliques;

	double lambda = 0;

	while (lambda_max - lambda_min >= precision)
	{
		lambda = (lambda_max + lambda_min) / 2;
		cpt = 0;
		bool ok = false;

		if (first)
			lambda = 0;
		


		while (cpt < maxRun && !ok)
		{
			auto [out_ok, res, interdit] = insert(lambda, sigma);
			ok = out_ok; 

			if (ok)
			{
				lambda_min = lambda;

				//======= on essaie d'ameliorer avec le PL "topologie"

				precedence.clear();
				cliques.clear();

				//1. calcul des precedences et cliques 
				calculPredClique(precedence, cliques);

				//2. lance le PL et maj  _debit, _dateDebut
				double resPL = reconstruireSolCourantePL_exact(precedence, cliques);

				if (!_ins->isSolution(_debutCour, _debitCour, resPL))
					stopProg("RCPSP_Algo_V2::GRASP : solution fausse a la fin de la phase d 'optimisation");

				if (bestRes._margeFinale < resPL)
				{
					bestRes._dateDebut = _debutCour;
					bestRes._debit = _debitCour;
					bestRes._margeFinale = resPL;
					bestRes._bestLambda = res._bestLambda;
					if (lambda < EPSILON)
						bestRes._marge0 = resPL;
				}

				//==== fin amelioration
			}
			else
			{
#ifdef _VERIF_
				if (interdit.first == _source)
					stopProg("RCPSP_Algo_V2::GRASP : on ne peut pas interdire la source avant");
#endif
				coupleInterdit[interdit.first].push_back(interdit.second);
				cpt++;
				bool sigmaOK = construireSigmaProba(sigma, coupleInterdit, proba_sigma);


				while (!sigmaOK)
				{
					//reset(coupleInterdit);
					menageAlea(coupleInterdit, 0.5);
					sigmaOK = construireSigmaProba(sigma, coupleInterdit, proba_sigma);

				}
			}
		}
		if (!ok && !first) //meme si lambda = 0 a echoue on va essayer avec une marge plus grande (car cela change la maniere de distribuer les flots)
			lambda_max = lambda;
		
		first = false;
		
	}//fin while

	
	cout << "marge GRASP = " << bestRes._margeFinale << endl;
	return bestRes;
}

//seconde version de l'heuristique : on ameliore a chaque iteration avec le PL derive du MIP flot dans lequel on fixe les x_ij
RCPSP_Algo_V2::Resultat RCPSP_Algo_V2::GRASP_V2(int maxRun, double precision)
{
	Resultat bestRes;


	int cpt = 0; //compte le nb de tours de boucle principal

	//  ============= INITIALISATION  =====================

	vector<vector<int>> coupleInterdit(_ins->_nbJob); //coupleInterdit[i] = liste des j tel que i ne doit pas etre avant j dans sigma

	//marges minimales et maximales (en supposant que l'instance est realisable)
	double lambda_min = 0, lambda_max = TIME_INFINITY;
	for (int i = 0; i < _ins->_nbJob; ++i)
		lambda_max = min(lambda_max, _ins->_LF[i] - (_ins->_ES[i] + (double)(_ins->_pop[i]) / _ins->_debitMax[i]));

	//liste ordonnee des jobs suivant leur priorite
	vector<int> sigma;
	construireSigmaProba(sigma, coupleInterdit, 1);

	//  ============= BOUCLE PRINCIPALE : construction d'une solution realisable =====================

	bool first = true;
	double proba_sigma = 0.5;

	double lambda = 0;

	while (lambda_max - lambda_min >= precision)
	{
		lambda = (lambda_max + lambda_min) / 2;
		cpt = 0;
		bool ok = false;

		if (first)
			lambda = 0;

		while (cpt < maxRun && !ok)
		{
			auto[out_ok, res, interdit] = construitSol(lambda, sigma);
			ok = out_ok;

			if (ok)
			{
				lambda_min = lambda;

				if (res._margeFinale > bestRes._margeFinale)
				{
					bestRes._dateDebut = res._dateDebut;
					bestRes._debit = res._debit;
					bestRes._margeFinale = res._margeFinale;
					bestRes._bestLambda = lambda;
					if (lambda < EPSILON)
						bestRes._marge0 = res._margeFinale;
				}

			}
			else
			{
#ifdef _VERIF_
				if (interdit.first == _source)
					stopProg("RCPSP_Algo_V2::GRASP : on ne peut pas interdire la source avant");
#endif
				coupleInterdit[interdit.first].push_back(interdit.second);
				cpt++;
				bool sigmaOK = construireSigmaProba(sigma, coupleInterdit, proba_sigma);


				while (!sigmaOK)
				{
					//reset(coupleInterdit);
					menageAlea(coupleInterdit, 0.5);
					sigmaOK = construireSigmaProba(sigma, coupleInterdit, proba_sigma);

				}
			}
		}
		if (!ok && !first) //meme si lambda = 0 a echoue on va essayer avec une marge plus grande (car cela change la maniere de distribuer les flots)
			lambda_max = lambda;

		first = false;

	}//fin while

	if (bestRes._margeFinale >= 0)
	{
#ifdef _VERIF_
		if (!_ins->isSolution(bestRes._dateDebut, bestRes._debit, bestRes._margeFinale))
			stopProg("GRASP_V2 : solution finale fausse");
#endif
	}
	else
		bestRes._margeFinale = -1;

	//cout << "GRASP_V2, best marge = " << bestRes._margeFinale << endl;
	return bestRes;

}


//troisieme version de l'heuristique : on ameliore a chaque iteration "a la main" puis a la fin avec le PL derive du MIP flot dans lequel on fixe les x_ij
RCPSP_Algo_V2::Resultat RCPSP_Algo_V2::GRASP_V3(int maxRun, double precision, bool PL_iter)
{
	Resultat bestRes;
	srand(777);//pour pouvoir rejouer ...

	int cpt = 0; //compte le nb de tours de boucle principal

	//  ============= INITIALISATION  =====================

	vector<vector<int>> coupleInterdit(_ins->_nbJob); //coupleInterdit[i] = liste des j tel que i ne doit pas etre avant j dans sigma

	//marges minimales et maximales (en supposant que l'instance est realisable)
	double lambda_min = 0, lambda_max = TIME_INFINITY;
	for (int i = 0; i < _ins->_nbJob; ++i)
		lambda_max = min(lambda_max, _ins->_LF[i] - (_ins->_ES[i] + (double)(_ins->_pop[i]) / _ins->_debitMax[i]));

	//liste ordonnee des jobs suivant leur priorite
	vector<int> sigma;
	construireSigmaProba(sigma, coupleInterdit, 1);

	//  ============= BOUCLE PRINCIPALE : construction d'une solution realisable =====================

	bool first = true;
	double proba_sigma = 0.5;

	double lambda = 0;

	while (lambda_max - lambda_min >= precision)
	{
		lambda = (lambda_max + lambda_min) / 2;
		cpt = 0;
		bool ok = false;

		if (first)
			lambda = 0;

		while (cpt < maxRun && !ok)
		{
			auto[out_ok, res, interdit] = construitSol_graspV3(lambda, sigma, PL_iter);
			ok = out_ok;

			if (ok)
			{
				lambda_min = lambda;

				if (res._margeFinale > bestRes._margeFinale)
				{
					bestRes._dateDebut = res._dateDebut;
					bestRes._debit = res._debit;
					bestRes._margeFinale = res._margeFinale;
					bestRes._bestLambda = lambda;
					if (lambda < EPSILON)
						bestRes._marge0 = res._margeFinale;
				}

			}
			else
			{
#ifdef _VERIF_
				if (interdit.first == _source)
					stopProg("RCPSP_Algo_V2::GRASP : on ne peut pas interdire la source avant");
#endif
				coupleInterdit[interdit.first].push_back(interdit.second);
				cpt++;
				bool sigmaOK = construireSigmaProba(sigma, coupleInterdit, proba_sigma);


				while (!sigmaOK)
				{
					reset(coupleInterdit);
					//menageAlea(coupleInterdit, 0.5);
					sigmaOK = construireSigmaProba(sigma, coupleInterdit, proba_sigma);

				}
			}
		}
		if (!ok && !first) //meme si lambda = 0 a echoue on va essayer avec une marge plus grande (car cela change la maniere de distribuer les flots)
			lambda_max = lambda;

		first = false;

	}//fin while

	if (bestRes._margeFinale >= 0)
	{
#ifdef _VERIF_
		if (!_ins->isSolution(bestRes._dateDebut, bestRes._debit, bestRes._margeFinale))
			stopProg("GRASP_V3 : solution finales fausse");
#endif
	}
	else
		bestRes._margeFinale = -1;

	cout << "GRASP_V3, best marge = " << bestRes._margeFinale << endl;
	return bestRes;

}


bool RCPSP_Algo_V2::construireSigmaAlea(vector<int> & sigma, const vector<vector<int>> & coupleInterdit)
{
	sigma.clear();
	int nbCandidatMax = 5;

	vector<bool> inserted(_ins->_nbJob, false);
	vector<pair<double, int>> candidat; //liste des jobs avec leur ES qu'on peut mettre a la fin de sigma sans viole les coupleInterdit

	while (sigma.size() < _ins->_nbJob)
	{
		// 1. on cherche les candidats
		candidat.clear();
		for (int i = 0; i < _ins->_nbJob; ++i)
		{
			if (!inserted[i])
			{
				int j = 0;
				while (j < coupleInterdit[i].size() && inserted[coupleInterdit[i][j]])
					j++;
				if (j == coupleInterdit[i].size())//pas de coupel interdit trouve, i est candidat
					candidat.push_back({ _ins->_ES[i], i });
			}

		}

		//si on ne trouve pas de candidat on renvoie un signal d echec
		if (candidat.size() == 0)
			return false;


		// 2. on tri sur les ES croissant si on a plus de nbCandidatMax candidats
		if (candidat.size() > nbCandidatMax)
			sort(candidat.begin(), candidat.end());

		int t = min(nbCandidatMax, (int)(candidat.size()));

		// 3. on choisit les 2 avec la plus petite valeur V = LF - ES - 2 * pop/(debMax-debMin)
		for (int k = 0; k < t; ++k)
		{
			int job = candidat[k].second;
			candidat[k].first = _ins->_LF[job] - _ins->_ES[job] - 2 * _ins->_pop[job] / (_ins->_debitMax[job] - _ins->_debitMin[job]);
		}
		//on tri les t premiers candidats suivant la valeur V
		sort(candidat.begin(), candidat.begin() + t);

		//4. on choisit alea parmi les 2 premiers 
		t = rand() % min((int)(candidat.size()), 2);



		//5. on insere le job choisi dans sigma
		sigma.push_back(candidat[t].second);
		inserted[candidat[t].second] = true;
	}
#ifdef _VERIF_
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		if (!inserted[i])
			stopProg("RCPSP_Algo_V2::construireSigma : des jobs non inseres");
	}
#endif

	return true;
}

//alpha \in [0,1] = probabilite de choisir le 1er candidat 
bool RCPSP_Algo_V2::construireSigmaProba(vector<int> & sigma, const vector<vector<int>> & coupleInterdit, double alpha)
{
	sigma.clear();
	
	vector<bool> inserted(_ins->_nbJob, false);
	vector<pair<double, int>> candidat; //liste des jobs avec leur ES qu'on peut mettre a la fin de sigma sans viole les coupleInterdit

	while (sigma.size() < _ins->_nbJob)
	{
		// 1. on cherche les candidats
		candidat.clear();
		for (int i = 0; i < _ins->_nbJob; ++i)
		{
			if (!inserted[i])
			{
				int j = 0;
				while (j < coupleInterdit[i].size() && inserted[coupleInterdit[i][j]])
					j++;
				if (j == coupleInterdit[i].size())//pas de coupel interdit trouve, i est candidat
				{
					//stress : 
					//candidat.push_back({ _ins->_LF[i] - _ins->_ES[i] - 2 * _ins->_pop[i] / (_ins->_debitMax[i] - _ins->_debitMin[i]), i });
					//fin : 
					//candidat.push_back({ _ins->_LF[i], i });
					//moyenne deb-fin : 
					candidat.push_back({ (_ins->_ES[i] + 2*_ins->_LF[i])/3, i });
					//moyenne debit min / max => MAUVAIS
					//candidat.push_back({ (_ins->_debitMin[i]+ _ins->_debitMax[i])/2, i });
					//pop => TRES MAUVAIS
					//candidat.push_back({ _ins->_pop[i], i });
				}
			}

		}

		//si on ne trouve pas de candidat on renvoie un signal d echec
		if (candidat.size() == 0)
			return false;


		//on tri les candidats suivant la valeur 
		std::sort(candidat.begin(), candidat.end());
		double alea = (double)(rand()) / RAND_MAX;
		int t = 1;
		if (alea < alpha || candidat.size() == 1)//chosit le 1er avec proba alpha, sinon choisit le 2eme
			t = 0;



		//5. on insere le job choisi dans sigma
		sigma.push_back(candidat[t].second);
		inserted[candidat[t].second] = true;
	}
#ifdef _VERIF_
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		if (!inserted[i])
			stopProg("RCPSP_Algo_V2::construireSigma : des jobs non inseres");
	}
#endif

	return true;
}



tuple<bool, RCPSP_Algo_V2::Resultat, pair<int, int>> RCPSP_Algo_V2::insert(double lambda, const vector<int> & sigma)
{
	bool use_trou = false;//mettre a vrai pour utiliser les "trous" i.e. inserer en priorite dans es trous entre 2 jobs
	bool ok = true;
	RCPSP_Algo_V2::Resultat res;
	pair<int, int> interdit = {-1,-1};//couple (i,j) tel que i << j interdit dans sigma en cas d echec

	// =========== INITIALISATION =================

	initInsert();
	vector<pair<double,int>> listePrec;      // listePrec = liste des jobs deja inseres et trie suivant leur fin croissante
	vector<pair<double, int>> listePrec_deb; // listePrec_deb = liste des jobs deja inseres et trie suivant leur debut decroissant


	// =========== BOUCLE PRINCIPALE =================

	for (int i = 0; i < sigma.size(); ++i)
	{
		int jobCour = sigma[i];
		double popJCour = (double)(_ins->_pop[jobCour]);
		double LFJCour = _ins->_LF[jobCour];

		//pretraitement : liste des jobs deja inseres suivant leur date de fin croissante
		listePrec.clear(); 
		listePrec_deb.clear();
		listePrec.push_back({0, _source}); //init avec la source
		listePrec_deb.push_back({ 0, _source });

		for (int j = 0; j < i; ++j)
		{
			listePrec.push_back({ _finCour[sigma[j]], sigma[j] });
			listePrec_deb.push_back({ -_debutCour[sigma[j]], sigma[j] }); //on met un "-" pour trier en decroissant
		}

		sort(listePrec.begin(), listePrec.end());
		sort(listePrec_deb.begin(), listePrec_deb.end());


		vector<double> delta(_ins->_nbArc);
		vector<double> flotCour(_ins->_nbArc);


		//   ===== etape 1 : on donne du flot a jobCour pour qu'il finisse avec une marge lambda  ============
		bool ok1; pair<int, int> interdit1; vector<double> detourne;
		if (!use_trou)
		{
			auto res1 = insert_etape1_fin(jobCour, popJCour, LFJCour, lambda, listePrec, flotCour, delta);
			ok1 = get<0>(res1); interdit1 = get<1>(res1);
		}
		else
		{
			auto res1 = insert_etape1_milieu(jobCour, popJCour, LFJCour, lambda, listePrec, listePrec_deb, flotCour, delta);
			ok1 = get<0>(res1); interdit1 = get<1>(res1); detourne = get<2>(res1);
		}
		
		if (!ok1)
			return { false, res, interdit1 };
		
		// ============= mise a jour des info pour jobCour ============

		_debitCour[jobCour] = 0;
		int e_max = -1; //arc qui donne la valeur du debit de jobCour
		for (int e : _ins->_jobToArc[jobCour])
		{
			if (flotCour[e] > _debitCour[jobCour])
			{
				_debitCour[jobCour] = flotCour[e];
				e_max = e;
			}
		}

		_debutCour[jobCour] = delta[e_max];
		_finCour[jobCour] = _debutCour[jobCour] + popJCour / _debitCour[jobCour];
		for (int e : _ins->_jobToArc[jobCour])
		{
			_ressATransmettre[jobCour][e] = _debitCour[jobCour];
			if (use_trou)
				_ressATransmettre[jobCour][e] -= detourne[e];
		}

		//   ===== etape 2 : jobCour doit recevoir _debitCour[jobCour] pour chaque arc (ce n'est pas forcement au terme de l'etape 1 : on essaie d'ajuster 
		bool ok2;
		pair<int, int> interdit2;

		if (!use_trou)
		{
			auto res2 = insert_etape2_fin(e_max, jobCour, listePrec, flotCour);
			ok2 = get<0>(res2); interdit2 = get<1>(res2);
		}
		else
		{
			auto res2 = insert_etape2_milieu(e_max, jobCour, popJCour, listePrec, listePrec_deb, flotCour);
			ok2 = get<0>(res2); interdit2 = get<1>(res2);
		}


		if (!ok2)
			return { false, res, interdit2 };

	}//fin for sigma

#ifdef _VERIF_
	if (!_ins->isSolution(_debutCour, _debitCour, lambda))
		stopProg("RCPSP_Algo_V2::insert : sol. non realisable");
	if (use_trou)
	{
		if (!isFlotOK()) //ATTENTION prend beaucoup de temps
			stopProg("insert: flot faux");
	}
#endif


	//si on na pas retourner un echec alors on a trouve une solution 
	res._dateDebut = _debutCour;
	res._debit = _debitCour;
	res._bestLambda = lambda;


	return { true, res, interdit };
}


// procedure principale du GRASP_V2 : construit une solution heuristique avec une marge au moins lambda
tuple<bool, RCPSP_Algo_V2::Resultat, pair<int, int>> RCPSP_Algo_V2::construitSol(double lambda, const vector<int> & sigma)
{
	bool ok = true;
	RCPSP_Algo_V2::Resultat res;
	pair<int, int> interdit = { -1,-1 };//couple (i,j) tel que i << j interdit dans sigma en cas d echec

	double obj_PL = -1;
	double bestSol = -1;
	vector<pair<int, int>> bestArcSupport;

	// =========== INITIALISATION =================

	initInsert();
	initConstruitSol();

	_inserted.assign(_ins->_nbJob, false);
	_arcSupport.clear();

	vector<pair<double, int>> listePrec;      // listePrec = liste des jobs deja inseres et trie suivant leur fin croissante


	// =========== BOUCLE PRINCIPALE =================

	for (int i = 0; i < sigma.size(); ++i)
	{
		int jobCour = sigma[i];
		double popJCour = (double)(_ins->_pop[jobCour]);
		double LFJCour = _ins->_LF[jobCour];

		//pretraitement : liste des jobs deja inseres suivant leur date de fin croissante
		listePrec.clear();
		listePrec.push_back({ 0, _source }); //init avec la source
		_dejaSupport.assign(_ins->_nbJob, false);


		for (int j = 0; j < i; ++j)
		{
			listePrec.push_back({ _finCour[sigma[j]], sigma[j] });
		}

		sort(listePrec.begin(), listePrec.end());
		
		vector<double> delta(_ins->_nbArc);
		vector<double> flotCour(_ins->_nbArc);


		//   ===== etape 1 : on donne du flot a jobCour pour qu'il finisse avec une marge lambda  ============

//#ifdef _VERIF_
//		if (!isFlotOK())
//			stopProg("RCPSP_Algo_V2::construitSol : flot faux avant etape 1");
//#endif
//

		auto[ok1, interdit1] = insert_etape1_fin(jobCour, popJCour, LFJCour, lambda, listePrec, flotCour, delta);
			
		if (!ok1)
			return { false, res, interdit1 };


		// ============= mise a jour des info pour jobCour ============

		_debitCour[jobCour] = 0;
		int e_max = -1; //arc qui donne la valeur du debit de jobCour
		for (int e : _ins->_jobToArc[jobCour])
		{
			if (flotCour[e] > _debitCour[jobCour])
			{
				_debitCour[jobCour] = flotCour[e];
				e_max = e;
			}
		}



		_debutCour[jobCour] = delta[e_max];
		_finCour[jobCour] = _debutCour[jobCour] + popJCour / _debitCour[jobCour];
		for (int e : _ins->_jobToArc[jobCour])
			_ressATransmettre[jobCour][e] = _debitCour[jobCour];

		//   ===== etape 2 : jobCour doit recevoir _debitCour[jobCour] pour chaque arc (ce n'est pas forcement au terme de l'etape 1 : on essaie d'ajuster 
		

		auto[ok2, interdit2] = insert_etape2_fin(e_max, jobCour, listePrec, flotCour);

		if (!ok2)
			return { false, res, interdit2 };

		//=======================================================================================
		// ========= etape 3 : on ameliore avec le MIP_flot_partiel l'insertion de jobCour

		_inserted[jobCour] = true;
//#ifdef _VERIF_
//		if (!isFlotOK())
//			stopProg("RCPSP_Algo_V2::construitSol : flot faux apres etape 2");
//#endif


		//3.1. on ameliore la solution courante avec le MIP flot : meilleure sol possible avec les arcs supports courant
		obj_PL = PL_flot_partiel();
		conserveSolPL();

#ifdef _VERIF_
		if (obj_PL < lambda)
			stopProg("RCPSP_Algo_V2::construitSol : sol moins bien ap MIP_flot_partiel");
#endif
		
		bestArcSupport = _arcSupport;
		bestSol = obj_PL;


		vector<pair<double, int>> listeSucc;
		for (int j = 0; j < i; ++j)
			listeSucc.push_back({ -_debutCour[sigma[j]], sigma[j] });
		
		//on trie suivnt les debuts decroissants
		sort(listeSucc.begin(), listeSucc.end());


		//--------------------------------------------------------------
		//3.2. on essaie d'ajouter /supprimer des arcs support sans creer de cycle

		for (auto s = 0; s < listeSucc.size(); ++s)
		{
			//choisir les jobs j dans l'ordre de listeSucc, ajouter jobCour->jsucc a arcsupport
			int succ = listeSucc[s].second;

			//pour etre sur de ne pas creer de cycle, on ajoute jobcour -> succ dans les arcs supports que si jobCour ne finit pas avant le debut de succ
			if (_debutCour[jobCour] < _finCour[succ])
			{
				// +++ on ajoute l'arc jobCour->succ tel que succ commence le plus tard possible
				_arcSupport.push_back({ jobCour,succ });

				//on regarde si on ameliore avec le PL MIP et ce nouvel arc support
				obj_PL = PL_flot_partiel();
				if (obj_PL > bestSol)
				{
					conserveSolPL();
					bestSol = obj_PL;
					bestArcSupport = _arcSupport;
				}
			}

			bool finMajSupport = false;
			while (!finMajSupport)
			{
				finMajSupport = true;

				// --- on supprime l'arc prec -> jobCour des arc support telque prec fini le plus tard possible (s'il en existe un avec prec != source)
				if (supprimeSupportfinMaxOrigine(jobCour))
				{
					finMajSupport = false;
					obj_PL = PL_flot_partiel();
					if (obj_PL > bestSol)
					{
						conserveSolPL();
						bestSol = obj_PL;
						bestArcSupport = _arcSupport;
					}
				}
			}
			

		}
		//on remet les arcs support qui correspondent a la meilleure solution (celle de laquelle on repart pour continuer les insertions)
		_arcSupport = bestArcSupport;
	}//fin for sigma


#ifdef _VERIF_
	if (!_ins->isSolution(_debutCour, _debitCour, bestSol))
		stopProg("RCPSP_Algo_V2::construitSol : sol. non realisable");
#endif


	//si on na pas retourner un echec alors on a trouve une solution 
	res._dateDebut = _debutCour;
	res._debit = _debitCour;
	res._bestLambda = lambda;
	res._margeFinale = bestSol;


	return { true, res, interdit };
}

// procedure principale du GRASP_V3 : construit une solution heuristique avec une marge au moins lambda
// correspond a A_Insert-MSM-RCPSP du doc Alain
tuple<bool, RCPSP_Algo_V2::Resultat, pair<int, int>> RCPSP_Algo_V2::construitSol_graspV3(double lambda, const vector<int> & sigma, bool PL_iter)
{

	RCPSP_Algo_V2::Resultat res;
	pair<int, int> interdit = { -1,-1 };//couple (i,j) tel que i << j interdit dans sigma en cas d echec
	double obj_PL = -1;


	// =========== INITIALISATION =================

	initInsert();
	initConstruitSol();

	_inserted.assign(_ins->_nbJob, false);
	vector<pair<double, int>> listePrec;       // listePrec = liste des jobs deja inseres et trie suivant leur fin croissante [J_Curr]
	vector<pair<double, int>> listePrec_deb;   // listePrec_deb = liste des jobs deja inseres et trie suivant leur debut decroissante [J*]

	_arcSupport.clear();


	vector<vector<double> > ressAtransmettreINIT;
	vector<vector<vector<double> > > flotCourINIT;
	vector<vector<double> > ressAtransmettreSAV;
	vector<vector<vector<double> > > flotCourSAV;

	// =========== BOUCLE PRINCIPALE =================


	for (int i = 0; i < sigma.size(); ++i)
	{
		int jobCour = sigma[i];
		double popJCour = (double)(_ins->_pop[jobCour]);
		double LFJCour = _ins->_LF[jobCour];

		//pretraitement : liste des jobs deja inseres suivant leur date de fin croissante
		listePrec.clear();
		listePrec.push_back({ 0, _source }); //init avec la source
	
		for (int j = 0; j < i; ++j)
		{
			listePrec.push_back({ _finCour[sigma[j]], sigma[j] });
		}

		sort(listePrec.begin(), listePrec.end());

		vector<double> delta(_ins->_nbArc);
		vector<double> flotCour(_ins->_nbArc);

		//on sauvegarde l'etat des flots avant d inserer jobCour (car on repart de cet etat dans la boucle d'amelioration)
		ressAtransmettreINIT = _ressATransmettre;
		flotCourINIT = _flotCour;

		//   ===== etape 1 : on donne du flot a jobCour pour qu'il finisse avec une marge lambda  ============

		auto[ok1, interdit1] = insert_etape1_fin(jobCour, popJCour, LFJCour, lambda, listePrec, flotCour, delta);

		if (!ok1)
		{
			interdit = interdit1;
			//return { false, res, interdit1 };
		}
		else
		{
			// ============= mise a jour des info pour jobCour ============

			_debitCour[jobCour] = 0;
			int e_max = -1; //arc qui donne la valeur du debit de jobCour
			for (int e : _ins->_jobToArc[jobCour])
			{
				if (flotCour[e] > _debitCour[jobCour])
				{
					_debitCour[jobCour] = flotCour[e];
					e_max = e;
				}
			}

			_debutCour[jobCour] = delta[e_max];
			_finCour[jobCour] = _debutCour[jobCour] + popJCour / _debitCour[jobCour];
			for (int e : _ins->_jobToArc[jobCour])
				_ressATransmettre[jobCour][e] = _debitCour[jobCour];

			//   ===== etape 2 : jobCour doit recevoir _debitCour[jobCour] pour chaque arc (ce n'est pas forcement au terme de l'etape 1 : on essaie d'ajuster 

			auto[ok2, interdit2] = insert_etape2_fin(e_max, jobCour, listePrec, flotCour);

			if (!ok2)
			{
				interdit = interdit2;
				//return { false, res, interdit2 };
			}
			else 
			{
				_inserted[jobCour] = true;
#ifdef _VERIF_
				if (!isFlotOK())
					stopProg("RCPSP_Algo_V2::construitSol_grapsV3 : flot faux avant bricolage");
#endif
			}
		}

		//=======================================================================================
		//  etape 3 : on essaie de mettre jobCour plus tot en fixant sa deadline au
		// debut du job qui demarre le plus tard, puis avant etc ...

		//RQ meme si on a echoue au tour precedent, on tente car ca peut permettre de debloquer la situation

		//3.1 on tri les jobs inseres (sauf jobCour) par debut decroissant
		double debitJCourSAV;
		double debutJCourSAV;
		double finJCourSAV;
		listePrec_deb.clear();

		for (int j = 0; j < i; ++j)
			listePrec_deb.push_back({ -_debutCour[sigma[j]], sigma[j] });

		sort(listePrec_deb.begin(), listePrec_deb.end());

		for (int k = 0; k < listePrec_deb.size(); ++k)
		{
			int jsucc = listePrec_deb[k].second; //[j rond]

			//on doit avoir la place de mettre jobCour entre sa deadline et le debut de jsucc, sinon on peut arreter la boucle (echec)
			if (_debutCour[jsucc] - _ins->_ES[jobCour] < popJCour / _ins->_debitMax[jobCour])//si longueur possible < longueur minimale
				break;

			// ---- on sauvegarde les flots en cas d'echec (si on a deja trouve une sol)--------------
			if (_inserted[jobCour])
			{
				ressAtransmettreSAV = _ressATransmettre;
				flotCourSAV = _flotCour;
				debitJCourSAV = _debitCour[jobCour];
				debutJCourSAV = _debutCour[jobCour];
				finJCourSAV = _finCour[jobCour];
			}


			//et on repart comme si jobCour n'avait pas ete insere
			_ressATransmettre = ressAtransmettreINIT;
			_flotCour = flotCourINIT; //attention, on ne peut pas juste remmettre les flots qui vont vers / depuis jobCour car _flotCour change pour d'autres jobs dans l'etape d 'amelioration
			
			
			//3.2 on fait comme si jsucc n etait plus dans le planing => il rend les ressources
			for (int e : _ins->_jobToArc[jsucc])
			{
				_ressATransmettre[_source][e] += _flotCour[_source][jsucc][e];
				for (int j = 0; j < i; ++j)
					_ressATransmettre[sigma[j]][e] += _flotCour[sigma[j]][jsucc][e];
			}

			//3.3 on veut mettre jobCour avnt jsucc => on adapte sa deadline
			LFJCour = min(-listePrec_deb[k].first, _ins->_LF[jobCour]);
			double newLambda = max(0.0, lambda - (_ins->_LF[jobCour] - LFJCour));

			// 3.4 on recommence l'insertino de jobCour
			//   ===== etape 1 : on donne du flot a jobCour pour qu'il finisse avec une marge lambda  ============

			auto[ok1, interdit1] = insert_etape1_fin(jobCour, popJCour, LFJCour, newLambda, listePrec, flotCour, delta);

			if (!ok1)
			{
				if (_inserted[jobCour])
				{
					_ressATransmettre = ressAtransmettreSAV;
					_flotCour = flotCourSAV;
				}
			}
			else
			{

				// ============= mise a jour des info pour jobCour ============

				_debitCour[jobCour] = 0;
				int e_max = -1; //arc qui donne la valeur du debit de jobCour
				for (int e : _ins->_jobToArc[jobCour])
				{
					if (flotCour[e] > _debitCour[jobCour])
					{
						_debitCour[jobCour] = flotCour[e];
						e_max = e;
					}
				}

				_debutCour[jobCour] = delta[e_max];
				_finCour[jobCour] = _debutCour[jobCour] + popJCour / _debitCour[jobCour];
				for (int e : _ins->_jobToArc[jobCour])
					_ressATransmettre[jobCour][e] = _debitCour[jobCour];

				//   ===== etape 2 : jobCour doit recevoir _debitCour[jobCour] pour chaque arc (ce n'est pas forcement au terme de l'etape 1 : on essaie d'ajuster 

				auto[ok2, interdit2] = insert_etape2_fin(e_max, jobCour, listePrec, flotCour);
				if (!ok2)
				{

					if (_inserted[jobCour])
					{
						_ressATransmettre = ressAtransmettreSAV;
						_flotCour = flotCourSAV;
						_debitCour[jobCour] = debitJCourSAV;
						_debutCour[jobCour] = debutJCourSAV;
						_finCour[jobCour] = finJCourSAV;
					}
				}
				else
				{
					//cout << "nouvelle sol !!" << endl;
					//on a une nouvelle solution 
					//il faut remettre correctement le flot qu'on a enleve vers jsucc
					_inserted[jobCour] = true;

					for (int e : _ins->_jobToArc[jsucc])
					{
						for (int j = 0; j < i; ++j)
						{
							int job = sigma[j]; //[j1] = jobs deja inseres
							if (job == jsucc)
								job = _source; //inutile de le faire pour jsucc, par contre il ne faut pas oublier la source
							
							double f_jsucc = _flotCour[job][jsucc][e];
						
							if (f_jsucc > 0)//si job ne transmettait rien a jsucc alors il ny a rien a faire
							{
								double f_jCour = _flotCour[job][jobCour][e];
								if (f_jCour >= f_jsucc)//si job trabsmet mainteant a jobCour au moins f => on les envoie vers jsucc
								{
									_flotCour[job][jsucc][e] -= f_jsucc;
									_flotCour[jobCour][jsucc][e] += f_jsucc;
									_ressATransmettre[jobCour][e] -= f_jsucc;
								}
								else//si job transmet f
								{
									_flotCour[job][jsucc][e] -= f_jCour;
									_flotCour[jobCour][jsucc][e] += f_jCour;
									_ressATransmettre[jobCour][e] -= f_jCour;
									_ressATransmettre[job][e] -= _flotCour[job][jsucc][e];
								}
							}
						}
					}//fin for e
#ifdef _VERIF_
					if (!isFlotOK())
						stopProg("RCPSP_Algo_V2::construitSol_grapsV3 : flot faux apres bricolage");
#endif
				}
			}



		}

		if (!_inserted[jobCour])
			return{ false, res, interdit };
		else
		{
			if (PL_iter)
			{
				recupereArcSupport();

				obj_PL = PL_flot_partiel();
				conserveSolPL();

#ifdef _VERIF_
				if (obj_PL < lambda)
					stopProg("RCPSP_Algo_V2::construitSol_graspV3 : sol moins bien ap MIP_flot_partiel");

#endif
			}
		}

	}//fin for sigma


	// ===================================================================================
	// AMELIORATION AVEC LE MIP (si pas fait a chaque it)
	// ===================================================================================
	
	//remplit _arcSupport en utilisant le flot courant
	if (!PL_iter)//on lance une seule fois le PL ameliorant a la fin
	{
		recupereArcSupport();

		obj_PL = PL_flot_partiel();
		conserveSolPL();

#ifdef _VERIF_
		if (obj_PL < lambda)
			stopProg("RCPSP_Algo_V2::construitSol_graspV3 : sol moins bien ap MIP_flot_partiel");
#endif
	}

	res._dateDebut = _debutCour;
	res._debit = _debitCour;
	res._bestLambda = lambda;
	res._margeFinale = obj_PL;


#ifdef _VERIF_
	if (!isFlotOK())
		stopProg("RCPSP_Algo_V2::construitSol_grapsV3 : flot faux sol finale");

	if (!_ins->isSolution(_debutCour, _debitCour, obj_PL))
		stopProg("RCPSP_Algo_V2::construitSol : sol. non realisable");
#endif

	return { true, res, interdit };
}





tuple<bool, pair<int, int>> RCPSP_Algo_V2::insert_etape1_fin(int jobCour, double popJCour, double LFJCour, double lambda,
	 const vector<pair<double, int>> & listePrec, vector<double> & flotCour, vector<double> & delta)
{
	bool fin = false;
	pair<int, int> interdit = { -1,-1 };

	for (int e : _ins->_jobToArc[jobCour])
	{
		flotCour[e] = 0; //qte de ressource e recu par jobCour
		bool ok_e = false; //on a donne assez de ressource e a jCour pour qu'il finisse avant sa deadline - marge lambda

		for (int k = 0; !ok_e && k < listePrec.size(); ++k)
		{
			int jprec = listePrec[k].second;


			//test ajoute le 30/03/2023 pour grasp V3 (il arrive qu'on donne plus de debit que le debit max possible => pourquoi on n'avait pas le pb avant ??)
			//si la duree de jobCour en commancant a LFJCour - lambda et finissant a _finCour[jprec] est plus petite que la duree min autorisee par le debit max => echec
			if (LFJCour - lambda - _finCour[jprec] < popJCour / _ins->_debitMax[jobCour])
			{
				interdit = { jprec, jobCour };
				break; //echec => on arrete la boucle
			}

			if (jprec == _source || _ins->_utilise[jprec][e])
			{
				double debutPossible = max(_finCour[jprec], _ins->_ES[jobCour]);

				//si jprec va fournir du flot a jobCour pour la premiere fois on ajoute jprec -> jobCour dans arcSupport (sert pour GRASP_V2)
				if (jprec != _source && _ressATransmettre[jprec][e] > EPSILON && !_dejaSupport[jprec])
				{
					_arcSupport.push_back({ jprec,jobCour });
					_dejaSupport[jprec] = true;
				}

				//si jprec a assez de ressources pour que jCour finisse a sa deadline
				if ( (_ressATransmettre[jprec][e] > EPSILON) && 
					 (debutPossible + popJCour / (flotCour[e] + _ressATransmettre[jprec][e]) <= LFJCour - lambda) )
				{
					double f = popJCour / (LFJCour - debutPossible - lambda) - flotCour[e];
					flotCour[e] += f;
					_ressATransmettre[jprec][e] -= f;
					_flotCour[jprec][jobCour][e] += f;
					delta[e] = debutPossible;
					ok_e = true;

					if (_ressATransmettre[jprec][e] < -EPSILON)
						stopProg("RCPSP_Algo_V2::insert_etape1_fin : _ressATransmettre < 0");
				}
				else
				{
					//jprec transmet tout ce qu'il a
					flotCour[e] += _ressATransmettre[jprec][e];
					_flotCour[jprec][jobCour][e] += _ressATransmettre[jprec][e];

					_ressATransmettre[jprec][e] = 0;
					//en cas d'echec on renvoie le dernier couple jprec, jobCour vu
					interdit = { jprec, jobCour };
			
				}

			}
		}//fin for (auto p : listePrec)

		//si on n'a pas reussi a fournir en ressource e => echec
		if (!ok_e)
			return { false, interdit };


	}// fin for (int e : _ins->_jobToArc[jobCour])

	return  { true, interdit };
}



tuple<bool, pair<int, int>, vector<double>> RCPSP_Algo_V2::insert_etape1_milieu(int jobCour, double popJCour, double LFJCour, double lambda,
	const vector<pair<double, int>> & listePrec, const vector<pair<double, int>> & listePrec_deb, vector<double> & flotCour, vector<double> & delta)
{
	bool fin = false;
	pair<int, int> interdit = { -1,-1 };
	vector<double> detourne (_ins->_nbArc,0);//quantite de flot detourne pour inserer jobCour (il faudra la deduire de ressAtransmettre dans insert)
	double debut_tmp = _ins->_ES[jobCour]; //debut de jobCour (va augmenter au fur et a mesure qu'on donne du flot)



	for (int e : _ins->_jobToArc[jobCour])
	{
		flotCour[e] = 0; //qte de ressource e recu par jobCour
		bool ok_e = false; //on a donne assez de ressource e a jCour pour qu'il finisse avant sa deadline - marge lambda

		// === 1. on regarde si on peut detourner du flot (on ne considere que les "trous" i.e. les i -> jsuiv tels que 
			// fin (i) <= ES(jobCour) et debut(jsuiv) >= LF(jobCour)-lambda )

		for (int k = 0; !ok_e && k < listePrec_deb.size(); ++k)
		{
			int jsuiv = listePrec_deb[k].second;

			// des qu'on tombe sur un job jsuiv qui commence trop tot par rapport a la fin de jobCour on stoppe la boucle : listePrec_deb est par debut decroissant
			if (_debutCour[jsuiv] + EPSILON < LFJCour - lambda)
				break;

			// si jsuiv utilise l'arc e, on regarde s'il y a un "trou" entre son predecesseur et lui pour mettre jobCour sans decaler les debuts et fins deja calculees
			if (_ins->_utilise[jsuiv][e])
			{
				for (int l = 0; !ok_e && l < listePrec.size(); ++l)
				{
					int prec_jsuiv = listePrec[l].second;

					//si jobCour tient dans le trou entre prec_jsuiv et jsuiv
					if (_finCour[prec_jsuiv] <= debut_tmp + EPSILON)
					{

						//si assez de ressources pour que jCour finisse a sa deadline
						if (debut_tmp + popJCour / (flotCour[e] + _flotCour[prec_jsuiv][jsuiv][e]) <= LFJCour - lambda)
						{
							double fderoute = popJCour / (LFJCour - debut_tmp - lambda) - flotCour[e];
							flotCour[e] += fderoute;
							_flotCour[prec_jsuiv][jsuiv][e] -= fderoute;
							_flotCour[prec_jsuiv][jobCour][e] += fderoute;
							_flotCour[jobCour][jsuiv][e] += fderoute;
							delta[e] = debut_tmp;
							ok_e = true;

						}
						else
						{
							// deroute tout le flot
							double fderoute = _flotCour[prec_jsuiv][jsuiv][e];
							flotCour[e] += fderoute;
							_flotCour[prec_jsuiv][jsuiv][e] = 0;
							_flotCour[prec_jsuiv][jobCour][e] += fderoute;
							_flotCour[jobCour][jsuiv][e] += fderoute;
						}
					}
					else
						break;//on arrete la boucle : les jobs qui restent dans listePrec commenceront tous trop trad pour donner du flot
				}
			}
		} //fin for listePrec_deb

		detourne[e] = flotCour[e];


		//2. on insere a la fin 
		for (int k = 0; !ok_e && k < listePrec.size(); ++k)
		{
			int jprec = listePrec[k].second;

			if (jprec == _source || _ins->_utilise[jprec][e])
			{
				double debutPossible = max(_finCour[jprec], _ins->_ES[jobCour]);

				//si jprec a assez de ressources pour que jCour finisse a sa deadline
				if (debutPossible + popJCour / (flotCour[e] + _ressATransmettre[jprec][e]) <= LFJCour - lambda)
				{
					double f = popJCour / (LFJCour - debutPossible - lambda) - flotCour[e];
					flotCour[e] += f;
					_ressATransmettre[jprec][e] -= f;
					_flotCour[jprec][jobCour][e] += f;
					delta[e] = debutPossible;
					debut_tmp = max(debut_tmp, delta[e]);
					ok_e = true;

				}
				else
				{
					//jprec transmet tout ce qu'il a
					flotCour[e] += _ressATransmettre[jprec][e];
					_flotCour[jprec][jobCour][e] += _ressATransmettre[jprec][e];

					_ressATransmettre[jprec][e] = 0;
					//en cas d'echec on renvoie le dernier couple jprec, jobCour vu
					interdit = { jprec, jobCour };

				}
			}
		}//fin for (auto p : listePrec)

		//si on n'a pas reussi a fournir en ressource e => echec
		if (!ok_e)
			return { false, interdit, detourne };


	}// fin for (int e : _ins->_jobToArc[jobCour])

	return  { true, interdit, detourne };
}

tuple<bool, pair<int, int>> RCPSP_Algo_V2::insert_etape2_fin(int e_max, int jobCour, 
	const vector<pair<double, int>> & listePrec, vector<double> & flotCour)
{
	pair<int, int> interdit = { -1,-1 };

	for (int e : _ins->_jobToArc[jobCour])
	{
		if (e != e_max)
		{

			bool ok_e = false; //on a donne assez de ressource e a jCour 


			for (int k = 0; !ok_e && k < listePrec.size(); ++k)
			{
				int jprec = listePrec[k].second;

				//si jprec va fournir du flot a jobCour pour la premiere fois on ajoute jprec -> jobCour dans arcSupport (sert pour GRASP_V2)
				if (jprec != _source && _ressATransmettre[jprec][e] > EPSILON && !_dejaSupport[jprec])
				{
					_arcSupport.push_back({ jprec,jobCour });
					_dejaSupport[jprec] = true;
				}

				if (_finCour[jprec] <= _debutCour[jobCour] + EPSILON)
				{
					//si jprec peut transferer assez a jcour pour qu'il finisse a sa deadline
					if (_ressATransmettre[jprec][e] >= _debitCour[jobCour] - flotCour[e])
					{
						ok_e = true;
						_ressATransmettre[jprec][e] = _ressATransmettre[jprec][e] - _debitCour[jobCour] + flotCour[e];
						_flotCour[jprec][jobCour][e] += _debitCour[jobCour] - flotCour[e];//supp??
						if (_ressATransmettre[jprec][e] < -EPSILON)
							stopProg("RCPSP_Algo_V2::insert_etape2_fin : _ressATransmettre < 0");

					}
					else
					{
						flotCour[e] += _ressATransmettre[jprec][e];
						_flotCour[jprec][jobCour][e] += _ressATransmettre[jprec][e];//??
						_ressATransmettre[jprec][e] = 0;
						//en cas d echec on retourne le dernier couple vu
						interdit = { jprec, jobCour };
						
					}
				}
				else
				{
					if (interdit.first == _source) //on ne peut pas interdire la source avant !
						interdit.first = jprec;
					return { false, interdit };//echec (on ne trouvera plus de jobs pour donner du flot car ils sont par ordre croissant des fins
				}
			}



			if (!ok_e)
				return { false, interdit };

		}

	}

	return { true, interdit };
}



tuple<bool, pair<int, int>> RCPSP_Algo_V2::insert_etape2_milieu(int e_max, int jobCour, double popCour,
	const vector<pair<double, int>> & listePrec, const vector<pair<double, int>> & listePrec_deb, vector<double> & flotCour)
{
	pair<int, int> interdit = { -1,-1 };

	for (int e : _ins->_jobToArc[jobCour])
	{
		if (e != e_max && _debitCour[jobCour] > flotCour[e] + EPSILON)
		{

			bool ok_e = false; //on a donne assez de ressource e a jCour 

			// === 1. on regarde si on peut detourner du flot (on ne considere que les "trous" i.e. les i -> j tels que 
			// fin (i) <= debut(jobCour) et fin(jobCour) <= debut(j) )

			for (int k = 0; !ok_e && k < listePrec_deb.size(); ++k)
			{
				int jsuiv = listePrec_deb[k].second;
				
				// des qu'on tombe sur un job jsuiv qui commence trop tot par rapport a la fin de jobCour on stoppe la boucle : listePrec_deb est par debut decroissant
				if (_debutCour[jsuiv] + EPSILON < _finCour[jobCour])
					break;

				// si jsuiv utilise l'arc e, on regarde s'il y a un "trou" entre son predecesseur et lui pour mettre jobCour sans decaler les debuts et fins deja calculees
				if (_ins->_utilise[jsuiv][e])
				{
					for (int l = 0; !ok_e && l < listePrec.size(); ++l)
					{
						int prec_jsuiv = listePrec[l].second;

						//si jobCour tient dans le trou entre prec_jsuiv et jsuiv
						if (_finCour[prec_jsuiv] <= _debutCour[jobCour] + EPSILON)
						{
							// quantite de flot a derouter
							double fderoute = min(_flotCour[prec_jsuiv][jsuiv][e], _debitCour[jobCour] - flotCour[e]);
							
							// on deroute
							_flotCour[prec_jsuiv][jsuiv][e] -= fderoute;
							_flotCour[prec_jsuiv][jobCour][e] += fderoute;
							_flotCour[jobCour][jsuiv][e] += fderoute;

							//on le deduit des ressources a transmettre puisqu'ici on le transmet direct a jsuiv
							_ressATransmettre[jobCour][e] -= fderoute;

							// on augmente la qte de flot recue pour e par jobCour
							flotCour[e] += fderoute;
							ok_e = (abs(_debitCour[jobCour] - flotCour[e]) < EPSILON);
						}
						else
							break;//on arrete la boucle : les jobs qui restent dans listePrec commenceront tous trop trad pour donner du flot
					}
				}
			} //fin for listePrec_deb

			// ==== 2. parcours des precs suivant leur fin croissante => on met a la fin
			for (int k = 0; !ok_e && k < listePrec.size(); ++k)
			{
				int jprec = listePrec[k].second;


				if (_finCour[jprec] <= _debutCour[jobCour] + EPSILON)
				{
					//si jprec peut transferer assez a jcour pour qu'il finisse a sa deadline
					if (_ressATransmettre[jprec][e] >= _debitCour[jobCour] - flotCour[e])
					{
						ok_e = true;
						_ressATransmettre[jprec][e] = _ressATransmettre[jprec][e] - _debitCour[jobCour] + flotCour[e];
						_flotCour[jprec][jobCour][e] += _debitCour[jobCour] - flotCour[e];
					}
					else
					{
						flotCour[e] += _ressATransmettre[jprec][e];
						_flotCour[jprec][jobCour][e] += _ressATransmettre[jprec][e];
						_ressATransmettre[jprec][e] = 0;
						//en cas d echec on retourne le dernier couple vu
						interdit = { jprec, jobCour };

					}
				}
				else
				{

					return { false, interdit };//echec (on ne trouvera plus de jobs pour donner du flot car ils sont par ordre croissant des fins
				}
			}
			//======



			if (!ok_e)
				return { false, interdit };

		}

	}

	return { true, interdit };
}


//init les SDD utilisee dans insert
void RCPSP_Algo_V2::initInsert()
{
	_debutCour.assign(_ins->_nbJob + 1, -1);
	_debitCour.assign(_ins->_nbJob + 1, -1);
	_finCour.assign(_ins->_nbJob + 1, -1);
	_ressATransmettre.assign(_ins->_nbJob + 1, vector<double>(_ins->_nbArc, 0));

	_debutCour[_source] = 0;
	_finCour[_source] = 0;

	for (int j = 0; j < _ins->_nbArc; ++j)
		_ressATransmettre[_source][j] = _ins->_cap[j];

	_flotCour.assign(_ins->_nbJob + 1, vector<vector<double>>(_ins->_nbJob + 1, vector<double>(_ins->_nbArc, 0)));

	_dejaSupport.assign(_ins->_nbJob, false);

}

//init les SDD utilisee dans construitSol
void RCPSP_Algo_V2::initConstruitSol()
{
	_ressATransmettre_PL.assign(_ins->_nbJob + 1, vector<double>(_ins->_nbArc, 0));
	_flotCour_PL.assign(_ins->_nbJob + 1, vector<vector<double>>(_ins->_nbJob + 1, vector<double>(_ins->_nbArc, 0)));
}


//==============================================================================
//----------- AUTRE -----------------------


//remet chaque element de v a 0
void reset(vector<vector<int>> & v)
{
	for (int i = 0; i < v.size(); ++i)
	{
		v[i].clear();
	}
}

//supprime les elements avec une proba = alpha
void menageAlea(vector<vector<int>> & v, double alpha)
{
	for (int i = 0; i < v.size(); ++i)
	{
		int j = 0;
		while (j < v[i].size())
		{
			double a = (double)(rand()) / RAND_MAX;
			if (a < alpha)
				v[i].erase(v[i].begin() + j);
			else
				j++;
		}

	}
}


//==============================================================================
// fonctions pour optimiser la solution obtenu par les flots : on garde la meme topologie et
//  on cherche la meilleure solution possible avec un PL

//ce PL est le même que BranchBoundDicho::reconstruireSolCourantePL_exact mais a ete adapte pour RCPSP_algo_v2
// le but est, a partir d'une solution obtenue heuristiquement et pour laquelle on a extrait les precedences et les cliques
// de construire une solution (i.e. date d'arrivee au safe node et debit) qui respecte toutes les contraintes avec en plus
// les memes precedence que la solution heuristique (imposer les precedences permet d'eviter les cycles)


double RCPSP_Algo_V2::reconstruireSolCourantePL_exact(const vector<pair<int, int> > & precedence, const vector<vector<int> > & cliques)
{
	double res = -1.0;

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	//	=============================================================================================
	// 0. INITIALISATION DES DONNEES

	
	vector<double> debitMin(_ins->_nbJob);
	vector<double> dead(_ins->_nbJob);

	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		dead[i] = _ins->_LF[i];
		debitMin[i] = _ins->_debitMin[i];
	}


	//==========================================================
	//variables

	IloNumVarArray debit(env, _ins->_nbJob); //debit[i] = debit du job i 
	IloNumVarArray debut(env, _ins->_nbJob); //debut[i] = debut du job i 
	IloNumVar margeMin = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, "M"); //marge minimum => quantite a maximiser

	for (int i = 0; i < _ins->_nbJob; ++i)
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


	for (int k = 2; k < nbPoint; ++k)
		S += 1.0 / k;



	for (int i = 0; i < _ins->_nbJob; ++i)
	{

		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->_debitMax[i] - debitMin[i]) / S;

		//1. debut au plus tot
		model.add(debut[i] >= _ins->_ES[i]);


		//2. la marge est plus petite ou egal que la fin des jobs : 
		//linearisation de la contrainte "dead(i) - debut(i) - pop(i)/debit(i) >= marge
		// => on approxime par la tangente

		double pas = (_ins->_debitMax[i] - debitMin[i]) / (nbPoint - 1);

		pas = y / (nbPoint - 1);

		double x = debitMin[i]; //abscisse auquel on calcul la tangente
		double x2 = x + pas;

		for (int k = 0; k < nbPoint; ++k)
		{

			double a = _ins->_pop[i] / x / x; //pop / x^2
			model.add(a * debit[i] - margeMin - debut[i] >= 2 * _ins->_pop[i] / x - dead[i]);


			//abscise pour la prochaine tangente
			x += pas;
			x2 += pas;

			if (k != nbPoint - 2)
				pas = y / (nbPoint - k - 2);
		}

	}

	//3. constraintes de precedence

	for (pair<int, int> p : precedence)
	{
		//p.first << p.second

		int i = p.first, j = p.second;

		double pas = (_ins->_debitMax[i] - debitMin[i]) / (nbPoint - 1);

		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->_debitMax[i] - debitMin[i]) / S;


		pas = y / (nbPoint - 1);

		double x = debitMin[i]; //abscisse a laquelle on calcul la tangente
		double x2 = x + pas;

		for (int k = 0; k < nbPoint; ++k)
		{


			double a = _ins->_pop[i] / x / x; //pop / x^2
			model.add(debut[j] - debut[i] + a * debit[i] >= 2 * _ins->_pop[i] / x);


			x += pas;
			x2 += pas;

			if (k != nbPoint - 2)
				pas = y / (nbPoint - k - 2);
		}

	}

	//==========================================================
	//contraintes de capacité

	//4. debits min et max
	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		model.add(debit[i] >= debitMin[i]);
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


	//cplex.exportModel("PL_cut.lp"); 


	cplex.setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	cplex.setOut(env.getNullStream());
	//cplex.setParam();

	cplex.solve();

	if (cplex.getStatus() == IloAlgorithm::Optimal)
	{
		res = cplex.getObjValue();


		//================================================================================
		//7. on genere des coupes dynamiquement

		//cout << "cutting..." << endl;

		bool fin = false;
		int cpt = 0;


		//vecteurs pour stocker la solution
		IloNumArray sol_debut(env, _ins->_nbJob);
		IloNumArray sol_debit(env, _ins->_nbJob);

		while (!fin)
		{
			bool ajout_ctr = false;

			cpt++;

			//au cas ou on boucle sans s'arreter => on renvoie qu'on ne trouve pas de solution (a priori n'arrive jamais)
			if (cpt > 1000)
				return -1.0;


			//on recupere la solution pour pouvoir ajouter plusieurs ctr d'un coup
			cplex.getValues(debut, sol_debut);
			cplex.getValues(debit, sol_debit);

			//7.1 test contrainte "marge <= fin"
			for (int i = 0; i < _ins->_nbJob; ++i)
			{


				//on teste la contrainte initiale (non lineaire)...
				if (dead[i] - sol_debut[i] - static_cast<double>(_ins->_pop[i]) / sol_debit[i] < res - EPSILON_CUT)
				{
					//... si elle n'est pas respectee alors on ajoute la tangente au point x = debit[i]

					double x = sol_debit[i];
					double a = _ins->_pop[i] / x / x; //pop / x^2
					model.add(a * debit[i] - margeMin - debut[i] >= 2 * _ins->_pop[i] / x - dead[i]);
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
				else
					stopProg("RCPSP_Algo_V2::reconstruireSolCourantePL_exact : pas de solution");
			}

		}//fin whlie (!fin)

		 //cout << " ==== nb tours cut = " << cpt << endl;


		//=============================================================================
		// on recupere la solution

		for (int i = 0; i < _ins->_nbJob; ++i)
		{
			_debutCour[i] = sol_debut[i];
			_debitCour[i] = sol_debit[i];
		}

	}
	else
		stopProg("RCPSP_Algo_V2::reconstruireSolCourantePL_exact : pas de solution");

	cplex.end();
	model.end();
	env.end();


	return res;
}


//calcul les jobs qui forment une clique (overlap) ou qui sont l'un avant l'autre en utilisant les dates de debut et debit des jobs de la solution courante
void RCPSP_Algo_V2::calculPredClique(vector <pair<int, int> > & precedence, vector <vector<int> > & clique)
{


	//1. on ordonne les dates de debut et fin

	vector<double> dates;

	//on a besoin de debut et fin d'arrivee 
	vector<double> fin(_ins->_nbJob, 0);

	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		fin[i] = _debutCour[i] + static_cast<double>(_ins->_pop[i]) / _debitCour[i];

		dates.push_back(_debutCour[i]);
		dates.push_back(fin[i]);
	}

	sort(dates.begin(), dates.end());

	auto it = unique(dates.begin(), dates.end());
	dates.erase(it, dates.end());


	//2. construction des cliques

	int K = static_cast<int> (dates.size());

	vector<vector<int>> jobsDebutent(K); //jobDebutent[k] = ensemble des jobs qui debutent a la date k
	vector<vector<int>> jobsFinissent(K); //jobsFinissent[k] = ensemble des jobs qui finissent a la date k
	clique.resize(K);



	//-----------------------------------------------------------------------
	// 2.1. on repere les jobs qui debutent / finissent en dates[k]


	for (int k = 0; k < K; ++k)
	{
		for (int i = 0; i < _ins->_nbJob; ++i)
		{
			//if (abs(debut[i] - dates[k]) < EPSILON_COMP)
			if (_debutCour[i] == dates[k])
				jobsDebutent[k].push_back(i);

			if (fin[i] == dates[k])
				jobsFinissent[k].push_back(i);
		}
	}

	//-------------------------------------------------------------------------
	// 2.2. on calcule les cliques "overlap"

	clique[0] = jobsDebutent[0];

	for (int k = 1; k < K; ++k)
	{
		clique[k] = clique[k - 1];

		//2.1 on enleve les jobs qui se finissent en k
		for (int j : jobsFinissent[k])
		{
			auto it = find(clique[k].begin(), clique[k].end(), j);

#ifdef _VERIF_
			if (it == clique[k].end())
			{
				stopProg("RCPSP_Algo : calculPredClique : on ne peut pas enlever un job qu'on n a pas mis");
			}
#endif
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


//verifie kirchoff sur les flots
bool RCPSP_Algo_V2::isFlotOK()
{
	int N = _ins->_nbJob;
	int M = _ins->_nbArc;


	// 0. on verifie que tout est positif
	for (int e = 0; e < M; ++e)
	{
		for (int i = 0; i < N; ++i) //vrai jobs (pas la source _source qui vaut N)
		{
			if (_ressATransmettre[i][e] < -EPSILON)
				return false;
			for (int j = 0; j < N; ++j)
			{
				if (_flotCour[i][j][e] < -EPSILON)
					return false;
			}
		}
	}


	// 1. pour chaque arc e, ce qui entre en job i = ce qui sort
	for (int e = 0; e < M; ++e)
	{
		for (int i = 0; i < N; ++i) //vrai jobs (pas la source _source qui vaut N)
		{
			double sumIn = _flotCour[_source][i][e]; //on commence par compter ce que donne la source
			double sumOut = _ressATransmettre[i][e]; //ce qui va au puits

			for (int j = 0; j < N; ++j)
			{
				sumIn += _flotCour[j][i][e];
				sumOut += _flotCour[i][j][e];
			}
			if (abs(sumIn - sumOut) > EPSILON_COMP)
				return false;
			if (_ins->_utilise[i][e] && _debitCour[i] > 0 && abs(sumIn - _debitCour[i]) > EPSILON_COMP)
				return false;
		}
		
	}  

	//ce qui sort de la source n excede pas la capa de l'arc

	for (int e = 0; e < M; ++e)
	{
		double sumOut = _ressATransmettre[_source][e]; //ce qui va au puits

		if (sumOut < -EPSILON)
			return false;

		for (int j = 0; j < N; ++j)
			sumOut += _flotCour[_source][j][e];

		if (sumOut > _ins->_cap[e] + EPSILON)
			return false;


	}

	return true;
}


double RCPSP_Algo_V2::PL_flot_partiel()
{
	MIP_flot_partiel MIPFlot(_ins, _inserted);
	int puits = _ins->_nbJob + 1; //indice du puits dans le MIP

	//lance l'optimisation
	double obj = MIPFlot.creeEtResout(_arcSupport);

	if (obj > 0)
	{
		//recupere la solution du MIP / PL
		_debutCour_PL.assign(_ins->_nbJob+1, -1);
		_debitCour_PL.assign(_ins->_nbJob+1, -1);
		_finCour_PL.assign(_ins->_nbJob+1, -1);
		_finCour_PL[_source] = 0;
		_debutCour_PL[_source] = 0;

		_ressATransmettre_PL.assign(_ins->_nbJob + 1, vector<double>(_ins->_nbArc, 0));
		_flotCour_PL.assign(_ins->_nbJob + 1, vector<vector<double>>(_ins->_nbJob + 1, vector<double>(_ins->_nbArc, 0)));

		for (int i = 0; i < _ins->_nbJob; ++i)
		{
			if (_inserted[i])
			{
				_debutCour_PL[i] = MIPFlot.getDebut(i + 1);//decalage de 1 a cause de la source
				_debitCour_PL[i] = MIPFlot.getDebit(i + 1);
				_finCour_PL[i] = MIPFlot.getFin(i + 1);
				for (int k = 0; k < _ins->_nbArc; ++k)
				{
					_ressATransmettre_PL[i][k] = MIPFlot.getFlot(i + 1, puits, k);
					_flotCour_PL[_source][i][k] = MIPFlot.getFlot(0, i + 1, k);
				}
			}
		}

		for (int k = 0; k < _ins->_nbArc; ++k)
			_ressATransmettre_PL[_source][k] = MIPFlot.getFlot(0, puits, k);

		for (auto p : _arcSupport)
		{
			int i = p.first, j = p.second;
			for (int k = 0; k < _ins->_nbArc; ++k)
			{
				_flotCour_PL[i][j][k] = MIPFlot.getFlot(i + 1, j + 1, k);
			}
		}
	}

	return obj;	
}



void RCPSP_Algo_V2::conserveSolPL()
{
	_debutCour = _debutCour_PL;
	_debitCour = _debitCour_PL;
	_finCour = _finCour_PL;
	_ressATransmettre = _ressATransmettre_PL;
	_flotCour = _flotCour_PL;

#ifdef _VERIF_
	if (!isFlotOK())
		stopProg("RCPSP_Algo_V2::conserveSolPL : pb flot");
#endif 

}

//supprime de la liste des arc support l'arc prec -> jobCour (prec != source) avec la date de fin la plus grande pour prec
bool RCPSP_Algo_V2::supprimeSupportfinMaxOrigine(int jobCour)
{
	bool supp = false;
	double finMax = 0;
	int bestI = -1;

	//_dejaSupport
	for (auto i = 0; i < _arcSupport.size(); ++i)
	{
		auto p = _arcSupport[i];
		if (p.second == jobCour)
		{
			int origine = p.first;
			if (finMax < _finCour[origine])
			{
				finMax = _finCour[origine];
				bestI = i;
				supp = true;
			}
		}
	}

	if (supp)
	{
		_arcSupport.erase(_arcSupport.begin() + bestI);
	}

	return supp;
}



//remplit _arcSupport en utilisant le flot courant
void RCPSP_Algo_V2::recupereArcSupport()
{

	_arcSupport.clear();

	for (int i = 0; i < _ins->_nbJob; ++i)
	{
		for (int j = 0; j < _ins->_nbJob; ++j)
		{
			for (int e : _ins->_jobToArc[i])
			{
				if (_flotCour[i][j][e] > EPSILON)
				{
					_arcSupport.push_back({ i,j });
					break;//on ne met l'arc qu'une fois (pas pour chaque e)
				}
			}
		}
	}

}