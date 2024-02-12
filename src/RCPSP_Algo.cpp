#include "RCPSP_Algo.h"
#include "PL_Optim_RCPSP.h"
#include "BSmarge.h"
#include <random>
#include <queue>
#include <algorithm>


//constructeur : lit le fichier en entree, appelle l'heuristique de RCPSP_Algo et remplit les vecteurs s,h,p avec
//la solution trouvee par RCPSP_Algo arrondie (donc degradee) (CPO ne prend que des entiers)
RCPSPHeuristic::RCPSPHeuristic(string name)
{

	Instance ins;


	//on lit le fichier d'instance (fichier genere par le programme python de Christian et Emmanuel, modifie par Alicia en juillet 2019)
	ins.readFileAlicia(name);


	//l instance peut contenir des noeuds a ala fois transit et evac ce qui n'est pas gérer par notre programme 
	//on duplique les noeuds et on adapte les distances et deadlines pour que les ensembles de noeuds evacuation et transit soit disjoints
	ins.dupliqueNodeEvacTransit();

	//dessine l'instance
	//ins.trace("arbre", false);

#ifdef _VERIF_
	if (!ins.verifArbreReduit() || !ins.traceChemin("InsChemins"))
	{
		stopProg("RCPSPHeuristic : instance non conforme");
	}
#endif


	RCPSP_Graphe graphRCPSP(&ins);


	//============================================================
	// BS heuritsique amelioree (raisonnement sur plusieurs arcs) 

	double estimationBorneSup = 0;
	BSmarge BSmarg2(&ins);
	vector<double> finOptimiste;
	vector<int> rang;

	//fournit une estimation de la meilleure marge minimum possible
	estimationBorneSup = BSmarg2.calculMargeArcAmelioree(finOptimiste, rang).first;

	//si pas de solution avec calculMargeArcAmelioree on essaie avec calculMargeArc (qui fournit une vraie borne sup)
	if (estimationBorneSup < -0.5)
	{
		BSmarge BSmarg(&ins);
		estimationBorneSup = BSmarg.calculMargeArc();
	}

	if (estimationBorneSup >= 0) //sinon pas de solution
	{

		RCPSP_Algo algo(&graphRCPSP, &ins);
		RCPSP_Algo::Resultat resRCPSP;

		//calcul d'une solution realisable en s'appuyant sur la marge estimationBorneSup
		resRCPSP = algo.GRASPstrategie2(estimationBorneSup, 5, 5);

		if (resRCPSP._solFinale > -0.5) //on a trouve une solution
		{

			//remplissage des vecteurs pour CPO
			algo.transformePourCPO();

			//on recalcule les marges en prenant en compte les dates de debut recalculees (s) avec les debits arrondis a l'entier inf (h)
			int margeCPO = algo.calculeMargeCPO();

			cout << "RCPSPHeuristic : minimal marge double = " << resRCPSP._solFinale << ", minimal marge int = " << margeCPO << endl;

			//on remplit les membres de RCPSPHeuristic qui serviront dans CPO
			isSolComputed = true;
			s = move(algo.s);
			h = move(algo.h);
			p = move(algo.p);
		}
	}
	if (!isSolComputed)
		cout << "pas de sol" << endl;
}



//=====================================================================
// algo principal


//appelle plusieurs fois la resolution (= init une solution + decscente avec le PL)
//en randomisant l'ordre initial
// maxReboot = nb de fois ou on recommence en enlevant 1 aux deadline pour essayer de generer des solutions initiales avec des marges plus petites
// maxRun = nb de fois max ou on execute l'algo (init + descente) avec des deadlines données
// maxIterInit = nb de tentatives maximal pour obtenir une solution initiale

RCPSP_Algo::Resultat RCPSP_Algo::GRASP(vector<double> & finOptimiste, int maxRun, int maxIterInit)
{
	Resultat bestRes;
	int nbChangement = 0;
	int cpt = 0; //compte le nb de tours de boucle principal


	//on va modifier les dates des deadlines des chemins, on sauvegarde pour les remettre a leur valeur initiale
	vector < vector< SomWithWindow > > cheminSav = _ins->_chemins;

	//on sotcke les fin max initiales pour les utiliser dans le PL (pour le calcul de la "vraie" fct objectif)

	_finMaxInit.assign(_ins->_lastEvacuationNode + 1, 0);


	pair<double, double> res;


	double alpha = 1;
	double alphaMax = 1.8;

	int maxSameAlpha = 2;

	while (alpha - EPSILON <= alphaMax)
	{
		cpt++;
		for (int k = 0; k < maxSameAlpha; ++k)
		{
			//-----------------------------------------------------------------------------
			// on ajuste la fin des chemins avec les fins optimistes donnees en parametre
			for (int i = 1; i <= _ins->_lastEvacuationNode; ++i)
			{
				//on commence par sauvegarder la fin max initiale
				_finMaxInit[i] = _ins->getFinMax(i);


				const int s = static_cast<int> (_ins->_chemins[i].size());

				int delta = 0;
				if (k == 0)//meme alpha pour tous les jobs
					delta = max(_ins->_chemins[i][s - 1]._LF - static_cast<int>(ceil(finOptimiste[i] * alpha)), 0);
				else
				{
					if (static_cast<double> (rand()) / RAND_MAX < 0.5)
						delta = max(_ins->_chemins[i][s - 1]._LF - static_cast<int>(ceil(finOptimiste[i] * alpha)), 0);
					else
						delta = max(_ins->_chemins[i][s - 1]._LF - static_cast<int>(ceil(finOptimiste[i] * (alpha + 0.1))), 0);
				}

				//on enleve delta partout
				for (int j = 0; j < s; ++j)
					_ins->_chemins[i][j]._LF -= delta;
			}

			//---------------------------------------------------------------------
			// test faisabilite : on doit avoir assez de ressource sur chaque arc de chaque chemin pour arriver a l'heure
			bool faisable = true;
			for (int i = 1; i <= _ins->_lastEvacuationNode; ++i)
			{
				int pred = _ins->_chemins[i][0]._som;
				for (int k = 1; k < _ins->_chemins[i].size(); ++k)
				{
					int cour = _ins->_chemins[i][k]._som;
					if (_ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / (_ins->_capArc[pred][cour]) - 1 > _ins->getFinMax(i))
					{
						//cout << _ins->getLongueurCh(i) << " + " << _ins->_pop[i] << "/" << _ins->_capArc[pred][cour] << " -1 = "
							//<< _ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / (_ins->_capArc[pred][cour]) - 1 << endl;
						//cout << "fin max = " << _ins->getFinMax(i) << endl;

						faisable = false;
					}
					pred = cour;
				}
			}

			if (faisable)
			{
				//-------------------------------------------
				//init generateur
				srand(7 * cpt);


				vector<int> ordre;

				// ordre commence par source (0) et finit par puits (nbSommet+1)
				// le premier ordre est calcule suivant une regle
				//ordre = genererOrdreInit2();
				ordre = genererOrdreInit();


				for (int i = 0; i < maxRun; ++i)
				{
					//genere une solution initial et tente d'ameliorer iterativement (i.e. augmente les debits a l'aide d'un PL)
					res = resolution(ordre, maxIterInit, nbChangement);


					if (res.second > bestRes._solFinale + EPSILON)
					{
						//cout << "AMELIORE !! : " << res.second << " remplace " << bestRes._solFinale << endl;
						bestRes = Resultat{ res.first, res.second, alpha, nbChangement };
					}

					//modif l'ordre (apres plusieurs tests pour randomiser genererOrdreInit on constate que shuffle est meilleur)
					std::random_device rd;
					std::mt19937 g(rd());

					shuffle(ordre.begin() + 1, ordre.begin() + ordre.size() - 2, g);



				}

			}


			//on restaure les chemins initiaux
			_ins->_chemins = cheminSav;
		}

		//on recommence avec une nouvelle valeur de alpha
		alpha += 0.05;
	}

	return bestRes;
}

RCPSP_Algo::Resultat RCPSP_Algo::GRASPstrategie2(double borneSup, int maxRun, int maxIterInit)
{
	Resultat bestRes;
	vector<double> bestDateDebut;
	vector<double> bestDebit;
	vector<vector<vector<double> > > bestflot;

	int nbChangement = 0;



	//on va modifier les dates des deadlines des chemins, on sauvegarde pour les remettre a leur valeur initiale
	vector < vector< SomWithWindow > > cheminSav = _ins->_chemins;

	//on sotcke les fin max initiales pour les utiliser dans le PL (pour le calcul de la "vraie" fct objectif)

	_finMaxInit.assign(_ins->_lastEvacuationNode + 1, 0);


	pair<double, double> res;

	double alpha = 1; //remettre 1 dans le cas general
	int iterAlpha = 20;

	double pas = 1.0 / iterAlpha;

	for (int ia = 0; alpha > -EPSILON && ia <= iterAlpha; ++ia)
	{
		//-----------------------------------------------------------------------------
		// on ajuste la fin des chemins avec les fins optimistes donnees en parametre

		//cout << "alpha = " << alpha << endl;
		for (int i = 1; i <= _ins->_lastEvacuationNode; ++i)
		{
			//on commence par sauvegarder la fin max initiale
			_finMaxInit[i] = _ins->getFinMax(i);

			const int s = static_cast<int> (_ins->_chemins[i].size());

			int delta = static_cast<int>(floor(alpha*borneSup + EPSILON));


			//si on a alpha = 0 alors il faut quand meme enlever 1 => une marge de 1 ici correspond a une marge de 0 dans l'article (il y
			// a un decalage d'un avec l'article)
			// dasn tous les cas, le resultat est donne avec _finMaxInit, donc sans compter le -1 ==> on fait -1 sur la marge qu'on trouve quand on affiche les resultats
			if (alpha < EPSILON)
				delta = 1;

			for (int j = 0; j < s; ++j)
				_ins->_chemins[i][j]._LF -= delta;
		}

		//maj alpha;
		alpha -= pas;

		//---------------------------------------------------------------------
		// test faisabilite : on doit avoir assez de ressource sur chaque arc de chaque chemin pour arriver a l'heure
		bool faisable = true;
		for (int i = 1; i <= _ins->_lastEvacuationNode; ++i)
		{
			int pred = _ins->_chemins[i][0]._som;
			for (int k = 1; k < _ins->_chemins[i].size(); ++k)
			{
				int cour = _ins->_chemins[i][k]._som;
				if (_ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / (_ins->_capArc[pred][cour]) - 1 > _ins->getFinMax(i))
				{
					//cout << _ins->getLongueurCh(i) << " + " << _ins->_pop[i] << "/" << _ins->_capArc[pred][cour] << " -1 = "
						//<< _ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / (_ins->_capArc[pred][cour]) - 1 << endl;
					//cout << "fin max = " << _ins->getFinMax(i) << endl;

					faisable = false;
				}
				pred = cour;
			}
		}

		if (faisable)
		{
			//-------------------------------------------
			//init generateur
			srand(7 * ia);


			vector<int> ordre;

			// ordre commence par source (0) et finit par puits (nbSommet+1)
			// le premier ordre est calcule suivant une regle
			//ordre = genererOrdreInit2();
			ordre = genererOrdreInit();


			for (int i = 0; i < maxRun; ++i)
			{
				//genere une solution initial et tente d'ameliorer iterativement (i.e. augmente les debits a l'aide d'un PL)
				//  res = resolution(ordre, maxIterInit, nbChangement);
				res = resolution_avec_PL_cut(ordre, maxIterInit, nbChangement);

				//attention : resolution_avec_PL_cut ne met pas à jour les flots (ils ne sont pas utilises dans la suite)

				//cout << res.second << " alpha = " << alpha + pas << endl;

				if (res.second > bestRes._solFinale + EPSILON)
				{
					//affiche();
					//afficheArcCritique();
					//cout << "AMELIORE !! : " << res.second << " remplace " << bestRes._solFinale << endl;
					bestRes = Resultat{ res.first, res.second, alpha + pas, nbChangement };//attention pour le alpha on a deja enleve le pas pour la prochaine it
					bestDateDebut = _dateDebut;
					bestDebit = _debit;
					bestflot = _flot;
				}

				//modif l'ordre (apres plusieurs tests pour randomiser genererOrdreInit on constate que shuffle est meilleur)
				std::random_device rd;
				std::mt19937 g(rd());
				shuffle(ordre.begin() + 1, ordre.begin() + ordre.size() - 2,g);
				//ordre = genererOrdreInit(rang); // teste et les resultats sont mauvais


			}

		}


		//on restaure les chemins initiaux
		_ins->_chemins = cheminSav;
	}


	_dateDebut = bestDateDebut;
	_debit = bestDebit;
	_flot = bestflot;
	//affiche();
	//cout << fixed << setprecision(2) << endl;

	//for (int i = 1; i <= _ins->_lastEvacuationNode; ++i)
	//{
	//	cout << i << endl;
	//	cout << "debit = " << _debit[i] << endl;
	//	cout << "debut depart " << _dateDebut[i] << endl;
	//	cout << "fin depart "  << _dateDebut[i] + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1 << endl;
	//	cout << "duree evac. " << static_cast<double>(_ins->_pop[i]) / _debit[i]  << endl;
	//	cout << "debut arrivee " << _dateDebut[i] + _ins->getLongueurCh(i) << endl;
	//	cout << "fin arrivee "  << _dateDebut[i] + static_cast<double>(_ins->_pop[i]) / _debit[i] /*- 1*/ + _ins->getLongueurCh(i) << endl;

	//	cout << " fin max " << _ins->getFinMax(i) << endl;

	//	cout << "marge = " << _ins->getFinMax(i) - _ins->getLongueurCh(i)
	//		- _dateDebut[i] - static_cast<double>(_ins->_pop[i]) / _debit[i] /*+ 1*/ << endl;

	//	cout << "--------------------------------------------------" << endl;
	//}

	return bestRes;
}

// initialise (si possible) une solution grace a une heuristique puis essaie de l'ameliorer 
//renvoie vrai si on a reussi a construire une solution, faux sinon
pair<double, double> RCPSP_Algo::resolution(vector<int> & ordre, int maxIterInit, int & nbChangement)
{
	double coutSolutionInit = -1;
	double coutSolutionFinal = -1;

	bool ok = false;


	vector<int> xMin, xMin2;

	int S = static_cast<int>(_graph->_arcs.size());

	vector<double> lambda(_graph->_nbSommet + 1);


	//====================================================================================
	//1. on calcule une solution initiale, ordre peut etre modifie si celui donne en parametre ne 
	//permet pas de trouver une solution
	ok = calculSolInit(maxIterInit, nbChangement, ordre);


	//cout << "solution initiale : " << endl;
	//traceCheminIndependant();
	//traceFlotArc(_graph->_arc2Ind[12][10]);
	//affiche();
	//afficheConsoRess();


	//traceFlotArc(11);


	//===================================================================
	// 2. on ameliore iterativement la solution initiale grace a un PL
	if (ok)
	{


#ifndef _MAX_MIN_
		coutSolutionInit = calculCoutSolSomme();
#else
		coutSolutionInit = calculCoutSolMin(xMin, xMin2);
#endif
		//affiche();
		if (!verifie())
			stopProg("RCPSP_Algo::resolution : sol init non ok");

		if (coutSolutionInit < EPSILON && coutSolutionInit > -EPSILON)
			coutSolutionInit = 0;

		double ancienCout = -1;
		int ancienSizexMin = -1;
		int ancienSizexMin2 = -1;
#ifndef _MAX_MIN_
		ancienCout = calculCoutSolSomme();
#else
		ancienCout = calculCoutSolMin(xMin, xMin2);
		ancienSizexMin = static_cast<int> (xMin.size());
		ancienSizexMin2 = static_cast<int> (xMin2.size());
#endif


		// fixer le pas (il sert pour multiplier le debit et le flot)
		double divPas = 1; // pas = 1 / divPas
		double pas = 1;


		//calcule de donnee d'entree pour le PL qui dependent de la solution initiale
		vector<pair<int, int>> arcActif; //couples qui portent du flot

		//tous les arcs actifs + ceux qu'on peut mettre en etant certain de ne pas creer de flot 
		//ni de decaler trop la date de debut
		arcActif = getArcOrdre(ordre);


		//indice de 1 a _nbSommet (debit de la sol initiale = debit pour etre realisable, on ne doit jamais repasser en dessous)
		vector<double> debitMin(_graph->_nbSommet + 1, 0);
		for (int i = 1; i <= _graph->_nbSommet; ++i)
		{
			debitMin[i] = _debit[i];
		}


		//---------------------------------------------------------
		//boucle principale d'amelioration : resol du PL (=calcul gradient projeté) 
		//puis application du gradient a la solution courante avec un pas multiplicatif

		bool stop = false;



		while (!stop && pas > 0.05)
		{

			//cout << "  ====================  boucle amelioratin pricipale  =======================  " << endl;

			//1. calcule des coefficients lamba
#ifndef _MAX_MIN_
			calculLambda(lambda);
#else
				// xMin et xMin2 sont mis a jour quand on calcule le cout de la nouvelle solution,
				// si la nouvelle solution n'est pas ameliorante alors xMin et xMin2 ne sont pas a jour
				//mais de toute facon la boucle s'arrete
			calculLambdaPourMin(lambda, xMin, xMin2);
#endif

			//3. resoudre le PL
			bool PLok;
			PL_Optim_RCPSP PL(_ins, this, _graph);
			PLok = PL.creeEtResout(debitMin, arcActif, lambda);

			//if (PLok)
				//PL.afficheSol();

			//4. si PL a une solution
			if (PLok)
			{
				//si on a une solution du PL on applique au debit et au flot la modification en consequence : 
				// on ajoute le delta avec un pas multiplicatif, si en faisant cela on degrade la solution, on reduit le pas
				// et on retente
				bool ameliore = false;

				while (!ameliore && pas > 0.05)
				{
					//repercuter les changements: 
					//4.1 mise a jour du debit avec la solution du PL
					for (int i = 1; i <= _graph->_nbSommet; ++i)
					{
						_debit[i] += PL.getDeltaDebit(i) * pas;
						//cout << "deltaD(" << i << ") * pas = " << PL.getDeltaDebit(i) * pas << endl;
					}

					//4.2 mise a jour du flot avec la solution du PL
					for (pair<int, int> a : arcActif)
					{
						int job1 = a.first;
						int job2 = a.second;

						for (int k = 0; k < _graph->_listeArcPartage[job1][job2].size(); ++k)
						{

							int o = _graph->_listeArcPartage[job1][job2][k]._orig;
							int d = _graph->_listeArcPartage[job1][job2][k]._dest;

							int e = _graph->_arc2Ind[o][d];

							_flot[job1][job2][e] += PL.getDeltaFlot(job1, job2, e) * pas; //maj flot pour e
							_flot[job1][job2][S] += PL.getDeltaFlot(job1, job2, e) * pas; //maj quantite total de flot de orig vers dest
						}

					}

					//traceFlotArc(_graph->_arc2Ind[10][9]);
					//traceFlotArc(11);

					//4.3 mise a jour des dates de debut : si deux jobs i et j se transmettent du flot (i -> j)
					//alors debut(j) >= debut(i) + TLcond(i,j,debit(i))
					auto dateDebutSav = _dateDebut; //on sauvegarde au cas où la nouvelle sol degrade, on retablira les anciens debuts
					calculeDebut();


					if (!verifie())
						stopProg("RCPSP_Algo::resolution : sol non ok ap PL");

					double nouveauCout = -1;
#ifndef _MAX_MIN_
					nouveauCout = calculCoutSolSomme();
#else
					nouveauCout = calculCoutSolMin(xMin, xMin2);
#endif

					// NB. on maximise => on veut que le nouveau cout soit strictement plus grand que l 'ancien 

					if (ancienCout > nouveauCout + EPSILON)//si on va dans le mauvais sens => on restaure la solution et on reduit le pas
					{

						// 1. restaure la sol
						// repercuter les changements: 

						for (int i = 1; i <= _graph->_nbSommet; ++i)
						{
							_debit[i] -= PL.getDeltaDebit(i) * pas;
						}


						for (pair<int, int> a : arcActif)
						{
							int job1 = a.first;
							int job2 = a.second;

							for (int k = 0; k < _graph->_listeArcPartage[job1][job2].size(); ++k)
							{

								int o = _graph->_listeArcPartage[job1][job2][k]._orig;
								int d = _graph->_listeArcPartage[job1][job2][k]._dest;

								int e = _graph->_arc2Ind[o][d];

								_flot[job1][job2][e] -= PL.getDeltaFlot(job1, job2, e) * pas; //maj flot pour e
								_flot[job1][job2][S] -= PL.getDeltaFlot(job1, job2, e) * pas; //maj quantite total de flot de orig vers dest
							}

						}

						//remettre les anciens debut
						_dateDebut = dateDebutSav;

					}
					else
					{
						ameliore = true; // on a ameliore la solution => on sort du while imbrique et on relance la PL avec des nouveaux lambas

						//si ancienCout = nouveauCout alors on fait du "sur place" => on a fini
						//remarque : on pourrait regarder si deltaDebit = 0 pour tout le monde mais c plus lent (*4)
#ifndef _MAX_MIN_
						if (abs(ancienCout - nouveauCout) < EPSILON)
#else
							//dans le cas où on maximise le min des marges, on peut avoir le même min mais moins de jobs qui atteignent ce min
							//alors on est meilleur
						if (abs(ancienCout - nouveauCout) < EPSILON
							&& (ancienSizexMin <= xMin.size()))//|| (ancienSizexMin == xMin.size() && ancienSizexMin2 == xMin2.size()) ) )
#endif
							stop = true;



						ancienCout = nouveauCout;
						ancienSizexMin = static_cast<int> (xMin.size());
						ancienSizexMin2 = static_cast<int> (xMin2.size());
					}

					// 2. reduit le pas

					divPas++;
					pas = 1.0 / (divPas - pas / 2);

					//cout << "pas = "<< pas << endl;

					//traceCheminIndependant();
					//traceFlotArc(_graph->_arc2Ind[10][9]);
					//affiche();
					//afficheConsoRess();


				}//fin while (!ameliore)
			}
			else
			{
				stop = true;
				stopProg("RCPSP_Algo::resolution : pourquoi le PL n'a pas de solution ? il suffit de prendre la solution initiale ?");
			}
		}//fin while principal : !stop
#ifndef _MAX_MIN_
		coutSolutionFinal = calculCoutSolSomme();
#else
		coutSolutionFinal = calculCoutSolMin(xMin, xMin2);
#endif

	}//fin if ok


	if (ok && !verifie())
		stopProg("RCPSP_Algo::resolution : sol init non ok");

	//traceCheminIndependant();
	return { coutSolutionInit, coutSolutionFinal };
}



// idem resolution mais on utilise le PL qui sert dans le "Branch&Bound Dicho" au moment ou on essaie de reconstruire
// une solution non preemptive en conservant l'ordre des actions

pair<double, double> RCPSP_Algo::resolution_avec_PL_cut(vector<int> & ordre, int maxIterInit, int & nbChangement)
{
	double coutSolutionInit = -1;
	double coutSolutionFinal = -1;

	bool ok = false;


	vector<int> xMin, xMin2; //necessaire pour appeler calculCoutSolMin mais non utilises

	int S = static_cast<int>(_graph->_arcs.size());

	vector<double> lambda(_graph->_nbSommet + 1);


	//====================================================================================
	//1. on calcule une solution initiale, ordre peut etre modifie si celui donne en parametre ne 
	//permet pas de trouver une solution
	ok = calculSolInit(maxIterInit, nbChangement, ordre);




	coutSolutionInit = calculCoutSolMin(xMin, xMin2);





	//cout << "solution initiale : " << endl;
	//traceCheminIndependant();
	//traceFlotArc(_graph->_arc2Ind[12][10]);
	//affiche();
	//afficheConsoRess();


	//traceFlotArc(11);


	//===================================================================
	// 2. on ameliore iterativement la solution initiale grace a un PL
	if (ok)
	{


#ifndef _MAX_MIN_
		stopProg("le PL cut / clique ne peut que maximiser le min des marges");
#endif

			//affiche();
			if (!verifie())
				stopProg("RCPSP_Algo::resolution_avec_PL_cut : sol init non ok");

		if (coutSolutionInit < EPSILON && coutSolutionInit > -EPSILON)
			coutSolutionInit = 0;

		double ancienCout = -1;

		calculeDebut();
		//affiche();


		//1. calcul des precedences et cliques 
		vector < pair<int, int> > precedence;
		vector < vector<int> > cliques;

		calculPredClique(precedence, cliques);


		//2. lance le PL et maj  _debit, _dateDebut
		double resPL = reconstruireSolCourantePL_exact( precedence,  cliques);

		//3. maj les flots ? a priori inutile car non utilises dans la suite


		coutSolutionFinal = calculCoutSolMin(xMin, xMin2);

	}//fin if ok


	//traceCheminIndependant();
	return { coutSolutionInit, coutSolutionFinal };
}


//calcul les jobs qui forment une clique (overlap) ou qui sont l'un avant l'autre en utilisant les dates de debut et debit des jobs de la solution courante
void RCPSP_Algo::calculPredClique(vector <pair<int, int> > & precedence, vector <vector<int> > & clique)
{

	int N = _graph->_nbSommet;


	//1. on ordonne les dates de debut et fin

	vector<double> dates;

	//on a besoin de debut et fin d'arrivee sur le dernier arc (ou au safe node ??? a voir...)
	vector<double> debut(N + 1, 0);
	vector<double> fin(N + 1, 0);

	for (int i = 1; i <= N; ++i)
	{
		debut[i] = _dateDebut[i] + _ins->getLongueurCh(i);
		fin[i] = debut[i] + static_cast<double>(_ins->_pop[i]) / _debit[i];

		dates.push_back(debut[i]);
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
		for (int i = 1; i <= N; ++i)
		{
			//if (abs(debut[i] - dates[k]) < EPSILON_COMP)
			if (debut[i] == dates[k]) 
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



//ce PL est le même que BranchBoundDicho::reconstruireSolCourantePL_exact mais a ete adapte aux SDD du RCPSP_Algo.
// le but est, a partir d'une solution obtenue heuristiquement et pour laquelle on a extrait les precedences et les cliques
// au niveau du dernier noeud de l'arbre d'evacuation ("safe node") 
// de construire une solution (i.e. date d'arrivee au safe node et debit) qui respecte toutes les contraintes avec en plus
// les memes precedence que la solution heuristique (imposer les precedences permet d'eviter les cycles)

//retourne la marge minimale 
//attention ce PL ne fonctionne que sur un graphe en forme d'arbre

double RCPSP_Algo::reconstruireSolCourantePL_exact(const vector<pair<int, int> > & precedence, const vector<vector<int> > & cliques)
{
	double res = -1.0;

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	//	=============================================================================================
	// 0. INITIALISATION DES DONNEES

	//dans la version initiale de ce PL ecrite pour le branch & bound (BranchBoundDicho::reconstruireSolCourantePL_exact)
	//on utilise une "InstanceReduite", il faut adapter => on recalcule ici toutes les donnees necessaires au PL

	int N = _graph->_nbSommet; //1..N

	vector<double> debitMin(N + 1);
	vector<double> dead(N + 1);

	for (int i = 1; i <= N; ++i)
	{
		dead[i] = _ins->getFinMax(i) + 1; //attention au plus 1 qui traine partout dans cette version RCPSP...
		debitMin[i] = _ins->_pop[i] / (dead[i] - _ins->getLongueurCh(i));
	}


	//==========================================================
	//variables

	IloNumVarArray debit(env, N + 1); //debit[i] = debit du job i 
	IloNumVarArray debut(env, N + 1); //debut[i] = debut du job i 
	IloNumVar margeMin = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, "M"); //marge minimum => quantite a maximiser

	for (int i = 0; i <= N; ++i)
	{
		string nomDebit = "D_" + to_string(i);
		string nomDebut = "T_" + to_string(i);

		debit[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, nomDebit.c_str());
		debut[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, nomDebut.c_str());
		//debit[0] et debut[0] ne sont pas utilisee mais on les declare pour eviter un segfault quand on recupere la sol automatiquement

	}//fin for i 

	model.add(debut[0] ==0);
	model.add(debit[0] == 0);

	//==========================================================
	//contraintes temporelles


	int nbPoint = 5; // 5*_ins->_nbJob;//nombre de tangentes par job (pour ctr 2. et 3.)
	double S = 1; // S = 1 + 1/2 + 1/3 +...+ 1/(nbPoint-1) utile si pas non constant


	for (int k = 2; k < nbPoint; ++k)
		S += 1.0 / k;



	for (int i = 1; i <= N; ++i)
	{

		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->getDebitMax(i) - debitMin[i]) / S;

		//1. debut au plus tot
		model.add(debut[i] >= _ins->getLongueurCh(i));


		//2. la marge est plus petite ou egal que la fin des jobs : 
		//linearisation de la contrainte "dead(i) - debut(i) - pop(i)/debit(i) >= marge
		// => on approxime par la tangente

		double pas = (_ins->getDebitMax(i) - debitMin[i]) / (nbPoint - 1);

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

		double pas = (_ins->getDebitMax(i) - debitMin[i]) / (nbPoint - 1);

		//on calcul y tel que debitMin = DebitMax - y * S
		// => si pas non constant alors pas = y/(nbPoint - 1 - k)
		double y = (_ins->getDebitMax(i) - debitMin[i]) / S;


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
	for (int i = 1; i <= N; ++i)
	{
		model.add(debit[i] >= debitMin[i]);
		model.add(debit[i] <= _ins->getDebitMax(i));
	}



	//5. ctr de capacité pour chaque clique

	int nbClique = static_cast<int>(cliques.size());
	bool nonVide = false;


	for (int a = 0; a < _nbArc; ++a)
	{
		int orig = _graph->_arcs[a]._orig;
		int dest = _graph->_arcs[a]._dest;


		for (int c = 0; c < nbClique; ++c)
		{
			IloExpr expr(env);

			for (int i : cliques[c])
			{
				if (_graph->_utilise[i][a])
				{
					expr += debit[i];
					nonVide = true;
				}
			}

			if (nonVide)
				model.add(expr <= _ins->_capArc[orig][dest]);

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
		//7. si on est dans le cas ou on genere des coupes dynamiquement

		//cout << "cutting..." << endl;

		bool fin = false;
		int cpt = 0;


		//vecteurs pour stocker la solution
		IloNumArray sol_debut(env, N + 1);
		IloNumArray sol_debit(env, N + 1);

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
			for (int i = 1; i <= N; ++i)
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
			}

		}//fin whlie (!fin)

		 //cout << " ==== nb tours cut = " << cpt << endl;


		//=============================================================================
		// on recupere la solution

		for (int i = 1; i <= N; ++i)
		{
			//attention dans le RCPSP le debut = le debut au noeud d'evacutaion (dans le PL c'est l'arrivee au safe node)
			_dateDebut[i] = sol_debut[i] - _ins->getLongueurCh(i);
			_debit[i] = sol_debit[i];
		}

	}

	cplex.end();
	model.end();
	env.end();


	return res;
}

//transforme la solution courante (dans _flot, _dateDebut, _debit) pour CPO 
//=> arrondir les debit a l'entier inf + recalculer les starting times avec le flot sachant que les 
//durees d'evacuation = pop / debit sont arrondies a l'entier superieur
//le resultat est stocke dans les vecteurs 
void RCPSP_Algo::transformePourCPO()
{

	//ATTENTION : entre les vecteurs h,s,p pour CPO et les vecteurs propre à RCPSP_Algo 
	//on a un decalage de 1 dans les indices (RCPSP_Algo considere une source, et pas CPO)

	h.resize(_ins->_lastEvacuationNode);
	s.resize(_ins->_lastEvacuationNode);
	p.resize(_ins->_lastEvacuationNode);

	// 1. les debits sont arrondis a l'entier inferieur
	for (int i = 0; i < _ins->_lastEvacuationNode; ++i)
	{
		h[i] = static_cast<int>(_debit[i + 1] + EPSILON); //attention au decalage de 1 sur les indices
	}

	// 2. calcule des durees d'evacuation
	for (int i = 0; i < _ins->_lastEvacuationNode; ++i)
	{
		p[i] = static_cast<int>(ceil((static_cast<double>(_ins->_pop[i + 1]) / h[i]) + EPSILON));//attention au decalage de 1 sur les indices
	}

	// 3. recalcule des starting time avec les debits entiers (h) et les precedences donnees par la solution flot
	//REMARQUE : si les noeuds a la fois transit et evac ont ete dupliques alors tous les chemins ont ete
	// allonge de 1 et aussi les deadlines => ca change les dates d'arrivee mais pas les dates de depart ni les marges

	calculeDebutPourCPO();

}



//mise a jour des dates de debut en fonction des valeurs de flots :
//si deux jobs i et j se transmettent du flot(i->j)
//alors debut(j) >= debut(i) + TLcond(i,j,debit(i))
void RCPSP_Algo::calculeDebut()
{
	int S = static_cast<int>(_graph->_arcs.size());

	queue<int> file;
	file.push(0);//source

	//1. reinit les debuts a 0
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		_dateDebut[i] = -1;
	}


	//2. boucle tant qu on augmente des dates 

	while (!file.empty())
	{
		int cour = file.front();
		file.pop();

		//on cherche les suivants : ceux a qui cour transmet du flot
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			if (_flot[cour][j][S] > EPSILON)//cour -> j
			{
				double deb = 0;
				if (cour != 0)
					deb = _dateDebut[cour] + _graph->getTLcond(cour, j, _debit[cour]);

				if (deb > _dateDebut[j] + EPSILON)
				{
					_dateDebut[j] = deb;
					file.push(j);
				}
			}
		}

	}//fin while

}



//on calcule les dates de debut pour CPO (on stocke le resultat dans le vecteur membre s)
//idem calculeDebut sauf que CPO necessite des donnees entieres => on degrade la solution en arrondissant 
//les debits a l'entier inf.
//si deux jobs i et j se transmettent du flot(i->j)
//alors debut(j) >= debut(i) + TLcond(i,j,debit(i))
//le vecteur membre h doit prealablement avoir ete calcule
//attention au decalage de 1 dans les indices entre CPO (pas de source) et RCPSP_Algo (0 = source)
void RCPSP_Algo::calculeDebutPourCPO()
{
	int S = static_cast<int>(_graph->_arcs.size());

	queue<int> file;
	file.push(0); //source

	//1. reinit les debuts a 0
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		s[i - 1] = -1;
	}

	//2. boucle tant qu on augmente des dates 
	while (!file.empty())
	{
		int cour = file.front();
		file.pop();

		//on cherche les suivants : ceux a qui cour transmet du flot
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			if (_flot[cour][j][S] > EPSILON)//cour -> j
			{
				int deb = 0;
				if (cour != 0)
					deb = s[cour - 1] + _graph->getTLcondArrondi(cour, j, h[cour - 1]);

				if (deb > s[j - 1] + EPSILON)
				{
					s[j - 1] = deb;
					file.push(j);
				}
			}
		}

	}//fin while
}

//on recalcule les marges en prenant en compte les dates de debut recalculees (s) avec les debits arrondis a l'entier inf (h)
// on retourne la marge min
int RCPSP_Algo::calculeMargeCPO()
{
	int minMarge = TIME_INFINITY;

	//Dans les calculs : attention au decalage d'un entre les vecteurs pour CPO (0..n-1)  et les autres (1..n)
	//remarque si on a duplique des noeuds on a ajoute 1 a la finMax et 1 a la longueur du chemin => ne change pas la marge
	for (int i = 0; i < _ins->_lastEvacuationNode; ++i)
	{
		//fin de l'action i : 
		int fin_i = s[i] + p[i] + _ins->getLongueurCh(i + 1);
		int marge_i = _ins->getFinMax(i + 1) - fin_i;

		if (marge_i < minMarge)
			minMarge = marge_i;
	}

	return minMarge;
}

//calcul les coeff lambda pour le PL PL_Optim_RCPSP
void RCPSP_Algo::calculLambda(vector<double> & lambda)
{
	double sum = 0;
	int S = static_cast<int>(_graph->_arcs.size());

	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		//i est dans son propre chemin critique
		sum = 1;// _ins->_pop[i];


		//on cherche les j  | i sur  chemin critique de j 
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			//le TL(i,j) doit etre > 0 pour que i soit sur le chemin critique et i,j doivent partager un arc (donc avoir une transmission de flot)
			if (i != j &&
				_graph->getTLcond(i, j, _debit[i]) > EPSILON &&
				_flot[i][j][S] > EPSILON
				)
			{
				double deb = _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]);

				//si i et j partage un/des arcs et debut(j) == debut(i) + TLcond(i,j) alors i est sur chemin critique de j
				if (abs(_dateDebut[j] - deb) <= EPSILON)
				{
					sum += 1;// _ins->_pop[j];
				}
			}
		}
		lambda[i] = sum;
		//cout << "lambda_" << i << " = " << sum << endl;
	}

}


//calcul les coeff lambda pour le PL PL_Optim_RCPSP dans le cas ou on veut maximiser le min des magres
void RCPSP_Algo::calculLambdaPourMin(vector<double> & lambda, vector<int> & xmin, vector<int> & xmin2)
{

	int S = static_cast<int>(_graph->_arcs.size());

	//1. on a un coeff 1 pour les jobs dans xmin et 1/2 pour ceux dans xmin2
	for (int k = 0; k < xmin.size(); ++k)
	{
		int i = xmin[k];
		lambda[i] = 1;
	}
	for (int k = 0; k < xmin2.size(); ++k)
	{
		int i = xmin2[k];
		lambda[i] = 0.5;
	}

	//2. on doit ajouter 1 à chaque fois que i est dans un chemin critique d'un job de xmin et 
	//   1/2  à chaque fois que i est dans un chemin critique d'un job de xmin2 
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{


		for (int k = 0; k < xmin.size(); ++k)
		{
			int j = xmin[k];

			//le TL(i,j) doit etre > 0 pour que i soit sur le chemin critique et i,j doivent partager un arc (donc avoir une transmission de flot)
			if (i != j &&
				_graph->getTLcond(i, j, _debit[i]) > EPSILON &&
				_flot[i][j][S] > EPSILON
				)
			{
				double deb = _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]);

				//si i et j partage un/des arcs et debut(j) == debut(i) + TLcond(i,j) alors i est sur chemin critique de j
				if (abs(_dateDebut[j] - deb) <= EPSILON)
				{
					lambda[i] += 1;
				}
			}
		}//fin k xmin
		for (int k = 0; k < xmin2.size(); ++k)
		{
			int j = xmin2[k];

			//le TL(i,j) doit etre > 0 pour que i soit sur le chemin critique et i,j doivent partager un arc (donc avoir une transmission de flot)
			if (i != j &&
				_graph->getTLcond(i, j, _debit[i]) > EPSILON &&
				_flot[i][j][S] > EPSILON
				)
			{
				double deb = _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]);

				//si i et j partage un/des arcs et debut(j) == debut(i) + TLcond(i,j) alors i est sur chemin critique de j
				if (abs(_dateDebut[j] - deb) <= EPSILON)
				{
					lambda[i] += 0.5;
				}
			}
		}//fin k xmin2

	}

}




//retourne la liste des arcs actifs (i.e. qui portent du flot dans la sol courante)
vector<pair<int, int>> RCPSP_Algo::getArcActif()
{
	vector<pair<int, int>> arcs;
	int S = static_cast<int> (_graph->_arcs.size());

	for (int i = 0; i <= _graph->_nbSommet; ++i)
	{
		for (int j = 1; j <= _graph->_nbSommet + 1; ++j)
		{
			if (_flot[i][j][S] > EPSILON)
			{
				arcs.push_back({ i,j });
			}
		}
	}

	return arcs;
}


//retourne tous les arcs (i,j) possible en suivant l'ordre tels que l'ajout de flot i -> j ne décale pas trop
//le debut de j
vector<pair<int, int>> RCPSP_Algo::getArcOrdre(vector<int> & ordre)
{
	vector<pair<int, int>> arcs;
	int S = static_cast<int> (ordre.size()) - 1;//-1 : puits = Cas particulier

	for (int i_i = 0; i_i < S; ++i_i)
	{
		int i = ordre[i_i];
		for (int i_j = i_i + 1; i_j < S; ++i_j)
		{

			int j = ordre[i_j];

			//on ne garde que les arcs qui ne dacale pas trop la date de debut de j
			/*double debut_j = _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]);
			double fin_j = dernier(j, debut_j);//date du dernier qui part
			*/
			//on ajoute i->j si on a deja debut(j) >= debut(i) + TL(i,j)

			if (_dateDebut[j] >= _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]))
				arcs.push_back({ i, j });

		}
		arcs.push_back({ i, S }); //on ajoute tous les i -> puits
	}

	return arcs;
}

//ordre : vecteur initialise dans la fonction, on le retourne car on en a besoin pour le PL
//d'amelioration des debits
bool RCPSP_Algo::calculSolInit(int maxIter, int & nbChangement, vector<int> & ordre)
{

	int actionEchec, actionCour;
	vector < pair<int, int> > coupleImpose;//liste des couples de sommets  dont l'ordre est impose : si (x,y) est dans coupleImpose alors x doit etre avant y dans ordre



	int cpt = 0;
	bool succes = false, echec = false;


	while (cpt < maxIter && !succes && !echec)
	{

		//cout << "cons gloutonne essai " << cpt << endl;
		cpt++;


		//algo principal
		succes = constructionGloutonne(ordre, actionEchec, actionCour);

		//modif de l'ordre random mais legere (si i < j alors il a une "bonne" proba de le rester)
		//if (!succes)//en cas de succes on doit garder l'ordre car on s'en sert dans le PL
		//{
		//	int deb = 1;
		//	if (static_cast<double> (rand()) / RAND_MAX < 0.5)
		//		deb = 2;

		//	for (int i = deb; i < ordre.size() - 2; ++i)
		//	{
		//		if (static_cast<double> (rand()) / RAND_MAX < 0.5)
		//			swap(ordre[i], ordre[i + 1]);
		//	}
		//}


		//modif totale random de l'ordre
		//if (!succes)//en cas de succes on doit garder l'ordre car on s'en sert dans le PL
		//	random_shuffle(ordre.begin() + 1, ordre.begin() + ordre.size() - 2);


		//modif "intelligente" de l'ordre
		if (!succes)
		{
			//on veut ajouter actionCour -> actionEchec
			//si ca cree un circuit
			if (circuit(coupleImpose, actionCour, actionEchec))
			{
				//on cherche x avant actionCour dans ordre, qui partage au moins un arc avec actionCour
				//et dont la date de debut est la plus grande possible
				int x = rechercheAction(ordre, actionCour); //-1 si non trouve

				if (x != -1 && !circuit(coupleImpose, actionCour, x))
				{
					coupleImpose.push_back({ actionCour, x });
				}
				else
					echec = true;
			}
			else
			{
				coupleImpose.push_back({ actionCour, actionEchec });
			}
		}
		//modifier l'ordre avec coupleImpose
		if (!echec && actionEchec != -1)
			modifierOrdre(ordre, coupleImpose);
	}


	nbChangement = cpt - 1;


	return succes;
}


//l'arc (prec,som) a ete ajoute a coupleImpose ce qui signifie que dans ordre som ne doit pas preceder prec
// => on fait les modifs necessaires
void RCPSP_Algo::modifierOrdre(vector<int> & ordre, const vector < pair<int, int> >  & coupleImpose)
{

	//dernier couple ajoute prec -> som
	int prec = coupleImpose[coupleImpose.size() - 1].first;
	int som = coupleImpose[coupleImpose.size() - 1].second;


#ifdef  _VERIF_
	//on verifie que som -> prec actuellement....
	bool trouve = false;
	for (int i = 0; i < ordre.size(); ++i)
	{
		if (ordre[i] == som)
			trouve = true;

		if (ordre[i] == prec && !trouve)
			stopProg("RCPSP_Algo::modifierOrdre : on a deja prec -> som");
	}
#endif //  _veri


	//on cherche som et tous ces successeurs et on les met apres prec

	//1. on fait un parcours en largeur a partir de som, comme ils sont deja dans ordre en respectant 
	//les relations de precendeces on peut parcourir ordre et regardaer pour chaque sommet si un de ses predesseurs
	// dasn coupleImpose a ete mis dans file, si oui on le met aussi


	vector<bool> isIn(_graph->_nbSommet + 2, false); //isIn[i] = vrai si i dans file
	vector<int> file;
	file.push_back(som);
	isIn[som] = true;

	int k = 0;


	//on positionne au nibeau de som dans ordre
	while (ordre[k] != som)
		++k;

	for (int j = k + 1; j < ordre.size(); ++j)
	{
		//si un des prec dans coupleImpose de ordre est dans la file
		if (appartientPrec(coupleImpose, isIn, ordre[j]))
		{
			file.push_back(ordre[j]);
			isIn[ordre[j]] = true;
		}
	}



	//2. on recopie ordre dans un autre vecteur (plus rapide que des decalages) en faisant attention que les 
	// sommets dans file soient recopie apres prec

	vector<int> ordre2(ordre.size());

	k = 0;
	int i = 0;
	while (ordre[i] != prec) //on sait qu on va trouver prec
	{
		if (!isIn[ordre[i]])
		{
			ordre2[k] = ordre[i];
			++k;
		}

		++i;
	}


	//arrivee la on est sur prec, on  recopie prec puis le contenu de file
	ordre2[k] = ordre[i];
	++k;
	++i;
	for (int j = 0; j < file.size(); ++j)
	{
		ordre2[k] = file[j];
		++k;
	}

	//on finit de recopier ce qui etait apres prec dans ordre2 en faisant attention de ne pas mettre le contenu de file en double


	while (i < ordre.size())
	{
		int cour = ordre[i];

		if (!isIn[cour])
		{
#ifdef _VERIF_
			if (k < 0 || k >= ordre2.size())
				stopProg("RCPSP_Algo::modifierOrdre : pb k");
#endif
			ordre2[k] = cour;
			++k;
		}
		++i;
	}


	//c ordre qui devait etre modifie...
	ordre = ordre2;
}


//renvoie vrai si un prec p (de coupleImpose) de som est tel que isIn[p] = vrai
bool RCPSP_Algo::appartientPrec(const vector<pair<int, int> > & coupleImpose, vector<bool> & isIn, int som)
{
	bool app = false;

	for (int i = 0; !app && i < coupleImpose.size(); ++i)
	{
		int prec = coupleImpose[i].first;
		int suiv = coupleImpose[i].second;

		if (suiv == som && isIn[prec])
			app = true;
	}

	return app;
}

//retourne vrai si l'ajout du couple (orig, dest) dans coupleImpose genereun cricuit
bool RCPSP_Algo::circuit(vector < pair<int, int> > & coupleImpose, int orig, int dest)
{
	bool circ = false;

	vector<int> file;
	file.push_back(dest);

	int k = 0;
	while (!circ && k < file.size())
	{
		int cour = file[k];

		//si on arrive a orig alors on a un circuit
		if (cour == orig)
			circ = true;
		else//sinon on met les successeurs dans file
		{
			for (int i = 0; i < coupleImpose.size(); ++i)
			{
				int a = coupleImpose[i].first;

				if (a == cour)
				{
					int b = coupleImpose[i].second;
					file.push_back(b);
				}
			}
			++k;
		}
	}//fin while

	return circ;
}

//on cherche x avant som dans ordre, qui partage au moins un arc avec som et dont la date de debut est la plus grande possible
//return -1 si x non trouve
int RCPSP_Algo::rechercheAction(const vector<int> & ordre, int som)
{
	int x = -1;

	double debMax = 0;

	int i = 0;
	while (ordre[i] != som)
	{
		int cour = ordre[i];
		if (_graph->_listeArcPartage[som][cour].size() >= 1)
		{
			if (_dateDebut[cour] > debMax)
			{
				debMax = _dateDebut[cour];
				x = cour;
			}
		}
		i++;
	}


	return x;

}

//alloue les vecteuirs
void RCPSP_Algo::alloc()
{
	_ressATransmettre.assign(_graph->_nbSommet + 2, vector<double>(_graph->_arcs.size(), 0));

	_flot.assign(_graph->_nbSommet + 2,
		vector < vector<double> >(_graph->_nbSommet + 2,
			vector<double>(_graph->_arcs.size() + 1, 0)));

	_dateDebut.assign(_graph->_nbSommet + 2, 0);

	_debit.assign(_graph->_nbSommet + 2, 0);

}


//renvoie vrai si construction ok
//faux sinon et place dans  (actionEchec,  actionCour) les actions qui ont echouees
bool RCPSP_Algo::constructionGloutonne(vector<int> & ordre, int & actionEchec, int & actionCour)
{
	bool ok = true;
	actionEchec = -1;
	actionCour = -1;


	//1. allocation memoire
	alloc();

	//2.algo general

	//2.1 init : la source peut transmettre tout le flot = capa des arcs

	_dateDebut[0] = 0;

	for (int k = 0; k < _graph->_arcs.size(); ++k)
	{
		int o = _graph->_arcs[k]._orig;
		int d = _graph->_arcs[k]._dest;

		_ressATransmettre[0][k] = _ins->_capArc[o][d];
	}

	//2.2 algo general

	bool stop = false;
	int cpt = 1; //ordre[0] = source, on commence donc a 1 avec les vraies actions
	int somCour = 0;


	while (!stop && cpt < ordre.size() - 1)
	{
		//recupere element courant 
		somCour = ordre[cpt];


		// 2.3 on essaie de placer somCour
		//retourne -1 si succes, sinon retourne un numero de somemet
		//cout << "assigner " << somCour << endl;
		actionEchec = assigner(ordre, somCour);

		if (actionEchec == -1)
		{
			//2.4  somCour peut donner son debit sur tous les arcs qui constituent son chemin

			int prec = _ins->_chemins[somCour][0]._som;
			for (int k = 1; k < _ins->_chemins[somCour].size(); ++k)
			{
				int cour = _ins->_chemins[somCour][k]._som;

				int ind = _graph->_arc2Ind[prec][cour];

#ifdef _VERIF_
				if (ind == -1)
					stopProg("RCPSP_Algo::constructionGloutonne indice = -1");
#endif

				_ressATransmettre[somCour][ind] = _debit[somCour];
				prec = cour;
			}

			// 2.5 pour tous les sommet p precedent somCour dans l'ordre en parametre, 
			//on retranche le flot transmis de p a somCour du flot que peut donner p

			for (int k = 0; k < cpt; ++k)
			{
				int som = ordre[k];

				for (int j = 0; j < _graph->_listeArcPartage[som][somCour].size(); ++j)
				{
					int o = _graph->_listeArcPartage[som][somCour][j]._orig;
					int d = _graph->_listeArcPartage[som][somCour][j]._dest;
					int ind = _graph->_arc2Ind[o][d];
					_ressATransmettre[som][ind] -= _flot[som][somCour][ind];
				}
			}

			//cout << "=============== cons gloutonne, cpt = " << cpt << "=================" << endl;
			//affiche();

#ifdef _VERIF_
			if (!verifieFlot(false))
				stopProg("RCPSP_Algo::constructionGloutonne : pb flot");
#endif

		}
		else
			stop = true;


		cpt++;

	}//fin while

	//valeur de retour 
	if (stop)
	{
		ok = false;
		actionCour = somCour;
	}
	else
	{
		for (int i = 0; i <= _graph->_nbSommet; ++i)
		{
			_flot[i][_puits][_nbArc] = 0;

			for (int e = 0; e < _nbArc; ++e)
			{
				// les ressources que som peut transmettre vont direct au puits
				_flot[i][_puits][e] = _ressATransmettre[i][e];
				_flot[i][_puits][_nbArc] += _flot[i][_puits][e];
			}
		}

#ifdef _VERIF_
		if (!verifieFlot())
		{
			affiche();
			stopProg("RCPSP_Algo::constructionGloutonne : pb flot");

		}
#endif
	}

	//
	//traceFlotArc(1);

	return ok;
	}

// 2.3 on essaie de placer somCour
//retourne -1 si succes, sinon retourne un numero de somemet
int  RCPSP_Algo::assigner(vector<int> & ordre, int somCour)
{
	int actionEchec = -1;
	int arcDebit = 0;

	actionEchec = assignerPhase1(ordre, somCour, arcDebit);



	if (actionEchec == -1) //pas d echec
	{



		actionEchec = assignerPhase2(ordre, somCour, arcDebit);



		if (actionEchec == -1)//pas d echec
		{
			assignerRegroupement(somCour);


		}
	}

	return actionEchec;//-1 si ok => verfier que actionEchec reste = -1 quand tout est ok
}

//pour chaque arc dans le chemin de somCour, on regarde si on peut l'alimenter en ressource "arc"
//de sorte que la fin de l'action d'evacuation de somCour ne depasse pas la date limite (LF(safenode))
//si echec renvoie un numero de noeud responsable de l echec, 
//si succes renvoie -1 et calcule 1) le debit _debit[somCour] 2) l'arc arcDebit qui donne ce debit 3) init les flots vers somCour

int RCPSP_Algo::assignerPhase1(vector<int> & ordre, int somCour, int & arcDebit)
{
	bool ok = true, stop = false;

	double debitMax = 0;
	int debitArcMax = -1;

	int actionEchec = -1;

	int prec = _ins->_chemins[somCour][0]._som;

	int k = 1;

	//while principal : parcourt le chemin qui part de somCour
	while (k < _ins->_chemins[somCour].size() && ok)
	{
		int cour = _ins->_chemins[somCour][k]._som;
		int indArcCour = _graph->_arc2Ind[prec][cour];

		//liste des sommets qui peuvent transmettre de la ressource de type [prec,cour] a somCour 
		vector<int> listeSom = construireListeSommetA1(somCour, ordre, prec, cour);

		stop = false;
		int cpt = 0;
		double debitTest = 0;

		while (!stop && cpt < listeSom.size())
		{
			int som = listeSom[cpt];
			if (_dateDebut[som] + _graph->getTLcond(som, somCour, _debit[som]) + _ins->getLongueurCh(somCour) +
				static_cast<double>(_ins->_pop[somCour]) / (debitTest + _ressATransmettre[som][indArcCour]) - 1 <= _ins->getFinMax(somCour))
			{
				//calculer augmentation max de v = w
				//som transmet suffisament de ressources a somCour pour avoir un debit suffisant pour que tout le monde soit sauf avant la deadline
				// => on augmente le flot et le debit en consequence et on passe a la suite (stop = true)
				double tmp = _ins->getFinMax(somCour) - _dateDebut[som] - _graph->getTLcond(som, somCour, _debit[som]) - _ins->getLongueurCh(somCour) + 1;
				double w = static_cast<double>(_ins->_pop[somCour]) / tmp - debitTest;

				debitTest += w;

				//attention il ne faut pas depasser le debit max
				if (debitTest <= _ins->_maxRateEvac[somCour] + EPSILON)
				{
					stop = true;

					_flot[som][somCour][_nbArc] += w - _flot[som][somCour][indArcCour];
					_flot[som][somCour][indArcCour] = w;
				}
				else //on depasse le debit max ==> echec!
					cpt = static_cast<int>(listeSom.size());//pour arreter la boucle
			}
			else
			{
				//som ne peut pas trabsmettre suffisament de ressource => il transmet ce qu'il peut et la boucle while continue...
				debitTest += _ressATransmettre[som][indArcCour];
				_flot[som][somCour][_nbArc] += _ressATransmettre[som][indArcCour] - _flot[som][somCour][indArcCour];
				_flot[som][somCour][indArcCour] = _ressATransmettre[som][indArcCour];
			}
			cpt++;

#ifdef _VERIF_
			if (!verifieFlotSum())
				stopProg("RCPSP_Algo::assignerPhase1 : pb flot ");
#endif

		}//fin while interne

		//si stop = false cela signifie qu'on n'a pas reussi a faire passer assez de ressource vers somCour
		if (!stop)
		{
			double valMax = -1;

			// actionEchec = action dasn listeSom avec le plus grand debut + TL
			for (int i = 0; i < listeSom.size(); ++i)
			{
				int som2 = listeSom[i];
				double val = _dateDebut[som2] + _graph->getTLcond(som2, somCour, _debit[som2]);
				if (val > valMax)
				{
					valMax = val;
					actionEchec = som2;
				}
			}

		}
		else//on a reussi a donner assez de flot a somCour
		{
			//on sauvegarde le plus grand debit qui a ete affecte a un arc du chemin => dans la seconde phase (assignerPhase2, on devra
			//avoir le meme debit pour tous les arcs du chemin)
			if (debitTest > debitMax)
			{
				debitMax = debitTest;
				debitArcMax = indArcCour;
			}
		}

		prec = cour;
		k++;
	}//fin while principal qui parcourt les arcs du chemin

	if (actionEchec == -1) //reussite
	{
		_debit[somCour] = debitMax;
		arcDebit = debitArcMax;
	}

	//cout << "ASSIGNER PHASE 1" << endl;
	//affiche();


#ifdef _VERIF_
	if (!verifieFlotSum())
		stopProg("RCPSP_Algo::assignerPhase1 : pb flot ");
#endif
#ifdef _VERIF_
	if (actionEchec == 0)
		stopProg("RCPSP_Algo::assignerPhase1 : source ne peut pas etre en echec");
#endif

	return actionEchec;
}

//construit la liste des sommets utilisee dans assignerPhase1 et assignerPhase2 : 
// ce sont les sommets susceptobles de transmettre de la ressource e = (orig,dest) à somCour
// les sommets x avant somCour dans ordre, tels que 
// 1. l'arc e =  (orig, dest) est dans le chemin partant de x 
// 2. _ressATransmettre[x][e] > 0

vector<int> RCPSP_Algo::construireListeSommetA1(int somCour, vector<int> & ordre, int orig, int dest)
{
	vector<pair<double, int>> vtmp;
	vector<int> v;
	int cpt = 0;

	int ind = _graph->_arc2Ind[orig][dest];

	while (ordre[cpt] != somCour)
	{
		int som = ordre[cpt];

		//si som a des ressources de type arc (orig,dest) a transmettre  
		//attention si arc (orig,dest) n appartient au chemin de som alors il faut que _ressATransmettre[som][ind] = 0
		if (_ressATransmettre[som][ind] > EPSILON)
		{
			double val = _dateDebut[som] + _graph->getTLcond(som, somCour, _debit[som]);
			vtmp.push_back(pair<double, int>(val, som));
		}
		cpt++;
	}

	sort(vtmp.begin(), vtmp.end());
	v.resize(vtmp.size());

	for (int i = 0; i < v.size(); ++i)
		v[i] = vtmp[i].second;

	return v;
}


int RCPSP_Algo::assignerPhase2(vector<int> & ordre, int somCour, int arcDebit)
{
	int actionEchec = -1;
#ifdef _VERIF_
	if (!verifieFlotSum())
	{
		stopProg("RCPSP_Algo::assignerPhase2 : pb flot");
	}
#endif


	//on parfcourt tous les arcs du chemin de somCour a lexception de arcDebit


	int prec = _ins->_chemins[somCour][0]._som;
	int k = 1;
	bool echec = false;

	//on parcourt tous les arcs dans le chemins de somCour, en evitant arcDebit : 
	//le but est que tous les arcs recoivent la quantite de flot _debit[somCour] (= plus grand debit possible sur le chemin de somCour)
	while (k < _ins->_chemins[somCour].size() && !echec)
	{
		int cour = _ins->_chemins[somCour][k]._som;
		int indArcCour = _graph->_arc2Ind[prec][cour];


		if (indArcCour != arcDebit)//on fait les modifs pour tous les arcs du chemin sauf arcDebit
		{
			//liste des sommets avavt somCour dans ordre et susceptibles de transmettre de la ressource (prec,cour)
			vector<int> listeSom = construireListeSommetA1(somCour, ordre, prec, cour);

			bool succes = false;
			bool echec2 = false;
			double v = 0;//qtantite de flot recu par somCour pour l'arc courant

			int cpt = 0;
			while (!succes && !echec2 && cpt < listeSom.size())
			{
				int som = listeSom[cpt];

				//si som ne peut pas transmetre de ressource a somCour (car decale trop la date de debut de somCour) => echec pour som
				if (_dateDebut[som] + _graph->getTLcond(som, somCour, _debit[som]) + _ins->getLongueurCh(somCour) - 1 +
					static_cast<double>(_ins->_pop[somCour]) / _debit[somCour] > _ins->getFinMax(somCour) + EPSILON)
				{
					echec2 = true;
					actionEchec = som;
					/*#ifdef _VERIF_
										cout << _dateDebut[som] << " + " << _graph->getTLcond(som, somCour, _debit[som]) << " + " <<
											_ins->getLongueurCh(somCour) << "- 1 + "
											<< _ins->_pop[somCour] << " / " << _debit[somCour] << " > " << _ins->getFinMax(somCour) << endl;
					#endif*/
				}
				else //som peut transmettre de la ressource a somCour
				{
					//si som peut transmettre tout ce qu'il manque, il le fait et succes devient vrai
					if (_debit[somCour] - v <= _ressATransmettre[som][indArcCour] + EPSILON)
					{
						_flot[som][somCour][_nbArc] += _debit[somCour] - v - _flot[som][somCour][indArcCour];
						_flot[som][somCour][indArcCour] = _debit[somCour] - v;
						succes = true;
					}
					else//sinon som transmet le max de ce qu'il peut pour cet arc 
					{
						v += _ressATransmettre[som][indArcCour];
						_flot[som][somCour][_nbArc] += _ressATransmettre[som][indArcCour] - _flot[som][somCour][indArcCour];
						_flot[som][somCour][indArcCour] = _ressATransmettre[som][indArcCour];
					}


				}
				cpt++;
			}//fin while parcourt listeSom

#ifdef _VERIF_
			if (!succes && !echec2)//on n'a pas reussi a donner assez de flot a somCour
			{

				stopProg("RCPSP_Algo::assignerPhase2 : si on n'a pas d echec on a forcmeent un succes car tous les jobs qui peuvent fournir du flot sont avant somCour dans ordre");
			}
#endif
			//=====================================
		}


		k++;
		prec = cour;
	}//fin while principal

#ifdef _VERIF_
	if (!verifieFlotSum())
	{
		stopProg("RCPSP_Algo::assignerPhase2 : pb flot");
	}
#endif
#ifdef _VERIF_
	if (actionEchec == 0)
		stopProg("RCPSP_Algo::assignerPhase2 : source ne peut pas etre en echec");
#endif

	return actionEchec;
}


//si pour un sommmet s donné s il n existe qu'un seul arc e tel que flot(s, somCour, e) > 0
// alors on essaie de vider cet arc en augmentant flot (y, somCour, e) pour un y qui partage 
// des arcs avec somCour et qui transmet deja de la ressource a somCour
void RCPSP_Algo::assignerRegroupement(int somCour)
{
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		//si i et somCour partage des arcs
		if (i != somCour && _graph->_listeArcPartage[i][somCour].size() >= 1) //
		{
			//on compte pour combien d arc on a du flot de i vers somCour
			int cpt = 0;
			int orig, dest, arc;

			for (int k = 0; k < _graph->_listeArcPartage[i][somCour].size(); ++k)
			{
				orig = _graph->_listeArcPartage[i][somCour][k]._orig;
				dest = _graph->_listeArcPartage[i][somCour][k]._dest;

				arc = _graph->_arc2Ind[orig][dest];

				if (_flot[i][somCour][arc] > EPSILON)
					cpt++;
			}

			//si cpt = 1 alors on a trouve un sommet i qui envoie une seule ressource arc vers somCour

			if (cpt == 1)
			{
				bool stop = false;
				for (int j = 1; !stop && j <= _graph->_nbSommet; ++j)
				{
					//si l'arc est a la fois dans le chemin de i et de j, et que j transmet cette ress a somCour
					if (j != i && j != somCour && _graph->isDansArcPartage(i, j, orig, dest)
						&& _flot[j][somCour][arc] > EPSILON
						&& _ressATransmettre[j][arc] >= _flot[i][somCour][arc] + _flot[j][somCour][arc])
					{
						_flot[j][somCour][_nbArc] += _flot[i][somCour][arc];
						_flot[j][somCour][arc] += _flot[i][somCour][arc];

						_flot[i][somCour][_nbArc] -= _flot[i][somCour][arc];
						_flot[i][somCour][arc] = 0;

						stop = true;
					}
				}

			}


		}
	}//fin for i


	//on cherche tous les i qui donne du flot a somCour
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		if (_flot[i][somCour][_nbArc] > EPSILON)
		{
			_dateDebut[somCour] = max(_dateDebut[somCour], _dateDebut[i] + _graph->getTLcond(i, somCour, _debit[i]));
		}
	}


}


//affiche la solution avec le detail des flots
void RCPSP_Algo::affiche()
{
	//on affiche uniquement le flot entre "vraies" actions pour eviter de surcharger l'affichage
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			for (int e = 0; e < _graph->_arcs.size(); ++e)
			{
				int o = _graph->_arcs[e]._orig;
				int d = _graph->_arcs[e]._dest;
				if (_flot[i][j][e] > EPSILON)
					cout << "flot[" << i << "][" << j << "][" << o << "->" << d << "] = " << _flot[i][j][e] << endl;
			}
		}
	}

	//traceFlotArc(5);

	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		cout << "debut " << i << " = " << _dateDebut[i] << " ET ";
		cout << "fin " << i << " = " << _dateDebut[i] + _ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1 << " ET ";
		cout << "debit " << i << " = " << _debit[i] << endl;
	}

	cout << "cout (somme) = " << calculCoutSolSomme() << endl;

	vector<int> x1, x2;
	cout << "cout (min) = " << calculCoutSolMin(x1, x2) << endl;
}

//affiche la solution avec le detail des flots
void RCPSP_Algo::afficheConsoRess()
{
	struct sav
	{
		int _numJob;
		double _debut;
		double _fin;
		double _debit;
	};

	vector<vector<sav>> vSav;
	vSav.resize(_graph->_arcs.size() + 1);

	//================================================================
	// etape 1 : calcul des infos 

	//on aprcourt chaque chemin et pour chaque arc du chemin on stocke l'arrive du job dessus et son depart et sa conso (=debit)
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		//on commence au deuxieme arc du chemin car le premier n'est pas partage avec un autre job, inutile de l'afficher
		int prec = _ins->_chemins[i][1]._som;


		double arrivee = _dateDebut[i] + _ins->_time[_ins->_chemins[i][0]._som][prec];//date d'arrivee en prec
		double depart = dernier(i, arrivee);

		for (int j = 2; j < _ins->_chemins[i].size(); ++j)
		{
			int cour = _ins->_chemins[i][j]._som;

			int ind = _graph->_arc2Ind[prec][cour];


			vSav[ind].push_back({ i, arrivee, depart, _debit[i] });

			//mise a jour pour les suivants
			arrivee += _ins->_time[prec][cour]; //arrivee du rpemier sur l'arc prec->cour
			depart = dernier(i, arrivee);//arrivee du dernier sur l'arc prec->cour

			prec = cour;
		}
	}//fin for


	//================================================================
	// etape 2 : affichage

	for (int k = 0; k < vSav.size(); ++k)
	{
		if (!vSav[k].empty())
		{
			cout << "arc ressource : " << _graph->_arcs[k]._orig << "->" << _graph->_arcs[k]._dest << "  ("
				<< _ins->_capArc[_graph->_arcs[k]._orig][_graph->_arcs[k]._dest] << ")" << endl;

			for (int j = 0; j < vSav[k].size(); ++j)
			{
				cout << vSav[k][j]._numJob << " : [" << vSav[k][j]._debut << "," << vSav[k][j]._fin << "], D = "
					<< vSav[k][j]._debit << endl;
			}
			cout << endl;
		}
	}
}

void RCPSP_Algo::traceFlotArc(int e)
{
	string name = "flotrcpsp";

	ofstream ficOut("flotrcpsp.txt", ios::out);// ouverture fichier de sortie

	//premiere ligne du fichier de sortie
	ficOut << "digraph G {" << endl;

	//graphe horizontal :
	ficOut << "rankdir=LR" << endl;

	for (int i = 0; i <= _graph->_nbSommet; ++i)
	{
		for (int j = 1; j <= _graph->_nbSommet + 1; ++j)
		{

			if (_flot[i][j][e] > EPSILON)
			{
				ficOut << i << " -> " << j << "[label=" << _flot[i][j][e];
				if (i == 0 || j == _graph->_nbSommet + 1)
					ficOut << ",style=dashed, color=grey";
				ficOut << "]" << endl;
			}
		}

	}//fin for i 



	ficOut << "}" << endl;
	ficOut.close();

	stringstream ligCom;


#if defined(_WIN32) || defined(_WIN64)
	ligCom << "dot.exe -Tpdf -o" << name << ".pdf " << name << ".txt";
	system(ligCom.str().c_str());
#endif


}


//dessine chaque chemin de manière independante 
void RCPSP_Algo::traceCheminIndependant()
{
	string name = "cheminIndpt";

	stringstream ss;
	ss << name << ".txt";
	ofstream ficOut(ss.str().c_str(), ios::out);// ouverture fichier de sortie

	ficOut << fixed << setprecision(2);


	ficOut << "digraph G {" << endl;

	int cpt = 0;

	//etape 1 : on est obligé de renumeroter les jobs pour avoir plusieurs graphes independants
	// => etape 1 donne un label au numero

	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{

		for (int j = 0; j < _ins->_chemins[i].size(); ++j)
		{
			int cour = _ins->_chemins[i][j]._som;

			ficOut << cpt << "[label=" << cour << "] " << endl;
			cpt++;
		}
	}

	// on definit des noeuds carre pour mettre les dates de depart et fin :
	// noeud numero 10000 + i pour le depart et noeud numero 20000 + i 
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		ficOut << 10000 + i << "[shape = square, label=" << _dateDebut[i] << "] " << endl;
		ficOut << 20000 + i << "[shape = square, label=<" << dernierSafe(i) << " <br /> <FONT COLOR=\"red\">"
			<< _ins->getFinMax(i) - dernierSafe(i) << "</FONT>>]" << endl;

		//cout << _ins->getFinMax(i) << " - " << dernierSafe(i) << endl;

	}

	//etape 2 : on est obligé de renumeroter les jobs pour avoir plusieurs graphes independants
	// => etape 2 donne un label au numero

	cpt = 0;
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{

		//A[label = "First Node" shape = "circle"]

		int prec = _ins->_chemins[i][0]._som;

		ficOut << 10000 + i << " -> " << cpt << "[style = dashed, color = grey]" << endl;

		for (int j = 1; j < _ins->_chemins[i].size(); ++j)
		{
			int cour = _ins->_chemins[i][j]._som;

			//on recalcul le debit sur cet arc a l'aide des flots car pdt l'algo les debits
			//ne sont pas forcement mis a jour mais les flots sont (théoriquement...) toujours à jour

			double conso = calculConsoFromFlot(i, _graph->_arc2Ind[prec][cour]);
			ficOut << cpt << " -> " << cpt + 1 << "[label=" << conso << "]" << endl;

			cpt++;

			prec = cour;
		}

		ficOut << cpt << " -> " << 20000 + i << "[style = dashed, color = grey]" << endl;

		cpt++;

	}


	ficOut << "}" << endl;
	ficOut.close();

	stringstream ligCom;


#if defined(_WIN32) || defined(_WIN64)
	ligCom << "dot.exe -Tpdf -o" << name << ".pdf " << name << ".txt";
#endif

	system(ligCom.str().c_str());
}


//renvoie le debit de i sur l'arc e ( recalcule avec les valeurs de flot)
double RCPSP_Algo::calculConsoFromFlot(int i, int e)
{
	double conso = 0;

	for (int j = 0; j <= _graph->_nbSommet; ++j)
	{
		conso += _flot[j][i][e];
	}

	return conso;
}


bool RCPSP_Algo::verifieFlotSum()
{
	bool ok = true;

	// 1. on verifie que le flot est coherent : on a bien mis la somme des flots sur les diffentes ressources de i vers j dans _flot[i][j][S]

	// 1.1 dans  _flot[i][j][S] on verifie qu on a bien la somme des flots de i a j sur tous les e
	int S = static_cast<int>(_graph->_arcs.size());

	for (int i = 0; ok && i <= _graph->_nbSommet + 1; ++i)
	{
		for (int j = 0; j <= _graph->_nbSommet + 1; ++j)
		{
			double sum = 0;
			for (int e = 0; e < _graph->_arcs.size(); ++e)
			{
				sum += _flot[i][j][e];
			}

			if (abs(sum - _flot[i][j][S]) > EPSILON)
				ok = false;
		}

	}
	return ok;
}

bool RCPSP_Algo::verifieFlot(bool verifKirch)
{


	// 1. on verifie que le flot est coherent : on a bien mis la somme des flots sur les diffentes ressources de i vers j dans _flot[i][j][S]

	// 1.1 dans  _flot[i][j][S] on verifie qu on a bien la somme des flots de i a j sur tous les e
	bool ok = verifieFlotSum();
	int S = static_cast<int>(_graph->_arcs.size());

	//

	if (verifKirch)
	{
		// 2. kirchhoff :pour chaque i lois de kirchhoff
		for (int i = 1; i <= _graph->_nbSommet; ++i)
		{
			for (int e = 0; e < _graph->_arcs.size(); ++e)
			{
				double sumIn = 0, sumOut = 0;

				for (int j = 0; j <= _graph->_nbSommet + 1; ++j)
				{

					sumOut += _flot[i][j][e];
					sumIn += _flot[j][i][e];
				}

				if (abs(sumIn - sumOut) > EPSILON)
					ok = false;
			}

		}


		// 3. ce qui sort de la source doit etre egal a la capa de l'arc
		//et ce qui entre au puits aussi
		for (int e = 0; e < _graph->_arcs.size(); ++e)
		{
			double sumOut = 0;
			double sumIn = 0;

			for (int j = 1; j <= _graph->_nbSommet + 1; ++j)
			{
				sumOut += _flot[0][j][e];
			}

			for (int j = 0; j <= _graph->_nbSommet; ++j)
			{
				sumIn += _flot[j][_graph->_nbSommet + 1][e];
			}

			int o = _graph->_arcs[e]._orig;
			int d = _graph->_arcs[e]._dest;

			if (abs(_ins->_capArc[o][d] - sumOut) > EPSILON)
				ok = false;

			if (abs(_ins->_capArc[o][d] - sumIn) > EPSILON)
				ok = false;
		}
	}

	// 4. si flot de i vers j alors debut(j) >= debut(i) + TL(i;j)
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			if (_flot[i][j][S] > EPSILON && _dateDebut[j] < _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]) - EPSILON)
			{
				ok = false;
				cout << _dateDebut[j] << "< ?" << _dateDebut[i] + _graph->getTLcond(i, j, _debit[i]) << endl;

			}
		}
	}

	// 5. [verif du TL] si flot de i vers j alors debut(j) sur e >= fin(i) sur e + 1 pour tous les arcs e partages par i et j

	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		for (int j = 1; j <= _graph->_nbSommet; ++j)
		{
			for (ArcCommun & A : _graph->_listeArcPartage[i][j])
			{
				int e = _graph->_arc2Ind[A._orig][A._dest];

				if (_flot[i][j][e] > EPSILON)
				{
					//debut de i sur e : 
					double deb_i = _dateDebut[i];
					double deb_j = _dateDebut[j];
					if (i < j)
					{
						deb_i += A._long1;
						deb_j += A._long2;
					}
					else
					{
						deb_i += A._long2;
						deb_j += A._long1;
					}

					double fin_i = deb_i + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1;

					//j ne peut passer que apres 1 unite de temps + fin de i (car i fini en 10 ==> il emprunte l'arc en 10 pendant une inité de temps)
					//il faut deb_j >= fin_i + 1
					if (fin_i + 1 - EPSILON > deb_j)
					{
						ok = false;
						cout << "RCPSP_Algo::verifieFlot : probleme TL : " << fin_i + 1 << "<=" << deb_j << endl;
						cout << "TL = " << _graph->getTLcond(i, j, _debit[i]) << endl;
					}

				}
			}
		}
	}


	return ok;
}

bool RCPSP_Algo::verifie()
{
	double minMarge = TIME_INFINITY;

	bool ok = true;
	vector< vector<debitTemporise> > arcDebit(_graph->_arcs.size()); //doit etre trie dans l ordre des debut croissant


	//1. on verifie que le flot est coherent

	verifieFlot();



	// si le premier part en t d'une origine o vers une dest d,
	// le dernier part en t  + (population / debit) - 1

	// le premier arrive en t + time[o][d] en d
	// le dernier arrive en t + time[o][d] + (population / debit) - 1


	//pour chaque chemin on verifie que le dernier arrive avant la deadline 

	for (int i = 1; ok && i <= _graph->_nbSommet; ++i)//pour chaque chemin i
	{
		double arrivePremierEnCour = _dateDebut[i];//arrivee du premier en cour
		int prec = _ins->_chemins[i][0]._som;

		for (int j = 1; ok && j < _ins->_chemins[i].size(); ++j)//pour chaque arc du chemin prec->cour
		{
			int cour = _ins->_chemins[i][j]._som;
			int indArc = _graph->_arc2Ind[prec][cour];

			//arrivee du premier au sommet  cour
			arrivePremierEnCour += _ins->_time[prec][cour];
			double arriveDernierEnCour = arrivePremierEnCour + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1;

			//on ajoute (_ins->_pop[i])/_debit[i] - 1 pour avoir la date du dernier arrivee
			if (arriveDernierEnCour > _ins->_chemins[i][j]._LF + EPSILON)
			{
				cout << " arriveDernierEnCour = " << arriveDernierEnCour << " LF = " << _ins->_chemins[i][j]._LF << endl;
				ok = false;
			}
			else
			{
				minMarge = min(minMarge, _ins->_chemins[i][j]._LF - arriveDernierEnCour);

				//ajouterArcDebit :  ajoute le debit et  verifie en meme temps qu on ne depasse pas la capacite de l'arc
				//attention on passe sur l'arc avant d'arriver en cour => on retranche donc la longuer de l'arc
				ok = ajouterArcDebit(arcDebit[indArc], _ins->_capArc[prec][cour], _debit[i],
					arrivePremierEnCour - _ins->_time[prec][cour], arriveDernierEnCour + 1 - _ins->_time[prec][cour]);// 23/10/2019 : ajout du +1 

				if (!ok)
				{
					cout << "pb debit" << endl;
					afficheUtilisation(indArc);
				}

			}
			prec = cour;
		}

	}

	//cout << minMarge << endl;

	return ok;

}

void RCPSP_Algo::afficheUtilisation(int indArc)
{

	for (int i = 1; i <= _graph->_nbSommet; ++i)//pour chaque chemin i
	{
		double arrivePremierEnCour = _dateDebut[i];//arrivee du premier en cour

		int prec = _ins->_chemins[i][0]._som;

		for (int j = 1; j < _ins->_chemins[i].size(); ++j)//pour chaque arc du chemin prec->cour
		{
			int cour = _ins->_chemins[i][j]._som;
			int indArcCour = _graph->_arc2Ind[prec][cour];

			if (indArcCour == indArc)
			{
				//arrivee du premier au sommet  cour
				arrivePremierEnCour += _ins->_time[prec][cour];
				double arriveDernierEnCour = arrivePremierEnCour + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1;

				cout << "chemin " << i << " debit = " << _debit[i]
					<< " -- 1er arrive sur arc en " << arrivePremierEnCour - _ins->_time[prec][cour]
					<< "  dernier en " << arriveDernierEnCour - _ins->_time[prec][cour] << endl;


			}
			prec = cour;
		}

	}
}

void RCPSP_Algo::afficheArcCritique()
{
	vector<int> x1, x2;
	calculCoutSolMin(x1, x2);


#ifdef _VERIF_
	//on verifie qu on a bien fait la somme des flots i->j avant de l'utiliser
	if (!verifieFlotSum())
		stopProg("RCPSP_Algo::afficheArcCritique: flot somme pas bon");
#endif

	//x1 contient les jobs critiques
	int S = static_cast<int>(_graph->_arcs.size());
	int cpt = 0;
	for (int i = 0; i < x1.size(); ++i)
	{
		int job1 = x1[i];
		for (int j = 0; j < x1.size(); ++j)
		{
			int job2 = x1[j];
			if (_flot[job1][job2][S] > EPSILON)
			{
				cout << job1 << " -> " << job2 << endl;
				cpt++;
			}
		}
	}

	if (cpt >= 1)
		cin >> cpt;
}


//ajoute un nouveau debit sur l'intervalle [debut,fin]
//le vecteir arcDebit contient deja une liste d'intervalle [deb,fin] avec une valeur de debit pour chaque intervalle
//il faut donc ajouter le nouvel interval[debut, fin] avec son debit en gardant les intervalles de arcDebit sans intersection et trie sur les deb croissants: 
// si l'intervalle est contenu dans un autre il faut faire la somme des debit, s'il est à cheval sur un autre il faut decouper en intervalle plus petits etc...
//retourne vrai si on peut ajouter le nouveau debit sans jamais depasser la capacité de l'arc (orig, dest), faux sinon
bool RCPSP_Algo::ajouterArcDebit(vector<debitTemporise> & arcDebit, const int CAP, double debit, double debut, double fin)
{
	bool ok = true;
	bool stop = false;//stop devient vrai quand on a fini de mettre la liste a jour




	if (CAP < debit - EPSILON)
		ok = false;

	//trivial : liste vide ou nouvel interval avant premier de la liste ou apres dernier de la liste
	if (arcDebit.size() == 0 || fin - EPSILON <= arcDebit[0]._deb || debut >= arcDebit[arcDebit.size() - 1]._fin - EPSILON)
	{
		arcDebit.push_back(debitTemporise(debut, fin, debit));
		stop = true;
	}


	// 2. sinon [debut, fin] intersecte un ou plusieurs intervalles : on va avancer dans la liste arcDebit
	//  et on va reduire [debut, fin] au fur et a mesure qu'on intersecte des intervalles dans arcDebit : 
	// [debut, fin] va donc etre decoupe au fur et a mesure, a chaque decoupage deb est augmente (on le renomme t)

	double t = debut; //debut du prochain intervalle
	vector<int> intervalleAsupp; //liste des intervalles a supprimer (on pourrait le faire en cours de creation des nouveaux mais ca complique l'algo)


	int i = 0; //indice de l'intervalle de arcDebit qu'on intersecte au tour de boucle courant

	int sizeInit = static_cast<int>(arcDebit.size());

	//2.1. on cherche le premier intervalle susceptible de s'intersecter avec [debut, fin]
	while (i < sizeInit && t >= arcDebit[i]._fin - EPSILON)
		i++;

	//2.2. tant qu'on a pas intergre le nouvel intervalle [debut, fin] avec son debit on continue
	while (ok && !stop && t + EPSILON < fin)
	{

		//--------------------------------------------------------------------------------------------------------
		//1er cas : t est au milieu d'un intervalle initial existant => l'intervalle existant va etre scinde en 
		// 2 ou 3 parties 
		if (i < sizeInit && t >= arcDebit[i]._deb - EPSILON)
		{
			//i sera supprime car on cree plusieurs intervalles a sa place
			intervalleAsupp.push_back(i);

			//1er morceau : de arcDebit[i]._deb a t si t != arcDebit[i]._deb 
			if (t > arcDebit[i]._deb + EPSILON)
				arcDebit.push_back(debitTemporise(arcDebit[i]._deb, t, arcDebit[i]._debit));

			//les deux fin tombent en meme temps = > on a juste un deuxieme intervalle a creer
			if (abs(fin - arcDebit[i]._fin) < EPSILON)
			{
				arcDebit.push_back(debitTemporise(t, fin, arcDebit[i]._debit + debit));
				ok = arcDebit[i]._debit + debit <= CAP + EPSILON;
				stop = true;
			}
			else //une fin precede l'autre 
			{
				//l intervalle [t, fin] est entierement dans l intervalle i => on va creer 2 autres morceaux
				if (fin < arcDebit[i]._fin)
				{
					//2sd morceau
					arcDebit.push_back(debitTemporise(t, fin, arcDebit[i]._debit + debit));
					ok = arcDebit[i]._debit + debit <= CAP + EPSILON;

					//3eme morceau
					arcDebit.push_back(debitTemporise(fin, arcDebit[i]._fin, arcDebit[i]._debit));

					//et on a fini : on a intergalement ajoute le nouvel intervalle
					stop = true;
				}
				else//l intervalle[t, fin] fini apres l intervalle i = > on va creer 1 autre morceau
				{
					// 2sd morceau
					arcDebit.push_back(debitTemporise(t, arcDebit[i]._fin, arcDebit[i]._debit + debit));
					ok = arcDebit[i]._debit + debit <= CAP + EPSILON;

					t = arcDebit[i]._fin;//on avance le nouveau debut courant
					i++;
				}
			}
		}
		//--------------------------------------------------------------------------------
		// 2sd cas t est avant intervalle i et apres intervalle i-1
		else
		{
			//debut du prochain intervalle deja dans la liste s'il existe
			double debutSuivant = TIME_INFINITY;
			if (i < sizeInit)
				debutSuivant = arcDebit[i]._deb;

			if (fin < debutSuivant) //on n intersecte pas d'intervalle
			{
				arcDebit.push_back(debitTemporise(t, fin, debit));
				stop = true;
			}
			else//on intersecte l'intervalle i+1 => on cree le l'intervalle qui s'arrete avant i+1 
			{
				arcDebit.push_back(debitTemporise(t, arcDebit[i]._deb, debit));
				t = arcDebit[i]._deb;

			}

		}


	}


	//====================================================================================
	//on fait le menage et on retrie en fonction des dates de debut

	for (int k = 0; k < intervalleAsupp.size(); ++k)
	{
		//attention au fur et a mesure qu'on enleve on decale les indices d ou le -k
		arcDebit.erase(arcDebit.begin() + intervalleAsupp[k] - k);
	}

	sort(arcDebit.begin(), arcDebit.end());

	return ok;
}

//renvoie vrai si [deb, fin] intersecte l'intervalle DT
bool RCPSP_Algo::intersection(debitTemporise & DT, double deb, double fin)
{
	double A = DT._deb, B = DT._fin;

	return (fin > A && fin < B) || (B > deb && B < fin);

}


vector<int> RCPSP_Algo::genererOrdreInit2()
{
	vector<int> ordre(_graph->_nbSommet + 2);

	//on met sigma dans l'ordre des plus petits chemins
	vector < pair<double, int> > ordre2;

	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		double l = _ins->getLongueurCh(i);
		ordre2.push_back({ l,i });
	}

	sort(ordre2.begin(), ordre2.end());

	for (int i = 0; i < _graph->_nbSommet; ++i)
	{
		ordre[i + 1] = ordre2[i].second;
	}
	ordre[0] = 0;//source
	ordre[_graph->_nbSommet + 1] = _graph->_nbSommet + 1;

	return ordre;
}

vector<int> RCPSP_Algo::genererOrdreInit()
{
	vector<int> ordre(_graph->_nbSommet + 2);

	vector<double> val(_graph->_nbSommet + 2);
	vector<pair<double, int>> v(_graph->_nbSommet);

	//1. on calcule la priorité (val) des actions qui correspond +/- a la marge entre la fin max pour ce chemin et une estimation de la fin au plus tot
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		val[i] = _ins->_maxRateEvac[i];

		int prec = _ins->_chemins[i][0]._som;

		for (int k = 1; k < _ins->_chemins[i].size(); ++k)
		{
			int cour = _ins->_chemins[i][k]._som;
			val[i] = min(val[i], static_cast<double>(_ins->_capArc[prec][cour]));

			prec = cour;
		}

		val[i] = _ins->getFinMax(i) - _ins->getLongueurCh(i) - static_cast<double>(_ins->_pop[i]) / 2 / val[i];
		v[i - 1] = { val[i], i };
	}


	//2. on tri par priorité croissante
	sort(v.begin(), v.end());

	//3. on remplit ordre
	ordre[0] = 0;//source
	ordre[ordre.size() - 1] = _graph->_nbSommet + 1;//puits

	for (int i = 0; i < _graph->_nbSommet; ++i)
	{
		ordre[i + 1] = v[i].second;
	}





	return ordre;
}

//dans cette version on utilise le rang des jobs (qui vient de la BS) pour generer un ordre : 
//les rangs les plus faibles en premier
vector<int> RCPSP_Algo::genererOrdreInit(vector<int> & rang)
{
	int r = 0;

	vector<int> ordre;
	vector<int> sousordre;
	ordre.reserve(_graph->_nbSommet + 2);

	ordre.push_back(0); //source


	while (r < _graph->_nbSommet)
	{
		//on ajoute tous les jobs de rang r de maniere random
		for (int i = 1; i <= _graph->_nbSommet; ++i)
			if (rang[i] == r)
				sousordre.push_back(i);

		if (!sousordre.empty())
		{
			std::random_device rd;
			std::mt19937 g(rd());
			shuffle(sousordre.begin(), sousordre.end(),g);
			ordre.insert(ordre.end(), sousordre.begin(), sousordre.end());
			sousordre.clear();
		}
		r++;
	}

	ordre.push_back(_graph->_nbSommet + 1);//puits

#ifdef _VERIF_
	if (ordre.size() != _graph->_nbSommet + 2)
		stopProg("RCPSP_Algo::genererOrdreInit(vector<int> & rang) : ordre incomplet");
#endif

	return ordre;
}

//renvoie la plus petite marge et remplit les vecteurs : 
//xMin avec tous les jobs qui donnent cette marge la plus faible
//xMin2 avec tous les jobs qui donnent la deuxième plus petite marge

double RCPSP_Algo::calculCoutSolMin(vector<int> & xMin, vector<int> & xMin2)
{
	double margeMin = TIME_INFINITY;
	double margeMin2 = TIME_INFINITY;
	xMin.clear(); xMin2.clear();

	//1. on cherche le min et le deuxieme min
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{

		//on calcule les marges avec la "vraie" fin max (qui vient des deadlines)
		//remarque : attention si on utilise_ins->getFinMax(i) car on a modifie les LF pour
		//la construction de la solution initiale
		double tmp = _finMaxInit[i] - dernierSafe(i);

		if (tmp < margeMin - EPSILON)
		{
			margeMin2 = margeMin;
			margeMin = tmp;
		}
		else
		{
			if (abs(tmp - margeMin) > EPSILON && tmp < margeMin2)
			{
				margeMin2 = tmp;
			}
		}
	}

	//2. on remplit les vecteurs
	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{


		//on calcule les marges avec la "vraie" fin max (qui vient des deadlines)
		//remarque : attention si on utilise_ins->getFinMax(i) car on a modifie les LF pour
		//la construction de la solution initiale
		double tmp = _finMaxInit[i] - dernierSafe(i);

		if (abs(tmp - margeMin) < EPSILON)
		{
			xMin.push_back(i);
		}
		else
		{
			if (abs(tmp - margeMin2) < EPSILON)
				xMin2.push_back(i);
		}
	}

	if (margeMin < EPSILON)
		margeMin = 0;

	return margeMin;
}



//renvoie la somme ponderee par la population des differences entre les deadlines et les derniers individus
double RCPSP_Algo::calculCoutSolSomme()
{
	double sum = 0;

	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		sum += /*_ins->_pop[i] * */ (_ins->getFinMax(i) - dernierSafe(i));
	}


	return sum;
}

//renvoie la somme ponderee par la population des dates de sorties des individus
double RCPSP_Algo::calculCoutSolSommeMilieu()
{
	double sum = 0;

	for (int i = 1; i <= _graph->_nbSommet; ++i)
	{
		//sum += _ins->_pop[i] * (ceil(dernierSafe(i)) - _ins->_pop[i] / 2 / _debit[i]);
		sum += _ins->_pop[i] * (ceil(dernierSafe(i)) + _dateDebut[i] + _ins->getLongueurCh(i)) / 2;

		//cout << _ins->_pop[i] << " * ( " << ceil(dernierSafe(i)) << " + " << _dateDebut[i] << " + " << _ins->getLongueurCh(i) << ") / 2" << endl;
	}


	return sum;
}

//renvoie la date a laquelle le dernier individus partant du noeud evac i arrive au safe node
double RCPSP_Algo::dernierSafe(int i)
{
	//sur un arc o->d,
	//si le premier part de o a la date t, le dernier arrive en d en t + distance(o,d) + (pop/debit) - 1
	return _dateDebut[i] + _ins->getLongueurCh(i) + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1;
}

//renvoie la date a laquelle le dernier individus de la popuplation i 
//arrive sur un noeud sachant la date d arrivee du premier
double RCPSP_Algo::dernier(int i, double arriveePremier)
{
	// le dernier arrive  (pop/debit) - 1 apres le premeir
	return arriveePremier + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1;
}


//on calcule le taux de parallélisme sur le dernier arc pour la best solution
//renvoie le taux de parallelisme moyen sur le dernier arc et le taux max
pair<double, int> RCPSP_Algo::statParallelisme()
{

	int N = _graph->_nbSommet;


	//1. on trie les dates par ordre croissant

	vector<double> date;
	date.reserve(2 * N);

	for (int i = 1; i <= N; ++i)
	{
		date.push_back(_dateDebut[i]);
		date.push_back(_dateDebut[i] + static_cast<double>(_ins->_pop[i]) / _debit[i]);

		//cout << _bestDebut[i] << " ====" << _bestDebut[i] + static_cast<double>(_ins->_pop[i]) / _bestDebit[i] << endl;
	}


	//on trie on ordre croissant et on enleve les doublons
	sort(date.begin(), date.end());
	auto it = unique(date.begin(), date.end());
	date.resize(distance(date.begin(), it));


	for (unsigned int k = 0; k < date.size() - 1; ++k)
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

	vector<bool> is_deb_fait(N + 1, false);
	vector<bool> is_fin_fait(N + 1, false);

	for (int i = 1; i <= N; ++i)
	{
		for (int k = 0; k < K; ++k)
		{
			if (!is_deb_fait[i] && abs(_dateDebut[i] - date[k]) < 5 * EPSILON_DATE) //attention : mettre un peu plus que EPSILON_DATE sinon on risque de rater des jobs
			{
				jobsDebutent[k].push_back(i);
				is_deb_fait[i] = true;
			}

			if (!is_fin_fait[i] && abs(_dateDebut[i] + static_cast<double>(_ins->_pop[i]) / _debit[i] - date[k]) < 5 * EPSILON_DATE)
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


	for (int k = 0; k < K - 1; ++k)//RQ : clique[K-1].size = 0, aucun job ne peut commencer a la fin
	{
		paralMax = max(paralMax, static_cast<int>(clique[k].size()));
		paralMoyen += clique[k].size() * (date[k + 1] - date[k]);

	}

	paralMoyen /= (date[K - 1] - date[0]);

	return { paralMoyen, paralMax };

}