#include "BSmarge.h"
#include <algorithm>
#include <algorithm>
#include "utilVector.h"

//on initialise les donnees d'entree pour calculMargeArc (jobs concernes par l'arc etc...)
//si on souhaite utiliser des donnees ES / LF / debitMax differente de celles de linstance 
//on passe un pointeur dataIns non nul
void BSmarge::initData(int orig, int dest, dataIns * data)
{
	_jobs.clear();
	_ES.clear();
	_LF.clear();
	_debitMax.clear();

	_nbJob = 0;

	_cap = _ins->_capArc[orig][dest];

	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{


		int prec = _ins->_chemins[i][0]._som;

		int k = 1;

		while (k < _ins->_chemins[i].size())
		{
			int cour = _ins->_chemins[i][k]._som;

			if (prec == orig && cour == dest)
			{
				_jobs.push_back(i);
				if (data)//data contient des ES / LF et debit max modifies, si on passe un pt non nul alos on les utilise a la place de ceux de l'instance
				{
					_ES.push_back(data->_ES[i]);
					_LF.push_back(data->_LF[i]);
					_debitMax.push_back(data->_debitMax[i]);
				}
				else
				{
					_ES.push_back(_ins->_chemins[i][k - 1]._ES);//date au plus tot d'arrivee a la porte
					_LF.push_back(_ins->_chemins[i][k - 1]._LF + 1);//date au plus tard de sortie de la porte (+1 car on va perdre une unite de temps en frachissant la porte imaginaire)
					_debitMax.push_back(_ins->getDebitMax(i));
				}

				_nbJob++;
			}

			prec = cour;
			++k;
		}
	}
}


//on initialise (renumerote) les arc d'entree pour calculMargeArcAmeliorer 
void BSmarge::initArc(int orig, int dest)
{
	_toArcs.assign(_ins->_nbNode + 1, vector<int>(_ins->_nbNode + 1, -1));
	_arc2jobs.resize(_ins->_nbNode);  //on ne connait pas le nb d'arc mais on est sur que c'est < au nb de noeuds
	_filsArc.resize(_ins->_nbNode);


	//1. init avec les arcs qui arrivent aux feuilles
	int numArc = 0;
	int h = 0; //hauteur max de l'arbre = longueur du plus long chemin
	for (int j = 0; j < _nbJob; ++j)
	{
		int i = _jobs[j];

		int prec = _ins->_chemins[i][0]._som;
		int cour = _ins->_chemins[i][1]._som;

		_toArcs[prec][cour] = numArc;
		
		_arc2jobs[numArc].push_back(j);
		_arcs.push_back({ prec, cour });
		_capaArc.push_back(_ins->_capArc[prec][cour]);

		numArc++ ;

		if (_ins->_chemins[i].size() > h)
			h = static_cast<int>(_ins->_chemins[i].size());
	}
	h--;


	//2. on avance dans l'arbre jusqu'à la racine en numerotant d'abord les arcs de hauteur h, puis h-1, h-2 etc...
	
	while (h >= 2) //pour h = 2 on va mettre l'arc racine, on pourrait arreter la boucle a h = 3 et mettre l'arc racine a la fin...
	{
		for (int j = 0; j < _nbJob; ++j)
		{
			int i = _jobs[j];

			int k = static_cast<int>(_ins->_chemins[i].size()) - h;//on regarde si le chemin arrive au moins    a la hauteur souhaitee
			if (k >= 1)
			{
				int prec = _ins->_chemins[i][k]._som;
				int cour = _ins->_chemins[i][k + 1]._som;
				int indArc = _toArcs[prec][cour];
				if (indArc == -1)
				{
					indArc = numArc; //pour remplir fils apres la boucle
					_toArcs[prec][cour] = numArc;
					_arcs.push_back({ prec, cour });
					_capaArc.push_back(_ins->_capArc[prec][cour]);

					numArc++;
				}
				_filsArc[indArc].push_back(_toArcs[_ins->_chemins[i][k - 1]._som][prec]);
				_arc2jobs[indArc].push_back(j);
			}
		}
		h--;
		
	}

	//on "nettoie" les fils en triant et en supprimant les doublons
	for (int a = 0; a < _filsArc.size(); ++a)
	{
		sort(_filsArc[a].begin(), _filsArc[a].end());
		supprimeDoublon(_filsArc[a]);

	}

}

void BSmarge::init()
{
	_tCour = TIME_INFINITY;
	_popRestante.resize(_nbJob);
	_fin.resize(_nbJob);
	_debut.assign(_nbJob, -1);
	_actif.resize(_nbJob);
	_debit.resize(_nbJob);
	_marge.assign(_nbJob, TIME_INFINITY);

	for (int k = 0; k < _nbJob; ++k)
	{
		int i = _jobs[k];

		_popRestante[k] = _ins->_pop[i];
		_actif[k] = true;
		_tCour = min(_tCour, _ES[k]);

		//il faut estimer le nb d'iterations max (on suppose < TMAX)
		_debit[k].resize(_ins->_TMAX);
	}

	_listeDate.push_back(_tCour);	
}
void BSmarge::cleanData()
{
	_arcs.clear();
	_capaArc.clear();
	_toArcs.clear();
	_arc2jobs.clear();
	_filsArc.clear();

	_popRestante.clear();
	_actif.clear();
	_listeDate.clear();
	_marge.clear();
	_debit.clear();
	_fin.clear();

}

//appelle calculMargeArc pour l'arc qui va vers le sommet "safe node"
//on suppose qu'il y en a qu'un seul (le même pour toutes les populations)
double BSmarge::calculMargeArc( )
{

	vector<SomWithWindow> & v = _ins->_chemins[1];

	return calculMargeArc(v[v.size() - 2]._som, v[v.size() - 1]._som);

}

//appelle calculMargeArcAmelioree pour l'arc qui va vers le sommet "safe node"
//on suppose qu'il y en a qu'un seul (le même pour toutes les populations)
pair<double, bool> BSmarge::calculMargeArcAmelioree(vector<double> & finOptimiste,  vector<int> & rang)
{
	double res = TIME_INFINITY;
	bool isTransfomable = true;

	vector<vector<bool>> avant(_ins->_lastEvacuationNode + 1, vector<bool>(_ins->_lastEvacuationNode + 1, false));
	rang.assign(_ins->_lastEvacuationNode + 1, 0);
	finOptimiste.resize(_ins->_lastEvacuationNode + 1);
	finOptimiste[0] = 0;
	vector<bool> test(_ins->_nbNode, false);//test(i) = vrai si la BS a ete appelee sur i -> safeNode
	

	//attention on suppose ici qu'il n'y a qu'un noeud d'evacuation
#ifdef _VERIF_
	if (_ins->_lastSafeNode != _ins->_firstSafeNode)
		stopProg("BSmarge::calculMargeArcAmelioree : il faut appeler la BS sur tous les arcs menant a un safe node");
#endif

	int safeNode = _ins->_lastSafeNode;


	//pour chaque arc qui va vers le noeud safe on appelle calculMargeArcAmelioree sur cet arc
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		vector<SomWithWindow> & v = _ins->_chemins[i];
		int precSafe = v[v.size() - 2]._som;

		if (!test[precSafe]) //on verifie qu'on n'a pas deja appele la methode principale sur cet arc
		{
			test[precSafe] = true;

			//===================== methode principale ==========================
			res = min(res,calculMargeArcAmelioree(precSafe, safeNode));
			//cout << "margeAmelioree : " << res << endl;
			isTransfomable = isTransfomable && isRealisable();

			//attention dans fin on a la date d'arrivee du dernier au début de l'arc, finOptimiste donne la fin a la destination de l'arc
			//attention on enleve 1 car on a ajoute 1 dans la borne sup pour franchir la "porte imaginaire"
			int longueur = _ins->_time[precSafe][safeNode] - 1;

			// ====================post traitement============================
			// on calcule des donnees qui nous servirons dans RCPSP_Algo pour construire une "vraie" solution
			for (int i = 0; i < _nbJob; ++i)
			{
				finOptimiste[_jobs[i]] = _fin[i] + longueur; //le franchissement de la porte imaginaire de duree 1 est deja compte dans longueur

				for (int j = i+1; j < _nbJob; ++j)
				{
					if (_debut[i] < _debut[j] && _fin[i] < _fin[j])
					{
						avant[_jobs[i]][_jobs[j]] = true;
						rang[_jobs[j]]++;
					}
					else
					{
						if (_debut[j] < _debut[i] && _fin[j] < _fin[i])
						{
							avant[_jobs[j]][_jobs[i]] = true;
							rang[_jobs[i]]++;
						}
					}
				}
			}

			//on clear les vectors pour la prochaine it
			cleanData();
		}
	}


	return { res, isTransfomable };
}

//calcul les marges optimiste (i.e. BS sur les marges) pour l'arc orig -> dest
double BSmarge::calculMargeArc(int orig, int dest)
{
	double res = -1;

	//on initialise les donnees d'entree (qui viennent de l'instance)   
	initData(orig, dest);

	//et les variables de travail :
	init();


	bool echec = false;
	int n = 0;//n = numero d'iteration
	vector<int> jobPossible;
	jobPossible.reserve(_nbJob);
	int cpt = 0;//compgte le nb de jobs finis

	while (!echec && cpt < _nbJob)
	{
		//=======================================================================================
		//1ere etape : on regarde quel jobs peuvent passer sur l'arc cible et avec quelle marge
		jobPossible.clear();
		_marge.assign(_nbJob, TIME_INFINITY);

		for (int j = 0; !echec && j < _nbJob; ++j)
		{
			//si le job peuvent passer sur l'arc cible, on calcule sa date
			if (_actif[j] && _ES[j] <= _tCour)
			{
				jobPossible.push_back(j);

				_marge[j] = _LF[j] - (_tCour + _popRestante[j] / _debitMax[j]);
				if (_marge[j] < -EPSILON)
					echec = true;
			}
			else//j ne peux pas encore commencer ou est fini => son debit est donc 0 pour cette iteration
				_debit[j][n] = 0;
		}

		if (!echec)
		{

			auto margeOrdonnees = ordonneMarge();
			auto listeJobOrdonnee = classeJobParMarge(margeOrdonnees);

			//============================================================
			//etape 2 : on fait passerle maximum le debit pour chaque job, avec priorite aux petites marges
			int nbMarge = static_cast<int>(margeOrdonnees.size());

			vector<double> alpha(nbMarge);

			bool stop = false;
			double capaRestante = _cap;

			int k = 0;//indice courant dans margeOrdonnees

			while (!stop && k < margeOrdonnees.size())
			{
				
					//retourne la somme des debits max pour les jobs dans listeJobOrdonnee[k]
					double beta = sommeDebitMax(listeJobOrdonnee[k]);

					//alpha = coefficient de reduction du debit pour les jobs dans listeJobOrdonnee[k]
					//de sorte de faire passer au plus vite tous les jobs quand la capa de l'arc est trop faible
					//pour que les jobs passent en même temps avec leur debit max
					alpha[k] = min(1.0, capaRestante / beta);

					//chaque job dans listeJobOrdonnee[k] recoit son debitMax * alpha
					for (int i = 0; i < listeJobOrdonnee[k].size(); ++i)
					{
						int j = listeJobOrdonnee[k][i];

						_debit[j][n] = alpha[k] * _debitMax[j];
					}

					if (beta >= capaRestante)//capa de l'arc completement utilisee a la date courante
					{
						stop = true;
						for (int l = k + 1; l < nbMarge; ++l)
						{
							alpha[l] = 0;
							fixeDebitNul(listeJobOrdonnee[l], n);
						}

					}
					else//encore de la place sur l'arc
					{
						capaRestante -= beta;
						k++;
					}
				
			}//fin while de l etape 2


			//========================================================================
			//etape 3 : on met à jour tCour = prochaine date "interessante" : fin d'un job ou
			// marge d'un job qui devient aussi petite qu'un autre et donc le job change de classe...

			//3.1 prochaine date où un job non actif devient actif :
			double t1 = TIME_INFINITY;
			for (int j = 0; j < _nbJob; ++j)
			{
				if (_actif[j] && _tCour + EPSILON < _ES[j])
					t1 = min(t1, _ES[j]);
			}


			//3.2 fin jobs en cours : 
			double t2 = finJobEnCours(n);

			//3.3 prochaine date ou la marge d'un job non en cours va devenir aussi petite que celle d'un job en cours
			double t3 = rattrapageMarge(alpha, margeOrdonnees);

			//c'est la plus petite des dates qui nous intéresse
			double t0 = min(t1, t2);
			t0 = min(t0, t3);

			//=====================================================================================
			//etape 4 : mise a jour des populations restantes... pour repartir sur une nouvelle iteration
			for (int j : jobPossible)
			{
				_popRestante[j] -= (t0-_tCour)*_debit[j][n];
				if (_popRestante[j] < EPSILON)
				{
					_actif[j] = false;
					_fin[j] = t0; 
					cpt++;
				}
				if (t0 - EPSILON > _LF[j])
					echec = true;
			}

			//==========================================================================
			//etape 5 : mise a jour donnes de l'iteration courante
			_tCour = t0;
			n++;
			_listeDate.push_back(_tCour);

		}//fin if !echec

	}//fin while
	
	if (!echec)
		res = calculMargeMin();
	
	return res;
}


//calcul les marges optimiste (i.e. BS sur les marges) en considerant le sous arbre d'evacuation 
// qui a sa racine en dest (version amelioree de calculMargeArc)
//ATTENTION cet algo donne souvent des marges optimistes (i.e. meilleure que la solution optimale)
//mais il peut se tromper et donner des marges pessimistes (moins bonnes que la solution optimale)
//cet algo n'est donc ni une borne sup, ni une borne inf mais il donne une "assez bonne" idée des marges optimales
double BSmarge::calculMargeArcAmelioree(int orig, int dest, dataIns * data)
{
	double res = -1;

	//on initialise les donnees d'entree (qui viennent de l'instance)   
	initData(orig, dest, data);

	//les arcs
	initArc(orig, dest);

	//et les variables de travail :
	init();

	int nbArc = static_cast<int>(_arcs.size());
	int arcRacine = nbArc - 1;

	bool echec = false;
	int n = 0;//n = numero d'iteration
	vector<int> jobPossible;
	jobPossible.reserve(_nbJob);
	int cpt = 0;//compgte le nb de jobs finis

	vector<double> capaRestante(nbArc);
	vector<double> capaInit(nbArc);

	for (int e = 0; e < nbArc; ++e)
	{
		pair<int, int> p = _arcs[e];
		capaInit[e] = _ins->_capArc[p.first][p.second];
	}


	while (!echec && cpt < _nbJob)
	{
		//=======================================================================================
		//1ere etape : on regarde quel jobs peuvent passer sur l'arc cible et avec quelle marge

		jobPossible.clear();
		_marge.assign(_nbJob, TIME_INFINITY);

		for (int j = 0; !echec && j < _nbJob; ++j)
		{
			//si le job peuvent passer sur l'arc cible, on calcule sa date
			if (_actif[j] && _ES[j] <= _tCour)
			{
				jobPossible.push_back(j);

				_marge[j] = _LF[j] - (_tCour + _popRestante[j] / _debitMax[j]);
				if (_marge[j] < -EPSILON)
					echec = true;
			}
			else//j ne peux pas encore commencer ou est fini => son debit est donc 0 pour cette iteration
				_debit[j][n] = 0;
		}

		if (!echec)
		{

			auto margeOrdonnees = ordonneMarge();
			auto listeJobOrdonnee = classeJobParMarge(margeOrdonnees);

			//============================================================
			//etape 2 : on fait passerle maximum le debit pour chaque job, avec priorite aux petites marges
			
			//2.1 init capa des arcs
			capaRestante = capaInit;

	
			int k = 0;
			vector<double> alpha(_nbJob);
			bool stop = false;
			//vector<int> jobTraiteOUenCours;//jobs dont le alpha a deja ete calcule ou est en cours de calcul
			
			while (k < margeOrdonnees.size() && !stop)
			{
				//init les jobs de marge = margeOrdonnees[k] avec alpha = 1
				for (int j : listeJobOrdonnee[k])
				{
					alpha[j] = 1;
					//jobTraiteOUenCours.push_back(j);
				}
				
				//----------------------
				//2.2 cherhce le coefficient reducteur du debit pour le groupe de jobs k

				//parcourt tous les arcs sauf les arcs dircetement connecte aux noeuds d 'evac
				//arc num _nbJob = premier arc de transit
				for (int e = _nbJob; e < nbArc; ++e)
				{
					//l'arc e nous interesse si des jobs dans listeJobOrdonnee[k] utilisent e : 
					//i.e. si l'intesection de listeJobOrdonnee[k] et _arc2jobs[e] est non vide

					vector<int> jobCourant = intersection(_arc2jobs[e], listeJobOrdonnee[k]);
					if ( !jobCourant.empty() )
					{
						//tri les valeurs de alpha du plus petit au plus grand en ajoutant 0 et en supprimant les doublons
						vector<double> alphaTri = triRestrictionAlpha(alpha, jobCourant);
						vector<vector<int> > jobOrdonneAlpha = ordonneJobSuivantAlpha(alpha, alphaTri, jobCourant);

						//indMaxAlpha = plus grand indice dans alpha tel que tout le monde passe sur l'arc e avec son alpha
						//indMaxAlpha = q* alain
						int indMaxAlpha = calculIndiceMaxAlpha(alphaTri, jobOrdonneAlpha, capaRestante[e]);//indMaxAlpha=q* dans doc Alain

						if (indMaxAlpha+1 < alphaTri.size())
						{
							double beta = calculBetaMargeAmelioree(indMaxAlpha, alphaTri, jobOrdonneAlpha, capaRestante[e]);

							for (int q = indMaxAlpha + 1; q < alphaTri.size(); ++q)
								for (int j : jobOrdonneAlpha[q])
									alpha[j] = min(beta, alpha[j]);//min(beta, 1.0);
						}

					}
				}//fin for e

				//----------------------
				//2.3 fixe les debits pour les jobs du groupe k
				for (int j : listeJobOrdonnee[k])
					_debit[j][n] = alpha[j] * _debitMax[j];

				//----------------------
				//2.4 recalcule la capacite qu'il reste sur l'arc
				

				for (int e = _nbJob; e < nbArc; ++e)
				{
					vector<int> jobCourant = intersection(_arc2jobs[e], listeJobOrdonnee[k]);

					double cap1 = 0, cap2 = 0;

					for (int f : _filsArc[e])
						cap1 += capaRestante[f];

					for (int j : jobCourant)
						cap2 -= alpha[j] * _debitMax[j];
					cap2 += capaRestante[e];

					capaRestante[e] = min(cap1, cap2);

#ifdef _VERIF_
					if (capaRestante[e] < -EPSILON)
						stopProg("BSmarge::calculMargeArcAmelioree : la capa devient negative :(");
#endif
				}

				//----------------------
				//2.5 critere d'arret : capa devenu nul sur l'arc racine
				if (capaRestante[arcRacine] < EPSILON)
				{
					stop = true;
					for (int l = k + 1; l < listeJobOrdonnee.size(); ++l)
						for (int j : listeJobOrdonnee[l])
							_debit[j][n] = 0;
				}
				else
					k++;//sinon on continue avec le groupe suivant


			}//fin while principal etape 2
			   

			//========================================================================
			//etape 3 : on met à jour tCour = prochaine date "interessante" : fin d'un job ou
			// marge d'un job qui devient aussi petite qu'un autre et donc le job change de classe...

			//3.1 prochaine date où un job non actif devient actif :
			double t1 = TIME_INFINITY;
			for (int j = 0; j < _nbJob; ++j)
			{
				if (_actif[j] && _tCour + EPSILON < _ES[j])
					t1 = min(t1, _ES[j]);
			}


			//3.2 fin jobs en cours : 
			double t2 = finJobEnCours(n);

			//3.3 prochaine date ou la marge d'un job non en cours va devenir aussi petite que celle d'un job en cours
			double t3 = rattrapageMargeAmelioree(jobPossible, alpha );

			//c'est la plus petite des dates qui nous intéresse
			double t0 = min(t1, t2);
			t0 = min(t0, t3);

#ifdef _VERIF_
			if (t0 < _tCour)
				stopProg("BSmarge::calculMargeArcAmelioree : t0 < _tCour");
#endif

			//=====================================================================================
			//etape 4 : mise a jour des populations restantes... pour repartir sur une nouvelle iteration
			for (int j : jobPossible)
			{

#ifdef _VERIF_
				if (_debit[j][n] < -EPSILON)
					stopProg("BSmarge::calculMargeArcAmelioree : _debit[j][n] < -EPSILON");
#endif
				if (_debit[j][n] > EPSILON && _debut[j] < -0.5) //_debut[j]=-1 si j non debute 
					_debut[j] = _tCour;

				_popRestante[j] -= (t0 - _tCour)*_debit[j][n];
				if (_popRestante[j] < EPSILON)
				{
					_actif[j] = false;
					_fin[j] = t0 ; 
					cpt++;
				}
				if (t0 - EPSILON > _LF[j])
					echec = true;
			}

			//==========================================================================
			//etape 5 : mise a jour donnes de l'iteration courante
			_tCour = t0;
			n++;
			_listeDate.push_back(_tCour);

		}//fin if !echec

	}//fin while

	if (!echec)
	{
		res = calculMargeMin();
#ifdef _VERIF_
		if (!verifieSolution())
			stopProg("BSmarge::calculMargeArcAmelioree");
#endif
	}

	return res;
}

bool BSmarge::verifieSolution() const
{
	bool ok = true;
	vector<double> popRestante;

	for (int j : _jobs)
		popRestante.push_back(_ins->_pop[j]);

	int nbDate = static_cast<int> (_listeDate.size());

	for (int i = 0; i < nbDate-1; ++i)
	{
		double somDebit = 0;
		double dateCour = _listeDate[i];
		double dateSuiv = _listeDate[i+1];

		for (int j = 0; j < _nbJob; ++j)
		{
			somDebit += _debit[j][i];

			//on verifie que le job ne commence pas avant sa release date et ne termine pas apres sa due date
			if (_debit[j][i] > EPSILON && (dateCour +EPSILON < _ES[j] || dateSuiv - EPSILON > _LF[j]) )
				ok = false;

			//on calcul la population restante
			if (_debit[j][i] > EPSILON)
				popRestante[j] -= (dateSuiv - dateCour) * _debit[j][i];

			//on verifie qu'on ne depasse pas le debit max
			if (_debit[j][i] > _debitMax[j] + EPSILON)
				ok = false;
		}

		//on verifie qu'on respecte la capa de l'arc racine
		if (somDebit > _cap + EPSILON)
			ok = false;
	}

	//on verifie que toutes les populations sont a 0
	for (int j = 0; j < _nbJob; ++j)
		if (popRestante[j] > EPSILON || popRestante[j] < -EPSILON)
			ok = false;

	return ok;
}






//retourne les marges triees de la plus petite a la plus grande et avec les doublons (a epsilon pres) supprimes
vector<double> BSmarge::ordonneMarge()
{
	//on ordonne les marges de la plus petite a la plus grande
	vector<double> margeOrdonnees = _marge;
	sort(margeOrdonnees.begin(), margeOrdonnees.end());

	double prec = margeOrdonnees[0];
	int k = 1;

	while (k < margeOrdonnees.size())
	{
		//attention si on a des marges infinies ca signifie que les jobs ne sont pas autorises a commencer
		if ( abs( margeOrdonnees[k] - prec ) < EPSILON || margeOrdonnees[k] >= TIME_INFINITY-EPSILON)
			margeOrdonnees.erase(margeOrdonnees.begin() + k);
		else
		{
			prec = margeOrdonnees[k];
			k++;
		}

	}

	return margeOrdonnees;
}

//retourne un vecteur v tel que v(i)  donne la liste des jobs avec une marge = margeOrdonnee(i)
vector<vector<int> > BSmarge::classeJobParMarge(vector<double> & margeOrdonnees)
{

	int nbM = static_cast<int>( margeOrdonnees.size());
	vector<vector<int> > v ;
	v.resize(nbM);

	for (int k = 0; k < nbM; ++k)
	{
		for (int j = 0; j < _jobs.size(); ++j)
		{
			if (margeOrdonnees[k] < TIME_INFINITY && abs(_marge[j] - margeOrdonnees[k]) < EPSILON)
			{
				v[k].push_back(j);

			}
		}
	}
	return v;
}

//renvoie la somme des debits max des jobs dans v
double  BSmarge::sommeDebitMax(const vector<int> & v)
{
	double s = 0;

	for (int k = 0; k < v.size(); ++k)
	{
		s += _debitMax[v[k]];

	}

	return s;
}

//fixe le debit a 0 pour tous les jobs dans v et pour l'iteration n
void BSmarge::fixeDebitNul(const vector<int> & v, int n)
{
	for (int j : v)
	{
		_debit[j][n] = 0;
	}
}

//retourne la date du job en cours qui fini le plus tot
double BSmarge::finJobEnCours(int n)
{
	double minFin = TIME_INFINITY;
	for (int j = 0; j < _nbJob; ++j)
	{
		if (_debit[j][n] > EPSILON)
			minFin = min(minFin, _popRestante[j] / _debit[j][n]);
	}


	return _tCour + minFin;
}


//renvoie la date min a laquelle la marge d'un job non en cours devient assez petit pour que le job passe ds les jobs en cours
double BSmarge::rattrapageMarge(const vector<double> & alpha, const vector<double> & margeOrdo)
{
	double t = TIME_INFINITY;

	for (int k = 0; k < margeOrdo.size()-1; ++k)
	{
		t = min(t, (margeOrdo[k + 1] - margeOrdo[k]) / (alpha[k] - alpha[k + 1]) );
	}


	return _tCour + t;
}

//idem rattrapageMarge mais adapte au cas particulier de marge amelioree où tous les alphas d'un meme groupe peuvent etre differents
double BSmarge::rattrapageMargeAmelioree(const vector<int> & jobPossible, const vector<double> & alpha)
{
	double t = TIME_INFINITY;//t = inf delta_k (notation alain)
	

	//on cherche le min (marge(j)-marge(i))/(alpha(i)-alpha(j)) pour tous les (i,j) dans jobPossible tels que 
	// m(i) < m(j) et alpha(i) > alpha(j)

	for (int i : jobPossible)
	{
		for (int j : jobPossible)
		{
			if (_marge[i] < _marge[j] - EPSILON && alpha[i] > alpha[j] + EPSILON)
				t = min(t, (_marge[j]-_marge[i]) / (alpha[i]-alpha[j]) );
		}
	}

#ifdef _VERIF_
	if (t < -EPSILON)
		stopProg("BSmarge::rattrapageMargeAmelioree : on recule");
#endif

	return _tCour + t;

}


//retourne la marge minimum (necissite que _fin et _LF soient initialises)
double BSmarge::calculMargeMin()
{
	double m = TIME_INFINITY;

	for (int i = 0; i < _nbJob; ++i)
	{
		m = min(m, _LF[i] - _fin[i]);
	}
	return m;
}

//renvoie un vecteur contenant les composantes de alpha triees en ordre croissant (en ajoutant 0 et en supprimant les doublons)
// on ne garde que les valeurs de alpha pour les jobs dans job
vector<double> BSmarge::triRestrictionAlpha(const vector<double> & alpha, const vector<int> & job)
{
	vector<double> res;

	for (int j : job)
	{
		res.push_back(alpha[j]);
	}

	//ajoute 0
	res.push_back(0);

	//tri
	sort(res.begin(), res.end());

	//supprime les doublons
	supprimeDoublon(res);

	return res;
}


vector<vector<int> > BSmarge::ordonneJobSuivantAlpha(const vector<double> & alpha, const vector<double> & alphaTri, const vector<int> & job)
{
	vector<vector<int> > res(alphaTri.size());
	
	for (int k = 0; k < alphaTri.size(); ++k)
	{
		double valk = alphaTri[k];

		for (int i = 0; i < job.size(); ++i)
		{
			int j = job[i];
			if (abs(alpha[j] - valk) < EPSILON)
				res[k].push_back(j);

		}
	}


	return res;
}



int BSmarge::calculIndiceMaxAlpha(const vector<double> & alphaTri, const vector<vector<int> > & jobOrdonneAlpha, double cap)
{
	int ind = -1;//q* dans doc alain
	int Q = static_cast<int>(alphaTri.size());

	double sum = 0, sum1 = 0, sum2 = 0;

	do {
		ind++;


		//sum1 contient la somme des sum_j _debitMax(j) * alphaTri(q) pour q <= ind 
		//(on garde sum1 en memoire et on ajoute juste la partie q = ind a chaque tour de boucle)
		for (int j : jobOrdonneAlpha[ind])
			sum1 += _debitMax[j] * alphaTri[ind];


		sum2 = 0;//recalculer a chaque fois, on pourrait aller plus vite en supprimant a chaque tour de boucle la sum correspondant a ind
		for (int q = ind + 1; q < Q; ++q)
		{
			for (int j : jobOrdonneAlpha[q])
				sum2 += _debitMax[j] * alphaTri[ind];
		}

		sum = sum1 + sum2;
	} while (ind < Q-1 && sum <= cap + EPSILON);

	
	if (sum > cap + EPSILON)
		ind--;

#ifdef _VERIF_
	if (ind == -1)
		stopProg("BSmarge::calculIndiceMaxAlpha : impossible de calculer un alpha max, l equation n est jamais verifiee");
#endif //_VERIF_


	return ind;
}

double BSmarge::calculBetaMargeAmelioree(int ind, const vector<double> & alphaTri, const vector<vector<int> > & jobOrdonneAlpha, double cap)
{
	double beta = 0;

	double sum1 = 0, sum2 = 0;

	for (int q = 0; q <= ind; ++q)
	{
		for (int j : jobOrdonneAlpha[q])
			sum1 += _debitMax[j] * alphaTri[q];
	}

	for (int q = ind+1; q < alphaTri.size(); ++q)
	{
		for (int j : jobOrdonneAlpha[q])
			sum2 += _debitMax[j] ;
	}


	//attention normalement cette fonction n'est appelee que si on est certain que sum2 > 0
#ifdef _VERIF_
	if (sum2 < EPSILON)
		stopProg("BSmarge::calculBetaMargeAmelioree : div par 0");
#endif

	beta = (cap - sum1) / sum2;


	return beta;
}

// renvoie vrai si la solution de la BS peut etre transformee en "vraie" solution : 
// on fait la moyenne des debits pour chaque job en ponderent par la duree des intervalles sur le quel elle s execute
//attention : si l'instance est divisee en sous arbre (pas d'arc commum a tous les jobs) alors il faut tester pour chaque sous arbre
bool BSmarge::isRealisable()
{
	bool ok = true;
	vector<double> debitMoyen(_nbJob, 0);

	const int nDate = static_cast<int> ( _listeDate.size() );

	//1. on calcule les debits moyens
	for (int i = 0; i < _nbJob; ++i)
	{
		for (int t = 0; t < nDate-1; ++t)
		{
			if (_debit[i][t] > EPSILON)
				debitMoyen[i] += _debit[i][t] * (_listeDate[t+1]-_listeDate[t]);
		}

		debitMoyen[i] /= (_fin[i] - _debut[i] + 1);
	}

	//2. on verifie que les capa des arcs sont ok avec ces debits

	int nbArc = static_cast<int>(_arcs.size());
	

	for (int e = 0; ok && e < nbArc; ++e)
	{
		pair<int, int> p = _arcs[e];
		int cap = _ins->_capArc[p.first][p.second];

		for (int t = 0; t < nDate - 1; ++t)
		{
			double sumDebit = 0;
			for (int i : _arc2jobs[e])
			{
				if (_debit[i][t] > EPSILON)
					sumDebit += debitMoyen[i];
			}
			if (sumDebit > cap + EPSILON)
				ok = false;
		}
	}



	return ok;
}

//on appelle iterativement l'algo calculMargeArcAmelioree qui calcul une BS,
// a chaque it on fixe le job le plus stresse a son debit moyen et on recommence,
// a la fin on obtient une "vraie" solution (tous les jobs ont un debit constant) 
double BSmarge::construireVraieSolution()
{

	double res = TIME_INFINITY;
	bool isTransfomable = true;
	vector<bool> test(_ins->_nbNode, false);//test(i) = vrai si la BS a ete appelee sur i -> safeNode


	dataIns data;
	allocDataIns(data);


	//attention on suppose ici qu'il n'y a qu'un noeud d'evacuation
#ifdef _VERIF_
	if (_ins->_lastSafeNode != _ins->_firstSafeNode)
		stopProg("BSmarge::calculMargeArcAmelioree : il faut appeler la BS sur tous les arcs menant a un safe node");
#endif

	int safeNode = _ins->_lastSafeNode;


	//pour chaque arc qui va vers le noeud safe on appelle calculMargeArcAmelioree sur cet arc
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		vector<SomWithWindow> & v = _ins->_chemins[i];
		int precSafe = v[v.size() - 2]._som;

		if (!test[precSafe]) //on verifie qu'on n'a pas deja appele la methode principale sur cet arc
		{
			test[precSafe] = true;
			res = min(res, construireVraieSolution(data, precSafe, safeNode));

			//on clear les vectors pour la prochaine it
			cleanData();
		}
	}

	return res;
}


//on appelle iterativement l'algo calculMargeArcAmelioree qui calcul une BS,
// a chaque it on fixe le job le plus stresse a son debit moyen et on recommence,
// a la fin on obtient une "vraie" solution (tous les jobs ont un debit constant) 
double BSmarge::construireVraieSolution(dataIns & data, int o, int d)
{
	double res = -1;
	double debitMoyen = 0;
	

	//===================== phase d'initialisation : on construit une premiere BS et on l'utilise pour remplir data
	res = calculMargeArcAmelioree(o, d);
	
	//on sauvagarde les LF initiales (pour le calcul de la "vraie" marge a la fin)
	vector<double> savLF = _LF;
	vector<bool> fix(_nbJob, false); //fix (i) = vrai si i a deja ete fixe

	for (int i = 0; i < _nbJob; ++i)
	{
		data._ES[_jobs[i]] = _ES[i];
		data._LF[_jobs[i]] = _LF[i];
		data._debitMax[_jobs[i]] = _debitMax[i];
	}		


	while (res > -EPSILON && !isRealisable())
	{


		//on cherche le job qui donne le min et on fixe son debit a son debit moyen, 
		//son ES a sa date de depart actuelle et son LF a sa date de fin actuelle
		const int nDate = static_cast<int> (_listeDate.size());
		double minMarge = TIME_INFINITY;
		int i; //job qu'on va fixer

		for (int k = 0; k < _nbJob; ++k)
		{
			if (!fix[k] && savLF[k] - _fin[k] < minMarge)
			{
				minMarge = savLF[k] - _fin[k];
				i = k;
			}
		}

		fix[i] = true;

		//on calcul son debit moyen
		debitMoyen = 0;
		for (int t = 0; t < nDate - 1; ++t)
		{
			if (_debit[i][t] > EPSILON)
				debitMoyen += _debit[i][t] * (_listeDate[t + 1] - _listeDate[t]);
		}

		debitMoyen /= (_fin[i] - _debut[i]);


		//on met a jour data
		data._ES[_jobs[i]] = _debut[i];
		data._LF[_jobs[i]] = _fin[i];
		data._debitMax[_jobs[i]] = debitMoyen;


		// on relance l'algo
		cleanData(); // 1. reinit sdd
		res = calculMargeArcAmelioree(o, d, &data); //2. algo
	}


	//on recalcule la marge min avec les vraies dates LF
	if (res != -1)
	{
		double minMarge = TIME_INFINITY;
		for (int k = 0; k < _nbJob; ++k)
			minMarge = min(minMarge, savLF[k] - _fin[k]);
	}

	return res;
}

//init data avec 
void BSmarge::allocDataIns(dataIns & data)
{
	data._ES.assign(_ins->_lastEvacuationNode + 1, -1);
	data._LF.assign(_ins->_lastEvacuationNode + 1, -1);
	data._debitMax.assign(_ins->_lastEvacuationNode + 1, -1);
}

//=========================================================================
// calcule de la solution optimale du probleme preemptif


// resout le probleme preemptif : i.e. debit variable le long d'un chemin
// en considerant le sous arbre d'evacuation qui a sa racine en dest 
double BSmarge::calculMargePLdicho()
{
	double margeMin = TIME_INFINITY;
	vector<bool> test(_ins->_nbNode, false);//test(i) = vrai si la BS a ete appelee sur i -> safeNode
	int safeNode = _ins->_lastSafeNode;

	//pour chaque arc qui va vers le noeud safe on appelle calculMargePLdichoArc sur cet arc
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		vector<SomWithWindow> & v = _ins->_chemins[i];
		int precSafe = v[v.size() - 2]._som;

		if (!test[precSafe]) //on verifie qu'on n'a pas deja appele la methode principale sur cet arc
		{
			test[precSafe] = true;


			double m = calculMargePLdichoArc(precSafe, safeNode);
			if (m < margeMin)
				margeMin = m;

			//on clear les vectors pour la prochaine it
			cleanData();
		}
	}

	return margeMin;
}

// resout le probleme preemptif : i.e. debit variable le long d'un chemin
// en considerant le sous arbre d'evacuation qui a sa racine en dest 
double BSmarge::calculMargePLdichoArc(int orig, int dest)
{
	double epsilon = 0.0001;
	double bestMargeCible = -1;

	//init des donnees

	//on initialise les donnees d'entree (qui viennent de l'instance)   
	initData(orig, dest);  

	//les arcs
	initArc(orig, dest);

	//reoslution



	double margeCibleMin = calculMargeCibleMin();
	double margeCibleMax = calculMargeCibleMax();


	
	vector<double> date;

	//boucle principale dicho
	while (abs(margeCibleMax - margeCibleMin) > epsilon)
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

	cleanData();

	return bestMargeCible;
}


//calcul la meilleure marge possible compte tenu des population
double BSmarge::calculMargeCibleMax()
{
	double margeMax = TIME_INFINITY;

	for (int i = 0; i < _nbJob; ++i)
	{
		int job = _jobs[i];

		double m = _LF[i] - (_ES[i] + static_cast<double>(_ins->_pop[job]) / _ins->_maxRateEvac[job]);
		if (margeMax > m)
			margeMax = m;
	}

	return margeMax;
}

//on donne la marge fournit par calculMargeAmelioree (donne une solution non optimale au probleme preemptif)
double BSmarge::calculMargeCibleMin()
{
	return 0;
}

//on ordonne l'ensemble des dates ESi, LFi-margeCible par ordre croissant
//renvoie le vecteur des valeurs triees
vector<double> BSmarge::trierDate(double margeCible)
{
	vector<double> date;
	date.reserve(_ins->_lastEvacuationNode * 2);

	//ajoute le vecteur _ES dans date
	copy(_ES.begin(), _ES.end(), back_inserter(date));
	transform(_LF.begin(), _LF.end(), back_inserter(date), [&](double lf) {return lf - margeCible; });

	sort(date.begin(), date.end());
	auto it = unique(date.begin(), date.end());
	
	date.resize(distance(date.begin(), it));

	return date;
}

//on lance un PL pour savoir si toutes les pop peuvent evacuer avec une marge >= margeCible  
bool BSmarge::isMargeRealisable(double margeCible, vector<double> & date)
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

	IloArray<IloNumVarArray> debit(env, _nbJob); //debit[i][k] = debit du job i entre les dates date[k] et date[k+1]

	for (int i = 0; i < _nbJob; ++i)
	{
		debit[i] = IloNumVarArray(env, nbDate-1);

		for (int k = 0; k < nbDate-1; ++k)
		{
			stringstream ss;
			ss << "d_" << i << "_" << k;
			debit[i][k] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, ss.str().c_str());
		}
	}//fin for i 

	//==========================================================
	//contraintes

	//1. le debit ne depasse pas le debit max
	for (int i = 0; i < _nbJob; ++i)
	{
		for (int k = 0; k < nbDate-1; ++k)
		{
			int job = _jobs[i];
			model.add(debit[i][k] <= _ins->_maxRateEvac[job]);
		}
		
	}//fin for i 


	//2. prise en compte des fenetres de temps [ES,LF-marge]
	for (int i = 0; i < _nbJob; ++i)
	{
		for (int k = 0; k < nbDate-1; ++k)
		{
			if (_ES[i] >= date[k + 1] - EPSILON || _LF[i] - margeCible <= date[k] + EPSILON)
				model.add(debit[i][k] == 0);
		}

	}//fin for i 

	//3. on evacue tout le monde
	for (int i = 0; i < _nbJob; ++i)
	{
		IloExpr expr(env);
		for (int k = 0; k < nbDate-1; ++k)
		{
			expr += debit[i][k] * (date[k+1]-date[k]);
		}

		int job = _jobs[i];

		model.add(expr == _ins->_pop[job]);
		expr.end();
	}//fin for i 


	//4. prise en compte des capacites des arcs
	int nbArc = static_cast<int> (_arcs.size());

	for (int a = 0; a < nbArc; ++a)
	{
		
		int nbJobSurA = static_cast<int> (_arc2jobs[a].size());

		for (int k = 0; k < nbDate-1; ++k)
		{
			
			IloExpr expr(env);

			for (int j = 0; j < nbJobSurA; ++j)
			{
				int i = _arc2jobs[a][j];
				expr += debit[i][k];
			}
			model.add(expr <= _capaArc[a]);
			expr.end();
		}


	}//fin for a

	//===========================================================================
	// resolution


	//_cplex->setParam(IloCplex::Param::Threads, NB_THREAD_CPLEX);
	//_cplex->setParam(IloCplex::Param::ClockType, CLOCK_TYPE_CPLEX);
	//_cplex->setParam(IloCplex::Param::TimeLimit, TIME_LIMIT_CPLEX);
	//_cplex->setParam(IloCplex::Param::MIP::Limits::TreeMemory, RAM_LIMIT_CPLEX);

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