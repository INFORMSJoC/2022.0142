#include "Instance.h"
#include <cmath>
#include <stack>

// display instance data
void Instance::display()
{
	cout << "nbNode = " << _nbNode << endl;
	cout << "TMAX = " << _TMAX << endl;
	cout << "evac nodes = " << _firstEvacuationNode << "..." << _lastEvacuationNode << endl;
	cout << "safe nodes = " << _firstSafeNode << "..." << _lastSafeNode << endl;
	cout << "transit nodes = " << _firstTransitNode << "..." << _lastTransitNode << endl;


	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
		cout << "pop in " << i << " = " << _pop[i] << endl;

	for (int i = 1; i <= _nbNode; ++i)
		cout << "capNode " << i << " = " << _capNode[i] << endl;

	cout << "TIME" << endl;

	for (int i = 1; i <= _nbNode; ++i)
	{
		for (int j = 1; j <= _nbNode; ++j)
		{
			if (_time[i][j] == TIME_INFINITY)
				cout << "inf";
			else
				cout << _time[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	cout << "CAP ARC" << endl;

	for (int i = 1; i <= _nbNode; ++i)
	{
		for (int j = 1; j <= _nbNode; ++j)
		{
			cout << _capArc[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	cout << "DEAD DATE ARC" << endl;

	for (int i = 1; i <= _nbNode; ++i)
	{
		for (int j = 1; j <= _nbNode; ++j)
		{
			cout << _dead[i][j] << " ";
		}
		cout << endl;
	}


}




//generate data randomly
void Instance::generateDataFile(string name, int nbEvacNode, int nbSafeNode, int nbCoucheInterne, int nbNodePerCouche)
{

	if (nbCoucheInterne < 1 || nbNodePerCouche < 1)
		stopProg("InstanceVS::generateDataFile: il faut au mois un noeud interne");

	ofstream fic(name);

	fic << nbEvacNode << " " << nbSafeNode << " " << nbCoucheInterne * nbNodePerCouche << endl;
	for (size_t i = 0; i < nbEvacNode; ++i)
	{
		fic << rand() % 10 + 1 << " " << 50 << endl;
	}

	size_t nbArc = nbEvacNode * nbNodePerCouche + nbSafeNode * nbNodePerCouche +
		nbNodePerCouche * nbNodePerCouche * (nbCoucheInterne - 1);



	fic << nbArc << endl;

	vector< vector<int> > couche(nbCoucheInterne + 2); //couche[i] = num des sommets dans la couche i

	int C = 0;

	//1ere couche = evac node
	for (int i = 1; i <= nbEvacNode; ++i)
		couche[C].emplace_back(i);
	C++;

	//middle couche
	int cpt = static_cast<int> (nbEvacNode + nbSafeNode + 1);
	for (size_t i = 0; i < nbCoucheInterne; ++i)
	{
		for (size_t j = 0; j < nbNodePerCouche; ++j)
		{
			couche[C].emplace_back(cpt);
			cpt++;
		}
		C++;
	}

	//derniere couche = safe node
	for (int i = nbEvacNode + 1; i <= nbEvacNode + nbSafeNode; ++i)
		couche[C].emplace_back(i);


	// pour chaque arc: origine, destination, deadline, longueur, capacite
	for (size_t c = 0; c < nbCoucheInterne + 1; ++c)
	{
		for (size_t j = 0; j < couche[c].size(); ++j)
		{
			for (size_t k = 0; k < couche[c + 1].size(); ++k)
			{
				fic << couche[c][j] << " " << couche[c + 1][k] << " " << (rand() % 11 + 15) * (c + 1)
					<< " " << rand() % 21 + 10 << " " << rand() % 4 + 2 << endl;
			}
		}

	}

}



//generate data randomly
//genere un graphe en arbre
void Instance::generateTree(int nbEvacNode, long seed)
{
	srand(seed);

	//on a deja tous le nodes evac + safe, on va ajouter les transit
	_nbNode = nbEvacNode + 1;


	_firstEvacuationNode = 1;
	_firstSafeNode = nbEvacNode + 1;
	_firstTransitNode = nbEvacNode + 2;
	_lastEvacuationNode = nbEvacNode;
	_lastSafeNode = nbEvacNode + 1;

	int numSom = nbEvacNode + 2;

	//----------------------------------------------
	// on construit un arbre de la racine vers les feuilles
	int nbMaxNode = 4 * nbEvacNode;
	vector<vector<int> > fils(nbMaxNode);//fils[i] = liste des fils de i

	int cour = _firstSafeNode, next = cour + 1;
	bool stop = false;
	int nbFeuille = 1;

	//on cree un premier arc du safe node vers un transit surlequel tout le monde passera
	fils[cour].push_back(next);
	next++;
	cour++;

	//arbre
	while (nbFeuille < nbEvacNode)
	{
		int nbFils = 2;

		if (nbEvacNode - nbFeuille >= 3)//s'il reste peu de feuilles manquantes on ne cree que 2 fils,
			//sinon on tire au hasard entre 2, 3 et 4 le nb de fils avec une proba p(2)> p(3) > p(4)
		{
			double alea = static_cast<double>(rand()) / RAND_MAX;
			if (alea > 0.5)
			{
				if (alea <= 0.8)
					nbFils = 3;
				else
					nbFils = 4;
			}
		}

		for (int i = 0; i < nbFils; ++i)
			fils[cour].push_back(next + i);


		cour++; //sommet courant
		next = next + nbFils; //prochain numero de sommet libre
		nbFeuille += nbFils - 1; //on a cree nbFils feuilles a partir d'une feuille

#ifdef _VERIF_
		if (next > nbMaxNode)
			stopProg("Instance::generateTree : il faut prevoir un vecteur fils plus grand");
#endif // _VERIF_
	}


#ifdef _VERIF_
	if (nbFeuille != nbEvacNode)
	{
		stopProg("Instance::generateTree : nb feuilles > NB EVAC");
	}
#endif

	//----------------------------------------------
	// 2. on construit les chemins a partir de l'arbre
	//a chaque fois qu'on arrive sur une feuille, chemin = contenu de la pile

	cour = _firstSafeNode;
	vector< pair<int, int> > pile; //empile sommet cour et numero du fils surlequel on est descendu
	int numFils = 0;

	//init pile
	pile.push_back({ cour, -1 });
	int last = 0;
	int nextChemin = 1;
	_chemins.resize(nbEvacNode + 1);


	while (!pile.empty())
	{
		cour = pile[last].first;
		numFils = pile[last].second + 1;//prochain fils a explorer
		pile[last].second = numFils; //enregistre dans la pile qu'on a avance dans les fils


		if (fils[cour].size() > numFils)
		{
			//on descend sur le prochain fils, on empile
			cour = fils[cour][numFils];
			pile.push_back({ cour,-1 });
			last++;
		}
		else
		{
			//on est sur une feuille => on remplit le chemin correspondant
			if (fils[cour].empty())
			{
				for (int j = last; j >= 0; --j)
				{
					_chemins[nextChemin].push_back(SomWithWindow(pile[j].first, 0, 0));
				}

				//chemin i commence en i => renumerote juste le premier du chemin
				_chemins[nextChemin][0]._som = nextChemin;

				_nbNode = max(_nbNode, _chemins[nextChemin][1]._som);
				nextChemin++;
			}
			//  on depile
			pile.erase(last + pile.begin());
			last--;

		}

	}


	//maj le nb de noeuds total
	_lastTransitNode = _nbNode;

	//----------------------------------------------
	// 3. allocations memoires
	_pop.resize(nbEvacNode + 1);
	_time.assign(_nbNode + 1, vector<int>(_nbNode + 1, TIME_INFINITY));
	_capArc.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_capNode.assign(_nbNode + 1, 0);
	_dead.assign(_nbNode + 1, vector<int>(_nbNode + 1, TIME_INFINITY));
	_maxRateEvac.resize(nbEvacNode + 1);
	_minDate.assign(_nbNode + 1, TIME_INFINITY);
	_succ.resize(_nbNode + 1);
	_TMAX = 0;

	//----------------------------------------------
	// 4. remplit le reste de l'instance (uniquement ce qui est utile dans le cas d'un arbre)

	int sumPop = 0;

	for (int i = 1; i <= nbEvacNode; ++i)
	{
		_pop[i] = 1;// rand() % 190 + 10;
		_capNode[i] = _pop[i];
		_maxRateEvac[i] = 5 + rand() % 46;
		_maxRateEvac[i] = 100000;//on va le reduire avec les capa min du chemin
		sumPop += _pop[i];
	}
	_capNode[_firstSafeNode] = sumPop;

	for (int i = _firstTransitNode; i <= _lastTransitNode; ++i)
	{
		_capNode[i] = 10;
	}

	//on parcourt les chemins en partant du safe node et on fixe des capa sur les arcs de plus en plus petites
	for (int i = 1; i <= nbEvacNode; ++i)
	{
		int s = static_cast<int>(_chemins[i].size() - 2);
		int cour = _firstSafeNode;
		int capPred = 11 * (s + 1) + 2; //cap suffisament grande pour enlever jusqu'à 5 unités à chaque fois

		int longueurCh = 0;

		for (int k = s; k >= 0; --k)
		{
			int prec = _chemins[i][k]._som;
			if (_capArc[prec][cour] == 0)
			{
				_capArc[prec][cour] = capPred - 5 - rand() % 6; //cap arc precedent entre 5 et 10 en moins
				//en meme temps on init time
				_time[prec][cour] = rand() % 21 + 2;

				_succ[prec].push_back(cour);

#ifdef _VERIF_
				if (_capArc[prec][cour] <= 0)
					stopProg(" Instance::generateTree : cap arc <= 0 ");
#endif
			}

			capPred = _capArc[prec][cour];
			longueurCh += _time[prec][cour];

			cour = prec;
		}


		//on fixe une deadline sur le chemin, les autres seront deduites
		//on suppose qu'on doit arriver au safe node avant longueurCh * r où r = alea entre 
		//on deduit la deadline pour le premier arc du chemin
		double a = static_cast<double>(rand()) / RAND_MAX;
		_dead[i][_chemins[i][1]._som] = static_cast<int>(longueurCh * (2 + a) - (longueurCh - _time[i][_chemins[i][1]._som]));
		_TMAX = max(_TMAX, static_cast<int>(longueurCh * (2 + a)) + 1);
	}


	computeMinDate();

	//dans le cas ou on a des chemins deja definis
	computeWindows();

	reduitMaxRateEvac();



}





void Instance::readFileLN(const string& name)
{
	int popTotal = 0;

	ifstream fic(name);
	if (!fic)
		stopProg("impossible to open " + name );

	string poub;

	int nbEvacNode, nbSafeNode, nbTransitNode, nbArc;

	// get nb nodes
	fic >> nbEvacNode >> nbSafeNode >> nbTransitNode;

	_nbNode = nbEvacNode + nbSafeNode + nbTransitNode;


	// ALLOC (on a une matrice pour les arcs mais le graphe n est pas forcement complet, dans ce cas on garde les valeurs par defaut
	_pop.resize(nbEvacNode + 1);
	_time.assign(_nbNode + 1, vector<int>(_nbNode + 1, TIME_INFINITY));
	_capArc.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_capNode.assign(_nbNode + 1, 0);
	_dead.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_maxRateEvac.resize(nbEvacNode + 1);
	_succ.resize(_nbNode + 1);
	_minDate.assign(_nbNode + 1, TIME_INFINITY);

	//
	_firstEvacuationNode = 1;
	_lastEvacuationNode = nbEvacNode;
	_firstSafeNode = nbEvacNode + 1;
	_lastSafeNode = nbEvacNode + nbSafeNode;
	_firstTransitNode = _lastSafeNode + 1;
	_lastTransitNode = _nbNode;

	// get info about evac node
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		fic >> _pop[i] >> _maxRateEvac[i];
		popTotal += _pop[i];
	}


	// get info about arc
	fic >> nbArc;

	cout << " =====>>> nbArc = " << nbArc << endl;

	int maxTime = 0;
	for (size_t i = 0; i < nbArc; ++i)
	{
		int o = 0, d = 0;
		fic >> o >> d;
		fic >> _dead[o][d] >> _time[o][d] >> _capArc[o][d];
		if (_dead[o][d] > maxTime)
			maxTime = _dead[o][d];

		_succ[o].push_back(d);
	}

	_TMAX = maxTime * 2;

	//chemins 
	_chemins.resize(_lastEvacuationNode + 1);
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int nb, som;
		fic >> nb;
		for (int k = 0; k < nb; ++k)
		{
			fic >> som;
			_chemins[i].push_back(SomWithWindow(som, -1, -1));
		}
	}



	//sur les nodes evac et safe on a une capa = a la pop
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		_capNode[i] = _pop[i];
	}

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		_capNode[i] = _pop[i];
	}

	for (int i = _firstSafeNode; i <= _lastSafeNode; ++i)
	{
		_capNode[i] = popTotal;
	}

	

	computeWindows();

	computeMinDate();
}


//read data from a file made by the software of C. Artigues
void Instance::readFileCArtigues(const string& name)
{

	int  popTotal = 0;
	_TMAX = 0;
	cout << "readFileCArtigues : attention on arrondi les distances ET les capa => on a juste une heuristiaue, OK ? Entrer un chiffre et entrer pour continuer" << endl;
	//int stop;cin >> stop;


	ifstream fic(name);
	if (!fic)
		stopProg("readFileCArtigues: Impossible d ouvrir le fichier de donnees");

	string poub;

	vector<int> CtoLN; //CtoLN[i] donne l indice dans les tableaux de la classe pour le noeud i dans le fichier de chritian
	vector<size_t> numEvac; //indice des node evac dans le fichier de christian


	//1. on recupere toutes les informatioms, les noeuds sont dans le desordre ( non numerotes evac puis safe puis transit)
	// => on les renumerote apres

	//=================================================================
	//1.1 il faut recuperer la taille pour les alloc

	getline(fic, poub);
	int nbEvacNode, idSafeNode, idNode = 0, nbArc;
	int nbSomCh;
	fic >> nbEvacNode >> idSafeNode;


	// on recupere les indices des node evac
	numEvac.resize(nbEvacNode);
	for (size_t i = 0; i < nbEvacNode; ++i)
	{
		fic >> numEvac[i] >> poub >> poub >> nbSomCh;
		for (int i = 0; i < nbSomCh; ++i)
		{
			fic >> poub;
		}
	}

	fic >> poub;

	if (poub != "c")
		stopProg("readFileCArtigues : pb lecture fichier, on ne recupere pas c");
	getline(fic, poub);

	fic >> _nbNode;

	if (_nbNode < 1)
		stopProg("readFileCArtigues : pb lecture fichier, nbNode < 1");

	cout << _nbNode << endl;

	//==============================================
	// allocation des donnees

	CtoLN.assign(_nbNode, 0);

	// ALLOC 1 (on a une matrice pour les arcs mais le graphe n est pas forcement complet, dans ce cas on garde les valeurs par defaut (TIME_INFINITY ou 0)
	_pop.resize(nbEvacNode + 1);
	_time.assign(_nbNode + 1, vector<int>(_nbNode + 1, TIME_INFINITY));
	_capArc.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_capNode.assign(_nbNode + 1, 0);
	_dead.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_maxRateEvac.resize(nbEvacNode + 1);
	_succ.resize(_nbNode + 1);
	_minDate.assign(_nbNode + 1, TIME_INFINITY);
	_chemins.resize(nbEvacNode + 1);


	_firstEvacuationNode = 1;
	_lastEvacuationNode = nbEvacNode;
	_firstSafeNode = nbEvacNode + 1;
	_lastSafeNode = nbEvacNode + 1; //un seul safe node chez christian
	_firstTransitNode = _lastSafeNode + 1;
	_lastTransitNode = _nbNode;

	// on remplit CtoLN
	for (int i = 0; i < _lastEvacuationNode; ++i)
	{
		CtoLN[numEvac[i]] = i + 1;
	}
	CtoLN[idSafeNode] = _firstSafeNode;

	// pour les cases non remplit (= 0) on remplit a partir de _firstSafeNode + 1 jusqu a _nbNode
	int cpt = _firstSafeNode + 1;
	for (int i = 0; i < _nbNode; ++i) // indice de christian vont de 0 a nbNode-1
	{
		if (CtoLN[i] == 0)
		{
			CtoLN[i] = cpt;
			cpt++;
		}
	}


	if (cpt != _nbNode + 1)
		stopProg("readFileCArtigues : manque t il des noeuds ??");

	//======================================================
	// retour au debut du fichier
	fic.clear();
	fic.seekg(0, ios::beg);

	getline(fic, poub);
	fic >> nbEvacNode >> idSafeNode;
	int idC1, idC2, id1, id2;

	// on recupere les indices des node evac
	for (size_t i = 0; i < nbEvacNode; ++i)
	{
		fic >> idNode;
		//cout << idNode << endl;
		fic >> _pop[CtoLN[idNode]] >> _maxRateEvac[CtoLN[idNode]];
		popTotal += _pop[CtoLN[idNode]];

		fic >> nbSomCh;

		//en premier on met le noeud d evac...
		_chemins[CtoLN[idNode]].push_back(SomWithWindow(CtoLN[idNode], -1, -1));

		//ensuite on met le chemin
		for (int i = 0; i < nbSomCh; ++i)
		{
			fic >> idC1;
			_chemins[CtoLN[idNode]].push_back(SomWithWindow(CtoLN[idC1], -1, -1));
		}


	}

	fic >> poub;

	if (poub != "c")
		stopProg("readFileCArtigues : pb lecture fichier, on ne recupere pas c (2eme passage)");
	getline(fic, poub);

	fic >> poub >> nbArc;


	double length, capacity;
	double duedate; //valeur ds fic trop grand pour un int

	for (size_t i = 0; i < nbArc; ++i)
	{
		fic >> idC1 >> idC2;
		id1 = static_cast<int>(CtoLN[idC1]);
		id2 = static_cast<int>(CtoLN[idC2]);
		fic >> duedate >> length >> capacity;

		length = max(length, 1.0);

		_time[id1][id2] = static_cast<int>(floor(length) + EPSILON); //notre modele traite des entiers, on arrondi a l entier inf ==> borne inf
		_capArc[id1][id2] = static_cast<int>(ceil(capacity) + EPSILON);

		_time[id2][id1] = static_cast<int>(floor(length) + EPSILON);
		_capArc[id2][id1] = static_cast<int>(ceil(capacity) + EPSILON);

		_succ[id1].push_back(id2);
		_succ[id2].push_back(id1);

		if (duedate > TIME_INFINITY)
		{
			_dead[id1][id2] = TIME_INFINITY;
			_dead[id2][id1] = TIME_INFINITY;
		}
		else
		{
			_dead[id1][id2] = static_cast<int> (duedate + EPSILON);
			_dead[id2][id1] = static_cast<int> (duedate + EPSILON);


			if (_dead[id1][id2] + _time[id1][id2] > _TMAX && _dead[id1][id2] + _time[id1][id2] < TIME_INFINITY)
				_TMAX = (_dead[id1][id2] + _time[id1][id2]);

			if (_dead[id1][id2] < 0)
			{
				stopProg("readFileCArtigues : une dead est < 0");
			}
		}

	}


	//sur les nodes evac et safe on a une capa = a la pop
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		_capNode[i] = _pop[i];
	}

	for (int i = _firstSafeNode; i <= _lastSafeNode; ++i)
	{
		_capNode[i] = popTotal;
	}


	//==================

	computeMinDate();

	//dans le cas ou on a des chemins deja definis
	computeWindows();

	reduitMaxRateEvac();
}


void Instance::readFileAlicia(const string & name)
{
	int  popTotal = 0;
	_TMAX = 0;
	//int stop;cin >> stop;


	ifstream fic(name);
	if (!fic)
		stopProg("readFileAlicia: Impossible d ouvrir le fichier de donnees : "+name);

	string poub;

	vector<int> CtoLN; //CtoLN[i] donne l indice dans les tableaux de la classe pour le noeud i dans le fichier de chritian
	vector<size_t> numEvac; //indice des node evac dans le fichier de christian


	//1. on recupere toutes les informatioms, les noeuds sont dans le desordre ( non numerotes evac puis safe puis transit)
	// => on les renumerote apres

	//=================================================================
	//1.1 il faut recuperer la taille pour les alloc

	getline(fic, poub);
	int nbEvacNode, idSafeNode, idNode = 0, nbArc;
	int nbSomCh;
	fic >> nbEvacNode >> idSafeNode;


	// on recupere les indices des node evac
	numEvac.resize(nbEvacNode);
	for (size_t i = 0; i < nbEvacNode; ++i)
	{
		fic >> numEvac[i] >> poub >> poub >> poub >> nbSomCh;
		for (int i = 0; i < nbSomCh; ++i)
		{
			fic >> poub;
		}
	}

	fic >> poub;
	if (poub != "c")
		stopProg("readFileAlicia : pb lecture fichier, on ne recupere pas c");

	getline(fic, poub);


	//nb de noeuds qu'on va effectivement utiliser = nb arcs de l'arbre reduit + 1
	fic >> _nbNode;//on recupere le nb arcs de l'arc = nb node  - 1
	_nbNode++;

	if (_nbNode < 2)
		stopProg("readFileAlicia : pb lecture fichier, nbNode < 2");

	while (poub != "c")
		fic >> poub;
	getline(fic, poub);

	int nbNodeTotal = 0;
	fic >> nbNodeTotal;

	if (nbNodeTotal < _nbNode)
		stopProg("readFileAlicia : pb lecture fichier, nbNodeTotal < _nbNode");

	//cout << _nbNode << endl;

	//==============================================
	// allocation des donnees

	CtoLN.assign(nbNodeTotal, 0);

	// ALLOC 1 (on a une matrice pour les arcs mais le graphe n est pas forcement complet, dans ce cas on garde les valeurs par defaut (TIME_INFINITY ou 0)
	_pop.resize(nbEvacNode + 1);
	_time.assign(_nbNode + 1, vector<int>(_nbNode + 1, TIME_INFINITY));
	_capArc.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_capNode.assign(_nbNode + 1, 0);
	_dead.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));
	_maxRateEvac.resize(nbEvacNode + 1);
	_succ.resize(_nbNode + 1);
	_minDate.assign(_nbNode + 1, TIME_INFINITY);
	_chemins.resize(nbEvacNode + 1);


	_firstEvacuationNode = 1;
	_lastEvacuationNode = nbEvacNode;
	_firstSafeNode = nbEvacNode + 1;
	_lastSafeNode = nbEvacNode + 1; //un seul safe node chez christian
	_firstTransitNode = _lastSafeNode + 1;
	_lastTransitNode = _nbNode;

	// on remplit CtoLN
	for (int i = 0; i < _lastEvacuationNode; ++i)
	{
		CtoLN[numEvac[i]] = i + 1;
	}
	CtoLN[idSafeNode] = _firstSafeNode;

	// on continuera de numeroter les noeuds en fcontion des noeuds dans l arbre reduit
	int nextNumNodeLN = _firstSafeNode + 1;




	//======================================================
	// retour au debut du fichier

	//=======================================================
	//lecture premiere partie :population , max rate, chemin non reduits (on ne stocke pas les chemins non reduits)


	fic.clear();
	fic.seekg(0, ios::beg);

	//derniere date possible pour le depart
	vector<int> due(nbEvacNode + 1);//date de depart du dernier au plus tard pour pop i 


	getline(fic, poub);
	fic >> nbEvacNode >> idSafeNode;
	int idC1;

	// on recupere les infos des node evac
	for (size_t i = 0; i < nbEvacNode; ++i)
	{
		fic >> idNode;
		//cout << idNode << endl;
		fic >> _pop[CtoLN[idNode]] >> _maxRateEvac[CtoLN[idNode]] >> due[CtoLN[idNode]];
		popTotal += _pop[CtoLN[idNode]];

		fic >> nbSomCh;


		//avance
		for (int i = 0; i < nbSomCh; ++i)
		{
			fic >> idC1;

		}

	}

	fic >> poub;

	if (poub != "c")
		stopProg("readFileAlicia : pb lecture fichier, on ne recupere pas c (2eme passage)");
	getline(fic, poub);



	//=======================================================
	//lecture deuxieme partie : arbre reduit
	fic >> nbArc;

	int fils, pere, longueur, cap;
	vector<int> peres(_nbNode + 1);

	// on recupere les indices des node evac
	for (int i = 0; i < nbArc; ++i)
	{
		fic >> fils >> pere >> longueur >> cap;
		if (CtoLN[fils] == 0)
			CtoLN[fils] = nextNumNodeLN++;

		if (CtoLN[pere] == 0)
			CtoLN[pere] = nextNumNodeLN++;

		_time[CtoLN[fils]][CtoLN[pere]] = longueur;
		_capArc[CtoLN[fils]][CtoLN[pere]] = cap;

		peres[CtoLN[fils]] = CtoLN[pere];

	}


	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		//en premier on met le noeud d evac...
		_chemins[i].push_back(SomWithWindow(i, -1, -1));
		int som = -1, prec = i;

		do
		{
			som = peres[prec];
			_chemins[i].push_back(SomWithWindow(som, -1, -1));
			prec = som;
		} while (som != _firstSafeNode);
	}


	//sur les nodes evac et safe on a une capa = a la pop
	_TMAX = 0;
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		_capNode[i] = _pop[i];
		_chemins[i][0]._LF = due[i];
		_chemins[i][0]._ES = 0;

		int prec = _chemins[i][0]._som;

		//on fixe les due dates des arcs sur le chemin et on en profite pour calculer ES et LF
		for (int k = 1; k < _chemins[i].size(); ++k)
		{
			int cour = _chemins[i][k]._som;
			_chemins[i][k]._LF = _chemins[i][k - 1]._LF + _time[prec][cour];
			_chemins[i][k]._ES = _chemins[i][k - 1]._ES + _time[prec][cour];
			_dead[prec][cour] = max(_dead[prec][cour], _chemins[i][k - 1]._LF);//attention dead(i,j) = derniere date pour s'engager sur (i,j)
			if (_succ[prec].empty())
				_succ[prec].push_back(cour);
			prec = cour;
		}
		_TMAX = max(_TMAX, _chemins[i][_chemins[i].size() - 1]._LF);
	}

	for (int i = _firstSafeNode; i <= _lastSafeNode; ++i)
	{
		_capNode[i] = popTotal;

	}


	//pour le PL
	computeMinDate();
	
	//on verfie apres duplication des noeuds evac / transit
	//if (!verifArbreReduit())
		//stopProg("readFileAlicia : arbre mal reduit");
}

//on verifie que l'arbre reduit a les bonnes caracteristiuqes (capa croissantes...)
bool Instance::verifArbreReduit()
{
	bool res = true;

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int prec = _chemins[i][0]._som;

		int capPred = 0;

		for (int k = 1; k < _chemins[i].size(); ++k)
		{
			int cour = _chemins[i][k]._som;
			int capCour = _capArc[prec][cour];
			if (capCour < capPred)
				res = false;

			prec = cour;
		}
		if (_maxRateEvac[i] != _capArc[_chemins[i][0]._som][_chemins[i][1]._som])
			res = false;

		if (_capArc[_chemins[i][0]._som][_chemins[i][1]._som] > _pop[i])
			res = false;
	}

	return res;
}

//utilise les plus courts chemins pour calculer les min dates
void Instance::computeMinDate()
{

	vector< vector<int> > ppcs;


	//1. pour chaque noeud d evac on doit calculer le plus court chemin de ce noeud vers tous les autres
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{

		vector<int> ppc;

		//init
		ppc.assign(_nbNode + 1, INT_INFINITY);
		ppc[i] = 0;

		int cpt = 0;

		vector<int> file;
		file.push_back(i); //on part de i

		while (cpt < toInt(file.size()))
		{

			int cour = file[cpt];
			cpt++;

			//on parcourt tous les succ
			for (size_t s = 0; s < _succ[cour].size(); ++s)
			{
				int succ = _succ[cour][s];

				int val = ppc[cour] + _time[cour][succ];

				if (val < ppc[succ])
				{
					ppc[succ] = val;
					file.push_back(succ);
				}

			}

		}//fin while

		ppcs.push_back(ppc);
	}

	//2. le plus court chemin depuis un noeud evac vers un noeud pas evac = min des ppc sur les noeuds avec
	for (int i = 1; i <= _nbNode; ++i)
	{
		_minDate[i] = ppcs[0][i];

		for (int j = 1; j < ppcs.size(); ++j)
		{
			_minDate[i] = min(_minDate[i], ppcs[j][i]);
		}

	}




}


void Instance::trace(string name, bool sym)
{

	stringstream ss;
	ss << name << ".txt";


	ofstream ficOut(ss.str(), ios::out);// ouverture fichier de sortie

	//premiere ligne du fichier de sortie
	ficOut << "digraph G {" << endl;



	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		ficOut << i << "[color = red];" << endl;

	}

	for (int i = _firstTransitNode; i <= _lastTransitNode; ++i)
	{
		ficOut << i << "[color = black];" << endl;

	}

	for (int i = _firstSafeNode; i <= _lastSafeNode; ++i)
	{
		ficOut << i << "[color = green];" << endl;

	}
	/*
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		ficOut << "0->" << i << "[label=" << _pop[i] << "]" << endl;

	}*/

	if (sym)
	{
		for (int i = 1; i <= _nbNode; ++i)
		{
			for (int j = i + 1; j <= _nbNode; ++j)
			{
				if (_capArc[i][j] > 0)
				{
					ficOut << i << " -> " << j << "[label=<" << _capArc[i][j] << "<br />" << _dead[i][j] << ">]" << endl;
				}
			}
		}
	}
	else
	{
		for (int i = 1; i <= _nbNode; ++i)
		{
			for (int j = 1; j <= _nbNode; ++j)
			{

				if (_capArc[i][j] > 0)
				{
					//ficOut << i << "->" << j << "[label=" << _capArc[i][j] /*<< "." << _dead[i][j]*/ << "]" << endl;
					ficOut << i << "->" << j << "[label=<" << _capArc[i][j] << "<br />"  << _time[i][j]  << ">]" << endl;
				}
				//if (_time[i][j] < TIME_INFINITY)
				//{
				//	ficOut << i << "->" << j << "[label=" << _time[i][j] /*<< "." << _dead[i][j]*/ << "]" << endl;
				//}
			}
		}
	}

	ficOut << "}" << endl;
	ficOut.close();

	stringstream ligCom;


#if defined(_WIN32) || defined(_WIN64)
	ligCom << "dot.exe -Tpdf -o" << name << ".pdf " << name << ".txt";
#endif

	system(ligCom.str().c_str());
}



//trace graphe
// si sym = true alors on trace des aretes sinon des arcs
bool Instance::traceChemin(string name)
{
	bool cheminOK = true;

	stringstream ss;
	ss << name << ".txt";

	vector<vector<bool> > exist(_nbNode + 1, vector<bool>(_nbNode + 1, false));

	ofstream ficOut(ss.str(), ios::out);// ouverture fichier de sortie

	//premiere ligne du fichier de sortie
	ficOut << "digraph G {" << endl;



	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		ficOut << i << "[color = red];" << endl;
	}
	/*
	for (int i = _firstTransitNode; i <= _lastTransitNode; ++i)
	{
		ficOut << i << "[color = black];" << endl;
	}*/

	for (int i = _firstSafeNode; i <= _lastSafeNode; ++i)
	{
		ficOut << i << "[color = green];" << endl;
	}


	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int prec = _chemins[i][0]._som;

		for (int j = 1; j < _chemins[i].size(); ++j)
		{
			int cour = _chemins[i][j]._som;

			if (cour <= _lastEvacuationNode)
				cheminOK = false;

			if (!exist[prec][cour])
			{
				ficOut << prec << " -> " << cour << "[label=< <font color = \"red\">" <<
					_capArc[prec][cour] << "</font><br />" << _time[prec][cour]
					<< ">]" << endl;
				exist[prec][cour] = true;
			}
			prec = cour;
		}
	}


	ficOut << "}" << endl;
	ficOut.close();

	stringstream ligCom;


#if defined(_WIN32) || defined(_WIN64)
	ligCom << "dot.exe -Tpdf -o" << name << ".pdf " << name << ".txt";
#endif

	system(ligCom.str().c_str());

	return cheminOK;
}




//dansle cas ou on a les chemins d evacuation predefini, on peut calculer facilement un Tmax et des fenetres de temps
//dans cette fonction on ajuste le Tmax avec les deadline
void Instance::computeWindowsEtTmax()
{

	// 1. calcul des dates au plus tot (chemin.ES)

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int cour = i; // = _chemins[i][0] puisque le premier noeud = i
		int date = 0;
		_chemins[i][0]._ES = 0;

		for (int j = 1; j < _chemins[i].size(); ++j)
		{

			int suiv = _chemins[i][j]._som;
			date = date + _time[cour][suiv];
			_chemins[i][j]._ES = date;

			cour = suiv;
		}
	}


	// 2. pour chaque chemin on calcule la distance entre le safe node et les noeuds du chemin
	// on les stocke dans la limite up de la fenetre
	// en meme temps on calcul la date max (TmaxPouri) a la quelle on peut arriver au safe node en partant de i
	// on stocke cette date dans _chemins[i][_chemins[i].size() - 1] ( +/- =  un TMax par chemin)
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int cour = _chemins[i][_chemins[i].size() - 1]._som;
		int l = 0;



		int TmaxPouri = TIME_INFINITY;

		for (int j = static_cast<int>(_chemins[i].size()) - 2; j >= 0; --j)
		{
			int prec = _chemins[i][j]._som;
			l = l + _time[prec][cour];
			_chemins[i][j]._LF = l;

			//date max a laquelle on peut arriver au safe node en partant de i
			TmaxPouri = min(TmaxPouri, _dead[prec][cour] + l);


			cour = prec;
		}

		_chemins[i][_chemins[i].size() - 1]._LF = TmaxPouri;

	}




	// 3. pour chaque chemin on remplace  _chemins[i][j]._LF qui nous a servi a stocker la longueur du chemin par la vrai borne sup de la fenetre associee au noeud
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int TmaxPouri = _chemins[i][_chemins[i].size() - 1]._LF;
		for (int j = 0; j < _chemins[i].size() - 1; ++j)
		{

			_chemins[i][j]._LF = TmaxPouri - _chemins[i][j]._LF;
		}
	}


	//4. Tmax = max des deadline de chaque chemin
	_TMAX = 0;
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		_TMAX = max(_TMAX, _chemins[i][_chemins[i].size() - 1]._LF);
	}

}

//dansle cas ou on a les chemins d evacuation predefini, on peut calculer facilement un Tmax et des fenetres de temps
//dans cette fonction on ne change pas le TMax, c'est une donnee du probleme
void Instance::computeWindows()
{

	// 1. calcul des dates au plus tot (chemin.ES)

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int cour = i; // = _chemins[i][0] puisque le premier noeud = i
		int date = 0;
		_chemins[i][0]._ES = 0;

		for (int j = 1; j < _chemins[i].size(); ++j)
		{

			int suiv = _chemins[i][j]._som;
			date = date + _time[cour][suiv];
			_chemins[i][j]._ES = date;

			cour = suiv;
		}
	}


	// 2. pour chaque chemin on calcule la distance entre le safe node et les noeuds du chemin
	// on les stocke dans la limite up de la fenetre
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int cour = _chemins[i][_chemins[i].size() - 1]._som;
		int l = 0;


		for (int j = static_cast<int>(_chemins[i].size()) - 2; j >= 0; --j)
		{
			int prec = _chemins[i][j]._som;
			l = l + _time[prec][cour];
			_chemins[i][j]._LF = l;
			cour = prec;
		}

	}


	// 3. pour chaque chemin on remplace  _chemins[i][j]._LF qui nous a servi a stocker la longueur du chemin 
	//par la vrai borne sup de la fenetre associee au noeud
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		// 3.1 on calcul la fin max du chemin en fonction du TMax et des duedates
		int finMax = _TMAX;

		for (int j = 0; j < _chemins[i].size() - 1; ++j)
		{
			int cour = _chemins[i][j]._som;
			int suiv = _chemins[i][j + 1]._som;
			finMax = min(finMax, _dead[cour][suiv] + _chemins[i][j]._LF); //dans LF strocke la longueur de j vers safe 
		}

		_chemins[i][_chemins[i].size() - 1]._LF = finMax;

		// 3.2 on remplit les LF
		for (int j = 0; j < _chemins[i].size() - 1; ++j)
		{
			int cour = _chemins[i][j]._som;
			int suiv = _chemins[i][j + 1]._som;
			_chemins[i][j]._LF = min(finMax - _chemins[i][j]._LF, _dead[cour][suiv]);

			//#ifdef _VERIF_
			//			if (_chemins[i][j]._LF < _chemins[i][j]._ES)
			//				stopProg("Instance::computeWindows : pb fenetre");
			//#endif
		}
	}



}




//retourne la liste des arcs a la fois dans le chemin partant de i et de j
vector< ArcCommun > Instance::arcEnCommun(int i, int j)
{
	vector< ArcCommun> v;

	//pour chaque arc sur le chemin de  i on regarde s il est dans le chemin de j

	for (int k = 0; k < _chemins[i].size() - 1; ++k)
	{
		int orig = _chemins[i][k]._som;
		int dest = _chemins[i][k + 1]._som;
		int l1 = _chemins[i][k]._ES;

		//retourne la longueur de j vers l arc (orig, dest) s'il existe dans le chemin de j, -1 sinon
		int l2 = appartientChemin(j, orig, dest);

		if (l2 != -1)
			v.push_back(ArcCommun(orig, dest, l1, l2));

	}


	return v;
}

//retourne la liste des arcs dans le chemin i sous forme de vector< ArcCommun >
vector< ArcCommun > Instance::arcCh_i(int i)
{
	vector< ArcCommun> v;



	for (int k = 0; k < _chemins[i].size() - 1; ++k)
	{
		int orig = _chemins[i][k]._som;
		int dest = _chemins[i][k + 1]._som;
		int l1 = _chemins[i][k]._ES;


		v.push_back(ArcCommun(orig, dest, l1, 0));

	}


	return v;
}




double Instance::getDebitMin(int i)
{
	return static_cast<double> (_pop[i]) / (this->getFinMax(i) - this->getLongueurCh(i) + 1);
}

double Instance::getDebitMax(int i) const
{
	return static_cast<double> (_maxRateEvac[i]);
}

//retourne la longueur de j vers l arc (orig, dest) s'il existe dans le chemin de j, -1 sinon
int Instance::appartientChemin(int j, int orig, int dest)
{

	int res = -1;


	int prec = _chemins[j][0]._som;

	int k = 1;

	while (res == -1 && k < _chemins[j].size())
	{
		int cour = _chemins[j][k]._som;

		if (prec == orig && cour == dest)
			res = _chemins[j][k - 1]._ES;

		prec = cour;
		k++;
	}

	return res;
}

//retourne la longueur du chemin de i vers le safe node
int Instance::getLongueurCh(int i)
{
	return _chemins[i][_chemins[i].size() - 1]._ES;
}

//retourne la date max pour arriver au safe node en partant de i (LF du dernier noeud du chemin)
int Instance::getFinMax(int i)
{
	return _chemins[i][_chemins[i].size() - 1]._LF;
}

//dans le cas ou des chemins sont definis, taux d'evac max <= CAP des arcs du chemin
void Instance::reduitMaxRateEvac()
{

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int maxRate = min(_maxRateEvac[i], _pop[i]);

		int prec = _chemins[i][0]._som;
		for (int k = 1; k < _chemins[i].size(); ++k)
		{
			int cour = _chemins[i][k]._som;

			maxRate = min(maxRate, _capArc[prec][cour]);
			prec = cour;
		}
		_maxRateEvac[i] = maxRate;
	}
}


//on recalcule les LF pour chaque job du chemin :
// la fin max d'un chemin devient ceil(finChemin[i]) * alpha
void Instance::recalculeLF(vector<double> & finChemin, double alpha)
{
	_TMAX = 0;
	_dead.assign(_nbNode + 1, vector<int>(_nbNode + 1, 0));

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int last = static_cast<int> (_chemins[i].size()) - 1;

		int finCh = static_cast<int> (alpha * ceil(finChemin[i]) + EPSILON);

		_TMAX = max(_TMAX, finCh);


		_chemins[i][last]._LF = finCh;
		int longCh = 0;

		for (int k = last; k >= 1; --k)
		{
			int orig = _chemins[i][k - 1]._som;
			int dest = _chemins[i][k]._som;

			longCh += _time[orig][dest];

			_chemins[i][k - 1]._LF = finCh - longCh;

			//on est oblige de prendre le max car plusieurs population 
			//avec des deadlines differentes peuvent emprunter le meme arc
			_dead[orig][dest] = max(_dead[orig][dest], _chemins[i][k - 1]._LF);


		}

	}

}


//on verifie que les fenetres de temps sont realisables
bool Instance::verifFen()
{
	bool ok = true;

	//1. verif fenetres
	for (int i = 1; ok && i <= _lastEvacuationNode; ++i)
	{
		for (int j = 0; ok && j < _chemins[i].size(); ++j)
		{
			if (_chemins[i][j]._ES > _chemins[i][j]._LF)
				ok = false;
		}
	}

	//2. verif deadline chemins
	for (int i = 1; ok && i <= _lastEvacuationNode; ++i)
	{
		if (getFinMax(i) < getLongueurCh(i) + static_cast<double>(_pop[i]) / _maxRateEvac[i] - 1)
			ok = false;
	}

	return ok;
}

//retourne l'ensemble des noeuds d'evacuation qui sont aussi noeuds de transit
// on memorise app | app[i] = vrai si i appartient a l'ensemble
// et prec | prec[i] = liste des precdent de i
vector<int> Instance::evacTransitNode(vector<bool> & app, vector < vector<int> > & vprec)
{
	vector<int> v;
	app.assign(_lastEvacuationNode + 1, false);
	vprec.resize(_lastEvacuationNode + 1);

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int prec = _chemins[i][0]._som;

		for (int j = 1; j < _chemins[i].size(); ++j)
		{
			int cour = _chemins[i][j]._som;

			if (cour <= _lastEvacuationNode)
			{
				vprec[cour].push_back(prec);
				if (!app[cour])
				{
					v.push_back(cour);
					app[cour] = true;
				}
			}
			prec = cour;
		}
	}

	return v;
}

//on duplique les noeuds qui sont a la fois evacuation et transit
//et on met a jour toutes les infos pour que la solution de l'instance dupliquée soit la meme que celle de l'instance d'origine
//ATTENTION : les longueurs des chemins et des deadlines sont augmentes de 1 avec cette manip, donc les marges restent les mêmes mais 
//les dates d'arrivee sont augmentee d'une unite
//(cette fonction ne fonctionne que si le graphe est un arbre)
void Instance::dupliqueNodeEvacTransit()
{

	// on duplique un noeud i qui est a la fois evac et transit en i_e = i et i_t : 
	// la partie evac i_e garde l indice initial et la partie transit i_t aura un nouvel indice

	// les nouveaux noeuds vont etre numerotes de _nbNode a _nbNode + nbNewNode - 1
	// ils seront du type transit : ils viennent a la suite de la liste des noeuds de transit
	// les SDD (pop, maxRate,...) qui concernent uniquement les noeuds evacuation ne sont pas impactees


	//===================================================
	//1. on identifie les jobs qui sont a la fois evac et transit
	//et on cree la SDD permettant de passer de l indice initial au nouvel indice
	vector <bool> app; //app[i] = vrai si i est dans node (alloue et init dans evacTransitNode)
	vector < vector<int> > prec;//prec[i] = liste des noeuds qui precede direct i (indice initiaux)

	vector<int> node = evacTransitNode(app, prec);
	int nbNewNode = static_cast<int> (node.size());

	if (nbNewNode > 0)
	{
		_hasDuplicateNode = true;

		
		vector<int> evac2transit(_lastEvacuationNode + 1);//evac2transit(i) donne l indice qu on attribue a la partie transit de i

		int cpt = _nbNode + 1;
		for (int i : node)
		{
			evac2transit[i] = cpt;
			cpt++;
		}


		//===================================================
		//2. on met a jour les SDD necessaires
		int nbAncienNode = _nbNode;

		_nbNode += nbNewNode;
		_lastTransitNode = _nbNode;//on indice a partir de 1
		_TMAX++;

		//--------------------------------------------------
		//2.1 modif de time

		//on alloue une nouvelle matrice car l'ancienne est trop petie
		vector< vector<int> > time2(_nbNode + 1, vector<int>(_nbNode + 1, TIME_INFINITY));

		//on commence par remplir time2 comme _time pour les noeuds qui existaient deja, ensuite on fera les modifications

		for (int i = 0; i <= nbAncienNode; ++i)
			for (int j = 0; j <= nbAncienNode; ++j)
				time2[i][j] = _time[i][j];

		//les distances des arcs qui arrivent sur les sommets dupliques deviennent +infini (arc supprime)
		//de meme que les capacites au successeur actuel (arc supprime)
		for (int i : node)
		{
			for (int j = 0; j < prec[i].size(); ++j)
				time2[prec[i][j]][i] = TIME_INFINITY;
			time2[i][_succ[i][0]] = TIME_INFINITY;
		}

		//on met une longueur 1 entre le noeud evac duplique et sont noeud transit correspondant...
		for (int i : node)
		{
			time2[i][evac2transit[i]] = 1;
		}

		//les nouveaux noeuds transit doivent être raccroche a leur succ
		for (int i : node)
		{
			int s = _succ[i][0];//on est dans un arbre : on n a qu un succ
			int s2 = s;
			if (s <= _lastEvacuationNode && app[s])
				s2 = evac2transit[s];//attention le successeur a peut etre change d'indice s'il fait parti des sommets dupliques

			time2[evac2transit[i]][s2] = _time[i][s];
		}

		for (int i : node)
		{
			for (int p : prec[i])
			{
				int p2 = p;
				if (p <= _lastEvacuationNode && app[p])
					p2 = evac2transit[p];//attention le successeur a peut etre change d'indice s'il fait parti des sommets dupliques

				time2[p2][evac2transit[i]] = _time[p][i];
			}
		}
		//du coup pour ne pas pénaliser les noeuds dupliques on ajoute 1 a la longueur de l'arc sortan d'un noeud evac non duplique
		// (on ajoutera aussi 1 aux deadline pour ne pas modifier les marges)

		for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
		{
			if (!app[i])
			{
				int s = _succ[i][0];//on est dans un arbre : on n a qu un succ
				int s2 = s;
				if (s <= _lastEvacuationNode && app[s])
					s2 = evac2transit[s];//attention le successeur a peut etre change d'indice s'il fait parti des sommets dupliques

				time2[i][s2] += 1;
			}
		}
		//
		_time = time2;

		//--------------------------------------------------
		//2.2 modif de capArc

		//on alloue une nouvelle matrice car l'ancienne est trop petie
		vector< vector<int> > cap2(_nbNode + 1, vector<int>(_nbNode + 1, 0));

		for (int i = 0; i <= nbAncienNode; ++i)
			for (int j = 0; j <= nbAncienNode; ++j)
				cap2[i][j] = _capArc[i][j];


		//les capacites des arcs qui arrivent sur les sommets dupliques deviennent nulles (arc supprime)
		//de meme que les capacites au successeur actuel (arc supprime)
		for (int i : node)
		{
			for (int j = 0; j < prec[i].size(); ++j)
				cap2[prec[i][j]][i] = 0;
			cap2[i][_succ[i][0]] = 0;

		}


		//capa a ajouter sur les arcs nouveaux : arcs nouveaux = arcs qui ont leur origine ou destination dans les nouveaux noeuds
		//

		for (int i : node)
			cap2[i][evac2transit[i]] = _maxRateEvac[i];


		//on doit ajouter les capacités sur les arcs dont l'origine ou la destination est un nouveau sommet
		for (int i : node)
		{
			int s = _succ[i][0];//on est dans un arbre : on n a qu un succ
			int s2 = s;
			if (s <= _lastEvacuationNode && app[s])
				s2 = evac2transit[s];//attention le successeur a peut etre change d'indice s'il fait parti des sommets dupliques

			cap2[evac2transit[i]][s2] = _capArc[i][s];
		}

		for (int i : node)
		{
			for (int p : prec[i])
			{
				int p2 = p;
				if (p <= _lastEvacuationNode && app[p])
					p2 = evac2transit[p];//attention le successeur a peut etre change d'indice s'il fait parti des sommets dupliques

				cap2[p2][evac2transit[i]] = _capArc[p][i];
			}
		}
		_capArc = cap2;

		//--------------------------------------------------------------
		//3.capNode, on ajoute les nouveaux sommets avec une capacite node = 0

		vector<int> capNode2(_nbNode + 1, 0);
		for (int i = 0; i <= nbAncienNode; ++i)
			capNode2[i] = _capNode[i];

		_capNode = capNode2;


		//-----------------------------------------------------------------------------
		//4. on remplit les successeurs a partir des chemins (on doit donc deja avoir adapte les chemins)

		vector<vector< SomWithWindow > > chemin2(_lastEvacuationNode + 1);

		for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
		{


			//le premier sommet du chemin est inchange : 
			chemin2[i].push_back(_chemins[i][0]);

			//si i a ete duplique il faut ajouter le noeud transit correspondant
			if (app[i])
			{
				chemin2[i].push_back({ evac2transit[i], 1, _chemins[i][0]._LF + 1 });
			}

			//pour les sommets suivent il faut verifier que leur indice n'a pas ete modifie
			//et augmente d'un les dates ES / LF
			for (int j = 1; j < _chemins[i].size(); ++j)
			{
				const SomWithWindow & S = _chemins[i][j];
				if (S._som > _lastEvacuationNode) //si le sommmet courant ne fait pas parti des sommets modifies
				{
					chemin2[i].push_back({ S._som, S._ES + 1, S._LF + 1 });

				}
				else //sinon mettre son nouvel indice
				{
					chemin2[i].push_back({ evac2transit[S._som], S._ES + 1, S._LF + 1 });

				}

			}
		}

		_chemins = chemin2;

		//-----------------------------------------------------------------------------
		//4.2 on verifie la capa du premier arc des chemins, elle doit etre = maxRate

		for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
		{
#ifdef _VERIF_
			if (_capArc[i][_chemins[i][1]._som] < _maxRateEvac[i])
				stopProg("Instance::" + static_cast<string>(__func__) + " : le maxRate est normalement <= a la capa du 1er arc");
#endif
			_capArc[i][_chemins[i][1]._som] = _maxRateEvac[i];
		}


		//-----------------------------------------------------------------------------
		//5. on remplit les successeurs a partir des chemins (on doit donc deja avoir adapte les chemins)
		vector<vector<int> > succ2(_nbNode + 1);

		for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
		{
			int prec = _chemins[i][0]._som;
			for (int j = 1; j < _chemins[i].size(); ++j)
			{
				int cour = _chemins[i][j]._som;
				if (succ2[prec].empty())
					succ2[prec].push_back(cour);
				prec = cour;
			}
		}
		_succ = succ2;


		for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
		{
			//on en profite pour mettre la capacite du premeir arc ) maxRate (la seule possibilite pour que ce ne soit pas deja le cas c'est
			// que la population est < a la capa de l'arc (la capa de l'arc ne peut que diminuer)
#ifdef _VERIF_
			if (_capArc[i][_chemins[i][1]._som] < _maxRateEvac[i])
				stopProg("Instance::" + static_cast<string>(__func__) + " : le maxRate est normalement <= a la capa du 1er arc");
#endif
			_capArc[i][_chemins[i][1]._som] = _maxRateEvac[i];
		}

		//--------------------------------------------------------------
		//6. deadline : on recalcule les deadline sur les arcs a partir des 
		//deadlines des chemins (la fonction realloue _dead si besoin)

		recalculeDead();

		//--------------------------------------------------------------------------
		//7. pour le PL ==> recalcule de _minDate

		_minDate.assign(_nbNode + 1, 0);
		computeMinDate();


	}
}


//a partir des deadlines des chemins on recalcule les deadlines sur les arcs 
void Instance::recalculeDead()
{
	_dead.assign(_nbNode+1, vector<int>(_nbNode+1, 0));

	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		int prec = _chemins[i][0]._som;
		for (int j = 1; j < _chemins[i].size(); ++j)
		{
			int cour = _chemins[i][j]._som;
			int LF = _chemins[i][j-1]._LF;
			_dead[prec][cour] = max (_dead[prec][cour], LF);//derniere date pour aller sur l'arc (i,j) = LF du prec
			prec = cour;
		}
	}

}



//attention cette fonction doit etre appelee quand tout est deja initialise
void Instance::calcul_stat()
{

	double sumRapportCap = 0;

	//1. stat longueur

	_stat_temps_total_min = 0;

	//ici on utilise la formule avec le -1 => qd on utilise cette longueur ^pour l'article il faut ajouter 1 
	for (int i = _firstEvacuationNode; i <= _lastEvacuationNode; ++i)
	{
		_stat_temps_total_min = max(_stat_temps_total_min,
			static_cast<double>(_pop[i]) / _maxRateEvac[i] -1 + getLongueurCh(i));
	}
	
	//2. stat capa

	for (int i = _firstTransitNode; i <= _lastTransitNode; ++i)
	{
		double capEntrante = 0;
		double capSortante = 0;

		for (int j = 0; j <= _nbNode; ++j)
		{
			//rq. si l'arc n'exsite pas alors capa  = 0

			capEntrante += _capArc[j][i];

			capSortante += _capArc[i][j];
		}
		//cout << capEntrante / capSortante << endl;

		sumRapportCap += capEntrante / capSortante;
	}

	_stat_capa = sumRapportCap / (_lastTransitNode - _firstTransitNode + 1);
}