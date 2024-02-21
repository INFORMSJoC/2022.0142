#include "InstanceReduite.h"
#include <cassert>



//initialise les attributs avec ins
void InstanceReduite::init(const Instance & ins, int orig, int dest)
{
	
	int Nmax =  ins._lastEvacuationNode - ins._firstEvacuationNode + 1;
	_nbArc = 0;

	_ESmin = TIME_INFINITY;


	_ES.clear(); _ES.reserve(Nmax);
	_LF.clear(); _LF.reserve(Nmax);
	_pop.clear(); _pop.reserve(Nmax);
	_cap.clear(); _cap.reserve(2* Nmax +1);
	_arcToJob.clear(); _arcToJob.resize(2* Nmax +1);
	_jobToArcInit.clear(); _jobToArcInit.reserve(Nmax);
	_jobToArc.clear(); _jobToArc.resize(Nmax);
	
	_debitMin.clear(); _debitMin.reserve(Nmax);
	_debitMax.clear(); _debitMax.reserve(Nmax);



	
	//1. on parcourt les arcs et on leur donne des numero

	vector<vector<int>> indArc(ins._nbNode+1, vector<int>(ins._nbNode+1, -1));//indArc[a][b] = indice de l'arc (a,b) si deja obetnu, -1 sinon
	_nbJob = 0;

	for (int i = ins._firstEvacuationNode; i <= ins._lastEvacuationNode; ++i)
	{
		const vector<SomWithWindow> & chemin_i = ins._chemins[i];
		
		//si le dernier arc de i est (orig, dest) on l'ajoute...
		if (chemin_i[chemin_i.size() - 2]._som == orig && chemin_i[chemin_i.size() - 1]._som == dest)
		{
			_pop.push_back(ins._pop[i]);
			_job.push_back(i);
			_jobToArcInit.push_back(_nbArc);
			double dureeMax = 0;

			//...et je parcours ses arcs
			int prec = chemin_i[0]._som;
			for (int k = 1; k < chemin_i.size(); ++k)
			{
				int cour = chemin_i[k]._som;

				if (prec == orig && cour == dest)
				{
					_ES.push_back(chemin_i[k - 1]._ES); //date au plus tot d'arrivee a la porte
					_LF.push_back(chemin_i[k - 1]._LF+1); //date au plus tard de sortie de la porte (mettre +1 si on considere la porte imaginaire) : +1 oblige pour etre compatible avec RCPCP_Algo
					_ESmin = min(_ESmin, *_ES.rbegin());
					dureeMax = *_LF.rbegin() - *_ES.rbegin();
				}
				if (indArc[prec][cour] == -1)
				{
					indArc[prec][cour] = _nbArc;
					_nbArc++;
					_cap.push_back(ins._capArc[prec][cour]);
				}
#ifdef _VERIF_
				if (_nbArc >= _arcToJob.size())
					stopProg("InstanceReduite::init : _arcToJob trop petit");
#endif
				_arcToJob[indArc[prec][cour]].push_back(_nbJob);
				_jobToArc[_nbJob].push_back(indArc[prec][cour]);
	

				prec = cour;
			}

#ifdef _VERIF_
			if (dureeMax < EPSILON)
				stopProg("InstanceReduite::init : dureeMax doit etre > 0");
#endif
			//maj debit min et max
			_debitMax.push_back(ins.getDebitMax(i));
			_debitMin.push_back(static_cast<double>(ins._pop[i]) / dureeMax);



			_nbJob++;
		}//fin if
	}


	//on met a jour utilise ici 
	_utilise.assign(_nbJob, vector<bool>(_nbArc, false));


	for (int i = 0; i < _nbJob; ++i)
	{
		for (int a : _jobToArc[i])
			_utilise[i][a] = true;
	}




#ifdef _VERIF_
	for (int a = 0; a < _nbArc; ++a)
	{
		if (_arcToJob[a].size() == 1 && _jobToArcInit[_arcToJob[a][0]] != a)
			stopProg("InstanceReduite::init : _jobToArcInit et _arcToJob incompatible");
	}
#endif

}
void InstanceReduite::alloc(int nbJob, int nbArc)
{
	_pop.resize(nbJob);
	_ES.resize(nbJob);
	_LF.resize(nbJob);

	_debitMax.resize(nbJob);
	_debitMin.resize(nbJob);

	_cap.resize(nbArc);

	_arcToJob.resize(nbArc);
	_jobToArc.resize(nbJob);
	_utilise.assign(nbJob, vector<bool>(nbArc, false));
}



//generation aleatoire d'une InstanceReduite : un job utilise un arc avec une proba alpha
void InstanceReduite::genereRandom(int nbJob, int nbArc, double alpha, int ESmax)
{

	_nbJob = nbJob;
	_nbArc = nbArc;

	// 1. alloc memoire
	alloc(nbJob, nbArc);


	// 2. generation 

	//2.1 donnees relatives aux jobs

	for (int i = 0; i < _nbJob; ++i)
	{
		_pop[i] = 1000 + rand() % 4000;
		_ES[i] = rand() % ESmax;
	}   

	//on laisse au moins un job commencer en 0 (le premier)
	_ES[0] = 0;
	_ESmin = 0;

	//2.2 donnees relatives aux arcs
	for (int j = 0; j < _nbArc; ++j)
	{
		_cap[j] = 50 + rand() % 100;
	}

	//2.3 utilisation des arcs par les jobs
	for (int i = 0; i < _nbJob; ++i)
	{
		int nbArcConso = 0;
		_debitMax[i] = INT_INFINITY;

		for (int j = 0; j < _nbArc; ++j)
		{
			//i utilise j avec une proba alpha
			if ((double)(rand()) / RAND_MAX < alpha)
			{
				nbArcConso++;
				_utilise[i][j] = true;
				_arcToJob[j].push_back(i);
				_jobToArc[i].push_back(j);
				_debitMax[i] = min(_debitMax[i], (double)(_cap[j]));
			}
		}
		//on verifie que le job consomme au moins un arc, sinon on en ajoute un random
		if (nbArcConso == 0)
		{
			int j = rand() % _nbArc;
			_utilise[i][j] = true;
			_arcToJob[j].push_back(i);
			_jobToArc[i].push_back(j);
			_debitMax[i] = min(_debitMax[i], (double)(_cap[j]));
		}
	}

	//2.4 on lance une heuristique pour fixer les dates de fin au plus tard
	fixeLF();
	//int gamma = _nbJob;
	//for (int i = 0; i < _nbJob; ++i)
	//	_LF[i] = _ES[i] + gamma*_pop[i] / _debitMax[i] + rand() % ESmax- ESmax; //+_ES[i]/2 

	for (int i = 0; i < _nbJob; ++i)
		_debitMin[i] = _pop[i] / (_LF[i] - _ES[i]);

	//si un arc n est utilise par aucun job, on l'ajoute a un job random
	for (int a = 0; a < _nbArc; ++a)
	{
		if (_arcToJob[a].size() == 0)
		{
			int irand = rand() % _nbJob;
			_utilise[irand][a] = true;
			_arcToJob[a].push_back(irand);
			_jobToArc[irand].push_back(a);
			_debitMax[irand] = min(_debitMax[irand], (double)(_cap[a]));
		}
	}
}


void InstanceReduite::fixeLF()
{
	//on va fixer le debit des jobs et on en deduit la longueur (durée) et hauteur (conso) des jobs 
	int longueur, hauteur;
	double maxLF = 0;


	//date de debut calcul par cette heuristique gloutonne (on prend les jobs un par un et on les met au plus tot)
	vector<double> debut(_nbJob);
	vector<double> debit(_nbJob);
	vector<pair<int, vector<int>>> conso; //conso donne les couple (t, conso) ordonnes par t croissant : consommation a la date t de chaque ressource (arc) (conso inchangee jusqu'au prochain t dans conso)

	//init conso pour les dates 0 et TMAX
	vector<int> v(_nbArc,0);
	conso.push_back({ 0,v });
	conso.push_back({ TIME_INFINITY,v });

	//pour chaque job on fixe son debit a min cap(e) / 3 et on construit un planning (heuristique gloutonne triviale)
	for (int i = 0; i < _nbJob; ++i)
	{
		int minCap = INT_INFINITY;

		// mars 2023 => generer des instances pas trop difficles
		for (int e : _jobToArc[i])
		{
			if (_arcToJob[e].size() >= 2)
				minCap = min(minCap, (int)(2 * _cap[e] / _arcToJob[e].size()));
			else
				minCap = min(minCap, _cap[e]);
		}
		int v = minCap; 

		// generer des instances difficiles (diviser minCap par 4)
		//for (int e : _jobToArc[i])
		//	minCap = min(minCap, _cap[e] );
		//int v = minCap / 4;

		longueur = (int)( ceil((double)(_pop[i]) / v));
		hauteur = v;
		debit[i] = v;

		int t = (int)(_ES[i]);
		
		//on cherche l'indice dans conso qui correspond a la plus grande date <= t
		int it = 0;
		while (conso[it].first < t)
			it++;
		if (conso[it].first > t)//si on a depasse t (cas ou t n'est pas dans conso), on revient en arriere
			it--;

		//on cherche la premiere date ou mettre le job i
		while (!consoOk(i, t, it, conso, longueur, hauteur))
		{
			it++;
			t = conso[it].first;
		}

		assert(t == max(conso[it].first, (int)(_ES[i])));
		debut[i] = t;
		insere(i, t, it, longueur, hauteur, conso);

		//on peut fixer LF
		_LF[i] = t + longueur; //MODIF mars 2023 :-100
	}

	//dessine(_nbJob, debut, debit, _pop);

	//if (!isSolution(debut, debit))
	//	stopProg("InstanceReduite::fixeLF : le glouton ne donne pas une sol. rea");
}


//consoOk renvoie vrai si on peut mettre le job a la date t telle que conso[it].first est la plus grande date <= t
bool InstanceReduite::consoOk(int job, int t, int it, const vector<pair<int, vector<int>>> & conso, int longueur, int hauteur)
{

		
	//pour chaque date dans conso a verifier
	while (conso[it].first < t + longueur)
	{
		const vector<int> & v = conso[it].second;

		//on verifie la capa de chaque job qui utilise e
		for (int e : _jobToArc[job])
		{
			if (v[e] + hauteur > _cap[e])
				return false;
		}
		it++;
	}

	return true;
}

//insere le job (longueur * hauteur) en t dans conso
	//it est l'indice du plus grand t dans conso <= t
void InstanceReduite::insere(int job, int t, int it, int longueur, int hauteur, vector<pair<int, vector<int>>> & conso)
{
	assert(conso[it].first <= t);

	//si conso[it].first != t alors la date t doit etre ajoutee dans conso
	if (conso[it].first != t)
	{
		vector<int> v = conso[it].second;
		conso.insert(conso.begin() + it + 1, { t,v });
		it++;
	}
	
	assert(conso[it].first == t);

	//pour toutes les dates \in [t, t+longueur[ on ajoute la conso du job

	while (conso[it].first < t + longueur)
	{
		for (int e : _jobToArc[job])
		{
			(conso[it].second)[e] += hauteur;
			assert((conso[it].second)[e] <= _cap[e]);
		}
		it++;
	}

	assert(conso[it].first >= t + longueur);

	//on ajoute la date t + longueur dans conso si elle n'y est pas deja
	if (conso[it].first > t + longueur)
	{
		//la conso en t + longueur est celle de la date precedente en enlevant la conso du job courant
		vector<int> v = conso[it - 1].second;
		for (int e : _jobToArc[job])
			v[e] -= hauteur;
		conso.insert(conso.begin() + it, { t + longueur, v });
	}

}

void InstanceReduite::ecrire(const string & ficName)
{
	ofstream fic(ficName);
	   
	fic << _nbJob << endl;
	fic << _nbArc << endl;

	for (int i = 0; i < _nbJob; ++i)
	{
		fic << _pop[i] << " " << _ES[i] << " " << _LF[i] << endl;
		fic << _jobToArc[i].size() << " ";
		for (int e : _jobToArc[i])
			fic << e << " ";
		fic << endl;
	}

	for (int j = 0; j < _nbArc; ++j)
		fic << _cap[j] << endl;

	fic.close();
}

void InstanceReduite::lire(const string & ficName)
{
	ifstream fic(ficName);
	if (!fic)
		stopProg("InstanceReduite::lire fichier non trouve : " + ficName);

	fic >> _nbJob;
	fic >> _nbArc;

	alloc(_nbJob, _nbArc);
	_ESmin = INT_INFINITY;

	for (int i = 0; i < _nbJob; ++i)
	{
		int m;
		fic >> _pop[i] >> _ES[i] >> _LF[i];
		fic >> m;
		_debitMax[i] = INT_INFINITY;

		for (int j = 0; j < m; ++j)
		{
			int e;
			fic >> e;
			_jobToArc[i].push_back(e);
			_arcToJob[e].push_back(i);
			_utilise[i][e] = true;
			
		}
		_ESmin = min(_ESmin, _ES[i]);
	}

	for (int j = 0; j < _nbArc; ++j)
		fic >> _cap[j];

	for (int i = 0; i < _nbJob; ++i)
	{
		_debitMin[i] = _pop[i] / (_LF[i] - _ES[i]);
		for (int e : _jobToArc[i])
			_debitMax[i] = min(_debitMax[i], (double)(_cap[e]));
	}


	for (int j = 0; j < _nbArc; ++j)
	{
		if (_arcToJob[j].size() == 0)
			stopProg("InstanceReduite::lire : des ressources non utilisees par les jobs");
	}



	fic.close();
}


bool InstanceReduite::isSolution(const vector<double> & dateDebut, const vector<double> & debit, double marge)
{
	vector<double> fin(_nbJob);

	vector<pair<double, int>> sigma;//vector de couples (date, job) ordonne sur les dates croissantes (les dates sont debut et fin)
	sigma.reserve(_nbJob * 2);

	for (int i = 0; i < _nbJob; ++i)
	{
		fin[i] = dateDebut[i] + _pop[i] / debit[i];
		if (fin[i] > _LF[i] + EPSILON_DATE)
		{
			cout << "fin_" << i << " = " << fin[i] << " et LF = " << _LF[i] << endl;
			return false;
		}
		if (dateDebut[i] + EPSILON_DATE < _ES[i])
		{
			cout << "deb_" << i << " = " << dateDebut[i] << "ES = " << _ES[i] << endl;
			return false;
		}
		if (debit[i] > _debitMax[i] + EPSILON_COMP || debit[i] < _debitMin[i] - EPSILON_COMP)
		{
			cout << "debit_" << i << " = " << debit[i] << "min = " << _debitMin[i] << "max = " << _debitMax[i] << endl;
			return false;
		}

		sigma.push_back({ dateDebut[i],i });
		sigma.push_back({ fin[i]-EPSILON_DATE,i }); //si une fin = un debut alors la fin doit etre avant dans sigma
	}

	//on verifie la conso de ressource 
	vector<double> conso(_nbArc,0);//conso courante pour chaque arc
	vector<bool> vu(_nbJob, false);//la premiere fois qu on recontre i dans sigma c est la date de debut (on augmente conso), sinon c est la fin (on diminue conso)
	sort(sigma.begin(), sigma.end());//on va parcourir toutes les dates debut et fin dans l'ordre croissant et mettre conso a jour en consequence

	for (auto p : sigma)
	{
		int job = p.second;
		int sens = 1;
		
		//si on a deja vu le job alors il finit => il faudra soustraire ses ressource a conso
		if (vu[job])
			sens = -1;
		vu[job] = true;

		for (int e : _jobToArc[job])
		{
			conso[e] += sens * debit[job];
			if (conso[e] > _cap[e] + EPSILON_COMP)
			{
				cout << "conso_" << e << " = " << conso[e] << "cap = " << _cap[e] << endl;
				return false;
			}
		}
	}

	//si on a passe la marge min correspondant a la solution, on la verifie
	if (marge > -EPSILON)
	{
		double marge_verif = INT_INFINITY;
		for (int i = 0; i < _nbJob; ++i)
			marge_verif = min(marge_verif, _LF[i] - fin[i]);

		if (abs(marge_verif - marge) > EPSILON_DATE)
		{
			cout << "marge = " << marge << "verif = " << marge_verif << endl;
			return false;
		}
	}


	return true;
}