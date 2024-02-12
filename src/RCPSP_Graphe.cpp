#include "RCPSP_Graphe.h"


void RCPSP_Graphe::constuireGraphe()
{
	vector<vector<bool> > exist(_ins->_nbNode + 1, vector<bool>(_ins->_nbNode + 1, false));
	_arc2Ind.assign(_ins->_nbNode + 2, vector<int>(_ins->_nbNode + 2, -1));

	//1. nb de sommets
	_nbSommet = _ins->_lastEvacuationNode - _ins->_firstEvacuationNode + 1; //les evec nodes sont numerotes de 1 a nbEvacNode dans instance donc on garde la mm numerotation
	
	//2. alloc des vecteurs (_TLcond = - infini si pas de TL)
	_listeArcPartage.resize(_nbSommet + 2, vector< vector<ArcCommun > > (_nbSommet+2, vector<ArcCommun> () ) );
	_TLcondFixe.assign(_nbSommet + 2, vector<int>(_nbSommet + 2, -TIME_INFINITY ) );



	//3.1. on remplit la liste des arcs partages
	for (int i = 1; i <= _nbSommet; ++i)
	{
		for (int j = i+1; j <= _nbSommet; ++j)
		{
			_listeArcPartage[i][j] = _ins->arcEnCommun(i, j);
			_listeArcPartage[j][i] = _listeArcPartage[i][j];
		}

		_listeArcPartage[0][i] = _ins->arcCh_i(i);//source => tous les arcs sont partages
		_listeArcPartage[i][_nbSommet+1] = _ins->arcCh_i(i);//puits => idem
		_listeArcPartage[i][0] = _listeArcPartage[0][i];
		_listeArcPartage[_nbSommet + 1][i] = _listeArcPartage[i][_nbSommet + 1];
	}
	//3.2 arcs partages entre source et puits => fait a la fin de la fonction car necessite _arcs initialise


	//4. on remplit les TL
	for (int i = 1; i <= _nbSommet; ++i) //<= car on doit remplie s -> i pour tout i
	{
		for (int j = i + 1; j <= _nbSommet; ++j)
		{

			//on veut le max de ES_i + time(arc) - ES_j pour tous les arcs partages entre i et j
			int maxTLij = - TIME_INFINITY;
			int maxTLji = - TIME_INFINITY;

			for (int k = 0; k < _listeArcPartage[i][j].size(); ++k)
			{
				ArcCommun A = _listeArcPartage[i][j][k];
				int o = A._orig;
				int d = A._dest;

				int TLij = A._long1  - A._long2;
				maxTLij = max(maxTLij, TLij);

				int TLji = A._long2  - A._long1;
				maxTLji = max(maxTLji, TLji);
			}

			_TLcondFixe[i][j] = maxTLij;
			_TLcondFixe[j][i] = maxTLji;
		}


		//5. cas particulier : source et puits
		_TLcondFixe[0][i] = 0;
		_TLcondFixe[i][_nbSommet+1] = _ins->getLongueurCh(i);
	}

	//6. rempilt les arcs

	for (int i = 1; i <= _nbSommet; ++i) 
	{
		int pred = _ins->_chemins[i][0]._som;
		for (int k = 1; k < _ins->_chemins[i].size(); ++k)
		{
			
			int cour = _ins->_chemins[i][k]._som;


			if (!exist[pred][cour])
			{
				_arc2Ind[pred][cour] = static_cast<int>(_arcs.size());
				_arcs.push_back(Arc(pred, cour));
				exist[pred][cour] = true;
			}

			pred = cour;
		}
	}


	//on remplit la SDD utilise (on ne peut pas le faire en meme temps car on  a besoin du nombre d'arcs total

	_nbArc = static_cast<int> (_arcs.size());
	_utilise.assign(_nbSommet + 1, vector<bool>(_nbArc, false));

	for (int i = 1; i <= _nbSommet; ++i)
	{
		int pred = _ins->_chemins[i][0]._som;
		for (int k = 1; k < _ins->_chemins[i].size(); ++k)
		{

			int cour = _ins->_chemins[i][k]._som;

			int a = _arc2Ind[pred][cour];
			_utilise[i][a] = true;

			pred = cour;
		}
	}



	// 3.2 pour eviter des cas particuliers dans les algos on suppose que source et puits partagent tous les arcs
	for (int e = 0; e < _arcs.size(); ++e)
	{
		int o = _arcs[e]._orig;
		int d = _arcs[e]._dest;

		_listeArcPartage[0][_nbSommet + 1].push_back({ o, d, 0, 0 });
	}
	_listeArcPartage[_nbSommet + 1][0] = _listeArcPartage[0][_nbSommet + 1];


}


// renvoie vrai si l'arc (o, d) est dans la liste d'arc partages entre i et j
bool RCPSP_Graphe::isDansArcPartage(int i, int j, int o, int d)
{
	bool isIn = false;

	for (int k = 0; !isIn && k < _listeArcPartage[i][j].size(); ++k)
	{
		isIn = _listeArcPartage[i][j][k]._orig == o && _listeArcPartage[i][j][k]._dest == d;
	}

	return isIn;
}