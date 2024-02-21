#ifndef RCPSP_GRAPHE_H
#define RCPSP_GRAPHE_H


#include "Instance.h"
#include "util.h"


class RCPSP_Graphe
{
	friend class RCPSP_Algo;
	friend class PL_Optim_RCPSP;
	friend class PL_RCPSP_fixeDead;

private:
	Instance * _ins;

	int _nbSommet; //0 = source, 1.._nbSommet = sommets(= noeuds evac dans l'instance) , _nbSommet+1 = puits
	int _nbArc;

	vector<vector<vector<ArcCommun > > > _listeArcPartage; // _listeArcParrtage[i][j] donne la liste des arcs dans ins qui sont a la fois dans le chemin de i et dans le chemin de j
													// si vide alors i et j non voisin dans RCPSP_Graphe

	vector<vector<int> >  _TLcondFixe; //_TLcondFixe[i][j] : partie fixe (qui depend que de i et j) du time lag cond entre i et j : debut(j) >= debut(i) + TLcond[i][j] (où TLcond est obtenu avec getTLcond)


	vector<Arc> _arcs; //liste des arcs existant (dans au moins un chemin)
	vector<vector<int> > _arc2Ind; //_arc2Ind[i][j] donne l indice de l'arc (i,j) dans _arcs
	vector<vector<bool> > _utilise; //_utilise[i][a] = vrai si job i = 1..N utilise a = 0..nbArc

public:


	RCPSP_Graphe(Instance * ins)
	{
		_ins = ins;
		
		constuireGraphe();
	}


	//retourne le TL conditionnelle compose d'un partie fixe et d'une partie variable qui depend du debit de i
	double getTLcond(int i, int j, double debit)
	{
		if (debit > EPSILON)
			return (_TLcondFixe[i][j] + static_cast<double>(_ins->_pop[i])/ debit);
		else
			return _TLcondFixe[i][j];
	}

	//retourne le TL conditionnelle compose d'un partie fixe et d'une partie variable qui depend du debit de i
	//on arrondi pop/debit a l'entier superieur (utilise pour calculer les donnees a passer a CPO, necessite des donnees entieres)
	int getTLcondArrondi(int i, int j, int debit)
	{
		if (debit > EPSILON)
			return (_TLcondFixe[i][j] + static_cast<int> ( ceil ( static_cast<double>(_ins->_pop[i]) / debit ) + EPSILON ) );
		else
			return _TLcondFixe[i][j];
	}

private:


	void constuireGraphe();

	// renvoie vrai si l'arc (o, d) est dans la liste d'arc partages entre i et j
	bool isDansArcPartage(int i, int j , int o, int d);

};



#endif // RCPSP_GRAPHE_H


