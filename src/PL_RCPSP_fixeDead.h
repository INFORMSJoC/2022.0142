#ifndef PL_RCPSP_FIXEDEAD_H
#define PL_RCPSP_FIXEDEAD_H



#include "Instance.h"
#include "RCPSP_Graphe.h"

#include "util.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using NumVarArray2 = IloArray<IloNumVarArray>;

using NumVarArray3 = IloArray<NumVarArray2>;

//=====================================================================
// 
// ce PL s'utilise sur une instance du probleme d'evacuation avec chemins fixes
// il sert a fixer des dealines aux differents chemins dans le but de creer des intances non triviales 
// pour tester la borne sup (RCPSP_Algo)


class PL_RCPSP_fixeDead
{

public:
	IloEnv * _env;
	IloModel * _mod;
	IloCplex * _cplex;

	Instance * _ins;
	RCPSP_Graphe * _graph;

	int _nbSommet;
	int _nbArc;

	vector<double> _debut; //date de debut des actions , on les calcul une fois qu on a le flot
	vector<double> _fin; //date de fin des actions, on les calcul une fois qu on a le flot

	NumVarArray3 _flot;  // _flot[i][j][e] : flot transmis du job i au job j pour la ressource "arc e" (i et j de 1 à nbSommet)
	IloNumVarArray _debit; // indexe de 1 a _nbSommet


public:

	PL_RCPSP_fixeDead(Instance * ins, RCPSP_Graphe * graph)
	{
		_ins = ins;
		_graph = graph;
		_nbSommet = _graph->_nbSommet;
		_nbArc = static_cast<int>(_graph->_arcs.size());

		_debut.resize(_nbSommet + 1);
		_fin.resize(_nbSommet + 1);

		_env = new IloEnv();
		_mod = new IloModel(*_env);
		_cplex = new IloCplex(*_mod);
	}


	~PL_RCPSP_fixeDead()
	{
		cout << "deleting PLNE..." << endl;
		delete _cplex;
		delete _mod;

		_env->end();
		delete _env;
	}



	//cree le modele et resout
	bool creeEtResout( );

	void afficheSol( );

	//dessine le graphe uniquement avec les arcs dont le flot pour la ressource arc e est non nulle
	void traceFlotArc( int e );

	void calculeDebutEtFin( );
	
	double getDebit(int i)
	{
		return _cplex->getValue(_debit[i]);
	}

private:


	//flot autorise de i vers j uniquement si i << j dans ordre
	//ordre ne doit pas contenir source et puits
	void creationModele(vector<int> & ordre);

	//allocation des variables
	void allocVar();

	//renvoie vrai si flot de i vers j pour au moins une ressource
	bool existeFlot(int i, int j);

};



#endif // 

