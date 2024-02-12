#ifndef PL_OPTIM_RCPSP_H
#define PL_OPTIM_RCPSP_H



#include "Instance.h"
#include "RCPSP_Algo.h"

#include "util.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using NumVarArray2 = IloArray<IloNumVarArray>;

using NumVarArray3 = IloArray<NumVarArray2>;


class PL_Optim_RCPSP
{

private:
	IloEnv * _env;
	IloModel * _mod;
	IloCplex * _cplex;

	Instance * _ins;
	RCPSP_Algo * _rcpspAlgo;
	RCPSP_Graphe * _graph;

	int _nbSommet;
	int _nbArc;

	NumVarArray3 _deltaFlot;  // _flot[i][j][e] : flot transmis du job i au job j pour la ressource "arc e" (i et j de 1 à nbSommet)
	IloNumVarArray _deltaDebit; // indexe de 1 a _nbSommet


public:

	PL_Optim_RCPSP(Instance * ins, RCPSP_Algo * rcpspAlgo, RCPSP_Graphe * graph)
	{
		_ins = ins;
		_rcpspAlgo = rcpspAlgo;
		_graph = graph;
		_nbSommet = _graph->_nbSommet;
		_nbArc = static_cast<int>(_graph->_arcs.size());

		_env = new IloEnv();
		_mod = new IloModel(*_env);
		_cplex = new IloCplex(*_mod);
	}


	~PL_Optim_RCPSP()
	{
		//cout << "deleting PLNE..." << endl;
		delete _cplex;
		delete _mod;

		_env->end();
		delete _env;
	}



	//cree le modele et resout
	bool creeEtResout(const vector<double> & debitMin, const vector<pair<int, int> > & arcActif,  const vector<double> & lambda);
	
	void afficheSol();


	double getDeltaDebit(int i)
	{
		double res = 0;

		if (_cplex->isExtracted(_deltaDebit[i]))
			res = _cplex->getValue(_deltaDebit[i]);

		return res;
	}

	double getDeltaFlot(int i, int j, int e)
	{
		double res = 0;

		if (_cplex->isExtracted(_deltaFlot[i][j][e]))
			res = _cplex->getValue(_deltaFlot[i][j][e]);

		return res;
	}

private:

	//on cree le modele sachant que les seuls arcs actifs sont ceux stockes dans arcActif : 
	//le flot sera nul pour les autres arcs (les variables sont crees par commodite mais non incluses dans les ctr)
	//debitMin : debit en dessous duquel il ne faut pas descendre 
	void creationModele(const vector<double> & debitMin, const vector<pair<int, int> > & arcActif, const vector<double> & lambda);

	//allocation des variables
	void allocVar(const vector<pair<int, int> > & arcActif);

};



#endif // 

