#ifndef MODELELINEAIRE_H
#define MODELELINEAIRE_H



#include "Instance.h"
#include "ExpandedNetwork.h"

#include "util.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using NumVarArray2 = IloArray<IloNumVarArray>;



class ModeleLineaire
{

private:
	IloEnv * _env;
	IloModel * _mod;
	IloCplex * _cplex;

	Instance * _ins;
	ExpandedNetwork * _expNet;

	IloNumVarArray _flot;  // _flot[i] : flot sur l arc i du ExpandedNetwork



public:


	ModeleLineaire(Instance * ins, ExpandedNetwork * expNet)
	{
		_ins = ins;
		_expNet = expNet;

		_env = new IloEnv();
		_mod = new IloModel(*_env);
		_cplex = new IloCplex(*_mod);
	}


	~ModeleLineaire()
	{

		cout << "deleting PLNE..." << endl;
		delete _cplex;
		delete _mod;

		_env->end();
		delete _env;
	}



	//creele modele avec fct obj si obj = true, sans sinon
	double creeEtResout(bool obj);




private:

	
	void creationModele( bool obj);

	//allocation des variables
	void allocVar();


	// si on veut mettre une fonction objectif on met obj a true,
	// si obj = false c juste un pb de faisabilite
	void creationCtr_Obj(bool obj);

	//on boucle jusqu a obtenir la solution qui permet d avoir les marges les plus grandes possibles : 
	//on tue iterativement les arcs qui minimise dead - t
	double resolutionIterative( );

	//on interdit tous les arcs transverses de la forme (i,t+k) -> (j,t+k+dist_ij), k >= 0
	void interdireArcs(int i, int j, int t);

	double resolutionSimple();
};



#endif // MODELELINEAIRE_H
