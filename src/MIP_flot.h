#pragma once

#include "InstanceReduite.h"
#include <queue>
#include <stack>


#include "util.h"

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

using NumVarArray2 = IloArray<IloNumVarArray>;
using NumVarArray3 = IloArray<NumVarArray2>;



//====================================================================================================
// ce MIP fait suite au mail de Christian du 9 mars 2021
// le but de ce MIP est de resoudre le probleme d'evacuation quand le graphe est un arbre
// en utilisant un MIP "flot RCPSP" + contrainte convexe "fin(j) = deb(j) + pop(j)/rate(j)"
// la ctr convexe est linearisee (par les tangentes) 
// et on ajoute les tangentes jusqu'a resolution exacte => Branch & Cut 
//=======================================================================================================

class MIP_flot
{

private:
	InstanceReduite * _ins;
	
	IloEnv * _env;
	IloModel * _mod;
	IloCplex * _cplex;

	//les indices de jobs vont de 0 (source) à nbJob+1 (puits)
	// ==> decalage d'un avec les donnees 

	NumVarArray3 _flot;  // _flot[i][j][k] : flot i -> j pour ressource k 
	NumVarArray2 _x; //_x[i][j] = 1 si i << j 
	IloNumVar _obj; //objectif = min des marges d'evacuation
	IloNumVarArray _start; //_start[i] = debut du job i (i = 0..N+1 où 0 = source, N+1 = puits)
	IloNumVarArray _end; //_end[i] = fin du job i (i = 0..N+1 où 0 = source, N+1 = puits)
	IloNumVarArray _rate; //_rate[i] = vitesse du job i (i = 0..N+1 où 0 = source, N+1 = puits)
	
public:

	long _nbCutLazy = 0; //nb coupes ajoutees dynamiquement (ie avec le callback lazy)
	long _nbCutUser = 0; //nb coupes ajoutees dynamiquement (ie avec le callback User)

	MIP_flot(InstanceReduite * ins)
	{
		_ins = ins;

		_env = new IloEnv();
		_mod = new IloModel(*_env);
		_cplex = new IloCplex(*_mod);
	}


	~MIP_flot()
	{

		//cout << "deleting PLNE..." << endl;
		delete _cplex;
		delete _mod;

		_env->end();
		delete _env;
	}


	//fonction principale : crée et résout le MIP
	//retourne borne inf, borne sup
	pair<double, double> creeEtResout(const vector<double> & debut, const vector<double> & debit, const vector<vector<vector<double>>> & flot);

	//retourne le nb de node explores par cplex
	long getExploredNode();



private:

	//fourni a cplex une solution initiale pour demarrer l'optimisation
	void warmStart1(const vector<double> & debut, const vector<double> & debit);

	//fourni a cplex une solution initiale pour demarrer l'optimisation
	void warmStart2(const vector<double> & debut, const vector<double> & debit);

	//fourni a cplex une solution initiale partielle pour demarrer l'optimisation : 
	void warmStart3(const vector<double> & debut, const vector<double> & debit, const vector<vector<vector<double>>> & flot);

	//fourni a cplex une solution initiale partielle pour demarrer l'optimisation : 
	void warmStart4(const vector<vector<vector<double>>> & flot);

	//allocation des variables
	void allocVar();

	//on cree les contraintes lineaires : toutes sauf "fin(j) = deb(j) + pop(j)/rate(j)"
	//BI_init :borne inf (obetnue par un heuristique par ex.), mettre -1 si pas de borne inf connue
	void creationCtrLineaire();

	// creation de la contrainte "fin(j) >= deb(j) + pop(j)/rate(j)"
	// ne fonctionne pas : cplex n accepte que des variables binaires dans les ctr quad convexe
	void creationCtrConvexe();

	//on cree nbPoints tangentes qui linearisent (partiellement)"fin(j) = deb(j) + pop(j)/rate(j)"
	//si uniforme = true alors on en crée nbPoint espacés régulièrement, sinon on en crée nbPoint non uniformément espacées
	//les autres tangentes seront générées dynamiquement (dans le lazy ctr callback)
	void creationCtrTangente(bool uniforme, int nbPoint);

	//idem creationCtrTangente mais on genere les points non uniformement 
	// avec une formule differente qui depend du carre du nb de points
	void creationCtrTangente_NonUniformeCarre(int nbPoint);

	//on cree la fonction objectif
	void creationFctObj();

	//affiche la sol. (suppose que _cplex->solve a ete appele !!)
	void afficheSol();

	//verif solution donnee par cplex
	bool verifSol();

	bool verifDebit(const vector<double> & dateDeb, const vector<double> & dateFin, const vector<double> & debit) const;
};
