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
// ce PL est une version relaxe du MIP flot (MIP_flot.h)
// ici on fixe les x_ij (il n'y a donc plus de variable entiere)
// et on travaille avec uniquement un sous-ensemble des jobs 
// => on obtient un planning partiel
// ce PL sert dans l'heuristique flot (ou il est appele de maniere iterative)
//=======================================================================================================

class MIP_flot_partiel
{

private:
	InstanceReduite * _ins;
	vector<bool> _existe; //_existe[i] = true si le job i de _ins doit etre considere (sinon les variables associees sont creees - pour ne pas introduire de decalage dans les indices - mais non utilisee)

	IloEnv * _env;
	IloModel * _mod;
	IloCplex * _cplex;

	//les indices de jobs vont de 0 (source) à nbJob+1 (puits)
	// ==> decalage d'un avec les donnees 

	NumVarArray3 _flot;  // _flot[i][j][k] : flot i -> j pour ressource k 
	IloNumVar _obj; //objectif = min des marges d'evacuation
	IloNumVarArray _start; //_start[i] = debut du job i (i = 0..N+1 où 0 = source, N+1 = puits)
	IloNumVarArray _end; //_end[i] = fin du job i (i = 0..N+1 où 0 = source, N+1 = puits)
	IloNumVarArray _rate; //_rate[i] = vitesse du job i (i = 0..N+1 où 0 = source, N+1 = puits)

public:


	MIP_flot_partiel(InstanceReduite * ins, const vector<bool> & existe)
	{
		_ins = ins;
		_existe = existe;

		_env = new IloEnv();
		_mod = new IloModel(*_env);
		_cplex = new IloCplex(*_mod);
	}


	~MIP_flot_partiel()
	{

		//cout << "deleting PLNE..." << endl;
		delete _cplex;
		delete _mod;

		_env->end();
		delete _env;
	}


	//fonction principale : crée et résout le PL
	// seul les i->j dans arcSupport peuvent porter du flot (utilise pour definir ctr de precedence)
	double creeEtResout(const vector<pair<int,int>> & arcSupport);

	double getDebut(int i) { return _cplex->getValue(_start[i]); }
	double getFin(int i) { return _cplex->getValue(_end[i]); }
	double getDebit(int i) { return _cplex->getValue(_rate[i]); }
	double getFlot(int i, int j, int k) 
	{
		double f = 0;
		if (_cplex->isExtracted(_flot[i][j][k])) 
			f = _cplex->getValue(_flot[i][j][k]);
		return f;
	}

private:



	//allocation des variables
	void allocVar();

	//on cree les contraintes lineaires : toutes sauf "fin(j) = deb(j) + pop(j)/rate(j)"
	//BI_init :borne inf (obetnue par un heuristique par ex.), mettre -1 si pas de borne inf connue
	void creationCtrLineaire(const vector<pair<int, int>> & arcSupport);

	//on cree nbPoints tangentes qui linearisent (partiellement)"fin(j) = deb(j) + pop(j)/rate(j)"
	//si uniforme = true alors on en crée nbPoint espacés régulièrement, sinon on en crée nbPoint non uniformément espacées
	//les autres tangentes seront générées dynamiquement (dans le lazy ctr callback)
	void creationCtrTangente(bool uniforme, int nbPoint);


	//on cree la fonction objectif
	void creationFctObj();

	//affiche la sol. (suppose que _cplex->solve a ete appele !!)
	void afficheSol();

	//verif solution donnee par cplex
	bool verifSol();

	bool verifDebit(const vector<double> & dateDeb, const vector<double> & dateFin, const vector<double> & debit) const;
};
