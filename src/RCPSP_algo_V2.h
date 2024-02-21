#pragma once


#include "InstanceReduite.h"

#include "util.h"



//=============================================================
// janvier 2023
// cette classe est faite pour resoudre le pb avec l'heuristique flot (celle de RCPSP_algo) mais adapte pour des instances generales
// (pas forcement un arbre d'evacutaion) comme demande par les reviewers


//==================================================
//classe principale
//=================================================

class RCPSP_Algo_V2
{
	friend class PL_Optim_RCPSP;
	friend class RCPSPHeuristic;

public:
	struct Resultat
	{
		double _bestLambda = -1; //lambda (= marge) qui a donne la meilleure solution initiale (donc de cout lambda)
		double _margeFinale = -1; //marge obtenu apres optimisation (etape 2)
		vector<double> _dateDebut; //_dateDebut[i] apres etape 2 si _margeFinale >= 0, sol. initiale sinon
		vector<double> _debit; //_debit[i] apres etape 2 si _margeFinale >= 0, sol. initiale sinon
		double _marge0 = -1;//solution obtenue avec le PL "topo" si on part de la solution du flot pour lambda = 0
		
		//remarque :  on ne connait pas le flot qui donne la best marge car elle vient du PL topologie
	}; 

private:


	InstanceReduite * _ins;
	int _source;


	//====================================================================
	// SDD pour les donnees internes aux algo


	vector<vector<double> > _ressATransmettre;//_ressATransmettre[i][e] : ressource que peut trabsmettre i pour ressource e
public:
	vector<vector<vector<double> > > _flotCour;//_flotCour[i][j][e] : ressource de type e que le job i transmet au job j
private:
	vector<pair<int, int>> _arcSupport;//contient la liste des i->j qui portent du flot


	vector<double> _debutCour;
	vector<double> _debitCour;
	vector<double> _finCour;

	// ====  SDD pour GRASP_V2
	vector<bool> _inserted; //_inserted[i] = true si le job i est insere dans la sol partielle
	vector<bool> _dejaSupport; //_dejaSupport[i] = true si l'arc i -> jobCour existe dans _inserted

	vector<vector<double> > _ressATransmettre_PL;//_ressATransmettre[i][e] : ressource que peut trabsmettre i pour ressource e
	vector<vector<vector<double> > > _flotCour_PL;//_flotCour[i][j][e] : ressource de type e que le job i transmet au job j

	vector<double> _debutCour_PL;
	vector<double> _debitCour_PL;
	vector<double> _finCour_PL;

	//=====


public:


	RCPSP_Algo_V2( InstanceReduite * ins)
	{
		_ins = ins;
		_source = _ins->_nbJob;//pour eviter de decaler les numero de jobs on donne a source le num nbJob
	}


	//=====================================================================
	// algo principal

	//premier version de l'heuristique (on ameliore avec le PL "topologie" a la fin)
	Resultat GRASP(int maxRun, double precision);

	//seconde version de l'heuristique : on ameliore a chaque iteration avec le PL derive du MIP flot dans lequel on fixe les x_ij
	Resultat GRASP_V2(int maxRun, double precision);

	//troisieme version de l'heuristique : on ameliore a chaque iteration "a la main" puis a la fin avec le PL derive du MIP flot dans lequel on fixe les x_ij
	//si PL_iter = vrai alors on lance le MIP relaxe a chaque iteration, sinon juste a la fin de construitSol_graspV3
	Resultat GRASP_V3(int maxRun, double precision, bool PL_iter);

private:

	bool construireSigmaAlea(vector<int> & sigma, const vector<vector<int>> & coupleInterdit);
	bool construireSigmaProba(vector<int> & sigma, const vector<vector<int>> & coupleInterdit, double alpha);

	//procedure principale du GRASP (version 1) => construit une solution heuristique avec une marge lambda
	//procedure simple_insert du document Alain (quand on ajoute le nouveau job en fin du planing)
	tuple<bool, Resultat, pair<int, int>> insert(double lambda, const vector<int> & sigma);

	//procedure principale du GRASP_V2 : construit une solution heuristique avec une marge au moins lambda
	tuple<bool, Resultat, pair<int, int>> construitSol(double lambda, const vector<int> & sigma);

	// procedure principale du GRASP_V3 : construit une solution heuristique avec une marge au moins lambda
	// correspond a A_Insert-MSM-RCPSP du doc Alain
	//si PL_iter = vrai alors on lance le MIP relaxe a chaque iteration, sinon juste a la fin
	tuple<bool, RCPSP_Algo_V2::Resultat, pair<int, int>> construitSol_graspV3(double lambda, const vector<int> & sigma, bool PL_iter);

	// les jobs deja inseres transmettent du flot a jobCour (si possible) pour qu'il finissent avec une marge lambda
	// pour chaque ressource la quantite de flot fourni depend de la fin du job qui transmet 
	// donc la quantite de flot fournie peut etre differente d'une ressource a l'autre
	tuple<bool, pair<int, int>>  insert_etape1_fin(int jobCour, double popJCour, double LFJCour, double lambda,
		 const vector<pair<double, int>> & listePrec, vector<double> & flotCour, vector<double> & delta);

	//idem insert_etape1_fin mais on peut profiter des "trous" pour inserer un job
	tuple<bool, pair<int, int>, vector<double>>  insert_etape1_milieu(int jobCour, double popJCour, double LFJCour, double lambda,
		const vector<pair<double, int>> & listePrec, const vector<pair<double, int>> & listePrec_deb, vector<double> & flotCour, vector<double> & delta);

	// e_max est la ressource pour laquelle jobCour a recu le plus de flot dans l'etape 1
	// son debit est donc = a cette quantite de flot
	// pour les autres ressources qu'il utilise il doit recevoir la meme quantite
	// on fait donc transmettre du flot depuis les jobs deja inseres vers jobCour 
	// jobCour est insere en fin du planning (on ne detourne pas de flot)
	tuple<bool, pair<int, int>>  insert_etape2_fin(int e_max, int jobCour,
		 const vector<pair<double, int>> & listePrec, vector<double> & flotCour);

	// idem insert_etape2_fin sauf qu'on peut detourner du flot pour donner a jobCour : 
	tuple<bool, pair<int, int>>  insert_etape2_milieu(int e_max, int jobCour, double popCour,
		const vector<pair<double, int>> & listePrec, const vector<pair<double, int>> & listePrec_deb, vector<double> & flotCour);

	//init les SDD utilisee dans insert
	void initInsert();

	//init les SDD utilisee dans construitSol
	void initConstruitSol();


	//======== optimisation post flots (PL "topologie)
	
	double reconstruireSolCourantePL_exact(const vector<pair<int, int> > & precedence, const vector<vector<int> > & cliques);
	void calculPredClique(vector <pair<int, int> > & precedence, vector <vector<int> > & clique);

	//======== optimisation a chaque iteration (PL extrait du MIP flot de Christian en fixant les x_ij => devient un PL)

	double PL_flot_partiel( );


	void conserveSolPL();

	//supprime de la liste des arc support l'arc prec -> jobCour (prec != source) avec la date de fin la plus grande pour prec
	//renvoie vrai si on a trouve un arc a supprimer
	bool supprimeSupportfinMaxOrigine(int jobCour);



	//remplit _arcSupport en utilisant le flot courant
	void recupereArcSupport();

	//============= verification

	//verifie kirchoff sur les flots et coherence avec le debit
	bool isFlotOK();
};


//remet chaque element de v a 0
void reset(vector<vector<int>> & v);

//supprime les elements avec une proba = alpha
void menageAlea(vector<vector<int>> & v, double alpha);