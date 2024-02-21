#pragma once


#include "Solution.h"

#include "InstanceReduite.h"
#include <queue>
#include <stack>

//=========================================================================
//  cette classe resout le problème d'evacuation non preemptif 
// avec un branch & bound fait "à la main"
//=========================================================================

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN


// type de contraintes pour le "PL_exact" : 
// TANGENTE : la ctr non lineaire est approximée par ses tangentes
// SECANTE : la ctr non lineaire est approximée par ses secante
// CUT : on part des tangente et on ajoute des coupes dynamiquement
enum TYPE_RELAX { TANGENTE, SECANTE, CUT };

class BranchBoundDicho
{
	Instance * _insInit; //instance intiale (juste utile pour faire des verifications)
	InstanceReduite * _ins;

	struct SignatureBranchement
	{
		int _type; //permet de connaitre le branchement effectue 0, 1 ou 2
		int _job = -1; //job surlequel on a branche
		double _alpha; //point de branchement 1
		double _beta; //point de branchment 2 (alpha < beta)
	};

	class Noeud
	{
	public:
		static unsigned int _nbNoeudTotal;//nb total de noeud crees (utile pour debug uniquement)
		unsigned int _id = 0; //identifiant du noeud (utile pour debug uniquement)
		unsigned int _idPere = -1; //identifiant du pere noeud (utile pour debug uniquement)

		//les vector suivant sont indexes sur 0...nbJob-1 et rel[j] <= _pointMil_1[j] < _pointMil_2[j] <= due[j]
		vector<double> _rel;
		vector<double> _due;
		vector<double> _pointMil_1; //_pointMil_1[i] = -1 si non defini
		vector<double> _pointMil_2; //_pointMil_2[i] = -1 si non defini

		vector<double> _LAdd; //liste de points additionnels
		double _marge; //meilleure marge qu'on peut espérer pour ce noeud (= borne sup)
		double _priorite; //priorite du noeud (plus grande priorite = plus prioritaire si on a une file a priorite)

		SignatureBranchement _sign; //contient les informations pour le branchement
		
		//un noeud est plus prioritaire si sa marge est plus grande
		bool operator < (const Noeud  & n2) const
		{
			return this->_priorite < n2._priorite;
		}
		

#ifdef _VERIF_BRANCHEMENT_
		vector<vector<double>> _sol; //_solCour[i][k] = debit du job i entre les dates _date[k] et _date[k+1] (on ajoute artificiellement 0 à la fin pour éviter les cas particuliers dans la recherche des param. de branchement)
		vector<double> _date; //vecteur de dates associes a solCour
#endif
	};



	//arbre du branch & bound => file à priorité pour un parcours en "Best first"
	
#ifdef _PARCOURS_PROF_
	stack<Noeud> _listeNoeud;
#else
	priority_queue<Noeud> _listeNoeud;
#endif

	//==================================================================
	// borne inf (= solution realisable car on maximise) 

	double _bestMargeRealisable; //meilleure marge (associee a une solution realisable) trouvee jusqu'a present = Borne inf (on maximise)
	double _coutRacine = -1; //cout de la solution a la racine (= cout de la sol. preemptive)


	//----------------------------------------------
	// "vraie" SOLUTION : dates de debut et debit
public:
	vector<double> _bestDebit;
	vector<double> _bestDebut; 
	bool _hasSol = false;
	//---------------------------------------------

	//====================================================================
	// donnees mises a jour apres l'evaluation d'un noeud (procedure eval)

private:
	vector<vector<double>> _solCour; //_solCour[i][k] = debit du job i entre les dates _date[k] et _date[k+1] (on ajoute artificiellement 0 à la fin pour éviter les cas particuliers dans la recherche des param. de branchement)
	vector<double> _dateCour; //vecteur de dates associes a solCour


	//=====================================================================
	// environnement cplex

	//on ne factorise pas l'environnement Cplex, il est beaucoup plus rapide (facteur 30)
	// de recréer et supprimer l'environnement à chaque fois => très bizarre ??

	//====================================================================
	// donnees extraites du noeud courant pour construire le PL (extraites par pretraitementPL())

	vector<double> _due_2; // = min(due, ins.LF-margeCour)
	vector<double> _date; //toutes les dates ordonnees en ordre croissant (_rel du noeud courant U _due_2 U _LAdd du noeud courant)
	vector<int> _indRel; //index sur 0...nbJob-1 : _indRel[j] = indice dans date de la date _rel[j] du  noeud courant
	vector<int> _indDue_2; //index sur 0...nbJob-1 : _indDue_2[j] = indice dans date de la date _due_2[j]
	vector<int> _indMil_1; //index sur 0...nbJob-1 : _indMil_1[j] = indice dans date tq _pointMil_1[j] in [date[k],date[k+1][ (= date.size() si n'existe pas)
	vector<int> _indMil_2; //index sur 0...nbJob-1 : _indMil_2[j] = indice dans date tq _pointMil_2[j] in ]date[k],date[k+1]] (= 0 si n'existe pas)

	//=========================================================================
	//infos sur l'arbre de recherche


	//nombre de fois que la branche 0 (repst. 1 et 2) n'a pas strictement diminue le coût de la solution pere
	unsigned int _nonAmelioreb0;
	unsigned int _nonAmelioreb1;
	unsigned int _nonAmelioreb2;

	unsigned int _nbNoeudExplore;//nb noeuds explores = nb appels a evalNoeud
	unsigned int _nbBranchType1;
	unsigned int _nbBranchType2ou3;

	//==========================================================================
	//infos sur le PL qui recalcule la meilleure solution quand l'ordre des jobs est connu
	unsigned int _nbCut = 0; // nombre d'iteration de la boucle qui ajoute les cuts (plusieurs cuts peuvent etre ajoutees a chaque tour de boucle) pour tous les appels a reconstruireSolCourantePL_exact
	int _nbAppelCut = 0; // nb appels a reconstruireSolCourantePL_exact => nbCut / nbAppelCut donne le nb moyen de cuts

	//===============================================================================
	//infos sur les heuristiques de reconstruction / sterilisation (solution du dernier noeud evalue)
public:
	double _steril1 = -1; //vrai si reconstruireSolutionCourante sterilise le noeud
	double _steril2 = -1; //vrai si reconstruireSolutionCouranteAmelioree sterilise le noeud
	double _res_DebitFixe = 999;
	double _res_heurCut = 999;
	double _res_tangente = 999;
	double _res_secante = 999;
	double _res_cut = 999;


public:
	//constructeur : BI est une marge realisable obtenue par une heuristique
	BranchBoundDicho(Instance * insInit, InstanceReduite * ins, double BI) : _insInit(insInit), _ins(ins), _bestMargeRealisable(BI) 
	{ init(); }

	//méthode principale : effectue le B&B et renvoie le coût de la meilleure solution (_bestMargeRealisable = BI) et de la borne sup
	// cpuMax = temps max en sec
	pair<double, double> run(double cpuMax);


	//getter pour recuperer les infos sur l'arbre
	unsigned int getNbNoeud() { return Noeud::_nbNoeudTotal; }
	unsigned int getNonAmelioreb0() { return _nonAmelioreb0; }
	unsigned int getNonAmelioreb1() { return _nonAmelioreb1; }
	unsigned int getNonAmelioreb2() { return _nonAmelioreb2; }
	unsigned int getNbNoeudExplore() { return _nbNoeudExplore; }
	unsigned int getNbBranchType1() { return _nbBranchType1; }
	unsigned int getnbBranchType2ou3() { return _nbBranchType2ou3; }
	unsigned int getNbCut() { return _nbCut; }
	int getNbAppelCut() { return _nbAppelCut; }


	//ecrit dans un fichier les sol des differentes heuristique au dernier noeud explore
	void ecrireHeurNoeud(const string & nomInst);


	//verifie que la sol best (bestDebit ; bestDebut) est ok : 
	//on verifie les capa, les marges et que la marge min = marge
	bool verifieSolBest(const double marge) const;

	//on calcule le taux de parallélisme sur le dernier arc pour la best solution
//renvoie le taux de parallelisme moyen sur le dernier arc et le taux max
	pair<double, int> statParallelisme();

private:

	//allocation / initialisation des attributs de classe
	void init();


	//evalue un noeud avec une procedure dicho basee sur un PL
	double evalNoeud(const Noeud & noeud, double marge);

	//calcule les donnees extraites du noeud  pour construire le PL
	void pretraitementPL(const Noeud & noeud, double marge);

	//on utilise un PL pour trouver une solution preemptive au noeud courant
	//retourne true si PL faisable
	bool resolutionPL(const Noeud & noeud, double marge);


	//modifie le noeud fils en fonction de sa signature et de b (numero de branche)
	void appliqueSignatureEtBranch(Noeud & fils, int b);

	//=================================================================================================
	//  fonctions pour reconstruire une solution premptive en solution non preemptive


	//calcul les infos liees a la solution courante : debut, fin, duree et debit moyen
	// si reduireDebit = true alors on  reduit les debits de 1% (et on allonge la duree en consequence) : 
	//cela pour eviter les problemes lies a la precision numerique : Cplex met parfois des debits qui violent "legerement" les ctr
	void calculInfoJob(vector<int> & i_deb, vector<int> & i_fin, vector<double> & duree, vector<double> & debit, bool reduireDebit);
	
	
	//calcul les jobs qui forment une clique (overlap) ou qui sont l'un avant l'autre en utilisant les debut et fin des jobs de la solution courante
	void calculPredClique(const vector <int> & i_deb, const vector<int> & i_fin, vector <pair<int, int> > & precedence, vector <vector<int> > & cliques);

	//a partir de la solution courante _solCour preemptve (obtenue par le PL) on tente de reconstruire une solution non preemptive
	//en faisant la moyenne des debits
	//si stocke =true alors on stocke la solution recalculee comme nouvelle meilleure solution 
	//on passe la marge correspond a la sol. courante (utile pour faire des stats uniquement)
	bool reconstruireSolutionCourante(bool stocke, const double margePere);

	//idem mais en modifiant les dates de debut et fin pour les jobs qui commencent avec un  petit débit ou finissent avec un petit débit
	//retourne vrai si la solution a ete modifiee par rapport a solMoy et qu'elle est realisable
	bool reconstruireSolCouranteAmeliore(const vector<int> & i_deb, const vector<int> & i_fin, const vector<vector<double>> & solMoy, double marge, bool stocke);

	//idem en utilisant un PL basé sur l'ordre des jobs donne par la solution preemptive
	//resolution heuristique
	double reconstruireSolCourantePL_heur(const vector<int> & i_deb, const vector<int> & i_fin, const vector<double> & duree, const vector<double> & debit);


	//idem reconstruireSolCourantePL_heur : utilise  l'ordre des jobs donne par la solution preemptive pour constrire une solution
	//mais on resout de maniere "quasi" exact en approximant la ctr non lineaire par ses tangentes (version optimiste) ou par ses secantes (version pessimiste)
	// ou avec une genration de coupes (version exacte mais risque de degenerescence)
	//precedence donne la liste des couples (i,j) avec i << j
	//cliques donne l'ensemble des cliques
	//pasConstant = true si o nutilise un pas constant, false sinon

	double reconstruireSolCourantePL_exact(const vector<pair<int, int> > & precedence, const vector<vector<int> > & cliques,
		const TYPE_RELAX & relax, bool pasConstant);

	//utilise les debits moyens courants pour essayer de reconstruire une solution avec un algo
	// glouton 
	//retourne la marge obtenue (-1 si echec)
	double reconstruireSolHeuristique(const vector<double> & debit, const vector<double> & longueur);

	//on choisit un job parmi l'ensmeble jobs (utilise dans reconstruireSolHeuristique)
	int choisir(double tCour, const vector<int> & jobs, const vector<double> & longueur,
		const vector<bool> & marq, const vector<double> & debit);

	//verifie que les debits sont coherents avec la capacité des arcs
	bool verifDebit(const vector<double> & dateDeb, const vector<double> & dateFin, const vector<double> & debit) const;

	//==============================================================
	//sous fonctions pour le calcul des paramètres de branchmeent

	//calcule la signature de branchement pour la solution courante associee au vecteur date courant
	SignatureBranchement calculParamBranchement();


	//calcul le creux le plus profond pour le job i s'il y existe, renvoie sa profondeur (-1 si n'existe pas)
	// et met a jour sign en consequence
	double calculCreux(int i, SignatureBranchement & sign);

	//calcul le branchement le plus large possible (on cherche le 1er et le dernier k | sol[k]!=0, si la solution n'est pas constante 
	//entre k1+1 et k2 on utilise les alpha beta associes pour faire un branchmeent de type 1
	// et met a jour sign en consequence
	//contient le branhcmeent creux dans le sens ou si creuxAmeliore n'existe pas alors creux non plus
	double calculCreuxAmeliore(int i, BranchBoundDicho::SignatureBranchement & sign);

	//calcul l'escalier  le plus haut pour le job i s'il y existe, renvoie sa hauteur (-1 si n'existe pas)
	// et met a jour sign en consequence
	double calculEscalierPositif(int i, SignatureBranchement & sign);

	//calcul l'escalier  le plus haut pour le job i s'il y existe, renvoie sa hauteur (-1 si n'existe pas)
	// et met a jour sign en consequence
	//rq : fonction appelee uniquement si pas de creux => l'escalier calcule est faux si la solution contient des creux
	double calculEscalierNegatif(int i, SignatureBranchement & sign);

	//on verifie que si on branche le job i suivant les indices k1 et k2 alors les 3 branches resultantes coupent la solution courante
	void verifieBranchementValide(int i, int k1, int k2);

	//calcul les paramètres de branchmeent pour les jobs de type ++- ou +--
	//on utilise les vecteurs sol et date en parametre dans lesquels les valeurs de debit identiques et consecutives ont ete regroupees
	double calculBranchType2ou3(int i, const vector<double> & sol, const vector<double> & date, BranchBoundDicho::SignatureBranchement & sign);

	//on regroupe les valeurs identiques consecutives dans _solCour[i] (evite de devoir traiter les cas = dans la recherche des param de branchement)
	void regroupe(int i, vector<double> & solGroup, vector<double> & dateGroup);

	//affiche la solution courante
	void afficheSolCour();

	//retourne le nombre de jobs non preemtif dans la solution courantes
	int nbJobResolu();

	//retourne la taille du plus petit intervalle dans dateCour
	double intervalMin();


};

