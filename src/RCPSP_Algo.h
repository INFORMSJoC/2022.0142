#ifndef RCPSP_ALGO_H
#define RCPSP_ALGO_H



#include "Instance.h"

#include "util.h"

#include "RCPSP_Graphe.h"


//=============================================================
// RCPSPHeuristic : classe utilisee uniquement pour passer la solution de RCPSP_Algo au code
// de Christian / Emmanuel qui utilise CPO
//=============================================================

class RCPSPHeuristic
{
public:

	bool isSolComputed = false;

	vector<int> s; // starting time (indice de 0 à nbNodeEvac-1) - entier
	vector<int> h; // debit (indice de 0 à nbNodeEvac-1) - entier
	vector<int> p; // duree de l'evacuation = val. entiere sup. de population/debit (indice de 0 à nbNodeEvac-1) - entier
	
	//constructeur : lit le fichier en entree, appelle l'heuristique de RCPSP_Algo et remplit les vecteurs s,h,p avec
	//la solution trouvee par RCPSP_Algo arrondie (donc degradee) (CPO ne prend que des entiers)
	RCPSPHeuristic(string nomFic);

};


//==================================================
//classe principale
//=================================================

class RCPSP_Algo
{
	friend class PL_Optim_RCPSP;
	friend class RCPSPHeuristic;
public:
	struct Resultat
	{
		double _solInit = -1;
		double _solFinale = -1;
		double _bestAlpha = -1; //alpha qui a donne la meilleure solution
		int _nbChangement = -1;//nb d'echanges effectues dans l ordre initial pour la best sol
	};

private:

	RCPSP_Graphe * _graph;
	Instance * _ins;

	// on stocke en pretraitement les "vraies" fins max des chemins obtenues avec les deadlines
	// (on s'en sert pour calculer les marges dans le PL PL_Optim_RCPSP)
	// en effet les LF des chemins ont ete modifies dans l'instance pour prendre en compte les marges
	// donnees par la BS
	// 
	vector<int> _finMaxInit; 

	int _nbArc;
	int _puits;

	//====================================================================
	// SDD pour les donnees internes aux algo


	vector<vector<vector<double> > > _flot; //_flot[i][j][e] = quantite de ressource e (e = un arc) transmise de i vers j
											// dans _flot[i][j][arc.size] on stocke la somme de flot de i vers j sur tout e

	vector<vector<double> >_ressATransmettre;//_ressATransmettre[i][e] : ressource que peut trabsmettre i pour ressource e


	//====================================================================
	// SDD pour la solution (on pourrait les mettre dans une classe a part)
public:
	vector<double> _dateDebut; //_dateDebut[i] = debuy de l evacuation du node evac i (0 = source, i = 1.. nbNodeEvac, puits)
	vector<double> _debit; //_debit[i] = _debit de l evacuation du node evac i (0 = source, i = 1.. nbNodeEvac, puits)

private:
	//====================================================================
	// SDD pour la solution qu'on donne a CPO (programme Christian / Emmanuel)
	// attention 1. dans CPO il n'y a pas la notion de source, donc l'indice 0 correspond au premier node evac
	//           2. CPO considere des variables entieres => il faut dégrader la solution obtenue par cette classe pour 
	//              avoir des variables s,h et p entieres
	vector<int> s; // starting time (indice de 0 à nbNodeEvac-1) - entier
	vector<int> h; // debit (indice de 0 à nbNodeEvac-1) - entier
	vector<int> p; // duree de l'evacuation = val. entiere sup. de population/debit (indice de 0 à nbNodeEvac-1) - entier

	//REMARQUE : calcul de la date de debut de i, i = node evac :
	// si i et j partage une ressource (= un arc) : debut(j) >= debut(i) + TLcond[i][j] (où TLcond est obtenu avec getTLcond)

public:


	RCPSP_Algo(RCPSP_Graphe * graph, Instance * ins)
	{
		_graph = graph;
		_ins = ins;
		_nbArc = static_cast<int>(_graph->_arcs.size());
		_puits = _graph->_nbSommet + 1;
	}


	//=====================================================================
	// algo principal

	//appelle plusieurs fois la resolution (= init une solution + decscente avec le PL)
	//en randomisant l'ordre initial
	// maxReboot = nb de fois ou on recommence en enlevant 1 aux deadline pour essayer de generer des solutions initiales avec des marges plus petites
	// maxRun = nb de fois max ou on execute l'algo (init + descente) avec des deadlines données
	// maxIterInit = nb de tentatives maximal pour obtenir une solution initiale
	// rang(i) donne le nb de job qui sont avant i (ie. debut < debut(i) et fin < fin(i), sert a calculer un ordre
	Resultat GRASP(vector<double> & finOptimiste,  int maxRun, int maxIterInit);
	
	//idem GRASP mais en initialisant les donnees differement
	Resultat GRASPstrategie2(double borneSup, int maxRun, int maxIterInit);


	// initialise (si possible) une solution grace a une heuristique puis essaie de l'ameliorer
	//renvoie le cout de la solution si on a reussi a construire une solution, -1 sinon
	pair<double, double> resolution(vector<int> & ordre, int maxIterInit, int & nbChangement);

	// idem resolution mais on utilise le PL qui sert dans le "Branch&Bound Dicho" au moment ou on essaie de reconstruire
	// une solution non preemptive en conservant l'ordre des actions
	pair<double, double> resolution_avec_PL_cut(vector<int> & ordre, int maxIterInit, int & nbChangement);

	//calcul les jobs qui forment une clique (overlap) ou qui sont l'un avant l'autre en utilisant les dates de debut et debit des jobs de la solution courante
	void calculPredClique(vector <pair<int, int> > & precedence, vector <vector<int> > & clique);

	//ce PL est le même que BranchBoundDicho::reconstruireSolCourantePL_exact mais a ete adapte aux SDD du RCPSP_Algo.
	// le but est, a partir d'une solution obtenue heuristiquement et pour laquelle on a extrait les precedences et les cliques
	// au niveau du dernier noeud de l'arbre d'evacuation ("safe node") 
	// de construire une solution (i.e. date d'arrivee au safe node et debit) qui respecte toutes les contraintes avec en plus
	// les memes precedence que la solution heuristique (imposer les precedences permet d'eviter les cycles)
	//retourne la marge minimale 
	//attention ce PL ne fonctionne que sur un graphe en forme d'arbre
	double reconstruireSolCourantePL_exact(const vector<pair<int, int> > & precedence, const vector<vector<int> > & cliques);

	//transforme la solution courante (dans _flot, _dateDebut, _debit) pour CPO 
	//=> arrondir les debit a l'entier inf + recalculer les starting times avec le flot sachant que les 
	//durees d'evacuation = pop / debit sont arrondies a l'entier superieur
	//le resultat est stocke dans les vecteurs 
	void transformePourCPO();

	//on recalcule les marges en prenant en compte les dates de debut recalculees (s) avec les debits arrondis a l'entier inf (h)
	// on retourne la marge min
	int calculeMargeCPO();


	//===============================================================
	// fonction pour apres la resolution

	//affiche la solution avec le detail des flots
	void affiche();

	//affiche quel action consomme quelle qte de quelle ressource et sur quelle periode
	void afficheConsoRess();

	//verifie que la solution respecte les deadline , le tmax, les capacites
	bool verifie();

	//dessine chaque chemin de manière independante 
	void traceCheminIndependant();

	//trace les flots pour l'arc e donne en parametre (dessine avec graphviz)
	void traceFlotArc(int e);

	//renvoie le debit de i sur l'arc e ( recalcule avec les valeurs de flot)
	double calculConsoFromFlot(int i, int e);

	//on calcule le taux de parallélisme sur le dernier arc pour la best solution
//renvoie le taux de parallelisme moyen sur le dernier arc et le taux max
	pair<double, int> statParallelisme();


private:


	//===========================================================================
	//sous-fonction algo general resolution 


	//construction initiale appelle constructionGloutonne et gere les echecs
	// en effectuant des modifs dans la liste "ordre" 
	//ordre : vecteur initialise dans la fonction, on le retourne car on en a besoin pour le PL
	//d'amelioration des debits
	bool calculSolInit(int maxIter, int & nbChangement, vector<int> & ordre);

	//essaie de construire un planing (date de debut + debit) a l aide d un flot de sorte
	//que les actions recoivent du flot en suivant "ordre"
	//renvoie vrai si construction ok
	//faux sinon et place dans  (actionEchec,  actionCour) les actions responsables de l'echec
	bool constructionGloutonne(vector<int> & ordre, int & actionEchec, int & actionCour);

	//mise a jour des dates de debut en fonction des valeurs de flots :
	//si deux jobs i et j se transmettent du flot(i->j)
	//alors debut(j) >= debut(i) + TLcond(i,j,debit(i))
	void calculeDebut();

	//mise a jour des dates de debut en fonction des valeurs de flots :
	//si deux jobs i et j se transmettent du flot(i->j)
	//alors debut(j) >= debut(i) + TLcond(i,j,debit(i))
	void calculeDebutPourCPO();




	//===========================================================================
	// sous fonctions pour la partie optimisation dans resolution
	
	//retourne la liste des arcs actifs (i.e. qui portent du flot dans la sol courante)
	vector<pair<int, int>> getArcActif();

	//retourne tous les arcs (i,j) possible en suivant l'ordre tels que l'ajout de flot i -> j ne décale pas trop
	//le debut de j
	vector<pair<int, int>> getArcOrdre(vector<int>& ordre);

	//calcul les coeff lambda pour le PL PL_Optim_RCPSP
	void calculLambda(vector<double> & lambda);

	//calcul les coeff lambda pour le PL PL_Optim_RCPSP dans le cas ou on veut maximiser le min des magres
	void calculLambdaPourMin(vector<double> & lambda, vector<int> & xmin, vector<int> & xmin2);

	//===========================================================================
	// sous fonctions pour constructionGloutonne


	//alloue les vecteurs
	void alloc();

	// essaie de placer somCour
	//retourne -1 si succes, sinon retourne un numero de somemet
	int assigner(vector<int> & ordre, int somCour);


	//pour chaque arc dans le chemin de somCour, on regarde si on peut l'alimenter en ressource "arc"
	//de sorte que la fin de l'action d'evacuation de somCour ne depasse pas la date limite (LF(safenode))
	//si echec renvoie un numero de noeud responsable de l echec, 
	//si succes renvoie -1 et calcule 1) le debit _debit[somCour] 2) l'arc arcDebit qui donne ce debit 3) init les flots vers somCour

	int assignerPhase1(vector<int> & ordre, int somCour, int & arcDebit);


	//construit la liste des sommets utilisee dans assignerPhase1 : 
	// les sommets x avant somCour dans ordre, tels que 
	// 1. l'arc e =  (orig, dest) est dans le chemin partant de x 
	// 2. _ressATransmettre[x][e] > 0

	vector<int> construireListeSommetA1(int somCour, vector<int> & ordre, int orig, int dest);


	//suite de assignerPhase1 si succes : ajuste les valeurs de flot pour les arcs != arcDebit
	//de sorte que le debit soit le meme partout sur le chemin de somCour
	//peut retourner un echec si on est oblige de trop decaler la date de debut de somCour
	int assignerPhase2(vector<int> & ordre, int somCour, int arcDebit);

	//si pour un sommmet s donné il n existe qu'un seul arc e tel que flot(s, somCour, e) > 0
	// alors on essaie de vider cet arc en augmentant flot (y, somCour, e) pour un y qui partage 
	// des arcs avec somCour et qui transmet deja de la ressource a somCour
	void assignerRegroupement(int somCour);

	//===========================================================================
	// sous fonctions pour algo general : init de l'ordre + resolution des echecs
	public:
	vector<int> genererOrdreInit();

	vector<int> genererOrdreInit2();

	vector<int> genererOrdreInit(vector<int> & rang);

	//l'arc (prec,som) a ete ajoute a coupleImpose ce qui signifie que dans ordre som ne doit pas preceder prec
	// => on fait les modifs necessaires
	void modifierOrdre(vector<int> & ordre, const vector < pair<int, int> >  & coupleImpose);


	//ranvoie vrai si un prec p (de coupleImpose) de som est tel que isIn[p] = vrai
	bool appartientPrec(const vector<pair<int, int> > & coupleImpose, vector<bool> & isIn, int som);

	//retourne vrai si l'ajout du couple (orig, dest) dans coupleImpose genereun cricuit
	bool circuit(vector < pair<int, int> > & coupleImpose, int orig, int dest);

	//on cherche x avant som dans ordre, qui partage au moins un arc avec som et dont la date de debut est la plus grande possible
	//return -1 si x non trouve
	int rechercheAction(const vector<int> & ordre, int som);

	//=====================================================================================================
	// SOLUTION

	//===========================================================================
	// sous fonctions pour la verification

	//ajoute un nouveau debit sur l'intervalle [debut,fin]
	//le vecteir arcDebit contient deja une liste d'intervalle [deb,fin] avec une valeur de debit pour chaque intervalle
	//il faut donc ajouter le nouvel interval[debut, fin] avec son debit en gardant les intervalles de arcDebit sans intersection et trie sur les deb croissants: 
	// si l'intervalle est contenu dans un autre il faut faire la somme des debit, s'il est à cheval sur un autre il faut decouper en intervalle plus petits etc...
	//retourne vrai si on peut ajouter le nouveau debit sans jamais depasser la capacité de l'arc (orig, dest), faux sinon

	bool ajouterArcDebit(vector<debitTemporise> & arcDebit, const int CAP, double debit, double debut, double fin);

	//renvoie vrai si [deb, fin] intersecte l'intervalle DT
	bool intersection(debitTemporise & DT, double deb, double fin);

	bool verifieFlot(bool verifKirch = true);

	bool verifieFlotSum();

	void afficheUtilisation(int indArc);

	void afficheArcCritique();

	//===========================================================================================
	// info sur la solution

	//renvoie la somme ponderee par la population des differences entre les deadlines et les derniers individus
	double calculCoutSolSomme();

	//renvoie la somme ponderee par la population des dates de sorties des individus
	double calculCoutSolSommeMilieu();

	//renvoie la plus petite marge et remplit les vecteurs : 
	//xMin avec tous les jobs qui donnent cette marge la plus faible
	//xMin2 avec tous les jobs qui donnent la deuxième plus petite marge
	double calculCoutSolMin(vector<int> & xMin, vector<int> & xMin2);

	//renvoie la date a laquelle le dernier individus partant du noeud evac i arrive au safe node
	double dernierSafe(int i);

	//renvoie la date a laquelle le dernier individus de la popuplation i 
	//arrive sur un noeud sachant la date d arrivee du premier
	double dernier(int i, double arriveePremier);


};



#endif // RCPSP_ALGO_H
