#ifndef INSTANCE_VS_
#define INSTANCE_VS_

#include <iostream>
#include <vector>
#include <utility>
#include <string>

#include <fstream>
#include <sstream>

#include <algorithm>



#include "util.h"

using namespace std;


//===================================================
// on suppose que les noeuds evac, transit, safe sont tous disjoints


class Instance
{
public:

    int _nbNode = 0; //in the vectors, index nodes go from 1 to _nbNode (they include evacuation nodes, safe nodes=shelter, transit nodes) / index 0 is for the source in the expanded network

    int _firstEvacuationNode = 0; //index of first evac node
    int _firstSafeNode = 0; //index of first safe node
    int _firstTransitNode = 0; //index of first transit node

    int _lastEvacuationNode = 0; //index of last evac node
    int _lastSafeNode = 0; //index of last safe node
    int _lastTransitNode = 0; //index of last transit node

	bool _hasDuplicateNode = false; //passe a vrai si on a des noeuds qui sont transit + evac (auquel cas ils sont dupliques)

    //REMARQUE: ON MET DES MATRICES POUR UN ACCES DIRECT MAIS LE GRAPHE N EST PAS COMPLTET
    // LES ARCS NON EXISTANT SONT REPERES PAR UNE CAPACITE = 0

    vector<int> _pop; // index from _firstEvacuationNode to _lastEvacuationNode: population

    vector<int> _maxRateEvac; // index from _firstEvacuationNode to _lastEvacuationNode: max evac rate at node

    vector< vector<int> > _time; //_time[i][j] = time to go from i to j

    vector< vector<int> > _capArc; //_capArc[i][j] = capacity of arc i,j (nb personn per minute)

    vector< int > _capNode; //_capNode[i] = capacity of node i (nb person that can stay in i during a while)

    vector< vector<int> > _dead; //_dead[i][j] = last time to go on arc ij

    int _TMAX; //total duration of the process

    vector <int> _minDate ; //_minDate[i] = date au plus tot a laquelle on peut arriver en i = plus court chemin depuis un noeud d evac vers i

    vector<vector<int> > _succ; //liste des successeurs de i

    vector<vector< SomWithWindow > > _chemins;//_chemins[i] = chemin d evacuation du node evac i (si instance vient de Christian)

	//============================================================================================
	//quelques donnees pour les stats
	double _stat_capa;//moyenne des capa entrantes / capa sortante pour chaque noeud interne
	double _stat_temps_total_min; //temps evacuation si tout le monde evacue au max

	void calcul_stat();
	//============================================================================================



    //generate data randomly
    //genere un graphe par couche : la 1ere = node evac, la derniere = node safe et entre les deux on met des couches de transit node
    void generateDataFile(string name, int nbEvacNode, int nbSafeNode, int nbCoucheInterne, int nbNodePerCouche);


	//generate data randomly
	//genere un graphe en arbre
	void generateTree(int nbEvacNode, long seed);

    //read data from a file (in our own format)
    void readFileLN(const string & name);

    //read data from a file made by the software of CArtigues
    void readFileCArtigues(const string & name);

	
	void readFileAlicia(const string & name);

	//on verifie que l'arbre reduit (donne sous la forme de ses chemins) a les bonnes caracteristiuqes (capa croissantes...)
	bool verifArbreReduit();

    // display instance data
    void display();



    //trace graphe
    // si sym = true alors on trace des aretes sinon des arcs
    void trace(string name, bool sym);

	//trace graphe
	// si sym = true alors on trace des aretes sinon des arcs
	bool traceChemin(string name);

    //utilise les plus courts chemins pour calculer les min dates
    void computeMinDate();

	//dansle cas ou on a les chemins d evacuation predefini, on peut calculer facilement un Tmax et des fenetres de temps
	//dans cette fonction on ajuste le Tmax avec les deadline
	void computeWindowsEtTmax();

	//dansle cas ou on a les chemins d evacuation predefini, on peut calculer facilement un Tmax et des fenetres de temps
	//dans cette fonction on ne change pas le TMax, c'est une donnee du probleme
	void computeWindows();

	//dans le cas ou des chemins sont definis, taux d'evac max <= CAP des arcs du chemin
	void reduitMaxRateEvac();

	//========================================================================
	// fonctions utilisees dans RCPSP_graphe


	//retourne la liste des arcs a la fois dans le chemin partant de i et de j
	vector< ArcCommun > arcEnCommun (int i, int j);

	//retourne la liste des arcs dans le chemin i sous forme de vector< ArcCommun >
	vector< ArcCommun > arcCh_i(int i);

	//retourne la longueur de j vers l arc(orig, dest) s'il existe dans le chemin de j, -1 sinon
	int appartientChemin(int j, int orig, int dest);

	//retourne la longueur du chemin de i vers le safe node
	int getLongueurCh(int i);

	//retourne la date max pour arriver au safe node en partant de i (LF du dernier noeud du chemin)
	int getFinMax(int i);

	double getDebitMin(int i);

	double getDebitMax(int i) const;

	//on recalcule les LF pour chaque job du chemin :
	// la fin max d'un chemin devient ceil(finChemin[i]) * alpha
	void recalculeLF(vector<double> & finChemin, double alpha);

	//on verifie que les fenetres de temps sont realisables
	bool verifFen();

	//=====================================================================================
	// fonction pour dupliquer les noeuds dans le cas où un meme noeud est a la fois
	// transit et evac

	//retourne l'ensemble des noeuds d'evacuation qui sont aussi noeuds de transit
	// on memorise app | app[i] = vrai si i appartient a l'ensemble
	// et prec | prec[i] = liste des precdent de i
	vector<int> evacTransitNode(vector<bool> & app, vector < vector<int> > & prec);

	// on duplique les noeuds qui sont a la fois evacuation et transit
	// et on met a jour toutes les infos pour que la solution de l'instance dupliquée soit la meme que celle de l'instance d'origine
	// ATTENTION : les longueurs des chemins et des deadlines sont augmentes de 1 avec cette manip, donc les marges restent les mêmes mais 
	// les dates d'arrivee sont augmentee d'une unite
	void dupliqueNodeEvacTransit();

	//a partir des deadlines des chemins on recalcule les deadlines sur les arcs 
	void recalculeDead();
};


#endif

















