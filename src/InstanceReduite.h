#ifndef INSTANCEREDUITE_H
#define INSTANCEREDUITE_H

#include "Instance.h"

//=========================================================
// on ne stocke dans cette classe que les donnees necessaires 
// a la resolution du problème d'evacuation non preemptif
// version 1 (avant janv. 2023) : les informations sont extraites de la classe Instance qui lit les instances generees par le programme de Christian / Emmanuel
// version 2 (janv. 2023) : on genere directement l'instance InstanceReduite de maniere aleatoire

//les ES et LF sont donnes pour entrer sur l'arc (orig, dest) passe en param du constructeur

class InstanceReduite 
{
public:

	int _nbJob; //nb jobs = nb de population a evacuer
	int _nbArc; //nb d'arcs du chemin

	vector<int> _job; //_job[i] = numero du job dans Istance

	vector<int> _pop;//_pop[j] = population (nombre de personnes a evacuer) pour j [0..._nbJob-1]
	vector<double> _ES;//_ES[j] = date au plus tot d'arrivee a la porte de j [0..._nbJob-1]
	vector<double> _LF;//_LF[j] = date au plus tard d'arrivee a la porte [0..._nbJob-1]
	double _ESmin; //ESmin = min des ES[j] = date au plus ou une population peut arriver a la porte
	vector<double> _debitMax; //debitMax[j] = debit max de j
	vector<double> _debitMin; //debitMin[j] = debit min de j (pour arriver avant sa deadline)

	vector<int> _cap; //_cap[i] = capacite de l arc i [0..._nbArc-1]
	vector<vector<int>> _arcToJob; //_arcToJob[i] = liste des jobs qui passent sur l'arc i
	vector<int> _jobToArcInit; //_jobToArcInit[j] = premier arc d'evacuation de j
	vector<vector<int>> _jobToArc; //_jobToArc[i] = liste des arcs empruntes par i

	vector<vector<bool> > _utilise; //utilise[i][a] = true si le job i utilise l'arc a

	//initialise les attributs avec ins
	void init(const Instance & ins, int orig, int dest);


	//=========== fonctions ajoutees en janv. 2023 pour les tests avec instances "generales"


	//================ fonctions pour generer une instance 

	//alloue les tableaux
	void alloc(int nbJob, int nbArc);


	//generation aleatoire d'une InstanceReduite : un job utilise un arc avec une proba alpha
	void genereRandom(int nbJob, int nbArc, double alpha, int ESmax);

	//fixe les LF a l'aide d'une heuristique gloutonnne pour etre sure que l'instance est faisable
	void fixeLF();

	//consoOk renvoie vrai si on peut mettre le job a la date t telle que conso[i].first est la plus grande date <= t
	bool consoOk(int job, int t, int it, const vector<pair<int, vector<int>>> & conso, int longueur, int hauteur);
	
	//insere le job i (longueur * hauteur) en t dans conso
	//it est l'indice du plus grand t dans conso <= t
	void insere(int job, int t, int it, int longueur, int hauteur, vector<pair<int, vector<int>>> & conso);


	//================ fonctions pour sauvegarder / lire une instance dans un fichier


	void lire(const string & ficName);

	void ecrire(const string & ficName);

	//================ fonctions de verification
	
	bool isSolution(const vector<double> & dateDebut, const vector<double> & debit, double marge = -1);



};






#endif // INSTANCEREDUITE_H
