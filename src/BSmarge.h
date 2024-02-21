
#ifndef BS_MARGE_H
#define BS_MARGE_H


#include "util.h"
#include "Instance.h"
#include <sstream>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN



//==============================================================================================
//cette classe implemente les fonctions pour calculer une estimation des marges qu'il est possible d'obtenir
// on a deux algos principaux : 
// - calculMargeArc : estimation optimiste = borne sup sur les marges 
// - calculMargeArcAmelioree : estimation parfois optimiste parfois pessimiste : ni borne sup, ni borne inf, mais donne une "bonne" estimation des marges 

// les solutions calculées par les algos pricipaux sont preemptives (debits non constants)

//la methode principale calculMargeArc(int orig, int dest) calcule les marges optimistes pour les jobs 
//qui empruntent cet arc, on init les attributs de classe avec toutes les donnees relatives a cet arc
//pour eviter d'aller les chercher dans l'instance a chaque fois qu'on en a besoin

//on travaille avec un sous ensemble des jobs "noeuds d'evacuation" : on considere uniquement les jobs qui empruntent
// larc orig->dest considere ==> toutes les variables de trvail de type vector sont numerotes sur le sous-ensemble 
// _jobs[j] donne le "vrai job" dans _instance qui correspond a j

//ATTENTION : dans cet modelisation on considere le franchissement d'une "porte imaginaire" qui correspond au fait 
//d'entrer sur l'arc cible. Le franchissement de cette porte est supposé de durée 1 (car le debit est en unité de temps), 
//on peut voir cette porte comme une ressource disponible en quantité la capacité de l'arc cible et cette ressource est 
//renouvelee toutes les unites de temps
// => pour passer d'une instance du probleme "Alicia / Ch. Artigues" a une instance de ce probleme il faut bien prendre en compte
//qu'on ajoute le franchissement d'une porte qui allonge le chemin du noeud evacuation vers safe node de 1 donc il faut reduire de 1
//la longueur de l'arc cible quand on calcule les LF dans ce modele et aussi quand on calcule les "fins optimistes"

class BSmarge
{
	//structure qu'on utilise dans l'heuristique construireVraieSolution : sauvegarde les modifs sur les ES / LF / debit
	// sans modifier l'instance
	struct dataIns 
	{
		vector<double> _ES;
		vector<double> _LF;
		vector<double> _debitMax;
	};

private:
	//====================================================================
	// donnees d entree communes a calculMargeArc et calculMargeArcAmelioree

	Instance * _ins;
	

	int _nbJob; //nombre de jobs = taille du vector _jobs
	vector<int> _jobs; //_jobs = ensemble des jobs sur lesuqles on travaille (indice a partir de 0) :
	//jobs(i) donne le numero du noeud d'evacuation dans _instance associé au job i

	double _cap; //capacite de l'arc considere
	vector<double> _ES; //_ES(i) = date au plus tot d'arrivee du job i sur l'arc cible (= arrivee a la porte)
	vector<double> _LF; //LF(i) = date au plus tard a laquelle on doit franchir la porte (= sortie de la porte, LF >= ES + 1 car il faut 1 pour franchir la porte)
	vector<double> _debitMax; //_debitMax(i) = debit max pour le job i

	//====================================================================
	// donnees d entree propres a calculMargeArcAmelioree et calculMargePLdicho

	vector< pair<int,int> > _arcs; //(*) si i < _nbJobs, _arcs[i] = arc de depart pour le job i (population correspondante dans _jobs [i]) 
								   //(*) si i >= _nbJobs, _arcs[i] = arc de transit numerote de sorte que les arcs peres (vers le safe node) > arcs fils
	vector<int> _capaArc;

	vector< vector<int> > _toArcs; //_toArcs[o][d] donne le numero dans _arcs de l'arc o -> d (-1 si n'existe pas)

	vector<vector<int> > _arc2jobs; //_arc2jobs[i] donne la liste des jobs qui utilise l 'arc i

	vector<vector<int> > _filsArc; //_filsArc[e] donne la liste des arcs fils de e
	
	//=============================================================
	//attributs "variables de travail" : tous les vecteurs sont indices de 0 a _jobs.size()-1
	vector<double> _popRestante; //_popRestante(i) = population restant a evacuer pour le job i 
	vector<bool> _actif; //_actif(i) = vrai si i non fini (il reste de la population a evacuer)
	vector<double> _listeDate; //liste des dates qui ont ete utilisees au cours de l'algo
	double _tCour; //date courante
	vector<double> _marge;//_marge(i) = estimation de la marge (optimiste) pour i à l'itération courante

	//=============================================================
	// variables de sortie

	vector<vector<double>> _debit; //_debit(j,n) donne le debit du job j a l'etape n
	vector<double> _fin; //_fin(j) = date de fin optimiste pour j (pour entrer sur l'arc considere)
	vector<double> _debut; //_debut(j) = date de _debut pour j (pour entrer sur l'arc considere)

public:

	BSmarge(Instance * ins) : _ins (ins) {}

	//appelle calculMargeArc pour l'arc qui va vers le sommet "safe node"
	//on suppose qu'il y en a qu'un seul (le même pour toutes les populations)
	double calculMargeArc();

	//appelle calculMargeArc pour l'arc qui va vers le sommet "safe node"
	//on suppose qu'il y en a qu'un seul (le même pour toutes les populations)
	//remplit la matrice avant de sorte avant (i,j) = 1 si i avant j (ie. deb(i) < deb(j) et fin(i) < fin(j) ) 
	//rang[i] donne le nombre de job avant i
	//renvoie la marge et un booleen qui vaut vrai si on peut transformer la solution BS en "vraie" solution (ie
	//en solution avec debit constant) juste en faisant la moyenne des debits pour chaque job
	pair<double, bool> calculMargeArcAmelioree(vector<double> & finOptimiste, vector<int> & rang);

	//calcul les marges optimiste (i.e. BS sur les marges) pour l'arc orig -> dest
	double calculMargeArc(int orig, int dest);

	//calcul les marges optimiste (i.e. BS sur les marges) en considerant le sous arbre d'evacuation 
	// qui a sa racine en dest (version amelioree de calculMargeArc)
	//ATTENTION cet algo donne souvent des marges optimistes (i.e. meilleure que la solution optimale)
	//mais il peut se tromper et donner des marges pessimistes (moins bonnes que la solution optimale)
	//cet algo n'est donc ni une borne sup, ni une borne inf mais il donne une "assez bonne" idée des marges optimales
	double calculMargeArcAmelioree(int orig, int dest, dataIns * data = 0);

	// resout le probleme preemptif : i.e. debit variable le long d'un chemin
	// en considerant le sous arbre d'evacuation qui a sa racine en dest 
	double calculMargePLdicho();

	//on appelle iterativement l'algo calculMargeArcAmelioree qui calcul une BS,
	// a chaque it on fixe le job le plus stresse a son debit moyen et on recommence,
	// a la fin on obtient une "vraie" solution (tous les jobs ont un debit constant) 
	double construireVraieSolution();

	//on appelle iterativement l'algo calculMargeArcAmelioree qui calcul une BS,
	// a chaque it on fixe le job le plus stresse a son debit moyen et on recommence,
	// a la fin on obtient une "vraie" solution (tous les jobs ont un debit constant) 
	double construireVraieSolution(dataIns & data, int o, int d);

	//on verifie que la solution definie par _listeDate, _debit verifie bien les contraintes suivantes : 
	// 1. la capacite de l'arc (_cap) est respectee a toute date  
	bool verifieSolution() const;



private:

	//on initialise les donnees d'entree pour calculMargeArc (jobs concernes par l'arc etc...)
	//si on souhaite utiliser des donnees ES / LF / debitMax differente de celles de linstance 
	//on passe un pointeur dataIns non nul
	void initData(int orig, int dest, dataIns * data = 0);

	//on initialise (renumerote) les arc d'entree pour calculMargeArcAmeliorer 
	//attention necessite que initData ait ete appele au prealable
	void initArc(int orig, int dest);

	//init les variables de travail
	void init();

	//init data avec 
	void allocDataIns(dataIns & data);

	//clear tous les vector de travail et sortie 
	void cleanData();

	//renvoie la somme des debits max des jobs dans v
	double sommeDebitMax(const vector<int>& v);

	//fixe le debit a 0 pour tous les jobs dans v et pour l'iteration n
	void fixeDebitNul(const vector<int>& v, int n);

	//retourne la date du job en cours qui fini le plus tot
	double finJobEnCours(int n);

	//renvoie la date min a laquelle la marge d'un job non en cours devient assez petit pour que le job passe ds les jobs en cours
	double rattrapageMarge(const vector<double>& alpha, const vector<double>& margeOrdo);

	//idem rattrapageMarge mais adapte au cas particulier de marge amelioree où tous les alphas d'un meme groupe peuvent etre differents
	double rattrapageMargeAmelioree(const vector<int> & jobPossible, const vector<double>& alpha);

	//retourne la marge minimum (necissite que _fin et _LF soient initialises)
	double calculMargeMin();


	//retourne les marges triees de la plus petite a la plus grande et avec les doublons (a epsilon pres) supprimes
	vector<double> ordonneMarge();
	
	//retourne un vecteur v tel que v(i)  donne la liste des jobs avec une marge = margeOrdonnee(i)
	vector<vector<int>> classeJobParMarge(vector<double>& margeOrdonnees);


	//renvoie un vecteur contenant les composantes de alpha triees en ordre croissant (en ajoutant 0 et en supprimant les doublons)
	// on ne garde que les valeurs de alpha pour les jobs dans job
	vector<double> triRestrictionAlpha(const vector<double> & alpha, const vector<int> & job);

	//tri les valeurs de alpha du plus petit au plus grand en ajoutant 0 et en supprimant les doublons
	vector<vector<int>> ordonneJobSuivantAlpha(const vector<double>& alpha, const vector<double>& alphaTri, const vector<int>& job);

	int calculIndiceMaxAlpha(const vector<double>& alphaTri, const vector<vector<int>>& jobOrdonneAlpha, double cap);

	double calculBetaMargeAmelioree(int ind, const vector<double>& alphaTri, const vector<vector<int>>& jobOrdonneAlpha, double cap);

	// renvoie vrai si la solution de la BS peut etre transformee en "vraie" solution : 
	// on fait la moyenne des debits pour chaque job en ponderent par la duree des intervalles sur le quel elle s execute
	//attention : si l'instance est divisee en sous arbre (pas d'arc commum a tous les jobs) alors il faudra appeler cette procedure
	//sur chaque sous arbre
	bool isRealisable();

	//===========================================================================================
	// sous fonction pour calculMargePLdicho
	
	//
	auto calculMargePLdichoArc(int orig, int dest) -> double;

	//on lance un PL pour savoir si toutes les pop peuvent evacuer avec une marge >= margeCible  
	bool isMargeRealisable(double margeCible, vector<double> & date);

	//on ordonne l'ensemble des dates ESi, LFi-margeCible par ordre croissant
	vector<double> trierDate(double margeCible);

	//on donne la marge fournit par calculMargeAmelioree (donne une solution non optimale au probleme preemptif)
	double calculMargeCibleMin();
	
	//calcul la meilleure marge possible compte tenu des population
	double calculMargeCibleMax();

};




#endif // BS_MARGE


