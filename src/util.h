#ifndef UTIL_
#define UTIL_

#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;


//#define _VERIF_ 


//#define _VERIF_BRANCHEMENT_
#define _MAX_MIN_ // si actif alors le RCPSP coulissant maximise le min des marges (utilise dans RCPSP_Algo::resolution);

//#define _PARCOURS_PROF_ //si _PARCOURS_PROF_ est defini on fait un parcours en profondeur (avec une pile) sinon un parcours en largeur (file à priorité)


const int TIME_INFINITY = 100000;
const int INT_INFINITY = 1000000000;
const double EPSILON = 0.000001;
const double EPSILON_DEBIT = 0.001; //deux debits consecutifs sont egaux s'ils sont egaux a epsilondebit pres => attention on autorise de violer la ctr de capacite a EPSILON_DEBIT pres
const double EPSILON_DEBIT_VERIF = EPSILON_DEBIT- EPSILON_DEBIT/2; //quand on verifie il faut un epsilon legerement plus petit a cause des problemes numeriques

const double EPSILON_DATE = 0.01; //deux dates consecutives sont egales si elles le sont a EPSILON_DATE pres
const double EPSILON_DATE_VERIF = EPSILON_DATE - EPSILON_DATE / 2;
const double EPSILON_CUT = 0.0001;
const double EPSILON_CUT_CPLEX = 0.00001;

const double EPSILON_COMP = 0.00001; //epsilon pour les comparaison (> ou <)


#ifdef _WIN64
const int NB_THREAD_CPLEX = 1; //pour PL dicho et MIP flot (
const int TIME_LIMIT_CPLEX = 3600;
const int RAM_LIMIT_CPLEX = 10000; //en Mo
const int CLOCK_TYPE_CPLEX = 2; //1 => Temps CPU (non supporte par WINDOWS !!) -- 2 => Temps horloge (temps physique total écoulé) ; par défaut
#else
const int NB_THREAD_CPLEX = 1; 
const int TIME_LIMIT_CPLEX = 3600;
const int RAM_LIMIT_CPLEX = 80000; 
const int CLOCK_TYPE_CPLEX = 1; 
#endif

inline void stopProg(string msg)
{
#ifdef WIN64
    int stop;
    cout << msg << endl;
    cin >> stop;
#else
	cout << "ERREUR : " << msg << endl;
	exit(-1);
#endif
}

inline size_t toInd(int i)
{
    return static_cast<size_t>(i);

}

inline int toInt(size_t i)
{
    return static_cast<int>(i);

}

//expanded node
struct ExpNode
{

    int _nodeIndex = 0; //index of real node in Instance
    int _time = -1; // in 0..TMAX EXCEPT source (0,-1) and sink (N+1, TMAX+1)

    ExpNode(int n, int t) : _nodeIndex(n), _time(t){}
    ExpNode(){}

	
};

inline ostream & operator <<(ostream &os, const ExpNode &node)
{
	os << "(" << node._nodeIndex << ", " << node._time << ")" ;

	return os;
}

struct ExpArcValue
{
	int _orig = 0;
	int _timeOrig = 0;
	int _dest = 0;
	int _timeDest = 0;
	double _value = 0;//on a besoin d un double car qd on reconstruit a partir du graphe reduit on divise par le coeff de reduc

	ExpArcValue(int o, int to, int d, int td, double v) : _orig(o), _timeOrig(to), _dest(d), _timeDest(td), _value(v) {}


	ExpArcValue() {}

	inline void afficheSiRoute()
	{
		if (_orig != _dest)
			cout << _orig << "," << _timeOrig << "->" << _dest << "," << _timeDest << "(" << _value << ")" << endl;
	}

};




struct ExpSomValue
{
    int _som;
    int _time;
    double _value;

    ExpSomValue(int s, int t, double v) : _som(s), _time(t), _value(v){}

};

struct SomWithWindow
{
    int _som;
    int _ES;
    int _LF;

    SomWithWindow(int s, int es, int lf) : _som(s), _ES(es), _LF(lf){}

};

struct Arc
{
	int _orig;
	int _dest;
	

	Arc(int o, int d) : _orig(o), _dest(d) {}
	Arc() {}
};


//sert a sauvegarder le debit sur un arc entres les dates debut et fin
struct debitTemporise
{
	double _deb;
	double _fin;
	double _debit;

	debitTemporise(double d, double f, double de) : _deb(d), _fin(f), _debit(de) {}
	debitTemporise() {}

	bool operator < (const debitTemporise & D)
	{
		return this->_deb < D._deb;
	}
};



//struct pour stocker un arc (_orig, _dest) qui se trouve a la fois dans deux chemins : 
// long1 = longueur pour arriver a l'origne de cet arc par le premier chemin
// long2 = -------------------------------------------------second ------

struct ArcCommun
{
	int _orig;
	int _dest;
	int _long1;
	int _long2;

	ArcCommun(int o, int d, int l1, int l2) : _orig(o), _dest(d), _long1(l1), _long2(l2) {}
	ArcCommun() {}
};


//renvoie vrai si e appartient a v (v par forcement trie)
inline bool appartient(int e, const vector<int> & v)
{
	for (int i : v)
		if (e == i)
			return true;

	return false;
}

//calcul l'intersction de v et w (vecteurs non necessairement tries)
inline vector<int> intersection(const vector<int> & v, const vector<int> & w)
{

	vector<int> res;

	for (int i : v)
	{
		for (int j : w)
		{
			if (i == j)
			{
				res.push_back(i);
				break;
			}
		}
	}

	return res;
}

inline string dateEtHeure()
{
	// date / heure actuelle basée sur le système actuel
	time_t tmm = time(0);

	// convertir en forme de chaîne
	char* dt = ctime(&tmm);
	

	return dt;

}

void dessine(int N, const vector<double> & debut, const vector<double> & debit, const vector<int> & pop);

#endif
