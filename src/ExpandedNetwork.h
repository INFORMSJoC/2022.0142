#ifndef EXPANDEDNETWORK_H
#define EXPANDEDNETWORK_H

#include "Instance.h"



class ExpandedNetwork
{
public:
    //=================================================================
    // some new definitions



    //expanded arc
    struct ExpArc
    {


        ExpNode _orig; // ExpNode origin
        ExpNode _dest; // ExpNode destination
        double _cost = 0.0; // cost of the arc in the objective function

        ExpArc(ExpNode o, ExpNode d) : _orig(o), _dest(d){}
    };


    //=================================================================
    // attributs de classe
    Instance * _ins;
    ExpNode _source = ExpNode(0, 0) ;
    ExpNode _sink ;
    vector<ExpArc> _arc;//all arcs in the order: _inputExpArc, _outputExpArc,_stationExpArc, _routeExpArc,

    vector<vector<vector<int>>> _arcIn;//_arcIn[i][t] = liste des arcs entrant dans (i,t)
    vector<vector<vector<int>>> _arcOut;//_arcIn[i][t] = liste des arcs sortant dans (i,t)

    int _firstInputExpArc; //index of first input arc
	int _firstOutputExpArc; //index of first output arc
	int _firstStationExpArc; //index of first station arc
	int _firstRouteExpArc; //index of first route arc

	int _lastInputExpArc; //index of last input arc
	int _lastOutputExpArc; //index of last output arc
	int _lastStationExpArc; //index of last station arc
	int _lastRouteExpArc; //index of last route arc

    //=================================================================
    // methodes
    ExpandedNetwork(Instance * ins) : _ins(ins)
    {
        if (_ins->_nbNode == 0)
            stopProg("ERREUR: ExpandedNetwork construteur ne peut etre init car instance non init");
        _sink = ExpNode(_ins->_nbNode+1, _ins->_TMAX+1);
    }


    //create the list of ExpNode
    //void createExpNode();

    virtual void createExpArc();

    // write a .lp file that cplex can read
    // si obj = true on met une focntion obj, sinon c juste un probleme de faisabilite
    void toFileLP(const string& name, bool obj = true);

    // return arc name as (x,t)(y,t')
    string getArcName(size_t i);

	// retourne les info sur l arc, range dans ExpArcValue (champs valeur inutil)
	ExpArcValue getArcInfo(size_t i);

    // read a sol file written by cplex
    void readSolFile(string name, vector<ExpArcValue> & solution);

    // trace la solution avec graphviz
    void traceSol(vector<ExpArcValue> & solution);

    // calcule et affiche la deadline la plus serre (i.e. la marge de securite pour l'arc utilise au plus proche de sa deadline)
    double calculMaxMin(vector<ExpArcValue> & solution);

	void afficheNetwork();

    //=======================================================================
    // fonctions pour extraire les chemins

    // extrait des chemins a aprtir d une solution
    // chaque chemins[i] est un chemin de la forme (som, time)... qui commence en (nodeEvac,0) ou nodeEvac est un des noeuds d evac
    void extraireChemin(vector<ExpArcValue> & solution, vector<vector<pair<int,int> > > chemins);

    //a partir de la solution ExpArcValue on cree solutionRange \ solutionRange[i][t] donne la liste des sommets successeurs dans la solution
    void rangeSolution(vector<ExpArcValue> & solution, vector<vector<vector<ExpSomValue> > > & solutionRange);

    //retourne l indice du meilleur suivant de (som,t)
    size_t  choixSuivant(size_t som, int t, vector<vector<vector<ExpSomValue> > > & solutionRange);

    //enleve delta unite de flot au chemin defini par cheminCour et stocke dans solutionRange
    void enleveDeltaChemin( vector<vector<vector<ExpSomValue> > > & solutionRange, vector<pair<int,int> > & cheminCour, vector<size_t> & indChemin, double delta );

    //=======================================================================
    // fonctions pour construire une vraie solution

    //voir GrapheReduit


};

#endif // EXPANDEDNETWORK_H


















