#ifndef GRAPHEREDUIT_H
#define GRAPHEREDUIT_H

#include "util.h"
#include "Instance.h"




class GrapheReduit
{
    Instance * _ins;
    int _coeffReduc;

public:

    //init avec une vraie instance
    GrapheReduit(Instance * ins, int coeffReduc) : _ins(ins), _coeffReduc(coeffReduc){}

    vector<ExpNode> _sommets;

    vector< vector<pair<int,double>> > _succ; // succ[i].first = indice dans _sommet du successeur du sommet i, succ[i].second = val flot i - succ i

    vector<int> _plusLongCh; // _plusLongCh[i] = plus long chemin de source vers i


    //==========================================================================
    //construit le graphe a partir d une solution (qui vient d une instance reduite)
    void constructionGraphe(vector<ExpArcValue> & solutionReduite);

    // pour chaque sommet calcule un plus long chemin
    void calculPlusLongCh();

    //a partir du graphe et des plus long chemins on construit la solution non reduite
    void getSolution(vector<ExpArcValue> & solutionFinale);

};







#endif // GRAPHEREDUIT_H
