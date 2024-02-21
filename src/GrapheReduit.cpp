#

#include "GrapheReduit.h"



//construit le graphe a partir d une solution
void GrapheReduit::constructionGraphe(vector<ExpArcValue> & solution)
{
    vector < vector<int> > ind;//ind[x][t] = indice du sommet x,t (-1 si n existe pas)

    ind.assign( _ins->_nbNode + 1, vector<int> (toInd(_ins->_TMAX + 1), -1)  );

    //premier sommet = source
    _sommets.push_back(ExpNode(0,0));
    ind[0][0] = -1;

    //========================================================================
    //on parcourt la solution une premiere fois pour connaitre les sommets
    int cpt = 1;

    for (size_t i = 0; i < solution.size() ; ++i)
    {
        int som_o = solution[i]._orig;
        int time_o = solution[i]._timeOrig;

        int som_d = solution[i]._dest;
        int time_d = solution[i]._timeDest;

        if (ind[toInd(som_o)][toInd(time_o)] == -1)
        {
            ind[toInd(som_o)][toInd(time_o)] = cpt;
            cpt++;
            _sommets.push_back(ExpNode(som_o,time_o));
        }


        if (ind[toInd(som_d)][toInd(time_d)] == -1)
        {
            ind[toInd(som_d)][toInd(time_d)] = cpt;
            cpt++;
            _sommets.push_back(ExpNode(som_d,time_d));
        }

    }


    //========================================================================
    //on parcourt la solution une deuxieme  fois pour connaitre les arcs

    _succ.resize(_sommets.size());

    for (size_t i = 0; i < solution.size() ; ++i)
    {
        size_t som_o = toInd(solution[i]._orig);
        size_t time_o = toInd(solution[i]._timeOrig);

        size_t som_d = toInd(solution[i]._dest);
        size_t time_d = toInd(solution[i]._timeDest);

        double val = solution[i]._value;

        _succ[ toInd(ind[som_o][time_o]) ].push_back( pair<int, double> ( ind[som_d][time_d], val) );



    }


}



// pour chaque sommet calcule un plus long chemin
void GrapheReduit::calculPlusLongCh()
{
    //init
    _plusLongCh.assign(_sommets.size(), TIME_INFINITY);
    _plusLongCh[0] = 0;

    size_t somExpCour = 0;

    vector<int> file;
    file.push_back(0);

    while( !file.empty() )
    {

        // sommet dans graphe init
        int somCour = _sommets[somExpCour]._nodeIndex;

        //on parcourt tous les succ
        for (size_t s = 0; s < _succ[somExpCour].size(); ++s)
        {
            int somExpSucc = _succ[somExpCour][s].first;
            int somSucc = _sommets[somExpSucc]._nodeIndex;

            int val = _plusLongCh[somExpCour] + _ins->_time[somCour][somSucc];

            if( val < _plusLongCh[somExpSucc]  )
            {
                _plusLongCh[somExpSucc] = val;
                file.push_back(somExpSucc);

            }

        }

    }//fin while


}


//a partir du graphe et des plus long chemins on construit la solution non reduite
void GrapheReduit::getSolution(vector<ExpArcValue> & solutionFinale)
{
    for (size_t i = 0; i < _sommets.size(); ++i)
    {
        for (size_t j = 0; j < _succ[i].size(); ++j)
        {
            int som_o = _sommets[i]._nodeIndex;
            int time_o = _plusLongCh[i];

            int suiv = _succ[i][j].first;
            double val = _succ[i][j].second;

            int som_d = _sommets[suiv]._nodeIndex;
            int time_d = time_o + _ins->_time[som_o][som_d];

            for (int c = 0; c < _coeffReduc; ++c)
            {

                solutionFinale.push_back( ExpArcValue( som_o, time_o+c, som_d, time_d+c, val / _coeffReduc ) );
            }
        }

    }


}





