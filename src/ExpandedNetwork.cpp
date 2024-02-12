#include "ExpandedNetwork.h"
#include <fstream>
#include <algorithm>
#include <sstream>



void ExpandedNetwork::createExpArc()
{

    size_t Nnodes = _ins->_nbNode;
    int T = _ins->_TMAX;

    int cpt = 0;

    _arcIn.resize(Nnodes+1, vector<vector<int>> (T+1 ) );
    _arcOut.resize(Nnodes+1, vector<vector<int>> (T+1 ) );

    // source -> evac nodes

    _firstInputExpArc = cpt;
    for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
    {
        _arc.emplace_back( _source, ExpNode(i,0));
        _arcIn[i][0].push_back(cpt);
        cpt++;
    }

    // safe nodes -> sink
    _lastInputExpArc = cpt-1;
    _firstOutputExpArc = cpt;
    for (int i =  _ins->_firstSafeNode; i <=  _ins->_lastSafeNode; ++i)
    {
        _arc.emplace_back( ExpNode(i,T), _sink);
        _arcOut[i][static_cast<size_t>(T)].push_back(cpt);
        cpt++;
    }


    // (i,t)-> (j,t+time(i,j)) i,j = any node, j != i and cap[i][j] > 0
    _lastOutputExpArc = cpt-1;
    _firstStationExpArc = cpt;
    for (int i =  1; i <=  _ins->_nbNode; ++i)
    {
        for (int t =  0; t < _ins->_TMAX; ++t)
        {
            _arc.emplace_back( ExpNode(i,t), ExpNode(i,t+1));

            _arcOut[i][static_cast<size_t>(t)].push_back(cpt);
            _arcIn[i][static_cast<size_t>(t+1)].push_back(cpt);

            cpt++;
        }
    }

    // (i , t) -> (i , t+1), i = any node
    _lastStationExpArc = cpt-1;
    _firstRouteExpArc = cpt;
    for (int i =  1; i <=  _ins->_nbNode; ++i)
    {
        for (int t =  _ins->_minDate[i]; t < _ins->_TMAX; ++t)
        {

            for (int j =  1; j <=  _ins->_nbNode; ++j)
            {
                if (_ins->_capArc[i][j] > 0 && i != j && t <= _ins->_dead[i][j] && t+_ins->_time[i][j] <= T)
                {

                    ExpArc arc ( ExpNode(i,t), ExpNode(j,t+_ins->_time[i][j]));

                    if (_ins->_dead[i][j] < TIME_INFINITY)
                        arc._cost = _ins->_dead[i][j]-t;

                    _arc.emplace_back(arc);

                    _arcOut[i][static_cast<size_t>(t)].push_back(cpt);
                    _arcIn[j][static_cast<size_t>(t+_ins->_time[i][j])].push_back(cpt);

                    cpt++;
                }

            }
        }
    }
    _lastRouteExpArc = cpt-1;

    cout << " ====>>>> arc exp = " << _arc.size() << endl;

    // cout << "cap = " << _ins->_capArc[2][242] << "  dead = " << _ins->_dead[2][242] << "  time = " << _ins->_time[2][242] << endl;
}


//write a .lp file that cplex can read
void ExpandedNetwork::toFileLP(const string& name, bool obj)
{
    ofstream fic (name);
    //size_t Narc = _arc.size();

    size_t Nnodes = _ins->_nbNode;
    size_t T = static_cast<size_t>(_ins->_TMAX);

    //================================
    //1. objective function

    fic << "Maximize" << endl;

    if (obj)
    {
        for (size_t i = _firstRouteExpArc; i <= _lastRouteExpArc; ++i)
        {

            fic << _arc[i]._cost <<  " F_" << getArcName(i) << " + ";
        }
        fic << "0" << endl; //+ to finish the line properly
    }


    //================================
    //2. constraint

    cout << "ecriture lp ... " << endl;


    fic << "Subject To" << endl;


    //================================
    //2.1. population
    for (size_t i = _firstInputExpArc; i <= _lastInputExpArc; ++i)
    {
        size_t nodeId = _arc[i]._dest._nodeIndex;
        int population = _ins->_pop[nodeId];
        fic << "pop_" << i << ": F_" << getArcName(i) << " = " << population << endl;
    }

    cout << "ecriture capa ... " << endl;

    //================================
    //2.2. capacity safe node

    for (size_t i = _firstOutputExpArc; i <= _lastOutputExpArc; ++i)
    {
        size_t nodeId = _arc[i]._orig._nodeIndex;
        int cap = _ins->_capNode[nodeId];
        fic << "capSafe_" << i << ": F_" << getArcName(i) << " <= " << cap  << endl;
    }


    //================================
    //2.3. capacity transit node
    for (size_t i = _firstStationExpArc; i <= _lastStationExpArc; ++i)
    {
        size_t nodeId = _arc[i]._orig._nodeIndex;
        int cap = _ins->_capNode[nodeId];
        fic << "capNode_" << i << ": F_" << getArcName(i) << " <= " << cap  << endl;
    }

    //================================
    //2.4. capacity route arc
    for (size_t i = _firstRouteExpArc; i <= _lastRouteExpArc; ++i)
    {
        size_t nodeIdOrig = _arc[i]._orig._nodeIndex;
        size_t nodeIdDest = _arc[i]._dest._nodeIndex;
        int cap = _ins->_capArc[nodeIdOrig][nodeIdDest];
        fic << "capRoute_" << i << ": F_" << getArcName(i) << " <= " << cap  << endl;
    }

    //================================
    //2.5. conservation du flot

    cout << "ecriture cons flot ... " << endl;

    for (size_t i = 1; i <= Nnodes; ++i )
    {


        for (size_t t = 0; t <= T; ++t )
        {
            size_t sizeIn = _arcIn[i][t].size();
            size_t sizeOut = _arcOut[i][t].size();

			if (sizeIn > 0 || sizeOut > 0)
			{

				fic << "Kirch_" << i << "_" << t << ": ";

				if (sizeIn > 0)
				{
					for (size_t k = 0; k < sizeIn - 1; ++k)
					{
						fic << "F_" << getArcName(_arcIn[i][t][k]) << " + ";
					}

					fic << "F_" << getArcName(_arcIn[i][t][sizeIn - 1]);
				}

				if (sizeOut > 0)
				{
					for (size_t k = 0; k < sizeOut; ++k)
					{
						fic << " - F_" << getArcName(_arcOut[i][t][k]);
					}

				}

				fic << " = 0 " << endl;

			}
        }
    }

    //=====================================
    // bounds
    /*
    fic << "Generals" << endl;
    for (size_t i = 0; i < Narc; ++i)
    {
        fic << "F_" << getArcName(i) << endl;
    }*/

    fic << "End" << endl;
    fic.close();
}

//return arc name as (x,t)(y,t')
string ExpandedNetwork::getArcName(size_t i)
{
    size_t nodeOrigId = _arc[i]._orig._nodeIndex;
    int nodeOrigTime = _arc[i]._orig._time;

    size_t nodeDestId = _arc[i]._dest._nodeIndex;
    int nodeDestTime = _arc[i]._dest._time;

    stringstream ss;


    ss << "(" << nodeOrigId << "," << nodeOrigTime << ")(" <<  nodeDestId << "," << nodeDestTime << ")";

    return ss.str();
}


ExpArcValue ExpandedNetwork::getArcInfo(size_t i)
{
	

	int nodeOrigId = _arc[i]._orig._nodeIndex;
	int nodeOrigTime = _arc[i]._orig._time;

	int nodeDestId = _arc[i]._dest._nodeIndex;
	int nodeDestTime = _arc[i]._dest._time;

	ExpArcValue E(nodeOrigId, nodeOrigTime, nodeDestId, nodeDestTime, 0);	

	return E;
}

// read a sol file written by cplex
void ExpandedNetwork::readSolFile(string name, vector<ExpArcValue> & solution)
{
    ifstream fic (name);
    string poub;

    ExpArcValue varc;


    while( poub != "<variables>" )
    {
        fic >> poub;
        //cout << poub << endl;
    }

    string varName, varValue, ligne, finLigne, finLigne2;

    while( !fic.eof() )
    {
        getline(fic, ligne, '"');
        size_t sizeLigne = ligne.length();

        if (sizeLigne > 6)
        {
            string savLigne = ligne;
            finLigne = ligne.substr(ligne.length() - 5, ligne.length());
            finLigne2 = savLigne.substr(savLigne.length() - 6, savLigne.length());
        }

        if ( finLigne == "name=" )
        {
            getline(fic, varName, '"');

            istringstream iss(varName);

            //le nom de variable ets de la forme F_(x,t)(y,s), on doit recuperer x t y et s

            getline(iss, poub, '(');
            getline(iss, poub, ',');
            varc._orig = atoi(poub.c_str());

            getline(iss, poub, ')');

            varc._timeOrig = atoi(poub.c_str());
            getline(iss, poub, '(');
            getline(iss, poub, ',');

            varc._dest = atoi(poub.c_str());

            getline(iss, poub, ')');
            varc._timeDest = atoi(poub.c_str());
        }

        if ( finLigne2 == "value=" )
        {
            getline(fic, varValue, '"');
            varc._value = atof(varValue.c_str());
            if ( varc._value > 0.000001 )
                solution.push_back(varc);
        }
    }
}


// trace la solution avec graphviz
void ExpandedNetwork::traceSol(vector<ExpArcValue> & solution)
{
    string name = "solExpNet";
    stringstream ss;
    ss << name << ".txt" ;


    ofstream ficOut (ss.str(), ios::out);// ouverture fichier de sortie

    //premiere ligne du fichier de sortie
    ficOut << "digraph G {" << endl;

    for (size_t i = 0; i < solution.size(); ++i)
    {

        ficOut << solution[i]._orig << "." << solution[i]._timeOrig << "->" << solution[i]._dest << "." << solution[i]._timeDest
               << "[label=" << solution[i]._value << "]" << endl;

    }

    ficOut << "}" << endl;
    ficOut.close();

    stringstream ligCom;

#ifdef _WIN64
    ligCom << "dot.exe -Tpdf -o" << name << ".pdf " << name << ".txt" ;
#endif

    system(ligCom.str().c_str());
}


// calcule et affiche la deadline la plus serre (i.e. la marge de securite pour l'arc utilise au plus proche de sa deadline)
double ExpandedNetwork::calculMaxMin(vector<ExpArcValue> & solution)
{
    double minMarge = TIME_INFINITY;
    double minMargeWithFlow = TIME_INFINITY; //marge multiplie par le nb de personnes qui empruntent l arc
    int fin = 0;

    for (size_t i = 0; i < solution.size(); ++i)
    {
        if (solution[i]._value < 0.000001)
        {
            stopProg("la solution contient des arcs nuls");
        }

        size_t o = static_cast<size_t> (solution[i]._orig);
        size_t d = static_cast<size_t> (solution[i]._dest);

        double marge = TIME_INFINITY;

        if ( d != _ins->_nbNode + 1 && o != d && o != 0) // pas de deadline avec le puits
        {
            marge = _ins->_dead[o][d] - solution[i]._timeOrig;
            fin = max(fin, solution[i]._timeDest);
        }
        if (marge < minMarge)
            minMarge = marge;

        if (marge * solution[i]._value < minMargeWithFlow)
            minMargeWithFlow = marge * solution[i]._value;
    }

    //cout << "  =====================  MARGES  =============================  " << endl;
    //cout << "minMarge = " << minMarge << " ET minMargeWithFlow = " << minMargeWithFlow << endl;


    //cout << "  =====================  TMAX  =============================  " << endl;
    //cout << "fin evac = " << fin  << " et TMax = " << _ins->_TMAX << endl;

	return minMarge;
}


// extrait des chemins a aprtir d une solution
// chaque chemins[i] est un chemin de la forme (som, time)... qui commence en (nodeEvac,0) ou nodeEvac est un des noeuds d evac
void ExpandedNetwork::extraireChemin(vector<ExpArcValue> & solution, vector<vector<pair<int,int> > > chemins)
{
    //1. on range la solution dans un format qui nous convient mieux pour l algo

    vector<vector<vector<ExpSomValue> > > solutionRange;
    rangeSolution( solution, solutionRange);


    //======================================================================
    //2. algo


    int popTotal = 0;

    vector<int> popRestante =  _ins->_pop;

    for (size_t i = _ins->_firstEvacuationNode; i < _ins->_lastEvacuationNode; ++i )
    {
        popTotal += _ins->_pop[i] ;
    }

    int nextEvacNode = _ins->_firstEvacuationNode; //prochain pour point de depart du chemin

    //======================================================================
    //  boucle principale : tant qu il reste du flot, on cherche des chemins

    while( popTotal > 0 )
    {
        while ( popRestante[nextEvacNode] == 0 )
            nextEvacNode++;

        double delta = popRestante[nextEvacNode];

        vector<pair<int,int> > cheminCour; // stocke la liste des sommets du chemin dans la forme (sommet,date) - le premier sommet est (evacNode,0)
        vector<size_t> indChemin; // indChemin[i] = indice du prochain sommet dans la liste solutionRange[cheminCour[i].first][cheminCour[i].second] => ca donne le prochain sommet de chemincour

        bool fin = false;
        pair<int,int> somExpCour (nextEvacNode, 0); // (sommet, date)

        cheminCour.push_back(somExpCour);

        //==================================================================
        // 2.1 on cherche un chemin
        while ( !fin )//on a finit quand on arrive au puits
        {
            //retourne le prochain sommet qui partage le plus de flot avec le courant
            size_t indSomSuiv = choixSuivant(somExpCour.first, somExpCour.second, solutionRange);
            ExpSomValue & expSomTmp =  solutionRange[somExpCour.first][toInd(somExpCour.second)][indSomSuiv];

            int somSuiv = expSomTmp._som;
            int tSuiv = expSomTmp._time;
            double val = expSomTmp._value;


            delta = min(delta, val);
            pair<int,int> somExpSuiv(somSuiv,tSuiv);

            cheminCour.push_back( somExpSuiv );
            indChemin.push_back(indSomSuiv);

            somExpCour = somExpSuiv;

            chemins.push_back(cheminCour);

        }//fin while recherche chemin

          // on enleve delta sur le chemin
        if (delta < 1 - EPSILON)
            stopProg("extraireChemin : chemin nul");

        enleveDeltaChemin( solutionRange, cheminCour, indChemin, delta );

        //==================================================================
        // 2.2 on essaie de reutiliser le chemin existant (memes stations mais dates differentes)

        // a faire

    }

}


void ExpandedNetwork::enleveDeltaChemin( vector<vector<vector<ExpSomValue> > > & solutionRange, vector<pair<int,int> > & cheminCour, vector<size_t> & indChemin, double delta )
{

    size_t somCour = toInd(cheminCour[0].first);
    size_t tCour =  toInd(cheminCour[0].second);

    for (size_t i = 0; i < indChemin.size(); ++i)
    {
        size_t k = indChemin[i];

        int somSuiv = solutionRange[somCour][tCour][k]._som;
        int tSuiv = solutionRange[somCour][tCour][k]._time;

        solutionRange[somCour][tCour][k]._value -= delta;

        somCour= toInd(somSuiv);
        tCour = toInd(tSuiv);

        if (toInd(cheminCour[i].first) != somCour || toInd(cheminCour[i].second) != tCour)
            stopProg(" ExpandedNetwork::enleveDeltaChemin: incoherence chemin");

    }

}

size_t ExpandedNetwork::choixSuivant(size_t som, int t, vector<vector<vector<ExpSomValue> > > & solutionRange)
{
    size_t bestInd = 0;
    double val, bestVal = 0;

    vector<ExpSomValue> & v = solutionRange[som][static_cast<size_t>(t)];

    for (size_t i = 0; i < v.size(); ++i)
    {

        val = v[i]._value;

        if (val > bestVal)
        {
            bestVal = val;
            bestInd = i;
        }

    }

    return bestInd;

}


void ExpandedNetwork::rangeSolution(vector<ExpArcValue> & solution, vector<vector<vector<ExpSomValue> > > & solutionRange)
{
    size_t T = static_cast<size_t> (_ins->_TMAX);

    solutionRange.resize ( _ins->_nbNode+1, vector<vector<ExpSomValue>>(T) );

    for ( size_t i = 0; i < solution.size(); ++i )
    {
        size_t orig = static_cast<size_t> (solution[i]._orig);
        int dest = solution[i]._dest;

        size_t timeOrig = static_cast<size_t> (solution[i]._timeOrig);
        int timeDest = solution[i]._timeDest;

        double val =  solution[i]._value;

        solutionRange[orig][timeOrig].push_back( ExpSomValue(dest, timeDest, val) );
    }

}

void ExpandedNetwork::afficheNetwork()
{
	for (ExpArc a : _arc)
	{
		cout << a._orig << " -> " << a._dest << ", cost " << a._cost << endl;
	}
}