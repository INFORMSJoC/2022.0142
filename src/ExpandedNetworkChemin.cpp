#include "ExpandedNetworkChemin.h"


//on redefnit createExpArc car on ne va utiliser que les chemins deja definis et les fenetres de temps associes aux noeuds
void ExpandedNetworkChemin::createExpArc()
{

	size_t Nnodes = _ins->_nbNode;
	int T = _ins->_TMAX;

	int cpt = 0;

	_arcIn.resize(Nnodes + 1, vector<vector<int>>(T + 1));
	_arcOut.resize(Nnodes + 1, vector<vector<int>>(T + 1));

	// source -> evac nodes

	_firstInputExpArc = cpt;
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		_arc.emplace_back(_source, ExpNode(i, 0));
		_arcIn[i][0].push_back(cpt);
		cpt++;
	}

	// safe nodes -> sink
	_lastInputExpArc = cpt - 1;
	_firstOutputExpArc = cpt;
	for (int i = _ins->_firstSafeNode; i <= _ins->_lastSafeNode; ++i)
	{
		_arc.emplace_back(ExpNode(i, T), _sink);
		_arcOut[i][static_cast<size_t>(T)].push_back(cpt);
		cpt++;
	}



	// arc station uniquement pour les evac node et safe

	_lastOutputExpArc = cpt - 1;
	_firstStationExpArc = cpt;


	// (i , t) -> (i , t+1), i = evac node
	for (int i = _ins->_firstEvacuationNode; i <= _ins->_lastEvacuationNode; ++i)
	{
		int deb = _ins->_chemins[i][0]._ES;
		int fin = _ins->_chemins[i][0]._LF;

		for (int t = deb; t < fin; ++t)
		{
			_arc.emplace_back(ExpNode(i, t), ExpNode(i, t + 1));

			_arcOut[i][static_cast<size_t>(t)].push_back(cpt);
			_arcIn[i][static_cast<size_t>(t + 1)].push_back(cpt);

			cpt++;
		}
	}

// (i , t) -> (i , t+1), i = safe node
	for (int i = _ins->_firstSafeNode; i <= _ins->_lastSafeNode; ++i)
	{

		for (int t = 0; t < _ins->_TMAX; ++t)
		{
			_arc.emplace_back(ExpNode(i, t), ExpNode(i, t + 1));

			_arcOut[i][static_cast<size_t>(t)].push_back(cpt);
			_arcIn[i][static_cast<size_t>(t + 1)].push_back(cpt);

			cpt++;
		}
	}

	
	_lastStationExpArc = cpt - 1;




	// arc route : (i,t)-> (j,t+time(i,j)) i,j = any node, j != i and cap[i][j] > 0
	
	_firstRouteExpArc = cpt;

	for (int i = 1; i < _ins->_chemins.size(); ++i)
	{
		int prec = _ins->_chemins[i][0]._som;
		

		for (int j = 1; j < _ins->_chemins[i].size(); ++j)
		{
			int cour = _ins->_chemins[i][j]._som;
			int deb = _ins->_chemins[i][j-1]._ES;
			int fin = _ins->_chemins[i][j-1]._LF;

			//les LF sont cacule de sorte d'etre suffisant => on n'utilise plus les deadlines, 
			//fin = min(_ins->_dead[prec][cour], fin);

			for (int t = deb; t <= fin; ++t)
			{
				// il se peut que l'arc (prec, t) -> (cour, t + _ins->_time[prec][cour]) existe deja (d un autre chemin), dans ce cas on ne l'ajoute pas

				if ( !existe( prec, t, cour, t + _ins->_time[prec][cour] ) )
				{
					ExpArc arc ( ExpNode(prec, t), ExpNode(cour, t + _ins->_time[prec][cour]) );


					_arc.emplace_back(arc);

					_arcOut[prec][static_cast<size_t>(t)].push_back(cpt);
					_arcIn[cour][static_cast<size_t>(t + _ins->_time[prec][cour])].push_back(cpt);

					cpt++;
				}
			}
			prec = cour;
		}
	}


	_lastRouteExpArc = cpt - 1;

	//cout << " ====>>>> arc exp = " << _arc.size() << endl;
}


//retourne vrai si l'arc (prec,t) -> (cour,s) existe deja
bool ExpandedNetworkChemin::existe(int prec, int t, int cour, int s )
{
	bool isIn = false;

	vector<int> & v = _arcOut[prec][t];

	for (int i = 0; !isIn && i < v.size(); ++i)
	{
		int suiv = _arc[v[i]]._dest._nodeIndex;
		int tsuiv = _arc[v[i]]._dest._time;

		if (suiv == cour && tsuiv == s)
			isIn = true;
	}

	return isIn;
}