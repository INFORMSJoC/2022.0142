#ifndef EXPANDEDNETWORK_CHEMIN_H
#define EXPANDEDNETWORK_CHEMIN_H

#include "Instance.h"
#include "ExpandedNetwork.h"

//===================================================================================================
// construit l expanded network dans le cas ou on conneit les chemins d 'evacuation (donnes dans l instance)
//===============================================================================


class ExpandedNetworkChemin : public ExpandedNetwork
{
public:

	//on redefnit createExpArc car on ne va utiliser que les chemins deja definis et les fenetres de temps associes aux noeuds
	void createExpArc() override;

	ExpandedNetworkChemin(Instance * ins) : ExpandedNetwork(ins) {}

	//retourne vrai si l'arc (prec,t) -> (cour,s) existe deja
	bool existe(int prec, int t, int cour, int s);

};

#endif // EXPANDEDNETWORK_CHEMIN_H














