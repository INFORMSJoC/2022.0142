#pragma once

#include "Instance.h"

//classe qui stocke une solution non preemptive
//donne pour chaque job on donne la date de début (départ du noeud d'evac), 
//et le débit d'evacuation 

class Solution
{
public:
	Instance * _ins;

	vector<double> _debit;
	vector<double> _debut; //_debut[i] = -1 si job non encore planifie (dans le cas d'une solution partielle)
	

	Solution(Instance * ins) : _ins(ins) {}

	bool verifieDebit();



	void init();

private:
	bool ajouterArcDebit(vector<debitTemporise> & arcDebit, const int CAP, double debit, double debut, double fin);

};