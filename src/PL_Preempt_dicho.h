#pragma once

#include "InstanceReduite.h"
//=============================================================================
// classe creee en janvier 2023 : reprend l'algorithme calculMargePLdicho de la classe BSMarge
// en l'adaptant sur une instance de type InstanceReduite
// dans le but de pouvoir fournir des resultats sur des instances "generales" en 
// reponse aux reviewers d'INFORMS
// ========================================================================


class PL_Preempt_dicho
{
private:
	InstanceReduite * _ins;

public:
	//methode principale : lance l'algo (recherche meilleure marge par dicho, resolution d'un PL pour determiner si une marge est realisable)
	double run(double precision = 0.0001);

	PL_Preempt_dicho(InstanceReduite * ins) : _ins(ins) {}

private:
	double calculMargeCibleMax();

	//on ordonne l'ensemble des dates ESi, LFi-margeCible par ordre croissant
	//renvoie le vecteur des valeurs triees
	vector<double> trierDate(double margeCible);

	//on lance un PL pour savoir si toutes les pop peuvent evacuer avec une marge >= margeCible  
	bool isMargeRealisable(double margeCible, vector<double> & date);
};