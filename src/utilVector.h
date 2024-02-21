#ifndef _UTIL_V_
#define _UTIL_V_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

//quelques fonctions utiles pour la manipulation de vecteurs


//supprime les doublons dans un vecteur d'entiers v deja trie
inline void supprimeDoublon(vector<int> & v)
{
	vector<int>::iterator it = unique(v.begin(), v.end());
	v.resize(distance(v.begin(), it));

}

//comparaison utilisee dans supprimeDoublon(vector<double> & v)
inline bool myComp(double i, double j)
{
	return (abs(i-j) < EPSILON);
}

//supprime les doublons ( a epsilon pres ) dans un vecteur de double v deja trie
inline void supprimeDoublon(vector<double> & v)
{

	vector<double>::iterator it = unique(v.begin(), v.end(), myComp);
	v.resize(distance(v.begin(), it));

}

#endif

