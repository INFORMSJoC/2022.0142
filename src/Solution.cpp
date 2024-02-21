#include "Solution.h"

void Solution::init()
{
	int N = _ins->_lastEvacuationNode;

	_debit.assign(N + 1, -1);
	_debut.assign(N + 1, -1);
}

bool Solution::verifieDebit()
{
	bool ok = true; 
	double minMarge;
	vector<vector<vector<debitTemporise>>> arcDebit (_ins->_nbNode+1, vector<vector<debitTemporise>>(_ins->_nbNode+1));

	for (int i = _ins->_firstEvacuationNode; ok && i <= _ins->_lastEvacuationNode; ++i)//pour chaque chemin i
	{
		double arrivePremierEnCour = _debut[i];//arrivee du premier en cour

		if (arrivePremierEnCour != -1) //pour les solutions partielles : seulement une partie des jobs est planifie
		{
			int prec = _ins->_chemins[i][0]._som;

			for (int j = 1; ok && j < _ins->_chemins[i].size(); ++j)//pour chaque arc du chemin prec->cour
			{
				int cour = _ins->_chemins[i][j]._som;


				//arrivee du premier au sommet  cour
				arrivePremierEnCour += _ins->_time[prec][cour];
				double arriveDernierEnCour = arrivePremierEnCour + static_cast<double>(_ins->_pop[i]) / _debit[i] - 1;

				//on ajoute (_ins->_pop[i])/_debit[i] - 1 pour avoir la date du dernier arrivee
				if (arriveDernierEnCour > _ins->_chemins[i][j]._LF + EPSILON)
				{
					cout << " arriveDernierEnCour = " << arriveDernierEnCour << " LF = " << _ins->_chemins[i][j]._LF << endl;
					ok = false;
				}
				else
				{
					minMarge = min(minMarge, _ins->_chemins[i][j]._LF - arriveDernierEnCour);

					//ajouterArcDebit :  ajoute le debit et  verifie en meme temps qu on ne depasse pas la capacite de l'arc
					//attention on passe sur l'arc avant d'arriver en cour => on retranche donc la longuer de l'arc
					ok = ajouterArcDebit(arcDebit[prec][cour], _ins->_capArc[prec][cour], _debit[i],
						arrivePremierEnCour - _ins->_time[prec][cour], arriveDernierEnCour + 1 - _ins->_time[prec][cour]);// 23/10/2019 : ajout du +1 
				}
				prec = cour;
			}
		}
	}
	return ok;
}



//ajoute un nouveau debit sur l'intervalle [debut,fin]
//le vecteir arcDebit contient deja une liste d'intervalle [deb,fin] avec une valeur de debit pour chaque intervalle
//il faut donc ajouter le nouvel interval[debut, fin] avec son debit en gardant les intervalles de arcDebit sans intersection et trie sur les deb croissants: 
// si l'intervalle est contenu dans un autre il faut faire la somme des debit, s'il est à cheval sur un autre il faut decouper en intervalle plus petits etc...
//retourne vrai si on peut ajouter le nouveau debit sans jamais depasser la capacité de l'arc (orig, dest), faux sinon
bool Solution::ajouterArcDebit(vector<debitTemporise> & arcDebit, const int CAP, double debit, double debut, double fin)
{
	bool ok = true;
	bool stop = false;//stop devient vrai quand on a fini de mettre la liste a jour




	if (CAP < debit - EPSILON_DEBIT)
		ok = false;

	//trivial : liste vide ou nouvel interval avant premier de la liste ou apres dernier de la liste
	if (arcDebit.size() == 0 || fin - EPSILON <= arcDebit[0]._deb || debut >= arcDebit[arcDebit.size() - 1]._fin - EPSILON)
	{
		arcDebit.push_back(debitTemporise(debut, fin, debit));
		stop = true;
	}


	// 2. sinon [debut, fin] intersecte un ou plusieurs intervalles : on va avancer dans la liste arcDebit
	//  et on va reduire [debut, fin] au fur et a mesure qu'on intersecte des intervalles dans arcDebit : 
	// [debut, fin] va donc etre decoupe au fur et a mesure, a chaque decoupage deb est augmente (on le renomme t)

	double t = debut; //debut du prochain intervalle
	vector<int> intervalleAsupp; //liste des intervalles a supprimer (on pourrait le faire en cours de creation des nouveaux mais ca complique l'algo)


	int i = 0; //indice de l'intervalle de arcDebit qu'on intersecte au tour de boucle courant

	int sizeInit = static_cast<int>(arcDebit.size());

	//2.1. on cherche le premier intervalle susceptible de s'intersecter avec [debut, fin]
	while (i < sizeInit && t >= arcDebit[i]._fin - EPSILON)
		i++;

	//2.2. tant qu'on a pas intergre le nouvel intervalle [debut, fin] avec son debit on continue
	while (ok && !stop && t + EPSILON < fin)
	{

		//--------------------------------------------------------------------------------------------------------
		//1er cas : t est au milieu d'un intervalle initial existant => l'intervalle existant va etre scinde en 
		// 2 ou 3 parties 
		if (i < sizeInit && t >= arcDebit[i]._deb - EPSILON)
		{
			//i sera supprime car on cree plusieurs intervalles a sa place
			intervalleAsupp.push_back(i);

			//1er morceau : de arcDebit[i]._deb a t si t != arcDebit[i]._deb 
			if (t > arcDebit[i]._deb + EPSILON)
				arcDebit.push_back(debitTemporise(arcDebit[i]._deb, t, arcDebit[i]._debit));

			//les deux fin tombent en meme temps = > on a juste un deuxieme intervalle a creer
			if (abs(fin - arcDebit[i]._fin) < EPSILON)
			{
				arcDebit.push_back(debitTemporise(t, fin, arcDebit[i]._debit + debit));
				ok = arcDebit[i]._debit + debit <= CAP + EPSILON_DEBIT;
				stop = true;
			}
			else //une fin precede l'autre 
			{
				//l intervalle [t, fin] est entierement dans l intervalle i => on va creer 2 autres morceaux
				if (fin < arcDebit[i]._fin)
				{
					//2sd morceau
					arcDebit.push_back(debitTemporise(t, fin, arcDebit[i]._debit + debit));
					ok = arcDebit[i]._debit + debit <= CAP + EPSILON_DEBIT;

					//3eme morceau
					arcDebit.push_back(debitTemporise(fin, arcDebit[i]._fin, arcDebit[i]._debit));

					//et on a fini : on a intergalement ajoute le nouvel intervalle
					stop = true;
				}
				else//l intervalle[t, fin] fini apres l intervalle i = > on va creer 1 autre morceau
				{
					// 2sd morceau
					arcDebit.push_back(debitTemporise(t, arcDebit[i]._fin, arcDebit[i]._debit + debit));
					ok = arcDebit[i]._debit + debit <= CAP + EPSILON_DEBIT;

					t = arcDebit[i]._fin;//on avance le nouveau debut courant
					i++;
				}
			}
		}
		//--------------------------------------------------------------------------------
		// 2sd cas t est avant intervalle i et apres intervalle i-1
		else
		{
			//debut du prochain intervalle deja dans la liste s'il existe
			double debutSuivant = TIME_INFINITY;
			if (i < sizeInit)
				debutSuivant = arcDebit[i]._deb;

			if (fin < debutSuivant) //on n intersecte pas d'intervalle
			{
				arcDebit.push_back(debitTemporise(t, fin, debit));
				stop = true;
			}
			else//on intersecte l'intervalle i+1 => on cree le l'intervalle qui s'arrete avant i+1 
			{
				arcDebit.push_back(debitTemporise(t, arcDebit[i]._deb, debit));
				t = arcDebit[i]._deb;

			}

		}


	}


	//====================================================================================
	//on fait le menage et on retrie en fonction des dates de debut

	for (int k = 0; k < intervalleAsupp.size(); ++k)
	{
		//attention au fur et a mesure qu'on enleve on decale les indices d ou le -k
		arcDebit.erase(arcDebit.begin() + intervalleAsupp[k] - k);
	}

	sort(arcDebit.begin(), arcDebit.end());

	return ok;
}
