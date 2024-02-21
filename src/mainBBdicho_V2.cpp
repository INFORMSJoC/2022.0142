# include <iostream>
#include <fstream>
#include "InstanceReduite.h"
#include "BranchBoundDicho.h"
#include "PL_Preempt_dicho.h"
#include "RCPSP_algo_V2.h"

#include "my_time.h"
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{

	ofstream ficRes;
	string nameCourt;

	if (argc != 4)
	{
		cout << "bad usage: exe nbJob nbArc proba" << endl;
		exit(-1);
	}
	else
	{
		ficRes.open("resBBgen1_" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) + ".csv");
	}

	ficRes << "(" << argv[1] << "," << argv[2] << "," << argv[3] << ");"
		<< "BI_init" << ";"//borne inf fourni au BBdicho
		<< "BI_BB" << ";"// exacte si BB fini sinon borne inf
		<< "BS_BB" << ";"// exacte si BB fini sinon borne sup
		<< "NbNoeudExplore" << ";" 
		<< "NbBranchType1" << ";"
		<< "nbBranchType2ou3" << ";"
		<< "cutMoy" << ";"
		<< "paralMoy" << ";"
		<< "paralMax" << ";"
		<< "compress" << ";"
		<< "finProcess" << ";"
		<< "cpu_BB" 
		<< endl;


	double moy_BI_init = 0;
	double moy_BI_BB = 0;
	double moy_BS_BB = 0;
	double moy_NbNoeudExplore = 0;
	double moy_NbBranchType1 = 0;
	double moy_NbBranchType2ou3 = 0;
	double moy_cutMoy = 0;
	double moy_cpu_BB = 0;
	double moy_paralMoy = 0;
	double moy_paralMax = 0;
	double moy_compress = 0;
	double moy_finProcess = 0;

	int nbSolObtenue = 0;


	//les instances sont par paquets de 10
	int nbIns = 10;
	for (int cpt = 0; cpt < nbIns; ++cpt)
	{
		nameCourt = "gen1_" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) + "_" + to_string(cpt) + ".txt";

		InstanceReduite ins;

		string name = "../../data/" + nameCourt;

		//on lit le fichier d'instance
		ins.lire(name);

		
		double solRCPSP = -1;
		double BI_BB = -1;
		double BS_BB = -1;

		double nbCut = 0, nbAppelCut = 0;
		my_time cpu_BB = 0;


		bool stat = false; // stats calcule ? 
		bool RCPSPOK = false; //vrai si RCPSP trouve une solution


		//============================================================
		// solution non preemptive heuristique : "RCPSP flot"

		double cpuRCPSP = 0;
		my_time debut = give_time();

		RCPSP_Algo_V2 algo(&ins);
		
		//auto resRCPSP = algo.GRASP(1000, 0.1); //1ere version du grasp pour les instances generales => pas terrible
		auto resRCPSP = algo.GRASP_V3(100, 0.1, true); //3 version du grasp pour les instances generales 

		cpuRCPSP = give_time() - debut;
		solRCPSP = resRCPSP._margeFinale;

		cout << "solution heuristique (Borne inf) = " << solRCPSP << endl;

		//============================================================
		// solution non preemptive exacte : on peut donner le cout de la sol heuristique pour aider le B&B dans BI_init

		my_time deb = give_time();
		double BI_init = solRCPSP;
		BranchBoundDicho BBdicho(0, &ins, BI_init);
		pair<double, double> sol = BBdicho.run(3600);
		cpu_BB = give_time() - deb;

		//si on a explore un seul noeud : on n'a pas de solution sauvegardee dans le BB : 
		//la solution heuristique du RCPSP était optimale => on la reprend ????


		//verif solution 
		//if (BBdicho._hasSol && sol.first > 0 && !ins.isSolution(BBdicho._bestDebut, BBdicho._bestDebit))
		//	stopProg("main : LE RESULTAT DE BB EST FAUX");


		//sauvegarde info sol pour affichage
		BI_BB = sol.first;
		BS_BB = sol.second;
		nbCut = BBdicho.getNbCut();
		nbAppelCut = BBdicho.getNbAppelCut();
		double paralMoy;
		int paralMax;

		//si le BB n'a pas de solution on reprend celle du RCPSP (s'il en a une)
		if (!BBdicho._hasSol && solRCPSP > 0)
		{
			BBdicho._bestDebut = resRCPSP._dateDebut;
			BBdicho._bestDebit = resRCPSP._debit;
		}

		auto p = BBdicho.statParallelisme();
		paralMoy = p.first;
		paralMax = p.second;

		double compress = 0;
		double finProcess = 0;
		for (int i = 0; i < ins._nbJob; ++i)
		{
			compress += ins._debitMax[i] / BBdicho._bestDebit[i];
			finProcess = max(finProcess, BBdicho._bestDebut[i] + static_cast<int>(ins._pop[i]) / BBdicho._bestDebit[i]);
		}

		


		cout << "solution exacte in " << BI_BB << " ; " << BS_BB << "en " << cpu_BB << "sec" << endl;


		double cutMoy = nbAppelCut = 0 ? 0 : (nbCut / nbAppelCut);
		compress /= ins._nbJob;// ATTENTION : premiers resultats faux car il manquait cette division (ajoutee 15/03/2023)

		ficRes << setprecision(3) << fixed;
		ficRes << nameCourt << ";"
			<< BI_init << ";"  // borne inf
			<< BI_BB << ";"    // exacte si BB fini sinon borne inf
			<< BS_BB << ";"    // exacte si BB fini sinon borne sup
			<< BBdicho.getNbNoeudExplore() << ";"
			<< BBdicho.getNbBranchType1() << ";"
			<< BBdicho.getnbBranchType2ou3() << ";"
			<< cutMoy << ";"
			<< paralMoy << ";"
			<< paralMax << ";"
			<< compress << ";"
			<< finProcess << ";" 
			<< cpu_BB << ";"
			<< endl;


		//==== calcul des moyennes ===========
		//
		if (BI_BB > -EPSILON)
		{
			moy_BI_init += BI_init;
			moy_BI_BB += BI_BB;
			moy_BS_BB += BS_BB;
			moy_NbNoeudExplore += BBdicho.getNbNoeudExplore();
			moy_NbBranchType1 += BBdicho.getNbBranchType1();
			moy_NbBranchType2ou3 += BBdicho.getnbBranchType2ou3();
			moy_cutMoy += cutMoy;
			moy_cpu_BB += cpu_BB;
			moy_paralMoy += paralMoy;
			moy_paralMax += paralMax;
			moy_compress += compress;
			moy_finProcess += finProcess;
			nbSolObtenue++;
		}

	}

	//=================  moyennes  ========================
	ficRes << nameCourt + "_moy" << ";"
		<< moy_BI_init / nbSolObtenue << ";"//borne inf
		<< moy_BI_BB / nbSolObtenue << ";"// exacte si BB fini sinon borne inf
		<< moy_BS_BB / nbSolObtenue << ";"// exacte si BB fini sinon borne sup
		<< moy_NbNoeudExplore / nbSolObtenue << ";"
		<< moy_NbBranchType1 / nbSolObtenue << ";"
		<< moy_NbBranchType2ou3 / nbSolObtenue << ";"
		<< moy_cutMoy / nbSolObtenue << ";"
		<< moy_paralMoy / nbSolObtenue << ";"
		<< moy_paralMax / nbSolObtenue << ";"
		<< moy_compress / nbSolObtenue << ";"
		<< moy_finProcess / nbSolObtenue << ";"
		<< moy_cpu_BB / nbSolObtenue << ";"
		<< endl;

	ficRes.close();

	return 0;
}

