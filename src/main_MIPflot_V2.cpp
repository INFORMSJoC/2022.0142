#include <iostream>
#include <fstream>
#include "Instance.h"

#include "MIP_flot.h"
#include "BranchBoundDicho.h"
#include "PL_Preempt_dicho.h"
#include "RCPSP_algo_V2.h"

#include "my_time.h"
#include <algorithm>

using namespace std;


//en entrée on a un fichier qui contient la liste des instances a resoudre :
//on fait la moyenne a la fin de ce main
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
		ficRes.open("resMIPFlotgen1_" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) + ".csv");
	}


	ficRes << "(" << argv[1] << "," << argv[2] << "," << argv[3] << ");"
		<< "BI_init" << ";"
		<< "BI_MIP" << ";"
		<< "BS_MIP" << ";"
		<< "nbCutLazy" << ";"
		<< "nbCutUser" << ";"
		<< "nbNode" << ";"
		<< "cpu(s)" << ";"
		<< "nbResolu"
		<< endl;




	double moy_BI_init = 0, moy_BI_MIP = 0, moy_BS_MIP = 0,
		moy_nbCutLazy = 0, moy_nbCutUser = 0, moy_nbNode = 0,
		moy_cpu = 0;

	int nbResolu = 0, nbMoy = 0;

	int nbIns = 10;

	for (int i = 0; i < nbIns; ++i)
	{

		nameCourt = "gen1_" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) + "_" + to_string(i) + ".txt";

		InstanceReduite ins;

		string name = "../../data/" + nameCourt;


		//on lit le fichier d'instance
		ins.lire(name);
		cout << " ============ " << nameCourt << " ============ " << endl;

		//============================================================
		// ====           heuristique flot pour avoir une borne initiale

		double BI_init = 0;
		
		RCPSP_Algo_V2 algo(&ins);

		//si jamais on veut une borne inf > 0 :
		//auto resRCPSP = algo.GRASP_V3(100, 0.1, true); //3 version du grasp pour les instances generales 
		BI_init = -1;// resRCPSP._margeFinale;

		

		//============================================================
		// ====          solution preemptive pour verif

		PL_Preempt_dicho PL_pre(&ins);
		double sol_pre = PL_pre.run();
		cout << " +++++  sol pre = " << sol_pre << endl;

		//============================================================
		// ====            MIP flot
		
		MIP_flot mipFlot(&ins);


		double BI_MIP = -1, BS_MIP = INT_INFINITY;

		my_time cpu_deb = give_time(); 
		
		//sans warm start : on passe des vecteurs vides
		auto resMip = mipFlot.creeEtResout({}, {}, {});

		//avec warm start
		//auto resMip = mipFlot.creeEtResout(resRCPSP._dateDebut, resRCPSP._debit, algo._flotCour);

		cout << "#nodes = " << mipFlot.getExploredNode() << endl;


		if (resMip.first > -EPSILON)//on a trouve une solution
			BI_MIP = resMip.first + 1; //+1 car dans l'algo on est oblige de faire "-1" pour les vieilles instances GEOSAFE
		if (resMip.second > -EPSILON)//on a trouve une borne sup
			BS_MIP = resMip.second + 1;
		
		my_time cpu_mip = give_time() - cpu_deb;

		ficRes << nameCourt << ";"
			<< BI_init << ";"
			<< BI_MIP << ";"
			<< BS_MIP << ";" 
			<< mipFlot._nbCutLazy << ";"
			<< mipFlot._nbCutUser << ";"
			<< mipFlot.getExploredNode() << ";"
			<< cpu_mip << ";"
			<< (abs(resMip.second - resMip.first) < EPSILON) << ";"
			<< endl;

		if (sol_pre < BI_MIP - EPSILON_DATE)
			stopProg("main: sol preemtive ou MIP flot faux");

		//--------- calcul moyennes (somme) --------------------
		
		if (BI_MIP > -0.5 && abs(BS_MIP - BI_MIP) < EPSILON)
			nbResolu++;

		//on ne fait la moyenne que sur les instances dont on obtient une BI et une BS
		if (BI_MIP > -EPSILON && BS_MIP > -EPSILON)
		{
			moy_BI_init += BI_init;
			moy_BI_MIP += BI_MIP;
			moy_BS_MIP += BS_MIP;
			moy_nbCutLazy += mipFlot._nbCutLazy;
			moy_nbCutUser += mipFlot._nbCutUser;
			moy_nbNode += mipFlot.getExploredNode();
			moy_cpu += cpu_mip;
			nbMoy++;
		}
	}

	//--------- calcul moyennes (division) et affichage --------------------

	ficRes << "(" << argv[1] << "," << argv[2] << "," << argv[3] << ");"
		<< moy_BI_init / nbMoy << ";"
		<< moy_BI_MIP / nbMoy << ";"
		<< moy_BS_MIP / nbMoy << ";"
		<< moy_nbCutLazy / nbMoy << ";"
		<< moy_nbCutUser / nbMoy << ";"
		<< moy_nbNode / nbMoy << ";"
		<< moy_cpu / nbMoy << ";"
		<< nbResolu << ";"
		<< endl;





	ficRes.close();


	return 0;
}


