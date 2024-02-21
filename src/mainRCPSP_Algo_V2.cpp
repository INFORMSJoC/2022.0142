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
		ficRes.open("resRCPSPgen1_" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) + ".csv");
	}

	ficRes << "(" << argv[1] << "," << argv[2] << "," << argv[3] << ");"
		<< "sol_pre" << ";"
		<< "cpu_pre" << ";"
		<< "solEtape1_rcpsp" << ";"
		<< "solFinale_rcpsp" << ";"
		<< "solLambda0_rcpsp" << ";"
		<< "cpu_rcpsp" << ";"
		<< "ok_rcpsp" << ";"
		<< "gap_pre_rcpsp" << ";"
		<< endl;


	double moy_sol_pre = 0;
	double moy_cpu_pre = 0;
	double moy_solEtape1_rcpsp = 0;
	double moy_solFinale_rcpsp = 0;
	double moy_solLambda0_rcpsp = 0;
	double moy_cpu_rcpsp = 0;
	double sum_ok_rcpsp = 0;
	double moy_gap_pre_rcpsp = 0;

	 
	//les instances sont par paquets de 10
	int nbIns = 10;

	for (int cpt = 0; cpt < nbIns; ++cpt)
	{
		nameCourt = "gen1_" + string(argv[1]) + "_" + string(argv[2]) + "_" + string(argv[3]) + "_" + to_string(cpt) + ".txt";

		InstanceReduite ins;

		string name = "../../data/" + nameCourt;


		//on lit le fichier d'instance
		ins.lire(name);

		//cout << cpt << " : " << endl;
		//srand(cpt * 17);
		//ins.genereRandom(20, 5, 0.4, 200);//4-> it. 207 

		//============================================================
		// ==== version preemptive (solution exacte)

		PL_Preempt_dicho PL_pre(&ins);
		my_time debut = give_time();
		double sol_pre = PL_pre.run();
		double cpu_pre = give_time() - debut;

		cout << "sol_pre = " << sol_pre << endl;

		//============================================================
		// solution non preemptive heuristique : "RCPSP flot"

		
		debut = give_time();

		RCPSP_Algo_V2 algo(&ins);
		
	
		//auto resRCPSP_1 = algo.GRASP(1000, 0.1);
		auto resRCPSP = algo.GRASP_V3(100, 0.1, true);
		//dessine(ins._nbJob, resRCPSP._dateDebut, resRCPSP._debit, ins._pop);

		//if (resRCPSP._margeFinale < resRCPSP_1._margeFinale - EPSILON_DATE -1)
		//	stopProg("main : pourquoi GRASP_V2 est moins bon que GRASP ?");
		//if (resRCPSP_1._margeFinale > 0 && resRCPSP._margeFinale < 0)
		//	stopProg("pas de sol ?");


		double cpu_rcpsp = give_time() - debut;
		bool ok_rcpsp = resRCPSP._margeFinale > -EPSILON;
		cout << " ========== GAP = " << max(0.0, 100*(sol_pre - resRCPSP._margeFinale) / sol_pre) << " =================" << endl;
		double gap_pre_rcpsp = 100 * (sol_pre - resRCPSP._margeFinale) / sol_pre;
		
		ficRes << setprecision(3) << fixed;

		ficRes << nameCourt << ";"
			<< sol_pre << ";"
			<< cpu_pre << ";"
			<< resRCPSP._bestLambda << ";"
			<< resRCPSP._margeFinale << ";"
			<< resRCPSP._marge0 << ";"
			<< cpu_rcpsp << ";"
			<< ok_rcpsp << ";"
			<< gap_pre_rcpsp << ";"
			<< endl;

		cout << "cpu = " << cpu_rcpsp << endl;

		//==== calcul des moyennes ===========
		if (ok_rcpsp)
		{
			moy_sol_pre += sol_pre;
			moy_cpu_pre += cpu_pre;
			moy_solEtape1_rcpsp += resRCPSP._bestLambda;
			moy_solFinale_rcpsp += resRCPSP._margeFinale;
			moy_solLambda0_rcpsp += resRCPSP._marge0;
			moy_cpu_rcpsp += cpu_rcpsp;
			sum_ok_rcpsp++;
			moy_gap_pre_rcpsp += gap_pre_rcpsp;
		}

	}

	//=================  moyennes  ========================
	ficRes << "moy_(" << argv[1] << "," << argv[2] << "," << argv[3] << ");" 
		<< moy_sol_pre / sum_ok_rcpsp << ";"
		<< moy_cpu_pre / sum_ok_rcpsp << ";"
		<< moy_solEtape1_rcpsp / sum_ok_rcpsp << ";"
		<< moy_solFinale_rcpsp / sum_ok_rcpsp << ";"
		<< moy_solLambda0_rcpsp / sum_ok_rcpsp << ";"
		<< moy_cpu_rcpsp / sum_ok_rcpsp << ";"
		<< sum_ok_rcpsp << ";"
		<< moy_gap_pre_rcpsp / sum_ok_rcpsp << ";"
		<< endl;


	ficRes.close();
	cout << " moy_solFinale_rcpsp = " << moy_solFinale_rcpsp/nbIns << endl;
	cout << " moy_gap_pre_rcpsp = " << moy_gap_pre_rcpsp / sum_ok_rcpsp << " % " << endl;

	return 0;
}


