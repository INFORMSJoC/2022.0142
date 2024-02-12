cmake_minimum_required(VERSION 3.12)


#variables globales (a mettre avant project)
SET(CMAKE_CXX_STANDARD 17)

#en ligne de commande : cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=/opt/rh/devtoolset-6/root/usr/bin/g++
SET(CMAKE_BUILD_TYPE Release) # RelWithDebInfo, Debug, RelWithDebInfo, Release
#SET(CMAKE_CXX_COMPILER /opt/rh/devtoolset-6/root/usr/bin/g++)


project(Geosafe)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Créer une variable avec les fichiers à compiler
set(SRC_RCPSP
	BSmarge.cpp
	ExpandedNetwork.cpp
	ExpandedNetworkChemin.cpp
	GrapheReduit.cpp
	Instance.cpp
	InstanceReduite.cpp
	ModeleLineaire.cpp
	PL_Optim_RCPSP.cpp
	PL_RCPSP_fixeDead.cpp
	RCPSP_Algo.cpp
	RCPSP_Graphe.cpp
	BranchBoundDicho.cpp
	Solution.cpp
	PL_Preempt_dicho.cpp
	RCPSP_algo_V2.cpp
	MIP_flot_partiel.cpp
	mainRCPSP_Algo_V2.cpp
#	mainRCPSP_Algo_V2_vielleIns.cpp
    )
	
#facultatif, permet de recompiler si on fait des chgent dans les .h
set(HDR_RCPSP
    BSmarge.h
	ExpandedNetwork.h
	ExpandedNetworkChemin.h
	GrapheReduit.h
	Instance.h
	InstanceReduite.h
	ModeleLineaire.h
	PL_Optim_RCPSP.h
	PL_RCPSP_fixeDead.h
	RCPSP_Algo.h
	RCPSP_Graphe.h
	util.h
	utilVector.h
	BranchBoundDicho.h
	PL_Preempt_dicho.h
	RCPSP_algo_V2.h
	MIP_flot_partiel.h
	Solution.h
    )	
	
	# Créer une variable avec les fichiers à compiler
set(SRC_BB
	BSmarge.cpp
	ExpandedNetwork.cpp
	ExpandedNetworkChemin.cpp
	GrapheReduit.cpp
	Instance.cpp
	InstanceReduite.cpp
	ModeleLineaire.cpp
	PL_Optim_RCPSP.cpp
	PL_RCPSP_fixeDead.cpp
	RCPSP_Algo.cpp
	RCPSP_Graphe.cpp
	BranchBoundDicho.cpp
	Solution.cpp
	RCPSP_algo_V2.cpp
	MIP_flot_partiel.cpp
	mainBBdicho_V2.cpp
#	mainBBdicho_V2_vieille.cpp
    )
	
#facultatif, permet de recompiler si on fait des chgent dans les .h
set(HDR_BB
    BSmarge.h
	ExpandedNetwork.h
	ExpandedNetworkChemin.h
	GrapheReduit.h
	Instance.h
	InstanceReduite.h
	ModeleLineaire.h
	PL_Optim_RCPSP.h
	PL_RCPSP_fixeDead.h
	RCPSP_Algo.h
	RCPSP_Graphe.h
	util.h
	utilVector.h
	BranchBoundDicho.h
	Solution.h
	MIP_flot_partiel.h
	RCPSP_algo_V2.h
    )	
	
	
# Créer une variable avec les fichiers à compiler
set(SRC_Flot
	Instance.cpp
	InstanceReduite.cpp
	MIP_flot.cpp
	RCPSP_algo_V2.cpp
	PL_Preempt_dicho.cpp
 	MIP_flot_partiel.cpp
	main_MIPflot_V2.cpp
    )
	
#facultatif, permet de recompiler si on fait des chgent dans les .h
set(HDR_Flot
	Instance.h
	InstanceReduite.h
	util.h
	utilVector.h
	MIP_flot.h
	RCPSP_algo_V2.h
 	MIP_flot_partiel.h	
	PL_Preempt_dicho.h
	callback_MIP_Flot.h
    )	
	


# creation d'un executable exe a partir des fichiers cpp 
add_executable(BB ${SRC_BB} ${HDR_BB})
add_executable(flowHeur ${SRC_RCPSP} ${HDR_RCPSP})
add_executable(flowMip ${SRC_Flot} ${HDR_Flot})

# Si ne fonctionne pas, il faut ajouter le chemin de ILOG à la ligne 19 de FindCPLEX.cmake 
#du dossier cmake
find_package(CPLEX)
if (TARGET cplex-library)
  target_link_libraries(flowMip cplex-library cplex-concert ilocplex cplex-cp)
  target_link_libraries(flowHeur cplex-library cplex-concert ilocplex cplex-cp)
  target_link_libraries(BB cplex-library cplex-concert ilocplex cplex-cp)
endif ()

