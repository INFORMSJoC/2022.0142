FORMAT DES FICHIERS

--------------------------------------------------
nombre de jobs
nombre de ressources
Pour chaque job i :
	population(i) ES(i) LF(i)
	nombre_de_ressources_utilisees_par_i liste_des_ressources_utilisees_par_i
Pour chaque ressource e :
	capacite(e)
	
------------------------------------------------------
Remarque : les ressources sont numerotees a partir de 0
-------------------------------------------------------


Exemple de fichier avec commentaires : 

3              ----> 3 jobs
2              ----> 2 ressources
500 0 50       ----> le premier job a une population = 500, une earliest starting time = 0 et une latest finish time = 50
2 0 1          ----> le premier job utilise 2 ressources : la 0 et la 1
280 30 120     
1 0            ----> le second job utilise une ressource : la 0
600 25 150     
1 1            ----> le troisieme job utilise une ressource : la 1
20             -----> la premiere ressource a une capacite 20
25             -----> la seconde ressource a une capacite 25