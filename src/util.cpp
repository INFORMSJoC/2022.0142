
#include "util.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;


void dessine(int N, const vector<double> & debut, const vector<double> & debit, const vector<int> & pop)
{
	vector<vector<bool>> occupe;
	vector<pair<double, int> > coord(N);

	//0. taille de occupe
	int maxX = 0, maxY = 0;

	// 1. on trie par date de debut croissante
	vector<pair<double, int> > v;
	for (int i = 0; i < N; ++i)
	{
		v.push_back({ debut[i],i });
		maxX = max(maxX, (int)(debut[i] + (double)(pop[i]) / debit[i])+1);
		maxY += debit[i];
	}
	sort(v.begin(), v.end());

	occupe.assign(maxX+N, vector<bool>(maxY+N, false));

	// 2. on trouve (x,y) coordonnee du point inf gauche pour chaque job
	for (int i = 0; i < N; ++i)
	{
		int job = v[i].second;
		bool ok = false;
		int h = 0;
		int x = (int)(debut[job]);
		int debit_job = (int)(debit[job]);
		while (!ok)
		{
			while (occupe[x][h])
				h++;
			int y = h;
			while (y < h + debit_job && !occupe[x][y])
				y++;
			if (y == h + debit_job)//on a suffisament de place en hauteur pout mettre le job en h
			{
				ok = true;
				coord[job] = { debut[job] ,h };

				//on met a jour la yable occupe
				for (auto p : v)
				{
					int xcour = (int)(p.first);
					if (xcour >= debut[job] -1 && xcour <= debut[job] +1 + (double)(pop[job]) / debit[job])
					{
						for (int ycour = h; ycour <= h + debit[job]; ++ycour)
							occupe[xcour][ycour] = true;
					}
				}
			}
			else
				h = y;
		}
	}

	//3. ecrire le fichier
	ofstream fic("data.txt");

	for (int i = 0; i < N; ++i)
		fic << debut[i] << " ";
	fic << endl;
	for (int i = 0; i < N; ++i)
		fic << coord[i].second << " ";
	fic << endl;
	for (int i = 0; i < N; ++i)
		fic << debit[i] << " ";
	fic << endl;
	for (int i = 0; i < N; ++i)
		fic <<(double)(pop[i])/debit[i] << " ";
	fic << endl;

	fic.close();


	system("python dessinRectangle.py data.txt");
}