//modello SIR con arrival time e path

#include <map>
#include <vector>
#include <string>
#include <utility>

#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
//#include "command_line.h"

int i,j,infettiW,infettiH;
int stampaInfetti = 0;
double wij,TRout,TRin,s;
int cityseed = 1; 		//città seed di partenza Ikk (cityseed = 1 -> (I)11 seed)

double R0 = 1.5;
double mu = 1/3.0;	// rate di guarigione [day^-1]           R0=beta/mu>1 per avere outbreak
double beta = R0*mu;	// trasmissibilità (probabilità di trasmissione dell'infezione)
int infettiIniz = 10;
int t;
int nrun = 100;
pair<int, int> wiw;	//path who infects who

map<int,vector<int> > neighbours;	//mappa i vicini di un nodo
map<pair<int,int>,double> weight;	//mappa il peso del link tra due nodi
map<int,int> Npop;			//mappa la popolazione di ogni link
map<int,vector<double> > trafficOUTmap;	//mappa il peso di tutti i link uscenti dalla città i-esima
map<int,vector<double> > trafficINmap;	//mappa il peso di tutti i link uscenti dalla città i-esima
map<int,int> trafficOUT;		//mappa il traffico uscente dalla città i-esima
map<int,int> trafficIN;			//mappa il traffico entrante nella città i-esima
map<int,int> diag;			//mappa della diagonale della matrice dei suscettibili diag[i]=(N)i-(TRout)i
map<int,int> NpopW;			//mappa il numero di persone presenti in i durante il work-time
map<int,int> NpopH;			//mappa il numero di persone presenti in i durante l'home-time
map<int,double> forceW;			//mappa la force of infection di ogni nodo nel WORK TIME
map<int,double> forceH;			//mappa la force of infection di ogni nodo nell'HOME TIME

vector<double> NULLrows;	//righe di arrayNULL
vector<double> Srows;		//righe di arrayS
vector<double> Irows;		//righe di arrayI
vector<double> Rrows;		//righe di arrayR
vector<double> deltaStot;	//variazione infetti alla fine della giornata
vector<double> Stemp;		//variabile temporanea
vector<pair<int,int> >path;	//path di infezione

vector<int>::size_type k; //istanzio oggetto di tipo (vector) size

vector<vector<double> > arrayNULL; //matrice identicamente nulla
vector<vector<double> > deltaS;
vector<vector<double> > deltaI;
vector<vector<double> > deltaR;
vector<vector<double> > arrayS;	//matrice S
vector<vector<double> > arrayI;	//matrice I
vector<vector<double> > arrayR;	//matrice R
vector<vector<double> > arraySr;//matrice Srun
vector<vector<double> > arrayIr;//matrice Irun
vector<vector<double> > arrayRr;//matrice Rrun


map<int,vector<int> >::iterator iter1;		//iteratore per la mappa dei vicini
map<pair<int,int>,double>::iterator iter2;	//iteratore per la mappa dei pesi
map<int,vector<double> >::iterator iter3;	//iteratore per la mappa trafficOUTmap e trafficINmap
map<int,int>::iterator iter4;			//iteratore per la mappa del traffico uscente/entrante e di Npop




int main (int argc, char * argv[]){


//	 if (!read_command_line(argc,argv)){

//		cout << "\t*Error: parameters must be set via command-line!" << endl;
//		cout << "\texample of usage: ./provaSIR -n100 -s268" << endl;
//		return -1;
//	}

	//cityseed = 1;

	/*APRO E LEGGO IL FILE DELLA POPOLAZIONE*/

	ifstream popFILE ("./pop_boro.txt");
	if ( popFILE.is_open() ){
		cout<<"Open pop file \t"<<endl;
		while ( popFILE.good() ){

			while( popFILE >> i >> j ){
				Npop[i] = j;
			}

			popFILE.close();

		}

	}else cout << "\nUnable to open the population file"<< endl;


	/*APRO E LEGGO IL FILE CON LA RETE*/

	ifstream networkFILE ("./commuting_boro.txt");
	if ( networkFILE.is_open() ){

		while ( networkFILE.good() ){

			while( networkFILE >> i >> j >> wij ){
				neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
				weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
				trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita
				trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata

			}

			networkFILE.close();

		}

	}else cout << "\nUnable to open the network file"<< endl;


	for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){

		TRout = 0;

		for( k = 0; k<(iter3->second).size(); k++ ){
			TRout = TRout + iter3->second[k];
		}

		trafficOUT[iter3->first] = TRout;				// ho la mappa del traffico in uscita di ogni nodo
	}




	for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

		TRin = 0;

		for( k = 0; k<(iter3->second).size(); k++ ){
			TRin = TRin + iter3->second[k];
		}

		trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
	}


	for( int i = 0; i < Npop.size(); i++ ){					// COSTRUISCO L'ARRAY IDENTICAMENTE NULLA

		for( int j = 0; j < Npop.size(); j++ ){
			NULLrows.push_back(0);
		}

		arrayNULL.push_back(NULLrows);
		NULLrows.clear();
	}

	double AT[Npop.size()];
	double ATtemp[Npop.size()];


	/***************************************INIZIALIZZO LA MATRICE DEI SUSCETTIBILI*************************************************/


	for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){
		diag[iter4->first] = Npop[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i
	}


	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

		for( int j = 0; j < Npop.size(); j++ ){

			if(i != j){
				Srows.push_back(weight[pair<int,int>(i,j)]);
			}else if(i==cityseed){
				Srows.push_back(diag[i]-infettiIniz);		// tolgo dalla città seed il numero di infetti iniziale
			  	}else Srows.push_back(diag[i]);			// infatti quei suscettibili sono già infetti
		}

		arrayS.push_back(Srows);					// inserisco le righe nella matrice
		Srows.clear();
	}


	// arrayS = {{N0-TRout0, w01,...},{w10, N1-TRout1,...},...}
	// arrayS[0] = riga0

	// vediamo se è giusto l'inserimento delle righe di S:

	/*cout << "\nmatrice S iniziale: " << endl;
	for( int i = 0; i < Npop.size(); i++ ){
		for( int j = 0; j < Npop.size(); j++){
			cout << arrayS[i][j] << " " ;
		}cout << endl;
	}cout << endl;*/


	/***************************************INIZIALIZZO LA MATRICE DEGLI INFETTI*************************************************/

	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

		for( int j = 0; j < Npop.size(); j++ ){
			if(i == j ){

				if(i==cityseed){
					Irows.push_back(infettiIniz);
				}else Irows.push_back(0);

			} else Irows.push_back(0);
		}

		arrayI.push_back(Irows);					// inserisco le righe nella matrice
		Irows.clear();
	}

	// vediamo se è giusto l'inserimento delle righe di I:

	/*cout << "\nmatrice I iniziale: " << endl;
	for( int i = 0; i < Npop.size(); i++ ){
		for( int j = 0; j < Npop.size(); j++){
			cout << arrayI[i][j] << " " ;
		}cout << endl;
	} cout << endl;	*/

	/***************************************INIZIALIZZO LA MATRICE DEI RECOVERED*************************************************/

	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

		for( int j = 0; j < Npop.size(); j++ ){
			Rrows.push_back(0);
		}

		arrayR.push_back(Rrows);					// inserisco le righe nella matrice
		Rrows.clear();
	}

	// vediamo se è giusto l'inserimento delle righe di R:

	/*cout << "\nmatrice R iniziale: " << endl;
	for( int i = 0; i < Npop.size(); i++ ){
		for( int j = 0; j < Npop.size(); j++){
			cout << arrayR[i][j] << " " ;
		}cout << endl;
	} cout << endl;	*/


	/********************************************        DINAMICA S I R         ****************************************************/

	unsigned long int seed = time(0);
	gsl_rng * r;
	const gsl_rng_type * T;

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,seed);

	//cout << "\ngenerator type: " << gsl_rng_name (r) << endl;


	for(int i = 0; i < Npop.size(); i++){
		NpopW[i] = diag[i] + trafficIN[i]; // persone presenti nel nodo i durante il worktime (COSTANTE FISSA)
	}


	for(int l = 0; l < nrun; l++){  /////////////////////////////////////////////// FOR SUI RUN
		cout << "run n.ro = " << l << endl;

		stringstream fn; //stream di stringhe
		fn.str("");
		fn << l;

		stringstream fm;
		fm.str("");
		fm << cityseed;

		string nomeOUTPUTprevalence; //crea directory outputmonile1
		nomeOUTPUTprevalence = "./OUTputMOBILE" + fm.str() + "/outputPREVALENCE" + fn.str() + ".txt";
		ofstream fileOUTPUTprevalence(nomeOUTPUTprevalence.c_str());


		string nomeOUTPUTincidence;
		nomeOUTPUTincidence = "./OUTputMOBILE" + fm.str() + "/outputINCIDENCE" + fn.str() + ".txt";
		ofstream fileOUTPUTincidence(nomeOUTPUTincidence.c_str());


		string nomeOUTPUTarrivalTime;
		nomeOUTPUTarrivalTime = "./OUTputMOBILE" + fm.str() + "/AT" + fn.str() + ".txt";
		ofstream fileOUTPUTarrivalTime(nomeOUTPUTarrivalTime.c_str());


		string nomeOUTPUTpath;
		nomeOUTPUTpath = "./OUTputMOBILE" + fm.str() + "/OUTPUTpath" + fn.str() + ".txt";
		ofstream fileOUTPUTpath(nomeOUTPUTpath.c_str());

		cout << "Output e inizializzo matrici SIR \t";
		//inizializzazione delle matrici S I R

		arraySr = arrayS;
		arrayIr = arrayI;
		arrayRr = arrayR;
		deltaS = arrayNULL;
		deltaR = arrayNULL;
		deltaI = arrayNULL;
		t = 0;
		infettiW = 0;
		infettiH = 0;


		//inizializzo la matrice dell' ARRIVAL TIME

		for(int i = 0; i < Npop.size(); i++){

			AT[i]=-1;
			cout << AT[i]<< endl;

		}


//mi fermo quando finiscono gli infetti
		while(arrayIr != arrayNULL){	//dinamica S I R
		cout << "infected";

		/******************************************************* WORK TIME ***************************************************************/


			//calcolo ARRIVAL TIME work:

			for(int i = 0; i < Npop.size(); i++){
				ATtemp[i]=0.0;
			}


			for(int j = 0; j < Npop.size(); j++){

				for(int i = 0; i < Npop.size(); i++){

					if(arrayIr[i][j]>0 && AT[i]==-1){

						ATtemp[i] = 1;
						wiw.first = j;
						wiw.second = i;
						path.push_back(wiw);
					}
				}
			}

			for(int m = 0; m < Npop.size(); m++){

				if(ATtemp[m]==1){

					AT[m] = t +0.1; //worktime
				}
			}


			//calcolo la force of infection	al lavoro

			for(int j = 0; j < Npop.size(); j++ ){

				infettiW = 0;

				for(int i = 0; i < Npop.size(); i++){
					//fisso j e sommo la colonna (I00+I10+I20..)
					infettiW = infettiW + arrayIr[i][j];
				}

				forceW[j]=0.5*infettiW/NpopW[j];		//mappatura force of infection during work time del nodo j
			}


			//aggiornamento matrici al lavoro

			for(int j = 0; j < Npop.size(); j++){

				for(int i = 0; i < Npop.size(); i++){

					deltaS[i][j] = gsl_ran_binomial(r,beta*forceW[j],arraySr[i][j]);
					deltaR[i][j] = gsl_ran_binomial(r,0.5*mu,arrayIr[i][j]);
					deltaI[i][j] = deltaS[i][j] - deltaR[i][j];

					arraySr[i][j] = arraySr[i][j] - deltaS[i][j];
					arrayIr[i][j] = arrayIr[i][j] + deltaS[i][j] - deltaR[i][j];
					arrayRr[i][j] = arrayRr[i][j] + deltaR[i][j];

				}

			}


			//variazione di infetti nel work time


			for(int i = 0; i < Npop.size(); i++){

				s = 0;
				for(int j = 0; j < Npop.size(); j++){

					s = s + deltaS[i][j]; //sommo sulle righe

				}

				Stemp.push_back(s);
			}




		/******************************************************* HOME TIME ***************************************************************/


			//calcolo ARRIVAL TIME home

			for(int i = 0; i < Npop.size(); i++){
				ATtemp[i]=0.0;
			}


			for(int i = 0; i < Npop.size(); i++){

				for(int j = 0; j < Npop.size(); j++){

					if(arrayIr[i][j]>0 && AT[j]==-1){

						ATtemp[j] = 1;
						wiw.first = i;
						wiw.second = j;
						path.push_back(wiw);

					}
				}
			}


			for(int m = 0; m < Npop.size(); m++){

				if(ATtemp[m]==1){

					AT[m] = t + 0.3; //hometime
				}
			}



			//calcolo la force of infection a casa

			for(int i = 0; i < Npop.size(); i++ ){

				infettiH = 0;

				for(int j = 0; j < Npop.size(); j++){

					//tengo fisso i e sommo sulle righe (I01+I02+I03)
					infettiH = infettiH + arrayIr[i][j];
				}

				forceH[i] = 0.5*infettiH/Npop[i];		//mappatura force of infection during home time del nodo i-esimo


			}

			//aggiorno le matrici a casa

			for(int i = 0; i < Npop.size(); i++){

				for(int j = 0; j < Npop.size(); j++){

					deltaS[i][j] = gsl_ran_binomial(r,beta*forceH[i],arraySr[i][j]);
					deltaR[i][j] = gsl_ran_binomial(r,0.5*mu,arrayIr[i][j]);
					deltaI[i][j] = deltaS[i][j] - deltaR[i][j];

					arraySr[i][j] = arraySr[i][j] - deltaS[i][j];
					arrayIr[i][j] = arrayIr[i][j] + deltaS[i][j] - deltaR[i][j];
					arrayRr[i][j] = arrayRr[i][j] + deltaR[i][j];

				}

			}


			//variazione di infetti nell'home time + variazione di infetti nel work time:

			for(int i = 0; i < Npop.size(); i++){

				s = 0;
				for(int j = 0; j < Npop.size(); j++){

					s = s + deltaS[i][j];		//sommo la riga

				}

				deltaStot.push_back(s+Stemp[i]);

			}

			Stemp.clear();


			t++;


			///////////OUTPUT fileOUTPUTprevalence
			fileOUTPUTprevalence << t << " " ;
			for(int i = 0; i < Npop.size(); i++){
					stampaInfetti = 0;
				for(int j = 0; j < Npop.size(); j++){
					stampaInfetti = stampaInfetti + arrayIr[i][j];
				}

				fileOUTPUTprevalence << stampaInfetti << " ";

			} fileOUTPUTprevalence << endl;


			///////////OUTPUT fileOUTPUTincidence means delta dei suscettibili

			fileOUTPUTincidence << t << " ";
			for(int i = 0; i < deltaStot.size(); i++){

				fileOUTPUTincidence << deltaStot[i] << " ";

			} fileOUTPUTincidence << endl;


			deltaStot.clear();



		}//chiudo while

		/*cout << "\nl'ultima matrice dei recovered R é: " << endl;
		for( int i = 0; i < Npop.size(); i++ ){
			for( int j = 0; j < Npop.size(); j++){
				cout << arrayRr[i][j] << " " ;
			}cout << endl;
		} cout << endl;*/



		for(int i = 0; i < Npop.size(); i++){
			fileOUTPUTarrivalTime << i << " " << AT[i] << endl;
		}

		for(int i = 0; i < path.size(); i++){
			fileOUTPUTpath << path[i].first << " " << path[i].second << endl;
		}

		path.clear();


		fileOUTPUTprevalence.close();
		fileOUTPUTincidence.close();
		fileOUTPUTarrivalTime.close();
		fileOUTPUTpath.close();


	}// chiudo for dei run



	gsl_rng_free (r);


return 0;
}
