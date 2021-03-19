//modello SIR con arrival time e path

#include <map>
#include <vector>
#include <string>
#include <utility>
#include <stdio.h>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <bits/stdc++.h>

using namespace std;

int i,j,infettiW,infettiH;
int stampaInfetti = 0;
double wij,TRout,TRin,s, p, c;
int cityseed = 1; 		//città seed di partenza Ikk (cityseed = 1 -> (I)11 seed)

double R0 = 2.5;
double mu = 1/3.0;	// rate di guarigione [day^-1]           R0=beta/mu>1 per avere outbreak
double beta = R0*mu;	// trasmissibilità (probabilità di trasmissione dell'infezione)
int infettiIniz = 10;
int t;//counter of the simulation
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

vector<int>::size_type k; //istanzio oggetto di tipo (vector) size

vector<double> infectedhome;
vector<double> recoveredhome;
vector<double> recoveredboro;
double infh, rech, rb; //support var for infectedhome and recoveredhome

vector<vector<double> > arrayNULL; //matrice identicamente nulla
vector<vector<double> > deltaS;
vector<vector<double> > deltaI;
vector<vector<double> > deltaR;
vector<vector<double> > arrayS;	//matrice Suscettibili
vector<vector<double> > arrayI;	//matrice I
vector<vector<double> > arrayR;	//matrice R
vector<vector<double> > arraySr;//matrice Srun
vector<vector<double> > arrayIr;//matrice Irun
vector<vector<double> > arrayRr;//matrice Rrun

vector<double> totpop;
vector<vector<double> > comm;//matrice commuters on the network
vector<double> commrows;//rows matrice commuters


map<int,vector<int> >::iterator iter1;		//iteratore per la mappa dei vicini
map<pair<int,int>,double>::iterator iter2;	//iteratore per la mappa dei pesi
map<int,vector<double> >::iterator iter3;	//iteratore per la mappa trafficOUTmap e trafficINmap
map<int,int>::iterator iter4;			//iteratore per la mappa del traffico uscente/entrante e di Npop


int main (int argc, char * argv[]){

	//APRO E LEGGO IL FILE DELLA POPOLAZIONE, I'll do this only once//
	ifstream popFILE ("./pop_boro.txt");
	if ( popFILE.is_open() ){
		cout<<"Open pop file \t"<<endl;
		while ( popFILE.good() ){

			while( popFILE >> i >> j ){//i node, 5 boroughs
				Npop[i] = j;
			}

			popFILE.close();
		}

	}else cout << "\nUnable to open the population file"<< endl;


	for( int i = 0; i < Npop.size(); i++ ){	// COSTRUISCO la matrice identicamente nulla

		for( int j = 0; j < Npop.size(); j++ ){
			NULLrows.push_back(0);
		}

		arrayNULL.push_back(NULLrows);
		NULLrows.clear();

	}


	/********************************************     DINAMICA S I R with MOBILITY   ****************************************************/

	unsigned long int seed = time(0); //funzione che restituisce il tempo macchina quando parte
	//unsigned long int seed = 123456789;
	gsl_rng * r;
	const gsl_rng_type * T; //sequenza dei numeri casuali ()

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,seed); //random generation of the seed

	cout << "\nrandomness done"<< endl;


	/*SIMULATIONS START*/
	for(int l = 0; l < nrun; l++){  /////////////////////////////////////////////// FOR SUI RUN, l will indicate the number of simulationS

		cout << "\nrun n.ro = " << l << endl;

		cout<<"popolazione network:\t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			cout << Npop[i] << " " ;
		} cout << endl;
		//APRO E LEGGO IL FILE DI COMMUTING_0, INIZIALIZZO LE MATRICI AL PRIMO GIORNO, le modificherò con un if nel while durante la singola simulazione
		ifstream networkFILE ("./commuting_0_4x4.txt");
		if ( networkFILE.is_open() ){

			trafficOUTmap.clear();
			trafficINmap.clear();
			weight.clear();
			neighbours.clear();

			while ( networkFILE.good() ){

				while( networkFILE >> i >> j >> wij ){
					neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
					weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
					if(i!=j){
					trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
					trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
					}		//attribuisce al nodo i esimo il traffico in entrata
				}
				networkFILE.close();

			}
		}else cout << "\nUnable to open the network file"<< endl;

			cout<<"FIRST WEEK TRAFFICOUTMAP è\t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << trafficOUTmap[i][j] << " " ;
				}cout << endl;
			} cout << endl;

			cout<<"FIRST WEEK TRAFFICINMAP è\t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << trafficINmap[i][j] << " " ;
				}cout << endl;
			} cout << endl;


		for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

			TRout = 0;
			for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
				TRout = TRout + iter3->second[k];
			}

			trafficOUT[iter3->first] = TRout;				// definisco quindi un array di array 1x1 di traffico in uscita da ogni nodo
		}

		cout<<"TRAFFIC OUT PER OGNI BORO:\t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			cout << trafficOUT[i] << " " ;
		} cout << endl;

		//stesso discorso per costruire array di array 1x1 della mappa trafficIN
		for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

			TRin = 0;
			for( k = 0; k<(iter3->second).size(); k++ ){
				TRin = TRin + iter3->second[k];
			}
			trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
		}

		cout<<"\nTRAFFIC IN per ogni BORO\t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			cout << trafficIN[i] << " " ;
		} cout << endl;

		arrayI.clear();
		arrayR.clear();
		arrayS.clear();

		//MATRIX OF SUSCEPTIBLE_week0//
		for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){
			diag[iter4->first] = Npop[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i
		}

		//cout<<"diag of susceptible matrix: \t"<< endl;
		//for( int j = 0; j < Npop.size(); j++){
		//	cout << diag[j] << " " ;
		//}cout << endl;

		for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

			for( int j = 0; j < Npop.size(); j++ ){

				if(i != j){
					Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
				}else if(i==cityseed){ //questo solo nella prima settimana
					Srows.push_back(diag[i]-infettiIniz);		// tolgo dalla città seed il numero di infetti iniziale
				}else Srows.push_back(diag[i]);
			}

			arrayS.push_back(Srows);					// inserisco le righe nella matrice
			Srows.clear();
		}

		cout<<"First week SUSCEPTIBLE matrix: \t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			for( int j = 0; j < Npop.size(); j++){
				cout << arrayS[i][j] << " " ;
			}cout << endl;
		} cout << endl;

		//MATRIX OF INFECTED AT WEEK0: at first I have only the infected in the cityseed at first day of simulation
		for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

			for( int j = 0; j < Npop.size(); j++ ){
				if(i == j){

					if(i==cityseed){
						Irows.push_back(infettiIniz);
					}else Irows.push_back(0);
				} else Irows.push_back(0);
			}

			arrayI.push_back(Irows);					// inserisco le righe nella matrice
			Irows.clear();
		}

		cout << "\nInfected matrix "<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			for( int j = 0; j < Npop.size(); j++){
				cout << arrayI[i][j] << " " ;
			}cout << endl;
		} cout << endl;

		//MATRIX of RECOVERED at WEEK0
		for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

			for( int j = 0; j < Npop.size(); j++ ){
				Rrows.push_back(0);
			}

			arrayR.push_back(Rrows);					// inserisco le righe nella matrice
			Rrows.clear();
		}

		cout << "\nrecovered matrix"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			for( int j = 0; j < Npop.size(); j++){
				cout << arrayR[i][j] << " " ;
			}cout << endl;
		} cout << endl;


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

		string nomeOUTPUTmatrices;
		nomeOUTPUTmatrices = "./OUTputMOBILE" + fm.str() + "/outputMATRICES" + fn.str() + ".txt";
		ofstream fileOUTPUTmatrices(nomeOUTPUTmatrices.c_str());

		cout << "Output e inizializzo matrici SIR \t" << endl;

		//ad ogni simulazione dovrò reiniziallizzare gli array alla network della prima settimana
		arraySr = arrayS;
		arrayIr = arrayI;
		arrayRr = arrayR;
		deltaS = arrayNULL;
		deltaR = arrayNULL;
		deltaI = arrayNULL;
		t = 0; //is the counter for the timesteps of a single simulation
		infettiW = 0;
		infettiH = 0;




		while(arrayIr != arrayNULL){	//dinamica S I R, arrayir verrà poi modificato dentro il while

			cout << "Time step: " << t <<endl;

			//entro in questo if solo se passo alla prima settimana
			if(t==7){
				/*CLEAR DELLE MATRICI, APRO E LEGGO IL FILE CON I DATI DI COMMUTING PER BORO ALLA SECONDA SETTIMANA*/
				cout << "Cambio settimana 1"<<endl;
				fileOUTPUTmatrices << "\n Passo alla prima settimana di variazione del commuting: " << endl;

				ifstream networkFILE ("./commuting_1_4x4.txt");
				if ( networkFILE.is_open() ){

					trafficINmap.clear();
					trafficOUTmap.clear();
					weight.clear();
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){
							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							if(i!=j){
							trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
							trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
							}		//attribuisce al nodo i esimo il traffico in entrata
						}
						networkFILE.close();
					}

				}else cout << "\nUnable to open the network file"<< endl;

				cout<<"second WEEK TRAFFICOUTMAP è\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++){
						cout << trafficOUTmap[i][j] << " " ;
					}cout << endl;
				} cout << endl;

				cout<<"second WEEK TRAFFICINMAP è\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++){
						cout << trafficINmap[i][j] << " " ;
					}cout << endl;
				} cout << endl;

				trafficIN.clear();
				trafficOUT.clear();

				for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

					TRout = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
						TRout = TRout + iter3->second[k];
					}
					trafficOUT[iter3->first] = TRout;	//ridefinisco la mappatura del traffico in base al nuovo file di commuting
				}

				cout<<"TRAFFIC OUT 2 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficOUT[i] << " " ;
				} cout << endl;

				//stesso discorso per costruire array di array 1x1 della mappa trafficIN
				for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

					TRin = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){
						TRin = TRin + iter3->second[k];
					}
					trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
				}

				cout<<"TRAFFIC In 2 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficIN[i] << " " ;
				} cout << endl;
				//initialize matrix of susceptible, infects and recovered  at second week
				/*****************************************************  MATRIX OF SUSCEPTIBLE *****************************************************/
				//compute the number of infected at each node i that come back to the home node

				infectedhome.clear();
				recoveredhome.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// QUI ATTENZIONE ALLA DIMENSIONE DELL'ARRAY

					infh = 0;
					rech = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						infectedhome.push_back(infh + arrayIr[i][j]);
						recoveredhome.push_back(rech + arrayRr[i][j]);
					}
				}

				cout<<"Infected and recovered per boro:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << infectedhome[i] << " " ;
				} cout << endl;

				for( int i = 0; i < Npop.size(); i++ ){
					cout << recoveredhome[i] << " " ;
				} cout << endl;

				arraySr.clear();
				arrayS.clear();

				for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){ //definisco la diag[i]

					diag[iter4->first] = Npop[iter4->first] - infectedhome[iter4->first] - recoveredhome[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i-INFETTI[i]-recovered[i]
				}

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){

						if(i != j){
							Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
						}else Srows.push_back(diag[i]);
					}
			  	arrayS.push_back(Srows);					// inserisco le righe nella matrice
			  	Srows.clear();
		  	}

				arraySr = arrayS; //utilizzando quindi i nuovi dati di commuting

				/************************************************** MATRIX OF INFECTED AT HOME at second week **************************************************/
				//la matrice degli infetti sarà una matrice diagonale di infetti
				arrayIr.clear();
				arrayI.clear();

		  	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
			 			if(i == j ){
							Irows.push_back(infectedhome[i]); ///salvo sulla diagonale gli infetti
						}else Irows.push_back(0);
					}

					arrayI.push_back(Irows);					// inserisco le righe nella matrice
					Irows.clear();
				}

				arrayIr = arrayI;

				/************************************************** MATRIX OF RECOVERED at home at second week ************************************************/
				arrayRr.clear();
				arrayR.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Rrows.push_back(recoveredhome[i]);
						}else Rrows.push_back(0);
					}
					arrayR.push_back(Rrows);					// inserisco le righe nella matrice
					Rrows.clear();
				}
				arrayRr = arrayR;

				//then define the people at worktime in each node: function of diag and trafficIN definite con il file di commuting

			}else if(t==1400000000000){
				/*CLEAR DELLE MATRICI, APRO E LEGGO IL FILE CON I DATI DI COMMUTING PER BORO ALLA SECONDA SETTIMANA*/
				cout << "Cambio settimana 2"<<endl;
				fileOUTPUTmatrices << "\n Passo alla seconda settimana di variazione del commuting: " << endl;
			 //in questo modo non entrerò più nell'if

				ifstream networkFILE ("./commuting_2_4x4.txt");
				if ( networkFILE.is_open() ){

					trafficINmap.clear();
					trafficOUTmap.clear();
					weight.clear();
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){
							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							if(i!=j){
							trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
							trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
							}		//attribuisce al nodo i esimo il traffico in entrata
						}
						networkFILE.close();
					}

				}else cout << "\nUnable to open the network file"<< endl;

				trafficIN.clear();
				trafficOUT.clear();

				for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

					TRout = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
						TRout = TRout + iter3->second[k];
					}
					trafficOUT[iter3->first] = TRout;	//ridefinisco la mappatura del traffico in base al nuovo file di commuting
				}

				cout<<"TRAFFIC OUT 3 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficOUT[i] << " " ;
				} cout << endl;

				//stesso discorso per costruire array di array 1x1 della mappa trafficIN
				for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

					TRin = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){
						TRin = TRin + iter3->second[k];
					}
					trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
				}

				cout<<"TRAFFIC In 3 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficIN[i] << " " ;
				} cout << endl;
				//initialize matrix of susceptible, infects and recovered  at second week
				/*****************************************************  MATRIX OF SUSCEPTIBLE *****************************************************/
				//compute the number of infected at each node i that come back to the home node

				infectedhome.clear();
				recoveredhome.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// QUI ATTENZIONE ALLA DIMENSIONE DELL'ARRAY

					infh = 0;
					rech = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						infectedhome.push_back( infh + arrayIr[i][j]);
						recoveredhome.push_back ( rech + arrayRr[i][j]);
					}
				}

				arraySr.clear();
				arrayS.clear();

				for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){ //definisco la diag[i]

					diag[iter4->first] = Npop[iter4->first] - infectedhome[iter4->first] - recoveredhome[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i-INFETTI[i]-recovered[i]
				}

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){

						if(i != j){
							Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
						}else Srows.push_back(diag[i]);
					}
			  	arrayS.push_back(Srows);					// inserisco le righe nella matrice
			  	Srows.clear();
		  	}

				arraySr = arrayS; //utilizzando quindi i nuovi dati di commuting

				/************************************************** MATRIX OF INFECTED AT HOME at second week **************************************************/
				//la matrice degli infetti sarà una matrice diagonale di infetti
				arrayIr.clear();
				arrayI.clear();

		  	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
			 			if(i == j ){
							Irows.push_back(infectedhome[i]); ///salvo sulla diagonale gli infetti
						}else Irows.push_back(0);
					}

					arrayI.push_back(Irows);					// inserisco le righe nella matrice
					Irows.clear();
				}

				arrayIr = arrayI;

				/************************************************** MATRIX OF RECOVERED at home at second week ************************************************/
				arrayRr.clear();
				arrayR.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Rrows.push_back(recoveredhome[i]);
						}else Rrows.push_back(0);
					}
					arrayR.push_back(Rrows);					// inserisco le righe nella matrice
					Rrows.clear();
				}
				arrayRr = arrayR;

				//then define the people at worktime in each node: function of diag and trafficIN definite con il file di commuting

			}else if(t==2100000000000){
				cout << "Cambio settimana 3"<<endl;
				fileOUTPUTmatrices << "\n Passo alla terza settimana di variazione del commuting: " << endl;

				ifstream networkFILE ("./commuting_3_4x4.txt");
				if ( networkFILE.is_open() ){

					trafficINmap.clear();
					trafficOUTmap.clear();
					weight.clear();
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){
							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							if(i!=j){
							trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
							trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
							}		//attribuisce al nodo i esimo il traffico in entrata
						}
						networkFILE.close();
					}

				}else cout << "\nUnable to open the network file"<< endl;

				trafficIN.clear();
				trafficOUT.clear();

				for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

					TRout = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
						TRout = TRout + iter3->second[k];
					}
					trafficOUT[iter3->first] = TRout;	//ridefinisco la mappatura del traffico in base al nuovo file di commuting
				}

				cout<<"TRAFFIC OUT 3 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficOUT[i] << " " ;
				} cout << endl;

				//stesso discorso per costruire array di array 1x1 della mappa trafficIN
				for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

					TRin = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){
						TRin = TRin + iter3->second[k];
					}
					trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
				}

				cout<<"TRAFFIC In 3 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficIN[i] << " " ;
				} cout << endl;
				//initialize matrix of susceptible, infects and recovered  at second week
				/*****************************************************  MATRIX OF SUSCEPTIBLE *****************************************************/
				//compute the number of infected at each node i that come back to the home node
				infectedhome.clear();
				recoveredhome.clear();

				for( int i = 0; i < Npop.size(); i++ ){

					infh = 0;
					rech = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						infectedhome.push_back( infh + arrayIr[i][j]);
						recoveredhome.push_back ( rech + arrayRr[i][j]);
					}
				}

				arraySr.clear();
				arrayS.clear();

				for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){ //definisco la diag[i]

					diag[iter4->first] = Npop[iter4->first] - infectedhome[iter4->first] - recoveredhome[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i-INFETTI[i]-recovered[i]
				}

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){

						if(i != j){
							Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
						}else Srows.push_back(diag[i]);
					}
			  	arrayS.push_back(Srows);					// inserisco le righe nella matrice
			  	Srows.clear();
		  	}

				arraySr = arrayS; //utilizzando quindi i nuovi dati di commuting

				/************************************************** MATRIX OF INFECTED  **************************************************/
				//la matrice degli infetti sarà una matrice diagonale di infetti
				arrayIr.clear();
				arrayI.clear();

		  	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
			 			if(i == j ){
							Irows.push_back(infectedhome[i]); ///salvo sulla diagonale gli infetti
						}else Irows.push_back(0);
					}

					arrayI.push_back(Irows);					// inserisco le righe nella matrice
					Irows.clear();
				}

				arrayIr = arrayI;

				/************************************************** MATRIX OF RECOVERED  ************************************************/
				arrayRr.clear();
				arrayR.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Rrows.push_back(recoveredhome[i]);
						}else Rrows.push_back(0);
					}
					arrayR.push_back(Rrows);					// inserisco le righe nella matrice
					Rrows.clear();
				}
				arrayRr = arrayR;

			}else if(t==2800000000000){
				cout << "Cambio settimana 5"<<endl;
				fileOUTPUTmatrices << "\n Passo alla quarta settimana di variazione del commuting: " << endl;

				ifstream networkFILE ("./commuting_4_4x4.txt");
				if ( networkFILE.is_open() ){

					trafficINmap.clear();
					trafficOUTmap.clear();
					weight.clear();
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){
							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							if(i!=j){
							trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
							trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
							}		//attribuisce al nodo i esimo il traffico in entrata
						}
						networkFILE.close();
					}

				}else cout << "\nUnable to open the network file"<< endl;

				trafficIN.clear();
				trafficOUT.clear();

				for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

					TRout = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
						TRout = TRout + iter3->second[k];
					}
					trafficOUT[iter3->first] = TRout;	//ridefinisco la mappatura del traffico in base al nuovo file di commuting
				}

				cout<<"TRAFFIC OUT 4 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficOUT[i] << " " ;
				} cout << endl;

				//stesso discorso per costruire array di array 1x1 della mappa trafficIN
				for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

					TRin = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){
						TRin = TRin + iter3->second[k];
					}
					trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
				}

				cout<<"TRAFFIC In 4 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficIN[i] << " " ;
				} cout << endl;
				//initialize matrix of susceptible, infects and recovered  at second week
				/*****************************************************  MATRIX OF SUSCEPTIBLE *****************************************************/
				//compute the number of infected at each node i that come back to the home node

				infectedhome.clear();
				recoveredhome.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// QUI ATTENZIONE ALLA DIMENSIONE DELL'ARRAY

					infh = 0;
					rech = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						infectedhome.push_back( infh + arrayIr[i][j]);
						recoveredhome.push_back ( rech + arrayRr[i][j]);
					}
				}

				arraySr.clear();
				arrayS.clear();

				for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){ //definisco la diag[i]

					diag[iter4->first] = Npop[iter4->first] - infectedhome[iter4->first] - recoveredhome[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i-INFETTI[i]-recovered[i]
				}

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){

						if(i != j){
							Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
						}else Srows.push_back(diag[i]);
					}
			  	arrayS.push_back(Srows);					// inserisco le righe nella matrice
			  	Srows.clear();
		  	}

				arraySr = arrayS; //utilizzando quindi i nuovi dati di commuting

				/************************************************** MATRIX OF INFECTED AT HOME at second week **************************************************/
				//la matrice degli infetti sarà una matrice diagonale di infetti
				arrayIr.clear();
				arrayI.clear();

		  	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
			 			if(i == j ){
							Irows.push_back(infectedhome[i]); ///salvo sulla diagonale gli infetti
						}else Irows.push_back(0);
					}

					arrayI.push_back(Irows);					// inserisco le righe nella matrice
					Irows.clear();
				}

				arrayIr = arrayI;

				/************************************************** MATRIX OF RECOVERED at home at second week ************************************************/
				arrayRr.clear();
				arrayR.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Rrows.push_back(recoveredhome[i]);
						}else Rrows.push_back(0);
					}
					arrayR.push_back(Rrows);					// inserisco le righe nella matrice
					Rrows.clear();
				}
				arrayRr = arrayR;
			}else if(t==3500000000000){
				cout << "Cambio settimana 5"<<endl;
				fileOUTPUTmatrices << "\n Passo alla quinta settimana di variazione del commuting: " << endl;

				ifstream networkFILE ("./commuting_5_4x4.txt");
				if ( networkFILE.is_open() ){

					trafficINmap.clear();
					trafficOUTmap.clear();
					weight.clear();
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){
							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							if(i!=j){
							trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
							trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
							}		//attribuisce al nodo i esimo il traffico in entrata
						}
						networkFILE.close();
					}

				}else cout << "\nUnable to open the network file"<< endl;

				trafficIN.clear();
				trafficOUT.clear();

				for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

					TRout = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
						TRout = TRout + iter3->second[k];
					}
					trafficOUT[iter3->first] = TRout;	//ridefinisco la mappatura del traffico in base al nuovo file di commuting
				}

				cout<<"TRAFFIC OUT 3 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficOUT[i] << " " ;
				} cout << endl;

				//stesso discorso per costruire array di array 1x1 della mappa trafficIN
				for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

					TRin = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){
						TRin = TRin + iter3->second[k];
					}
					trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
				}

				cout<<"TRAFFIC In 3 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficIN[i] << " " ;
				} cout << endl;
				//initialize matrix of susceptible, infects and recovered  at second week
				/*****************************************************  MATRIX OF SUSCEPTIBLE *****************************************************/
				//compute the number of infected at each node i that come back to the home node

				infectedhome.clear();
				recoveredhome.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// QUI ATTENZIONE ALLA DIMENSIONE DELL'ARRAY

					infh = 0;
					rech = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						infectedhome.push_back( infh + arrayIr[i][j]);
						recoveredhome.push_back ( rech + arrayRr[i][j]);
					}
				}

				arraySr.clear();
				arrayS.clear();

				for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){ //definisco la diag[i]

					diag[iter4->first] = Npop[iter4->first] - infectedhome[iter4->first] - recoveredhome[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i-INFETTI[i]-recovered[i]
				}

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){

						if(i != j){
							Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
						}else Srows.push_back(diag[i]);
					}
			  	arrayS.push_back(Srows);					// inserisco le righe nella matrice
			  	Srows.clear();
		  	}

				arraySr = arrayS; //utilizzando quindi i nuovi dati di commuting

				/************************************************** MATRIX OF INFECTED AT HOME at second week **************************************************/
				//la matrice degli infetti sarà una matrice diagonale di infetti
				arrayIr.clear();
				arrayI.clear();

		  	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
			 			if(i == j ){
							Irows.push_back(infectedhome[i]); ///salvo sulla diagonale gli infetti
						}else Irows.push_back(0);
					}

					arrayI.push_back(Irows);					// inserisco le righe nella matrice
					Irows.clear();
				}

				arrayIr = arrayI;

				/************************************************** MATRIX OF RECOVERED at home at second week ************************************************/
				arrayRr.clear();
				arrayR.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Rrows.push_back(recoveredhome[i]);
						}else Rrows.push_back(0);
					}
					arrayR.push_back(Rrows);					// inserisco le righe nella matrice
					Rrows.clear();
				}
				arrayRr = arrayR;
			}else if(t==4200000000000){
				cout << "Cambio settimana 6"<<endl;
				fileOUTPUTmatrices << "\n Passo alla sesta settimana di variazione del commuting: " << endl;

				ifstream networkFILE ("./commuting_6_4x4.txt");
				if ( networkFILE.is_open() ){

					trafficINmap.clear();
					trafficOUTmap.clear();
					weight.clear();
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){
							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							if(i!=j){
							trafficOUTmap[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)
							trafficINmap[j].push_back(wij);		//attribuisce al nodo i esimo il traffico in entrata
							}		//attribuisce al nodo i esimo il traffico in entrata
						}
						networkFILE.close();
					}

				}else cout << "\nUnable to open the network file"<< endl;

				trafficIN.clear();
				trafficOUT.clear();

				for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){ //ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo

					TRout = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){ //infatti ogni nodo avrà 5 nodi in uscita (size di iter3)
						TRout = TRout + iter3->second[k];
					}
					trafficOUT[iter3->first] = TRout;	//ridefinisco la mappatura del traffico in base al nuovo file di commuting
				}

				cout<<"TRAFFIC OUT 6 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficOUT[i] << " " ;
				} cout << endl;

				//stesso discorso per costruire array di array 1x1 della mappa trafficIN
				for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

					TRin = 0;
					for( k = 0; k<(iter3->second).size(); k++ ){
						TRin = TRin + iter3->second[k];
					}
					trafficIN[iter3->first] = TRin;					// ho la mappa del traffico in entrata di ogni nodo
				}

				cout<<"TRAFFIC In 6 settimana PER OGNI BORO:\t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					cout << trafficIN[i] << " " ;
				} cout << endl;
				//initialize matrix of susceptible, infects and recovered  at second week
				/*****************************************************  MATRIX OF SUSCEPTIBLE *****************************************************/
				//compute the number of infected at each node i that come back to the home node

				infectedhome.clear();
				recoveredhome.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// QUI ATTENZIONE ALLA DIMENSIONE DELL'ARRAY

					infh = 0;
					rech = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						infectedhome.push_back( infh + arrayIr[i][j]);
						recoveredhome.push_back ( rech + arrayRr[i][j]);
					}
				}

				arraySr.clear();
				arrayS.clear();

				for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){ //definisco la diag[i]

					diag[iter4->first] = Npop[iter4->first] - infectedhome[iter4->first] - recoveredhome[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i-INFETTI[i]-recovered[i]
				}

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){

						if(i != j){
							Srows.push_back(weight[pair<int,int>(i,j)]);//riempio i link del network di suscettibili
						}else Srows.push_back(diag[i]);
					}
			  	arrayS.push_back(Srows);					// inserisco le righe nella matrice
			  	Srows.clear();
		  	}

				arraySr = arrayS; //utilizzando quindi i nuovi dati di commuting

				/************************************************** MATRIX OF INFECTED AT HOME at second week **************************************************/
				//la matrice degli infetti sarà una matrice diagonale di infetti
				arrayIr.clear();
				arrayI.clear();

		  	for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
			 			if(i == j ){
							Irows.push_back(infectedhome[i]); ///salvo sulla diagonale gli infetti
						}else Irows.push_back(0);
					}

					arrayI.push_back(Irows);					// inserisco le righe nella matrice
					Irows.clear();
				}

				arrayIr = arrayI;

				/************************************************** MATRIX OF RECOVERED at home at second week ************************************************/
				arrayRr.clear();
				arrayR.clear();

				for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe

					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Rrows.push_back(recoveredhome[i]);
						}else Rrows.push_back(0);
					}
					arrayR.push_back(Rrows);					// inserisco le righe nella matrice
					Rrows.clear();
				}
				arrayRr = arrayR;
			}

			/******************************************************* WORK TIME ***************************************************************/
			fileOUTPUTmatrices << "\nTimestep of the simulation:" << t << endl;

			fileOUTPUTmatrices<<"Starting ArrayIr: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayIr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices<< endl;

			fileOUTPUTmatrices<<"Starting ArraySr:\t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arraySr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"Starting ArrayRr: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayRr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;


			//then define the people at worktime in each node: function of diag and trafficIN che sono stati definiti con la lettura del file di commuting_0
			for(int i = 0; i < Npop.size(); i++){
				NpopW[i] = diag[i] + trafficIN[i];
			}// persone presenti nel nodo i durante il worktime


			//calcolo la force of infection	durante il worktime
			for(int j = 0; j < Npop.size(); j++ ){

				infettiW = 0;

				for(int i = 0; i < Npop.size(); i++){
					//fisso j e sommo la colonna (I00+I10+I20..)
					infettiW = infettiW + arrayIr[i][j]; //arrayIr will be update correctly from deltaS / deltaR
				}

				forceW[j]=0.5*infettiW/NpopW[j];		//mappatura force of infection during work time del nodo j
			}



			//aggiornamento matrici durante il worktime
			for(int j = 0; j < Npop.size(); j++){

				for(int i = 0; i < Npop.size(); i++){

					deltaS[i][j] = gsl_ran_binomial(r,beta*forceW[j],arraySr[i][j]);//is ok to use arraySr of the last time step
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

			fileOUTPUTmatrices<<"Worktime matrices: \t"<< endl;

			fileOUTPUTmatrices << "Npopw worktime:" << endl;
			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices<< NpopW[i] << " " ;
			}fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"Susceptible worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arraySr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"DeltaS worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << deltaS[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"Infected worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayIr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"Recovered worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayRr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

		/******************************************************* HOME TIME ***************************************************************/
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


			cout<<"ArrayIr è\t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << arrayIr[i][j] << " " ;
				}cout << endl;
			} cout << endl;

			cout<<"ArrayRr è\t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << arrayRr[i][j] << " " ;
				}cout << endl;
			} cout << endl;


			//check matrice commuters su ij costante
			for( int i = 0; i < Npop.size(); i++ ){					// costruisco le righe
				for( int j = 0; j < Npop.size(); j++ ){
					c=0;
					c = arraySr[i][j] + arrayRr[i][j] + arrayIr[i][j];
					commrows.push_back(c);
				}
				comm.push_back(commrows);					// inserisco le righe nella matrice
				commrows.clear();
			}

			cout<<"commuters on ij node: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << comm[i][j] << " " ;
				}cout << endl;
			} cout << endl;

			comm.clear();


			for(int i = 0; i < Npop.size(); i++){

				p = 0;
				for(int j = 0; j < Npop.size(); j++){

					p = p + arraySr[i][j] + arrayRr[i][j] + arrayIr[i][j]; //sommo sulle righe

				}

				totpop.push_back(p);
			}

			cout << " tot pop per boro:" << endl;
			for(int i = 0; i < Npop.size(); i++){
				cout << totpop[i] << " " ;
			}cout << endl;

			t++;//end of a single step in the simulation

			///////////OUTPUT fileOUTPUTprevalence
			fileOUTPUTprevalence << t << " \t" ;
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

			/////////OUTPUT MATRICES SIR///////////////
			fileOUTPUTmatrices<<"Hometime matrices: \t"<< endl;
			fileOUTPUTmatrices << "\n Susceptible hometime: \t" << endl;
			for( int i = 0; i < Npop.size(); i++ ){

				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arraySr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;

			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices << "\n DeltaS hometime: \t" << endl;
			for( int i = 0; i < Npop.size(); i++ ){

				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << deltaS[i][j] << " " ;
				}fileOUTPUTmatrices << endl;

			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices << "\n Infected hometime: \t" << endl;
			for( int i = 0; i < Npop.size(); i++ ){

				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayIr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;

			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices << "\n Recovered hometime: \t" << endl;
			for( int i = 0; i < Npop.size(); i++ ){

				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayRr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;

			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices << "\n Recovered per boro" <<endl;
			for( int i = 0; i < Npop.size(); i++ ){

				rb = 0;
				for( int j = 0; j < Npop.size(); j++ ){
					rb = rb + arrayRr[i][j];
				} recoveredboro.push_back(rb);
			}

			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices << recoveredboro[i] << " " ;
			}fileOUTPUTmatrices << endl;



			fileOUTPUTmatrices << "\nTot pop per node of the network:" << endl;
			for(int i = 0; i < Npop.size(); i++){

				fileOUTPUTmatrices << totpop[i] << " " ;
			}fileOUTPUTmatrices << endl;


			///end of writing of output files
			recoveredboro.clear();
			totpop.clear();
			deltaStot.clear();

		}//chiudo while, i.e. il numero degli infetti è tornato a zero

		fileOUTPUTprevalence.close();
		fileOUTPUTincidence.close();
		fileOUTPUTmatrices.close();


	}// chiudo for dei run (di tutte le simulazioni)


	gsl_rng_free (r);

return 0;

}
