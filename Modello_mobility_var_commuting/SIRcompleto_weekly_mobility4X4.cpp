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

double R0 = 3;
double mu = 1/3.0;	// rate di guarigione [day^-1]           R0=beta/mu>1 per avere outbreak
double beta = R0*mu;	// trasmissibilità (probabilità di trasmissione dell'infezione)
int infettiIniz = 10;
int t;//counter of the simulation
int nrun = 1;
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

//file dove salvare temporaneamente le variazioni relative
map<int,vector<double> > rel_var_matrix;	//mappa delle variazioni sui link del numero di commuters
map<pair<int,int>,double> rel_var_weight;


vector<double> NULLrows;	//righe di arrayNULL
vector<double> Srows;		//righe di arrayS
vector<double> Irows;		//righe di arrayI
vector<double> Rrows;		//righe di arrayR
vector<double> deltaStot;	//variazione infetti alla fine della giornata
vector<double> Stemp;		//variabile temporanea
vector<double> rel_var_rows;

vector<int>::size_type k; //istanzio oggetto di tipo (vector) size

vector<double> susceptiblehome; //array di suscettibili da riportare a casa
vector<double> infectedhome;
vector<double> recoveredhome;
vector<double> infectedboro;
vector<double> recoveredboro;
double infh, rech, susch, infb, rb; //support var for infectedhome and recoveredhome

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

//creo supporti temporanei per il riaggiornamento delle matrici temporanee di variazione del commuting
vector<double> tempSr_rows;
vector<double> tempIr_rows;
vector<double> tempRr_rows;
vector<vector<double> > Stemp_matrix;//matrice variazioni commuters suscettibili
vector<vector<double> > Itemp_matrix;//matrice variazioni commuters infetti
vector<vector<double> > Rtemp_matrix;//matrice variazioni commuters recovered

double news, sr_perboro, stemp_perboro, infr_perboro, itemp_perboro, recr_perboro, rtemp_perboro;
int newsint;//per forzare ad intero l'operazione di aggiornamento

vector<double> totpop;
vector<vector<double> > comm;//matrice commuters on the network
vector<double> commrows;//rows matrice commuters


map<int,vector<int> >::iterator iter1;		//iteratore per la mappa dei vicini
map<pair<int,int>,double>::iterator iter2;	//iteratore per la mappa dei pesi
map<int,vector<double> >::iterator iter3;	//iteratore per la mappa trafficOUTmap e trafficINmap
map<int,int>::iterator iter4;			//iteratore per la mappa del traffico uscente/entrante e di Npop


int main (int argc, char * argv[]){
	//APRO E LEGGO IL FILE DELLA POPOLAZIONE
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

		cout << "\n numero dela simulazione: " << l << endl;

		cout<<"Popolazione network:\t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			cout << Npop[i] << " " ;
		} cout << endl;

		ifstream networkFILE ("./commuting_0_4x4.txt");
		if ( networkFILE.is_open() ){ //open network file of initial commuting

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

		//ad ogni iter3 selezionerò il nodo di partenza e l'array che indica i nodi di arrivo
		for( iter3 = trafficOUTmap.begin(); iter3 != trafficOUTmap.end(); iter3 ++ ){

			TRout = 0;
			for( k = 0; k<(iter3->second).size(); k++ ){
				TRout = TRout + iter3->second[k];
			}

			trafficOUT[iter3->first] = TRout;
		}

		for( iter3 = trafficINmap.begin(); iter3 != trafficINmap.end(); iter3 ++ ){

			TRin = 0;
			for( k = 0; k<(iter3->second).size(); k++ ){
				TRin = TRin + iter3->second[k];
			}
			trafficIN[iter3->first] = TRin;
		}


		cout<<"\nTraffic out per boro:\t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			cout << trafficOUT[i] << " " ;
		} cout << endl;

		cout<<"\nTraffic in per boro\t"<< endl;
		for( int i = 0; i < Npop.size(); i++ ){
			cout << trafficIN[i] << " " ;
		} cout << endl;


		// MATRIX OF SUSCEPTIBLE AT FIRST WEEK OF SIMULATION//

		//define first the diagonal of the S matrix//
		for( iter4 = trafficOUT.begin(); iter4 != trafficOUT.end(); iter4 ++ ){
			diag[iter4->first] = Npop[iter4->first] - iter4->second;	// ho la mappa della diagonale di S: diag[i]=(N)i-(TRout)i
		}

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

		cout<<"\n Initial SUSCEPTIBLE matrix: \t"<< endl;
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

		cout << "\n Initial INFECTED matrix "<< endl;
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

		cout << "\n Initial RECOVERED matrix"<< endl;
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

		arraySr = arrayS;
		arrayIr = arrayI;
		arrayRr = arrayR;
		deltaS = arrayNULL;
		deltaR = arrayNULL;
		deltaI = arrayNULL;
		t = 0; //is the counter for the timesteps of a single simulation
		infettiW = 0;
		infettiH = 0;


		//*********************************************** DINAMICA SIR ********************************//
		while(arrayIr != arrayNULL){

			fileOUTPUTmatrices << "\nTimestep of the simulation:" << t << endl;
			cout << "\nTimestep of the simulation:" << t << endl;

			//if per aggiornare le matrici SIR con la variazione relativa di commuting
			if(t==7){
				//apro il file di variazione relativa del commuting alla prima settimana
				cout << "\n Commuting variation at first week of simulation: \t" << endl;
				fileOUTPUTmatrices << "\n Aggiorno le matrici SIR alla prima settimana: \t" <<endl;

				ifstream networkFILE ("./relative_var_1_4x4.txt");
				if ( networkFILE.is_open() ){ //open network file of initial commuting

					rel_var_matrix.clear(); //map
					neighbours.clear();

					while ( networkFILE.good() ){

						while( networkFILE >> i >> j >> wij ){

							neighbours[i].push_back(j);		//aggiunge ai vicini del nodo i il nodo j
							weight[pair<int,int>(i,j)]= wij;	//attribuisce al link tra il nodo i e il nodo j il peso wij
							rel_var_matrix[i].push_back(wij);	//attribuisce al nodo i esimo il traffico in uscita (tiene conto anche di i --> i)

						}
						networkFILE.close();

					}
				}else cout << "\nUnable to open the network file"<< endl;

				cout << "\n Variazione relativa commuting prima settimana: \t" <<endl;

				fileOUTPUTmatrices<<"\n Variazione relativa commuting prima settimana: \t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++){
						fileOUTPUTmatrices << rel_var_matrix[i][j] << " " ;
					}fileOUTPUTmatrices << endl;
				}fileOUTPUTmatrices<< endl;

				/***********************************MATRICE S, I, R  VAR PRIMA SETTIMANA***********************************/

				susceptiblehome.clear();
				Stemp_matrix.clear();

				for( int i = 0; i < Npop.size(); i++ ){
					sr_perboro = 0;
					stemp_perboro = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						news = 0;
						newsint = 0;
						if(i != j){ //off diagonal terms con variazione dei commuters suscettibili
							news = arraySr[i][j] * rel_var_matrix[i][j];
							newsint = news;
							tempSr_rows.push_back(newsint);//riempio i link del network di suscettibili

							sr_perboro = sr_perboro + arraySr[i][j]; //off diagonal terms in susceptible matrix
							stemp_perboro = stemp_perboro + newsint;

						}else tempSr_rows.push_back(0);
					}
					susceptiblehome.push_back(sr_perboro - stemp_perboro); //potrebbe venire negativo quando il numero di commutes aumenta
					Stemp_matrix.push_back(tempSr_rows);// inserisco le righe nella matrice
					tempSr_rows.clear();
				}

				fileOUTPUTmatrices << "Suscettibili at home per boro: \t" << endl;
				for ( int i = 0; i < Npop.size(); i++){
					fileOUTPUTmatrices << susceptiblehome[i] << endl;
				} fileOUTPUTmatrices << endl;

				infectedhome.clear();
				Itemp_matrix.clear();

				for( int i = 0; i < Npop.size(); i++ ){
					infr_perboro = 0;
					itemp_perboro = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						news = 0;
						newsint = 0;
						if(i != j){ //off diagonal terms con variazione dei commuters infetti
							news = arrayIr[i][j] * rel_var_matrix[i][j];
							newsint = news;
							tempIr_rows.push_back(newsint);//riempio i link del network di suscettibili

							infr_perboro = infr_perboro + arrayIr[i][j]; //off diagonal terms in susceptible matrix
							itemp_perboro = itemp_perboro + newsint;

						}else tempIr_rows.push_back(0);
					}
					infectedhome.push_back(infr_perboro - itemp_perboro); //potrebbe venire negativo quando il numero di commutes aumenta
					Itemp_matrix.push_back(tempIr_rows);// inserisco le righe nella matrice
					tempIr_rows.clear();
				}

				fileOUTPUTmatrices << "Infetti at home per boro: \t" << endl;
				for ( int i = 0; i < Npop.size(); i++){
					fileOUTPUTmatrices << infectedhome[i] << endl;
				} fileOUTPUTmatrices << endl;


				recoveredhome.clear();
				Rtemp_matrix.clear();

				for( int i = 0; i < Npop.size(); i++ ){
					recr_perboro = 0;
					rtemp_perboro = 0;

					for( int j = 0; j < Npop.size(); j++ ){

						news = 0;
						newsint = 0;
						if(i != j){ //off diagonal terms con variazione dei commuters infetti
							news = arrayRr[i][j] * rel_var_matrix[i][j];
							newsint = news;
							tempRr_rows.push_back(newsint);//riempio i link del network di suscettibili

							recr_perboro = recr_perboro + arrayRr[i][j]; //off diagonal terms in susceptible matrix
							rtemp_perboro = rtemp_perboro + newsint;

						}else tempRr_rows.push_back(0);
					}
					recoveredhome.push_back(recr_perboro - rtemp_perboro); //potrebbe venire negativo quando il numero di commutes aumenta
					Rtemp_matrix.push_back(tempRr_rows);// inserisco le righe nella matrice
					tempRr_rows.clear();
				}

				fileOUTPUTmatrices << "Recovered at home per boro: \t" << endl;
				for ( int i = 0; i < Npop.size(); i++){
					fileOUTPUTmatrices << recoveredhome[i] << endl;
				} fileOUTPUTmatrices << endl;

				// AGGIORNO LE MATRICI RUN RIPORTANDO SULLA DIAGONALE I NUOVI SUSCETTIBILI, INFETTI E RECOVERED CHE NON VIAGGIANO

				//DIAGONALI delle matrici Stemp Itemp Rtemp

				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++ ){
						if(i==j){
							Stemp_matrix[i][j] = arraySr[i][j] + susceptiblehome[i];
							Itemp_matrix[i][j] = arrayIr[i][j] + infectedhome[i];
							Rtemp_matrix[i][j] = arrayRr[i][j] + recoveredhome[i];
							trafficIN[i] = susceptiblehome[i] + infectedhome[i] + recoveredhome[i]; //traffico in ingresso
						}
					}
				}

				arraySr = Stemp_matrix;
				arrayIr = Itemp_matrix;
				arrayRr = Rtemp_matrix;

				fileOUTPUTmatrices<<"Matrice suscettibili aggiornata first week: \t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++){
						fileOUTPUTmatrices << arraySr[i][j] << " " ;
					}fileOUTPUTmatrices << endl;
				} fileOUTPUTmatrices << endl;

				fileOUTPUTmatrices<<"Matrice infetti aggiornata first week: \t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++){
						fileOUTPUTmatrices << arrayIr[i][j] << " " ;
					}fileOUTPUTmatrices << endl;
				} fileOUTPUTmatrices << endl;

				fileOUTPUTmatrices<<"Matrice recovered aggiornata first week: \t"<< endl;
				for( int i = 0; i < Npop.size(); i++ ){
					for( int j = 0; j < Npop.size(); j++){
						fileOUTPUTmatrices << arrayRr[i][j] << " " ;
					}fileOUTPUTmatrices << endl;
				} fileOUTPUTmatrices << endl;


			}//end dell'IF di aggiornamento della prima settimana


			/******************************************************* WORK TIME ***************************************************************/

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


			//compute total pop per boro
			for(int i = 0; i < Npop.size(); i++){
				p = 0;
				for(int j = 0; j < Npop.size(); j++){
					p = p + arraySr[i][j] + arrayRr[i][j] + arrayIr[i][j];
				}
				totpop.push_back(p);
			}

			fileOUTPUTmatrices << "\n Tot population at starting of the timestep: \t" << endl;
			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices << totpop[i] << " " ;
			}fileOUTPUTmatrices << endl;
			totpop.clear();

			//diagonale della matrice dei suscettibili
			for(int i = 0; i < Npop.size(); i++){
				for(int j = 0; j < Npop.size(); j++){
					if(i==j){
						diag[i] = arraySr[i][j];
					}
				}
			}

			fileOUTPUTmatrices << "\n Diagonale matrice dei suscettibili: \t" << endl;
			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices << diag[i] << " " ;
			}fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices << "\n TrafficIN per boro: \t" << endl;
			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices << trafficIN[i] << " " ;
			}fileOUTPUTmatrices << endl;

			//ho bisogno di calcolare anche la diagonale della matrice dei suscettibili
			for(int i = 0; i < Npop.size(); i++){
				NpopW[i] = diag[i] + trafficIN[i];
			}// persone presenti nel nodo i durante il worktime


			//calcolo la force of infection	durante il worktime
			for(int j = 0; j < Npop.size(); j++ ){

				infettiW = 0;

				for(int i = 0; i < Npop.size(); i++){

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

			//Variazione numero di infetti durante il worktime (variazione dei suscettibili)
			for(int i = 0; i < Npop.size(); i++){
				s = 0;
				for(int j = 0; j < Npop.size(); j++){

					s = s + deltaS[i][j]; //sommo sulle righe

				}
				Stemp.push_back(s);
			}

			fileOUTPUTmatrices<<"\n Worktime matrices: \t"<< endl;

			fileOUTPUTmatrices << "\n Npopw worktime:" << endl;
			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices<< NpopW[i] << " " ;
			}fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"\n Susceptible worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arraySr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"\n DeltaS worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << deltaS[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"\n Infected worktime: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					fileOUTPUTmatrices << arrayIr[i][j] << " " ;
				}fileOUTPUTmatrices << endl;
			} fileOUTPUTmatrices << endl;

			fileOUTPUTmatrices<<"\n Recovered worktime: \t"<< endl;
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

			cout<<"\n INFECTED matrix: \t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << arrayIr[i][j] << " " ;
				}cout << endl;
			} cout << endl;

			cout<<"\n RECOVERED matrix:\t"<< endl;
			for( int i = 0; i < Npop.size(); i++ ){
				for( int j = 0; j < Npop.size(); j++){
					cout << arrayRr[i][j] << " " ;
				}cout << endl;
			} cout << endl;

			//Verify the total pop per boro
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

			infectedboro.clear();
			recoveredboro.clear();

			fileOUTPUTmatrices << "\n Infected per boro \t" << endl;
			for( int i = 0; i < Npop.size(); i++){

				infb = 0;
				for (int j = 0; j < Npop.size(); j++){
					infb = infb + arrayIr[i][j];
				}infectedboro.push_back(infb);
			}

			for(int i = 0; i < Npop.size(); i++){
				fileOUTPUTmatrices << infectedboro[i] << " " ;
			}fileOUTPUTmatrices << endl;


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
