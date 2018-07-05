//#include"stdafx.h"
#include<boost/lexical_cast.hpp>
#include <vector>
#include<stdio.h>
#include<iostream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include"search_for.h"
using namespace std;

	vector<double> rx = { 0.125 ,0.125 ,-0.125 ,-0.125 };
	vector<double> ry = { 0.125 ,-0.125 ,0.125 ,-0.125 };
	vector<double> rz = { 0.125 ,-0.125 ,-0.125 ,0.125 };
double ourrandom()
{
	return rand() / (RAND_MAX + 1.);
}

int randint(int nsites)
{
	return rand() % nsites;
}

vector<int> create_fm_config(int L)
{
	vector<int> config(4*L*L*L);
	for (int i = 0; i<L*L*L; i++) config[i] = -1;
	for (int i = L*L*L; i<2*L*L*L; i++) config[i] =1;
	for (int i = 2*L*L*L; i<3*L*L*L; i++) config[i] = -1;
	for (int i = 3*L*L*L; i<4*L*L*L; i++) config[i] = 1;
	return config;
}

void make_pyrochlore4(int &L, vector< vector<int> > &neighbors)
{
	vector<int>    nnentry;
	
	//////////////////////////////////////////////
	//             Pyrochlore lattice
	//////////////////////////////////////////////
	// Make sites on cube of dimensions L,L,L
	// Associate 16 sites with every n1,n2,n3
	//
	// Tetrahedron coords (used in Lucile/Ross paper)
	// r0 = 1/8 (+1,+1,+1)
	// r1 = 1/8 (+1,-1,-1)
	// r2 = 1/8 (-1,+1,-1)
	// r3 = 1/8 (-1,-1,+1)
	for (int i = 0 ; i <  L*L*L; i++) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		nnentry.push_back(i + L*L*L); nnentry.push_back(i +3* L*L*L); nnentry.push_back(i + 2 * L*L*L);
		nnentry.push_back((i % (L*L) - L + L*L) % (L*L)  + (i / (L*L))*L*L +L*L*L);
		nnentry.push_back((i % L - 1 + L) % L + (i / L)*L +2* L*L*L);
		nnentry.push_back((i % (L*L*L) - L*L + L*L*L) % (L*L*L) +3 * L*L*L);
		neighbors.push_back(nnentry);
		nnentry.clear();
	}
	for (int i = L*L*L; i < 2*L*L*L; i++) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		nnentry.push_back(i + L*L*L); nnentry.push_back(i - L*L*L); nnentry.push_back(i + 2 * L*L*L);
		nnentry.push_back((i % (L*L) + L ) % (L*L) - i%L + (i % L - 1 + L) % L + (i / (L*L)*L*L) +  L*L*L);
		nnentry.push_back((i % (L*L*L) - L*L + L*L*L) % (L*L*L) + (i % (L*L) + L ) % (L*L) - i % (L*L) + (i / (L*L*L)*L*L*L) + 2 * L*L*L);
		nnentry.push_back((i % (L*L) + L ) % (L*L) + (i / (L*L))*L*L - L*L*L);
		neighbors.push_back(nnentry);
		nnentry.clear();
	}
	for (int i = 2 * L*L*L; i < 3 * L*L*L; i++) {
		
		///////////////////////////////////////////////////////////////////////////////////////////////
		nnentry.push_back(i + L*L*L); nnentry.push_back(i - L*L*L); nnentry.push_back(i - 2 * L*L*L);
		nnentry.push_back((i % (L*L*L) - L*L + L*L*L) % (L*L*L) + (i % L + 1 ) % L + (i / (L*L*L))*L*L*L - i%L+L*L*L);
		nnentry.push_back((i % L + 1 ) % L + (i / L)*L - 2 * L*L*L);
		nnentry.push_back((i % (L*L) - L + L*L) % (L*L) - i%L + (i % L + 1 ) % L +(i/(L*L)*L*L)-  L*L*L);
		neighbors.push_back(nnentry);
		nnentry.clear();
	}
	for (int i=3*L*L*L;i<4*L*L*L;i++){
		nnentry.push_back(i - 2*L*L*L); nnentry.push_back(i - L*L*L); nnentry.push_back(i - 3 * L*L*L);
		nnentry.push_back((i % (L*L*L) +L*L ) % (L*L*L)+ (i % (L*L) - L + L*L) % (L*L)-i%(L*L)+(i/(L*L*L)*L*L*L)-2*L*L*L);
		nnentry.push_back((i % (L*L*L) + L*L ) % (L*L*L)+(i/(L*L*L)*L*L*L)-3*L*L*L);
		nnentry.push_back((i % (L*L*L) + L*L ) % (L*L*L) + (i % L - 1 + L) % L-i%L+ (i / (L*L*L)*L*L*L)-L*L*L);
		neighbors.push_back(nnentry);
		nnentry.clear();
	}

}
 


double energy(int&nsites,int&L,double &Jnn, double &Dnn,double &H, vector<int> &config, vector< std::vector<int> > &neighbors,
std::vector< std::vector<double> > &fullcoords,
std::vector< std::vector<int> > &ijkt,
std::vector< std::vector<double> > &deltaJ)
{
	
	//make_pyrochlore4(L,  fullcoords, ijkt, neighbors);
	double energycalc = 0.0;

	for (int i = 0; i<nsites; i++)
	{
	
for (int k = 0; k < neighbors[i].size(); k++)
		{
			energycalc = energycalc + (deltaJ[i][k])*double(config[i]) * double(config[neighbors[i][k]]);
			//cout<< double(config[neighbors[i][k]])<<endl;		
}
	}
	return energycalc/(2.0);
}

double diffenergy_nnmodel(int& nsites,double &Jnn, double &Dnn,int &L, double &H,int &chosensite, int oldvalue, int newvalue, vector<int> &config,vector< std::vector<int> > &neighbors,
	std::vector< std::vector<double> > &fullcoords,
	std::vector< std::vector<int> > &ijkt, std::vector< std::vector<double> > &deltaJ)
{
	
//make_pyrochlore4(L, fullcoords, ijkt, neighbors);
//cout<<"3"<<endl;
	double energycalc = 0.0;
	//energycalc = energycalc + (8/3)*H*double(newvalue-oldvalue)*(rx[ijkt[chosensite][3]]+ry[ijkt[chosensite][3]]+rz[ijkt[chosensite][3]]);
	for (int k = 0; k < neighbors[chosensite].size(); k++)
	{
		energycalc = energycalc + (deltaJ[chosensite][k])*double(newvalue - oldvalue) * double(config[neighbors[chosensite][k]]);
	}
	return energycalc;
}


int str_to_int(string str)
{
	return boost::lexical_cast<int>(str);
}
long long int str_to_long_long(string str)
{
	return boost::lexical_cast<long long>(str);
}
double str_to_d(string str)
{
	return boost::lexical_cast<long double>(str);
}

void read_file(string infilename, int &L, double &Jnn, double &Dnn, double& dj,double &T1,double &T2, double &H,long long int &nsamples, int &nwait, int &seed,string &outfilename)
{
	string str_ret; 
	bool found; 

	search_for(string("L"),infilename,str_ret,found);
	if (found) {L=str_to_int(str_ret);} else {L=4;}

	search_for(string("Jnn"),infilename,str_ret,found);
	if (found) {Jnn=str_to_d(str_ret);} else {Jnn=0;}
	
	search_for(string("Dnn"),infilename,str_ret,found);
	if (found) {Dnn=str_to_d(str_ret);} else {Dnn=0;}
	
	search_for(string("dj"), infilename, str_ret, found);
	if (found) {dj=str_to_d(str_ret); }else {dj=0;}

        search_for(string("T1"),infilename,str_ret,found);
	if (found) {T1=str_to_d(str_ret);} else {T1=1.00;}
         search_for(string("T2"),infilename,str_ret,found);
	if (found) {T2=str_to_d(str_ret);} else {T2=1.00;}
       
        search_for(string("H"),infilename,str_ret,found);
	if (found) {H=str_to_d(str_ret);} else {H=0;}
	
        search_for(string("nsamples"),infilename,str_ret,found);
	if (found) {nsamples=str_to_long_long(str_ret);} else {nsamples=10000000;}
        
        search_for(string("nwait"),infilename,str_ret,found);
	if (found) {nwait=str_to_int(str_ret);} else {nwait=0;}
        
	search_for(string("seed"),infilename,str_ret,found);
	if (found) {seed=str_to_int(str_ret);} else {seed=1;}

        search_for(string("outfilename"),infilename,outfilename,found);
	if (not found) {outfilename="C.txt";}
}


int main(int argc, char *argv[])
{
	double Jnn=0.0;
	double Dnn=0.0;
	double dj = 0.0;
	double H=0.0;
	long long int nsamples=10;
	int nwait;
	double T1;
	double T2;
	double Tem;
	int seed;
	int L; 
	double q=2.0;
	
   string outfilename;
	string infilename; 

	infilename=argv[1];

	read_file(infilename,L,Jnn,Dnn,dj,T1,T2,H,nsamples,nwait,seed,outfilename);
	std::srand(seed);
	int nsites=4*L*L*L;
//cout<<"L   = "<<L<<endl;
	//cout<<"Jnn = "<<Jnn<<endl;
        //int nsamples = 1000000;
	//int nwait = 900000;
	//double Jnn = -1.24; //Jnn=J/3.
	//double Dnn =2.35 *pow(1 / 0.35530, 3);//Dnn=5D/3.
	//int L = 3;
	//double H = 0;
	////double T = 100;q
	//nsamples=40000000000LL;
	
	//for (int i=0;i<100;i++) cout<<"i = "<<i<<endl;
	
	std::vector< std::vector<int> > neighbors;
	std::vector< std::vector<double> > fullcoords;
	std::vector< std::vector<int> > ijkt;
	make_pyrochlore4(L,  neighbors);
	std::vector <vector <double> > deltaJ;
	for (int i = 0; i < nsites; i++)
	{
		vector<double> x;
			for (int j = 0; j < neighbors[i].size(); j++)
			{
				double ran = ourrandom();
				double y = dj*(ran - 0.5)+Jnn;
				//cout << y << endl;
				x.push_back(y);
			}
		deltaJ.push_back(x);
		x.clear();
	}
	
	//cout <<"nsites="<< nsites << endl;


	//string outfilename = "C.txt";
	ofstream outfile;
	const char *cstr = (outfilename).c_str();
	outfile.open(cstr);
	
	//for (int t = 0; t < 100; t++) 
 //{
//	double T = double(t)/10.0;
	int accept = 0;
	int reject = 0;
	int nmeas = 0;
	
	//// Create configuration and measure energy

	vector<int> config = create_fm_config(L);
	//for (int i=0;i<nsites;i++) cout<<config[i];

//cout<<config.size()<<endl;
//	cout << "start Metropolis" << endl;
//	cout<<nsamples<<endl;
	// Run Metropolis
double energycurrent=0.0;	
//double energycurrent = energy(nsites, L, Jnn, Dnn,H, config, neighbors, fullcoords, ijkt);
	//cout<<energycurrent/double(nsites)<<endl;
	double totalenergy = 0.0;
	double energysquare = 0.0;
        double totalmag = 0.0; 	
        double magcurrent=0.0;
       
//cout<<"1"<<endl;
 double sqnow= 0.0;
//for(int i=0;i<(config.size());i++){magcurrent+=4.6188*double(config[i])*(rx[ijkt[i][3]]+ry[ijkt[i][3]]+rz[ijkt[i][3]]);};                
//vector<double> sq(nsites);
/*for(int i=0;i<nsites;i++)
{sq[i]=0;}
;       		      for (int m=0;m<(config.size());m++)
         {
for (int i=0;i<(config.size());i++){
         sq[m]+=double(config[m])*double(config[i])*cos(2*q*M_PI*(fullcoords[i][2]-fullcoords[m][2]));
         };;
sqnow+=sq[m];
}*/
//cout<<sqnow/pow(double(nsites),2)<<endl;
//cout<<magcurrent*magcurrent/pow(double(nsites),2)<<endl;
//for (int i = 0; i < 4; i++){cout<<ry[i]<<"	"<<endl;};

//cout<<nsamples<<endl;
for (long long int n = 0; n < nsamples; n++)

{Tem=T1;
//cout<<nsites<<endl;
if(n>(nsamples/5)*2){Tem=T2;}
// else H=0.335; }
;		int chosensite = randint(nsites); 
//cout<<chosensite<<endl;		
double energydiff = diffenergy_nnmodel(nsites, Jnn, Dnn, L, H, chosensite, config[chosensite], -config[chosensite], config, neighbors, fullcoords, ijkt,deltaJ);
//cout<<config[chosensite]<<config[neighbors[chosensite][0]]<<config[neighbors[chosensite][1]]<<config[neighbors[chosensite][2]]<<config[neighbors[chosensite][3]]<<config[neighbors[chosensite][4]]<<config[neighbors[chosensite][5]]<<energydiff<<endl;	
	double boltzmann = exp(-energydiff / Tem);
		double r = ourrandom(); //cout << r<<"	" ; cout << boltzmann << "	";
		//outfile << "T, E  =" << T << "   " << energycurrent << endl;
		//cout << energycurrent<<","<<endl;
		//cout<<"Random number generated is = "<<r<<endl;
//if (abs(energydiff)<0.01)
//{
//	if (r < x[chosensite])
//	{
//cout<<"3"<<endl;	
//	energycurrent = energycurrent + energydiff;
//		config[chosensite] = -config[chosensite];
//		accept += 1;
		//  magcurrent=magcurrent+2*4.6188*double(config[chosensite])*(rx[ijkt[chosensite][3]]+ry[ijkt[chosensite][3]]+rz[ijkt[chosensite][3]]);;
		//outfile<<rx[ijkt[chosensite][3]]<<"	"<<ry[ijkt[chosensite][3]]<<"	"<<rz[ijkt[chosensite][3]]<<endl;
		/*for(int i=0;i<nsites;i++){
		if(i!=chosensite){ sqnow+=4*double(config[chosensite])*double(config[i])*cos(2*q*M_PI*(fullcoords[i][2]-fullcoords[chosensite][2]));}
		;


		}*/
		
//	}
//	else
//	{
//		reject += 1;
//	}
//}
//else 
//{

	if (r < boltzmann)
	{
		energycurrent = energycurrent + energydiff;
		config[chosensite] = -config[chosensite];
		accept += 1;
		//  magcurrent=magcurrent+2*4.6188*double(config[chosensite])*(rx[ijkt[chosensite][3]]+ry[ijkt[chosensite][3]]+rz[ijkt[chosensite][3]]);;
//outfile<<rx[ijkt[chosensite][3]]<<"	"<<ry[ijkt[chosensite][3]]<<"	"<<rz[ijkt[chosensite][3]]<<endl;
/*for(int i=0;i<nsites;i++){
		if(i!=chosensite){ sqnow+=4*double(config[chosensite])*double(config[i])*cos(2*q*M_PI*(fullcoords[i][2]-fullcoords[chosensite][2]));}
;


	}*/
		
	}
	else
	{
		reject += 1;
	}
//}

//}
		if ( n%nsites==0)
{


		//	totalenergy += (energycurrent); nmeas += 1; //cout << totalenergy<<"	";
		//	energysquare += pow(energycurrent, 2.0); //cout << sqrt(energysquare)<< "	";
                    //    totalmag+=magcurrent;
	//	}
	//	if (n%(2*nsites)==0) {
	outfile<<n/nsites<<" "<<energycurrent/double(nsites)/*<<" "<<sqnow/pow(double(nsites),2)*/<<endl;
	}
	}
//	totalenergy = totalenergy / double(nmeas);
  //      totalmag=totalmag/(double(nmeas)*double(nsites));
//	outfile << "Accept =" << accept << endl;
  //      outfile << "<M>="<<totalmag<<endl;
//	outfile << "Reject =" << reject << endl;
        //outfile << "<E>=" << totalenergy << endl;
//	double spheat=((energysquare / double(nmeas)) - pow(totalenergy, 2))/pow(T,2);
	//double spheatpersite=spheat/double(nsites);
	//cout << "T=" << T << "c=" << energysquare / double(nmeas) - pow(totalenergy, 2) << endl;
//	outfile << "T="<< T <<"spheat=" << spheat << "spheat/site="<< spheatpersite << endl;
	//outfile << T << config << endl;
	outfile.close();
	return 0;
  //      cout <<"program finished"<<endl;
}
