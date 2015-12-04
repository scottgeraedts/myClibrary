#ifndef SCOTTUTILS_H
#define SCOTTUTILS_H

#include<complex>
#include<string>
#include<fstream>
#include<algorithm>
#include<iostream>
#include<sstream>
#include <vector>
#include<bitset>
#include<Eigen/Core>

using namespace std;

//given a parameters file, it reads a value and returns it. but it also returns a default value if it can't find a value to read
template<class T>
T value_from_file(ifstream &infile, T def){
	string line,data;
	stringstream convert;
	T out;
	
	if(getline(infile,line)){
		//get line from file, put it in line
		//put line in convert
		convert.str("");
		convert<<line;
		//put part before pound sign into data
		getline(convert,data,'#');
		//put data into convert
		convert.str("");
		convert<<line;
		
		//put convert into out
		convert>>out;
		return out;
	}
	else return def;
}	

//given a vector, computes its level spacing ratio

void density_of_states(const vector<double> &x, vector<double> &p, const vector<double> &energy_grid, int start=-1, int end=-1);
double level_spacings(const vector<double> &x, const vector<double> &p, const vector<double> &energy_grid, int start=-1, int end=-1);
double stupid_spacings(const vector<double> &x, int label, int start=-1, int end=-1);

//computes kullback-leibler divergences (see arxiv 1411.0660)
template<class ART>
double kullback_leibler(const vector<ART> &x, const vector<ART> &y, int label){
	double out=0,p,q,temp;
	stringstream filename;
	filename<<"klhist"<<label;
	ofstream outfile;
	outfile.open(filename.str().c_str(),ios::app);
	for(unsigned int i=0;i<x.size();i++){
		p=abs(x[i])*abs(x[i]);
		q=abs(y[i])*abs(y[i]);
		temp=(p-q)*log(p/q);
		if(i%100==0) outfile<<temp<<endl; //no need to save all points for each ROD, a few will be enough for a nice histogram
		out+=temp;
	}
	outfile.close();
	return out;
}

//writes a vector to a provided file, optionally divides the vector by something first
void write_vector(Eigen::Matrix<double,-1,1> &data, string filename, double C=1.);

//counts the number of set bits in an integer
int count_bits(int x);

///***Some functions related to calculating Clebsch-Gordan coefficients, which uses a lot of factorials, etc
//computes the products of numbers from start to stop (a partial factorial)
double factorial(int start,int stop);
double factorial(int n);
//x choose p
long int comb(int x,int p);
//a function for computing powers of integers
int intpow(int x, int p);

double ClebschGordan(int a,int b,int L, int NPhi);

//***given an integer (which can be thought of as a bitstring) and a set of integers (bits to flip) and a vector of integers (possible end states)
//flips the bits and finds their positions in the vector
int lookup_flipped(int state, int a, const vector<int> &states);
int lookup_flipped(int state,int a, int b, const vector<int> &states);
#endif
