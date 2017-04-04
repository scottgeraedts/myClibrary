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
#include<cstdarg>

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

template<class T>
void Eigen_To_Std(const Eigen::Matrix<T,-1,1> &invec, vector<T> &outvec){
	//turn Eigen::Vector into array
	T *w=new T[invec.size()];
	Eigen::Map <Eigen::Matrix<T,-1,1> > (w,invec.size(),1)=invec; //using just out.data() fails for an unknown reason

	//turn array into std vector
	outvec=vector<T>(w,w+invec.size());
	delete [] w;
}	

template<class T>
Eigen::Matrix<T,-1,1> Std_To_Eigen(vector<T> &invec){
	T *v=&invec[0];

	//turn invec into an Eigen::Vector
	Eigen::Map <Eigen::Matrix<T,-1,1> > mapped_v(v,invec.size());
	return mapped_v;

}
//given a vector, computes its level spacing ratio
void density_of_states(const vector<double> &x, vector<double> &p, const vector<double> &energy_grid, int start=-1, int end=-1);
double level_spacings(const vector<double> &x, const vector<double> &p, const vector<double> &energy_grid, int start=-1, int end=-1);
double stupid_spacings(const vector<double> &x, int label=0, int start=-1, int end=-1);

vector<double> unfoldE(const vector<double> &x, int mesh=200);
double compute_r(const vector<double> &s, int start=-1, int end=-1);
vector<double> spacings(const vector<double> &x, int start=-1, int end=-1);

vector<double> make_grid(const vector<double> &x, int mesh);
vector<double> make_DOS(const vector<double> &x, const vector<double> &energy_grid);
vector<double> make_S(const vector<double> &x, const vector<double> &energy_grid, const vector<double> &integrated_DOS);

//computes kullback-leibler divergences (see arxiv 1411.0660)
template<class ART>
double kullback_leibler(const vector<ART> &x, const vector<ART> &y){
	double out=0,p,q,temp;
	for(unsigned int i=0;i<x.size();i++){
		p=abs(x[i])*abs(x[i]);
		q=abs(y[i])*abs(y[i]);
		temp=(p-q)*log(p/q);
		temp=(p)*log(p/q);
		out+=temp;
	}
	return out;
}

//writes a vector to a provided file, optionally divides the vector by something first
void write_vector(Eigen::Matrix<double,-1,1> &data, string filename, double C=1.);

//prints a std:vector
template<class T> void print_vector(const vector<T> &in){
	cout<<"--------------"<<endl;
	for(int i=0; i<in.size(); i++){
		cout<<in[i]<<endl;
	}
}

//counts the number of set bits in an integer
int count_bits(unsigned int x);
//returns the locations of flipped bits in a bitset
vector<int> bitset_to_pos(unsigned int x,int NPhi);
//sees if a certain bit is set
int bittest(unsigned int state, int bit);
//advances all the positions of bits in an integer, but only mod NPhi
template<class T>
T cycle_bits(T in, int NPhi){
	T out=0,one=1;
	for(int i=0;i<NPhi;i++){
		if( in & one<<i){
			if(i==NPhi-1) out+=1;
			else out+=one<<(i+1);
		}
	}
//	cout<<"cycle in: "<<(bitset<12>)in<<" out: "<<(bitset<12>)out<<endl;
	return out;
}
unsigned int invert_bits(unsigned int in, int NPhi);
//computes how different two bitstrings are
int distance_calc(const vector<int> &a, const vector<int> &b);

template<class T>
T move_bit(T in, int NPhi, int x, int dx){
	T one=1;
	if (!(in & one<<x) || (in & one<<((x+dx)%NPhi) && dx!=0) ){
		throw 0;
	}
	//cout<<(bitset<6>)in<<" "<<(bitset<6>)(in ^ 1<<x)<<" "<<(bitset<6>)((in^1<<x) ^ (1<<(x+dx)%NPhi))<<endl;
	return (in ^ one<<x) ^ (one<<(x+dx)%NPhi);
}

///***Some functions related to calculating Clebsch-Gordan coefficients, which uses a lot of factorials, etc
//computes the products of numbers from start to stop (a partial factorial)
double factorial(int start,int stop);
double factorial(int n);
//x choose p
long int comb(int x,int p);
//a function for computing powers of integers
unsigned long int intpow(int x, int p);

//useful when computing bloch symmetries, moves all bits by M
template <class T>
T cycle_M(T in, int NPhi, int M, int &sign){
	T out=0,old=in;
	for(int i=0;i<M;i++){
		out=cycle_bits(old,NPhi);
		if (out<old && NPhi%4==0) sign*=-1;
		old=out;
	}
	return out;
}

double ClebschGordan(int dj1, int dj2, int dm1, int dm2, int dJ);
double Wigner3j(int dj1, int dj2, int dj3, int dm1, int dm2, int dm3);
double Wigner6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6);

//***given an integer (which can be thought of as a bitstring) and a set of integers (bits to flip) and a vector of integers (possible end states)
//flips the bits and finds their positions in the vector
//numbits: the number of bits to flip
//remaining arguments: which bits to flip
template<class T>
T lookup_flipped(int i, const vector<T> &states, int numbits, ...){
	T compare=states[i];
	va_list ap;
	T temp;
	T one=1;
	va_start(ap,numbits);
	for(int j=0;j<numbits;j++){
		temp=va_arg(ap,T);
		compare=compare ^  ((T) (one<<temp) );
	}
	va_end(ap);	
	
	typename vector<T>::const_iterator low;
//	auto low=states.begin();
	low=lower_bound(states.begin(),states.end(),compare);
	if(low!=states.end()) return (low-states.begin());
	else{
		cout<<"error in lookup_flipped: "<<(bitset<40>)states[i]<<" "<<(bitset<40>)compare<<" "<<states[i]<<" "<<compare<<endl;
		va_list ap;
		va_start(ap,numbits);
		for(int j=0;j<numbits;j++){
			cout<<va_arg(ap,int)<<endl;
		}
		va_end(ap);
		exit(0);	
		return 0;
	}

}
//int lookup_flipped(int state,int a, int b, const vector<int> &states);
int permute_sign(int n, ...);
int lil_sign(int x);


//puts in minus signs for correct Fermi statistics
//for four-body terms, we have something like c^dagger_a c_d. This gives a (-1) for every electron between a and d
//this also works fine when you consider that we are performing two hops
template<class T>
int adjust_sign(int a,int b,T state){
	int sign=1;
	T one=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state & one<<i) sign*=-1;
	return sign;	
}
//puts in minus signs for correct Fermi statistics
//for four-body terms, we have something like c^dagger_a c_d. This gives a (-1) for every electron between a and d
//this also works fine when you consider that we are performing two hops
template<class T>
int adjust_sign(int a,int b,int c,int d,T state){
	int sign=1;
	T one=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state & one<<i && i!=c && i!=d) sign*=-1;
	start=c, end=d;
	if (c>d){ start=d; end=c;}
	for(int i=start+1;i<end;i++)
		if(state & one<<i && i!=a && i!=b) sign*=-1;
	return sign;	
}
template<class T>
int adjust_sign(int a, int b, int c, int d, int e, int f, T state){
	int sign=1;
	T one=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state & one<<i && i!=c && i!=d && i!=e && i!=f) sign*=-1;
	start=c, end=d;
	if (c>d){ start=d; end=c;}
	for(int i=start+1;i<end;i++)
		if(state & one<<i && i!=a && i!=b && i!=e && i!=f) sign*=-1;
	start=e, end=f;
	if (e>f){ start=f; end=e;}
	for(int i=start+1;i<end;i++)
		if(state & one<<i && i!=a && i!=b && i!=c && i!=d) sign*=-1;

	return sign;	
}

vector<int> sort_indexes(const vector<double> &x);
int supermod(int,int);
#endif
