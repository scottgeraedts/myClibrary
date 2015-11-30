#include<complex>
#include<string>
#include<fstream>
#include<algorithm>
#include<iostream>
#include<sstream>
#include <vector>

using namespace std;

//given a parameters file, it reads a value and returns it. but it also returns a default value if it can't find a value to read
template<class T>
T value_from_file(ifstream &infile, T def){
	string line;
	stringstream convert;
	T out;
	if(getline(infile,line)){
		convert.str("");
		convert<<line;
		convert>>out;
		return out;
	}
	else return def;
}	

//given a vector, computes its level spacing ratio

void density_of_states(const vector<double> &x, vector<double> &p, const vector<double> &energy_grid){
	int count, mark=0;
	double dE=energy_grid[1]-energy_grid[0];
	for(int i=0;i<(signed)energy_grid.size();i++){
		//for density of states, count how many states are within a certain energy window
		count=0;
		while(true){
			if(x[mark]<energy_grid[i]-0.5*dE) cout<<"something didn't make sense in level_spacings "<<x[mark]<<" "<<energy_grid[i]<<endl;
			if(x[mark]<energy_grid[i]+0.5*dE){
				count++;
//				cout<<"added energy "<<x[mark]<<" to grid point "<<energy_grid[i]<<endl;
				mark++;
				if(mark==(signed)x.size()) break;
			}
			else{
//				cout<<"done with grid point "<<energy_grid[i]<<", energy at "<<x[mark]<<endl;
				break;
			}
		}
		p[i]+=count/(1.*x.size()*dE);
		if(mark==(signed)x.size()) break;
	}

}
double level_spacings(const vector<double> &x, const vector<double> &p, const vector<double> &energy_grid){

	vector<double> s(x.size(),0);
//	vector<double> ps(ngrid,0);
	
	vector<double> integrated_DOS(p.size(),0);
	//compute integrated density of states with trapezoid rule
	for(unsigned int i=1;i<p.size();i++)
		integrated_DOS[i]=0.5*(p[i]+p[i-1])*(energy_grid[i]-energy_grid[i-1])+integrated_DOS[i-1];
	
	//use integrated density of states to get S, do a linear approximation between the different vales of integrated_DOS
	vector<double>::const_iterator low;
	int pos;
	for(unsigned int i=0;i<s.size();i++){
		low=lower_bound(energy_grid.begin(),energy_grid.end(),x[i]);
		pos=low-energy_grid.begin();
//		s[i]=integrated_DOS[pos];
		s[i]=integrated_DOS[pos-1]+(x[i]-energy_grid[pos-1])*(integrated_DOS[pos]-integrated_DOS[pos-1])/(energy_grid[pos]-energy_grid[pos-1]);
	}
	
//	//integrate density to get s
//	for(int i=0;i<(signed)s.size();i++){
//		for(int j=0;j<=(signed)p.size();j++){
//			if(energy_grid[j]>x[i]) break;
//			else s[i]+=p[j];
//		}
//	}
//	for(int i=0;i<s.size();i++) cout<<s[i]<<endl;
//	//for testing purposes, compute density of s
//	for(int i=0;i<ngrid;i++){
//		//for density of states, count how many states are within a certain energy window
//		count=0;
//		while(true){
//			if(s[mark]<energy_grid[i]-0.5*dE) cout<<"something didn't make sense in level_spacings test"<<endl;
//			if(s[mark]<energy_grid[i]+0.5*dE){
//				count++;
//				mark++;
//			}
//			else break;
//		}
//		ps[i]+=count;
//	}
//	for(int i=0;i<p.size();i++) cout<<energy_grid[i]<<" "<<p[i]<<" "<<s[i]<<" "<<ps[i]<<endl;
	
	//compute r
	ofstream csout;
	csout.open("tempS",ios::app);
	double r=0;
	int count=0;
	for(unsigned int i=s.size()/6;i<5*s.size()/6;i++){
		if(s[i+1]-s[i]>s[i]-s[i-1]) r+=(s[i]-s[i-1])/(s[i+1]-s[i]);
		else r+=(s[i+1]-s[i])/(s[i]-s[i-1]);
			count++;
//			if(x[i+1]-x[i]>x[i]-x[i-1]) r+=(x[i]-x[i-1])/(x[i+1]-x[i]);
//			else r+=(x[i+1]-x[i])/(x[i]-x[i-1]);
			csout<<x[i]<<" "<<s[i]<<" ";
			if(x[i+1]-x[i]>x[i]-x[i-1]) csout<<(x[i]-x[i-1])/(x[i+1]-x[i])<<" ";
			else csout<<(x[i+1]-x[i])/(x[i]-x[i-1])<<" ";
			if(s[i+1]-s[i]>s[i]-s[i-1]) csout<<(s[i]-s[i-1])/(s[i+1]-s[i])<<" "<<endl;
			else csout<<(s[i+1]-s[i])/(s[i]-s[i-1])<<endl;
//		}
	}
	csout.close();
//	for(int i=0;i<p.size();i++) cout<<energy_grid[i]<<" "<<p[i]<<endl;
	return r/(1.*count);
}

//counts the number of set bits in an integer
int count_bits(int x){
	int out=0, i=0,found_bits=0;
	while(x!=found_bits){
		if(1<<i & x){
			out++;
			found_bits+=1<<i;
		}
		i++;
	}
	return out;
}

///***Some functions related to calculating Clebsch-Gordan coefficients, which uses a lot of factorials, etc
//computes the products of numbers from start to stop (a partial factorial)
double factorial(int start,int stop){
	if (stop==start) return 1;
	else return start*factorial(start-1,stop);
}	
double factorial(int n){
	if (n==1 || n==0) return 1;    
	return n*factorial(n-1);
}	
//x choose p
long int comb(int x,int p){
	return factorial(x,x-p)/factorial(p);	
}
//a function for computing powers of integers
int intpow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = intpow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}
double ClebschGordan(int a,int b,int L, int NPhi){
//calculate CG coefficients from the formula on wikipedia
//m1=a-Q,m2=b-Q, 2Q=NPhi-1,these are stored this way because sometimes they are half-integer
//may eventually need to tabulate these since they are pretty slow to calculate
	double prefactor1=sqrt((2*L+1)*factorial(L)*factorial(L)*factorial((NPhi-1)-L)/(1.*factorial((NPhi-1)+L+1)) );
	double prefactor2=sqrt(factorial(L+a+b-(NPhi-1))*factorial(L-a-b+(NPhi-1))*factorial((NPhi-1)-a)*factorial((NPhi-1)-b)*factorial(a)*factorial(b));
	double sum=0.;
	int sign=1;
	for(int k=0;k<=b;k++){
		if (k%2==1) sign=-1;
		else sign=1;
		if(L-b+k<0) continue;
		if((NPhi-1)-L-k<0 || L-(NPhi-1)+a+k<0 || (NPhi-1)-a-k<0) continue;
		sum+=(sign*1.)/(1.*factorial((NPhi-1)-L-k)*factorial((NPhi-1)-a-k)*factorial(b-k)*factorial(L-(NPhi-1)+a+k)*factorial(L-b+k)*factorial(k));
	}
	return prefactor1*prefactor2*sum;
}

