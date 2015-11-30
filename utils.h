#include<complex>
#include<string>
#include<fstream>
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
		p[i]+=count/(1.*x.size());
		if(mark==(signed)x.size()) break;
	}

}
double level_spacings(const vector<double> &x, const vector<double> &p, const vector<double> &energy_grid){

	vector<double> s(x.size(),0);
//	vector<double> ps(ngrid,0);
	
	//integrate density to get s
	for(int i=0;i<(signed)s.size();i++){
		for(int j=0;j<=(signed)p.size();j++){
			if(energy_grid[j]>x[i]) break;
			else s[i]+=p[j];
		}
	}
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
	for(int i=1;i<s.size()-1;i++){
//		if(s[i+1]-s[i]>s[i]-s[i-1]) r+=(s[i]-s[i-1])/(s[i+1]-s[i]);
//		else r+=(s[i+1]-s[i])/(s[i]-s[i-1]);
		if(x[i]>-4 && x[i]<4){
			count++;
			if(x[i+1]-x[i]>x[i]-x[i-1]) r+=(x[i]-x[i-1])/(x[i+1]-x[i]);
			else r+=(x[i+1]-x[i])/(x[i]-x[i-1]);
			csout<<x[i]<<" "<<s[i]<<" ";
			if(x[i+1]-x[i]>x[i]-x[i-1]) csout<<(x[i]-x[i-1])/(x[i+1]-x[i])<<" ";
			else csout<<(x[i+1]-x[i])/(x[i]-x[i-1])<<" ";
			if(s[i+1]-s[i]>s[i]-s[i-1]) csout<<(s[i]-s[i-1])/(s[i+1]-s[i])<<" "<<endl;
			else csout<<(s[i+1]-s[i])/(s[i]-s[i-1])<<endl;
		}
	}
	csout.close();
//	for(int i=0;i<p.size();i++) cout<<energy_grid[i]<<" "<<p[i]<<endl;
	return r/(1.*count);
}
