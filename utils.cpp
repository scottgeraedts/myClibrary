#include "utils.h"

////given a parameters file, it reads a value and returns it. but it also returns a default value if it can't find a value to read
//template<class T>
//T value_from_file(ifstream &infile, T def){
//	string line;
//	stringstream convert;
//	T out;
//	if(getline(infile,line)){
//		convert.str("");
//		convert<<line;
//		convert>>out;
//		return out;
//	}
//	else return def;
//}	

double compute_r(const vector<double> &s, int start, int end){
	double r=0;
	int count=0;
	if(start==-1) start=0;
	if(end==-1) end=s.size();
	for(int i=start+1;i<end-start-1;i++){
		if(s[i+1]-s[i]>s[i]-s[i-1]) r+=(s[i]-s[i-1])/(s[i+1]-s[i]);
		else r+=(s[i+1]-s[i])/(s[i]-s[i-1]);
		count++;
	}
	return r/(1.*count);
}
//computes raw level spacings (not r)
vector<double> spacings(const vector<double> &x, int start, int end){
	if(start==-1) start=0;
	if(end==-1) end=x.size();
	vector<double> out(end-start-1);
	for(int i=start;i<end-1;i++){
		out[i-start]=(x[i+1]-x[i])/(x[end-1]-x[start])*(end-start);
	}
	return out;
}
vector<double> unfoldE(const vector<double> &x, int mesh){
	vector<double> s(x.size());
	vector<double> energy_grid(mesh);
	vector<double> p(mesh); //density of states
	double dE=(x.back()-x[0])/(1.*mesh-5.);
	for(int i=0;i<mesh;i++) energy_grid[i]=x[0]+dE*(i-1);
	int mark=0,count;

	for(int i=0;i<mesh;i++){
		//for density of states, count how many states are within a certain energy window
		count=0;
		while(true){
//			if(x[mark]<energy_grid[i]) cout<<"something didn't make sense in level_spacings "<<x[mark]<<" "<<energy_grid[i]<<" "<<mark<<" "<<i<<endl;
			if(x[mark]<energy_grid[i]+dE){
				count++;
				mark++;
				if(mark==x.size()) break;
			}
			else{
				//cout<<"done with grid point "<<energy_grid[i]<<", energy at "<<x[mark]<<endl;
				break;
			}
		}
		p[i]+=count/(1.*x.size()*dE);
		if(mark==x.size()) break;
	}
	//compute integrated density of states with trapezoid rule
	vector<double> integrated_DOS(mesh,0);
	for(int i=1;i<mesh;i++)
		integrated_DOS[i]=0.5*(p[i]+p[i-1])*(energy_grid[i]-energy_grid[i-1])+integrated_DOS[i-1];

	ofstream dos;
	dos.open("dos");
	for(int i=0;i<mesh;i++) dos<<energy_grid[i]<<" "<<p[i]<<" "<<integrated_DOS[i]<<endl;
	dos.close();
	//use integrated density of states to get S, do a linear approximation between the different vales of integrated_DOS
	vector<double>::const_iterator low;
	int pos;
	for(int i=0;i<(signed)x.size();i++){
		low=lower_bound(energy_grid.begin(),energy_grid.end(),x[i]);
		pos=low-energy_grid.begin();
		s[i]=integrated_DOS[pos-1]+(x[i]-energy_grid[pos-1])*(integrated_DOS[pos]-integrated_DOS[pos-1])/(energy_grid[pos]-energy_grid[pos-1]);
	}
	return s;
}
void density_of_states(const vector<double> &x, vector<double> &p, const vector<double> &energy_grid, int start, int end){
	if(start==-1) start=0;
	if(end==-1) end=x.size();
	int count, mark=start;
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
				if(mark==end) break;
			}
			else{
				break;
			}
		}
		p[i]+=count/(1.*(end-start)*dE);
		if(mark==end) break;
	}

}
/*x: data to take spacings of
  p: density of states
  energy_grid: grid to measure denstiy of states on
*/
double level_spacings(const vector<double> &x, const vector<double> &p, const vector<double> &energy_grid, int start, int end){

	if(start==-1) start=0;
	if(end==-1) end=x.size();	
	vector<double> s(end-start,0);
//	vector<double> ps(ngrid,0);
	
	vector<double> integrated_DOS(p.size(),0);
	//compute integrated density of states with trapezoid rule
	for(unsigned int i=1;i<p.size();i++)
		integrated_DOS[i]=0.5*(p[i]+p[i-1])*(energy_grid[i]-energy_grid[i-1])+integrated_DOS[i-1];
	
	//use integrated density of states to get S, do a linear approximation between the different vales of integrated_DOS
	vector<double>::const_iterator low;
	int pos;
	for(int i=0;i<end-start;i++){
		low=lower_bound(energy_grid.begin(),energy_grid.end(),x[i+start]);
		pos=low-energy_grid.begin();
//		s[i]=integrated_DOS[pos];
		s[i]=integrated_DOS[pos-1]+(x[i+start]-energy_grid[pos-1])*(integrated_DOS[pos]-integrated_DOS[pos-1])/(energy_grid[pos]-energy_grid[pos-1]);
	}
	
	//compute r
//	ofstream csout;
//	csout.open("tempS",ios::app);
	double r=0;
	int count=0;
	for(int i=1;i<end-start-1;i++){
		if(s[i+1]-s[i]>s[i]-s[i-1]) r+=(s[i]-s[i-1])/(s[i+1]-s[i]);
		else r+=(s[i+1]-s[i])/(s[i]-s[i-1]);
			count++;

	}
//	csout.close();
//	for(int i=0;i<p.size();i++) cout<<energy_grid[i]<<" "<<p[i]<<endl;
	return r/(1.*count);
}

double stupid_spacings(const vector<double> &x, int label, int start, int end){
	if(start==-1) start=0;
	if(end==-1) end=x.size();
	int count=0;
	double r=0,temp;
	stringstream filename;
	filename<<"stupidhist"<<label;
	ofstream outfile;
	outfile.open(filename.str().c_str(),ios::app);
	for(int i=start+1;i<end-1;i++){
		if(x[i+1]-x[i]>x[i]-x[i-1]) temp=(x[i]-x[i-1])/(x[i+1]-x[i]);
		else temp=(x[i+1]-x[i])/(x[i]-x[i-1]);
		r+=temp;
		outfile<<temp<<endl;
		count++;
	}
	outfile.close();
	return r/(1.*count);
}
//computes kullback-leibler divergences (see arxiv 1411.0660)
//template<class ART>
//double kullback_leibler(const vector<ART> &x, const vector<ART> &y){
//	double out=0,p,q;
//	for(unsigned int i=0;i<x.size();i++){
//		p=abs(x[i])*abs(x[i]);
//		q=abs(y[i])*abs(y[i]);
//		out+=(p-q)*log(p/q);
//	}
//	return out;
//}

//writes a vector to a provided file, optionally divides the vector by something first
void write_vector(Eigen::Matrix<double,-1,1> &data, string filename, double C){
	ofstream out;
	out.open(filename.c_str());
	data/=C;
	for(int i=0;i<data.size();i++) out<<data[i]<<" ";
	out<<endl;
	out.close();
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
double ClebschGordan(int dj1,int dj2, int dm1, int dm2, int dJ){
//calculate CG coefficients from the formula on wikipedia
//annoyingly, the inputs to CG coefficients can be half-integer. Therefore this function takes TWICE the actual number as input
//it uses this to calculate a bunch of factorials, the arguments to these factorials should always be integer so we just have to divide by 2
//may eventually need to tabulate these since they are pretty slow to calculate
	if( abs(dm1+dm2) > dJ || ( abs(dm1)%2!=dj1%2) || ( abs(dm2)%2!=dj2%2) || abs(dm1+dm2)%2!= dJ%2 || dj1+dj2<dJ ){
		cout<<"bad arguments to clebsch gordan calculator"<<endl;
		cout<<dj1<<" "<<dj2<<" "<<dm1<<" "<<dm2<<" "<<dJ<<endl;
		return 0;
	}
	double prefactor1=sqrt((dJ+1)*factorial( (dJ+dj1-dj2)/2 )*factorial( (dJ-dj1+dj2)/2 )*factorial( (dj1+dj2-dJ)/2 )/(1.*factorial( (dj1+dj2+dJ+2)/2 ) ) );
	double prefactor2=sqrt( factorial( (dJ+dm1+dm2)/2 )*factorial( (dJ-dm1-dm2)/2 )*factorial( (dj1-dm1)/2 )*factorial( (dj1+dm1)/2 )*factorial( (dj2-dm2)/2 )*factorial( (dj2+dm2)/2 ));
	double sum=0.,temp;
	int sign=1,farg;
	for(int k=0;k<=100;k++){
		if (k%2==1) sign=-1;
		else sign=1;
		temp=1/(1.*factorial(k));

		farg=(dj1+dj2-dJ-2*k)/2;
		if(farg<0) continue; else temp/=(1.*factorial(farg));

		farg=(dj1-dm1-2*k)/2;
		if(farg<0) continue; else temp/=(1.*factorial(farg));

		farg=(dj2+dm2-2*k)/2;
		if(farg<0) continue; else temp/=(1.*factorial(farg));

		farg=(dJ-dj2+dm1+2*k)/2;
		if(farg<0) continue; else temp/=(1.*factorial(farg));

		farg=(dJ-dj1-dm2+2*k)/2;
		if(farg<0) continue; else temp/=(1.*factorial(farg));

		sum+=(sign*1.)*temp;
	}
//	cout<<prefactor1<<" "<<prefactor2<<" "<<sum<<endl;
	return prefactor1*prefactor2*sum;
}

double Wigner3j(int dj1,int dj2,int dj3,int dm1,int dm2,int dm3){
	if(dm1+dm2+dm3!=0){
		cout<<"arguments in Wigner3j don't make sense "<<dm1<<" "<<dm2<<" "<<dm3<<endl;
		exit(0);
	}
	double out=1;
	if( (dj1+dj2-dm3)%4==2) out=-1;
	out/=sqrt(dj3+1.);
//	cout<<"3j: "<<dj1<<" "<<dj2<<" "<<dj3<<" "<<dm1<<" "<<dm2<<" "<<dm3<<" "<<out*ClebschGordan(dj1,dj2,dm1,dm2,dj3)<<endl;
	return out*ClebschGordan(dj1,dj2,dm1,dm2,dj3);
}

double Wigner6j(int dj1, int dj2, int dj3, int dj4, int dj5, int dj6){
	double out=0,temp1,temp2,temp3;
	int sign=1;
	for(int dm1=-dj1;dm1<=dj1;dm1+=2){
		for(int dm2=-dj2;dm2<=dj2;dm2+=2){
			for(int dm3=-dj3;dm3<=dj3;dm3+=2){
				//cout<<dm1<<" "<<dm2<<" "<<dm3<<endl;
				if(dm1+dm2-dm3!=0) continue;
				temp1=Wigner3j(dj1,dj2,dj3,dm1,dm2,-dm3);
//				cout<<"----"<<endl;
				for(int dm4=-dj4;dm4<=dj4;dm4+=2){
					for(int dm5=-dj5;dm5<=dj5;dm5+=2){
						//cout<<"    "<<dm4<<" "<<dm5<<endl;
						if(dm4-dm5+dm3!=0) continue;
						temp2=temp1*Wigner3j(dj4,dj5,dj3,dm4,-dm5,dm3);
//						cout<<"***"<<endl;
						for(int dm6=-dj6;dm6<=dj6;dm6+=2){
							//cout<<"         "<<dm6<<endl;
							if(-dm1+dm5+dm6!=0 || -dm4-dm2-dm6!=0) continue;
							if( (dj1+dj2+dj3+dj4+dj5+dj6-dm1-dm2-dm3-dm4-dm5-dm6)%4==0) sign=1;
							else sign=-1;
							temp3=sign*temp2*Wigner3j(dj1,dj5,dj6,-dm1,dm5,dm6)*Wigner3j(dj4,dj2,dj6,-dm4,-dm2,-dm6);
							out+=temp3;
//							cout<<dm1<<" "<<dm2<<" "<<dm3<<" "<<dm4<<" "<<dm5<<" "<<dm6<<" "<<temp3<<endl;
						}
					}
				}
			}
		}
	}
	return out;
}
//***given an integer (which can be thought of as a bitstring) and a set of integers (bits to flip) and a vector of integers (possible end states)
//flips the bits and finds their positions in the vector
int lookup_flipped(int i, const vector<int> &states, int numbits, ...){
	int compare=states[i];
	va_list ap;
	va_start(ap,numbits);
	for(int j=0;j<numbits;j++){
		compare=compare ^ 1<<va_arg(ap,int);
	}
	va_end(ap);
		
	vector<int>::const_iterator low;
	low=lower_bound(states.begin(),states.end(),compare);
	if(low!=states.end()) return (low-states.begin());
	else{
		cout<<"error in lookup_flipped: "<<(bitset<30>)states[i]<<" "<<(bitset<30>)compare<<endl;
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
int permute_sign(int n, ...){
	int sign=1;
	va_list ap;
	va_start(ap,n);
	vector<int> elements(n,0);
	
	for(int j=0; j<n; j++){
		elements[j]=va_arg(ap,int);
	}
	va_end(ap);
	
	int temp;
	for(int j=n;j>0;j--){
		for(int i=0;i<j-1;i++){
			if(elements[i]>elements[i+1]){
				 temp=elements[i];
				 elements[i]=elements[i+1];
				 elements[i+1]=temp;
				 sign*=-1;
			}
		}
	}
	return sign;
}		
int lil_sign(int x){
	if(x%2==0) return 1;
	else return -1;
}	
//int lookup_flipped(int state,int a, int b, const vector<int> &states){
//	int compare=state ^ 1<<a;
//	compare=compare ^ 1<<b;
//	vector<int>::const_iterator low;
//	low=lower_bound(states.begin(),states.end(),compare);
//	if(low!=states.end()) return (low-states.begin());
//	else{
//		cout<<"error in lookup_flipped: "<<(bitset<30>)state<<" "<<(bitset<30>)compare<<" "<<a<<" "<<b<<endl;
//		exit(0);	
//		return 0;
//	}
//}

