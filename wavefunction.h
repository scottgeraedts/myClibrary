//library of functions that can be applied to any wavefunction stores as a collection of spins
//right now the only thing here is entanglement entropy calculations, but more could be added
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include<Eigen/Dense>
#include<iostream>
#include<bitset>
#include<algorithm>

using namespace std;

template <class ART>
class Wavefunction{
private:
	int nBits; //number of physical bits (e.g. spins or flux quanta)
	Eigen::Matrix<ART,-1,-1> rho;

public:
	void init_wavefunction(int);

	//entanglement entropy stuff
	int ee_setup(int trunc_start, int trunc_end, const vector<int> &statep);
	void ee_compute_rho(const vector<ART> &evec, Eigen::Matrix<ART,-1,-1> &rho2,  const vector<int> &statep, double coeff);
	double ee_eval_rho(Eigen::Matrix<ART,-1,-1> &rho2);
	double entanglement_entropy(const vector< vector<ART> > &evec, const vector <int> &statep, int start, int end);
	double safe_mult(double,double);
	complex<double> safe_mult( complex<double>, complex<double>);
	vector<int> trunc_states;
	int trunc_part;
	
	//entanglement entropy stuff using SVD
	vector<double> entanglement_spectrum_SVD(const vector<ART> &evec, const vector<int> &statep, int to_trace);
	
	int rangeToBitstring(int start,int end);
};

template<class ART>
void Wavefunction<ART>::init_wavefunction(int n){ nBits=n; }

template<class ART>
int Wavefunction<ART>::rangeToBitstring(int start, int end){
	int out=0; //bitwise-AND with this gives just the component in the traced over states
	for(int i=0;i<nBits;i++){
		if ((i>=start && i <end && end>start) || (end<start && (i < end || i >= start) ) )
			out=out | 1<<i;
	}
	return out;
}

template<class ART>
vector<double> Wavefunction<ART>::entanglement_spectrum_SVD(const vector<ART> &evec, const vector<int> &statep, int to_trace){
	//get  mapping of original states to truncated/untruncated states
	vector<int> traced_states,untraced_states;
	bool found;
	for(int i=0;i<(signed)statep.size();i++){
		found=false;
		for(unsigned int j=0;j<traced_states.size();j++){
			if( ( statep[i] & to_trace) ==traced_states[j]){
				found=true;
				break;
			}
		}
		if(!found) traced_states.push_back( statep[i] & to_trace);
		found=false;
		for(unsigned int j=0;j<untraced_states.size();j++){
			if( ( statep[i] & ~to_trace) ==untraced_states[j]){
				found=true;
				break;
			}
		}
		if(!found) untraced_states.push_back( statep[i] & ~to_trace);
	}
		
	//loop over evec, putting elements in the right place
	Eigen::Matrix<ART,-1,-1> square_wf=Eigen::Matrix<ART,-1,-1>::Zero(traced_states.size(),untraced_states.size());
	vector<int>::iterator it;
	int traced_index, untraced_index;
	for(int i=0;i<(signed)evec.size();i++){
		it=find(traced_states.begin(),traced_states.end(),statep[i] & to_trace);
		traced_index=it-traced_states.begin();
		it=find(untraced_states.begin(),untraced_states.end(),statep[i] & ~to_trace);
		untraced_index=it-untraced_states.begin();
		square_wf(traced_index,untraced_index)=evec[i];
	}
	//SVD
	Eigen::JacobiSVD< Eigen::Matrix<ART,-1,-1> > svd(square_wf);
	//return square of output
	int outsize=traced_states.size();
	if(untraced_states.size()<traced_states.size()) outsize=untraced_states.size();
	vector<double> output(outsize);
	for(int i=0;i<outsize;i++){
		output[i]=pow(svd.singularValues()(i),2);
	}
	return output;
}
	
template<class ART>
double Wavefunction<ART>::entanglement_entropy(const vector< vector<ART> > &evec, const vector<int> &statep, int start, int end=-1){
	double out=0;
	if(end==-1) end=start+1;
	for(int i=0;i<nBits;i+=2){
		ee_setup(i,(i+nBits/2)%nBits,statep);
		rho=Eigen::Matrix<ART,-1,-1>::Zero(trunc_states.size(), trunc_states.size() );
		for(int j=start;j<end;j++)
			ee_compute_rho(evec[j],rho,statep,1/(1.*(end-start)) );
		out+=ee_eval_rho(rho)/(1.*nBits);
	}
	return out;
}
//sets up a vector of bitstrings, for states with parts traced over
template<class ART>
int Wavefunction<ART>::ee_setup(int trace_start, int trace_end, const vector<int> &statep){
	bool found;

	//compute trunc_part
	trunc_part=0; //bitwise-AND with this gives just the component in the traced over states
	for(int i=0;i<nBits;i++){
		if ((i>=trace_start && i <trace_end && trace_end>trace_start) || (trace_end<trace_start && (i < trace_end || i >= trace_start) ) )
			trunc_part=trunc_part | 1<<i;
	}

	//figure out how many reduced states there are		
	trunc_states.clear();
	for(unsigned int i=0;i<statep.size();i++){
		found=false;
		for(unsigned int j=0;j<trunc_states.size();j++)
			if( ( statep[i] & ~trunc_part) ==trunc_states[j]) found=true;
		if(!found) trunc_states.push_back( statep[i] & ~trunc_part);
	}
	return trunc_states.size();
}

template<class ART>
void Wavefunction<ART>::ee_compute_rho(const vector<ART> &evec, Eigen::Matrix<ART,-1,-1> &rho2, const vector<int> &statep, double coeff=1){
	//make matrix by looping over all states that aren't traced over
	unsigned int ti,tj;
	for(unsigned int i=0;i<evec.size();i++){
		for(unsigned int j=i;j<evec.size();j++){
			if( ( statep[i] & trunc_part) == ( statep[j] & trunc_part) ){
				for(ti=0;ti<trunc_states.size();ti++)
					if( ( statep[i] & ~trunc_part) == trunc_states[ti]) break;
				for(tj=0;tj<trunc_states.size();tj++)
					if( ( statep[j] & ~trunc_part) == trunc_states[tj]) break;
				if(ti==trunc_states.size() || tj==trunc_states.size()){
					cout<<ti<<" "<<tj<<endl;
					cout<<(bitset<16>)statep[i]<<" "<<(bitset<16>)statep[j]<<" "<<(bitset<16>)trunc_part<<endl;
				}
				rho2(ti,tj)+=coeff*safe_mult(evec[i],evec[j]);
				if(ti!=tj) rho2(tj,ti)+=coeff*safe_mult(evec[j],evec[i]);
			}
		}
	}
}
template<class ART>
double Wavefunction<ART>::ee_eval_rho(Eigen::Matrix<ART,-1,-1> &rho2){
	double out=0;
	//diagonalize matrix
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART,-1,-1> > rs(rho2);
	//output sum
	for(unsigned int i=0;i<trunc_states.size();i++){
		if(rs.eigenvalues()(i)>0){
			out-=rs.eigenvalues()(i)*log(rs.eigenvalues()(i));
			//cout<<log(rs.eigenvalues()(i))<<endl;
		}
	}
	return out;
}
template< >
inline double Wavefunction<double>::safe_mult(double x,double y){return x*y;}
template< >
inline complex<double> Wavefunction<complex<double> >::safe_mult(complex<double> x, complex<double> y){ return x*conj(y); }
#endif
