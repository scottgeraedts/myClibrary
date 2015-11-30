//library of functions that can be applied to any wavefunction stores as a collection of spins
//right now the only thing here is entanglement entropy calculations, but more could be added
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include<Eigen/Dense>
#include<iostream>
#include<bitset>

using namespace std;

template <class ART>
class Wavefunction{
private:
	int nBits; //number of physical bits (e.g. spins or flux quanta)
	Eigen::Matrix<ART,-1,-1> rho;

public:
	void init_wavefunction(int);

	//entanglement entropy stuff
	void ee_setup(int trunc_start, int trunc_end, const vector<int> &statep);
	void ee_compute_rho(const vector<ART> &evec, Eigen::Matrix<ART,-1,-1> &rho2,  const vector<int> &statep, double coeff);
	double ee_eval_rho(Eigen::Matrix<ART,-1,-1> &rho2);
	double entanglement_entropy(const vector<ART> &evec, const vector <int> &statep);
	vector<int> trunc_states;
	int trunc_part;
};

template<class ART>
void Wavefunction<ART>::init_wavefunction(int n){ nBits=n; }

template<class ART>
double Wavefunction<ART>::entanglement_entropy(const vector<ART> &evec, const vector<int> &statep){
	double out=0;
	for(int i=0;i<nBits;i++){
		ee_setup(i,(i+nBits/2)%nBits,statep);
		rho=Eigen::Matrix<ART,-1,-1>::Zero(trunc_states.size(), trunc_states.size() );
		ee_compute_rho(evec,rho,statep,1.);
		out+=ee_eval_rho(rho)/(1.*nBits);
	}
	return out;
}
template<class ART>
void Wavefunction<ART>::ee_setup(int trace_start, int trace_end, const vector<int> &statep){
	bool found;

	//compute trunc_part
	trunc_part=0; //bitwise-AND with this gives just the component in the traced over states
	for(int i=0;i<nBits;i++){
		if ((i>=trace_start && i <trace_end && trace_end>trace_start) || (trace_end<trace_start && (i < trace_end || i >= trace_start) ) )
			trunc_part=trunc_part | 1<<i;
	}

	//figure out how many reduced states there are		
	trunc_states.clear();
	for(int i=0;i<statep.size();i++){
		found=false;
		for(int j=0;j<trunc_states.size();j++)
			if( ( statep[i] & ~trunc_part) ==trunc_states[j]) found=true;
		if(!found) trunc_states.push_back( statep[i] & ~trunc_part);
	}
}

template<class ART>
void Wavefunction<ART>::ee_compute_rho(const vector<ART> &evec, Eigen::Matrix<ART,-1,-1> &rho2, const vector<int> &statep, double coeff=1){
	//make matrix by looping over all states that aren't traced over
	int ti,tj;
	for(int i=0;i<evec.size();i++){
		for(int j=0;j<evec.size();j++){
			if( ( statep[i] & trunc_part) == ( statep[j] & trunc_part) ){
				for(ti=0;ti<trunc_states.size();ti++)
					if( ( statep[i] & ~trunc_part) == trunc_states[ti]) break;
				for(tj=0;tj<trunc_states.size();tj++)
					if( ( statep[j] & ~trunc_part) == trunc_states[tj]) break;
				if(ti==trunc_states.size() || tj==trunc_states.size()){
					cout<<ti<<" "<<tj<<endl;
					cout<<(bitset<16>)statep[i]<<" "<<(bitset<16>)statep[j]<<" "<<(bitset<16>)trunc_part<<endl;
				}
				rho2(ti,tj)+=coeff*(evec[i]*conj(evec[j]));
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
	for(int i=0;i<trunc_states.size();i++) 
		if(rs.eigenvalues()(i)>0) out-=rs.eigenvalues()(i)*log(rs.eigenvalues()(i));
	return out;
}
	
#endif
