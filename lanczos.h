#ifndef LANCZOS_H
#define LANCZOS_H

#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <iomanip>

typedef Eigen::VectorXd Vec;
using namespace std;

template<class ARFLOAT, class ARFOP>
class Lanczos{

public:
	typedef Vec (ARFOP::* TypeOPx)(Vec);

	Lanczos(ARFOP* objOPp, Vec (ARFOP::* MultOPxp)(Vec) );
	void find_lowest();
	void polish(double gs, Vec start);
	void polishtest(double gs, Vec start);

private:

	ARFOP   *objOP;     // Object that has MultOPx as a member function.
	TypeOPx MultOPx;    // Function that evaluates the product OP*x.
	vector<double> t0,t1;
	vector<double> es;
	
};

template<class ARFLOAT, class ARFOP>
Lanczos<ARFLOAT, ARFOP>::Lanczos(ARFOP* objOPp, Vec (ARFOP::* MultOPxp)(Vec) ){
	objOP=objOPp;
	MultOPx=MultOPxp;
}
template<class ARFLOAT, class ARFOP>
void Lanczos<ARFLOAT, ARFOP>::find_lowest(){
	int maxstep=100;
	Vec start=Eigen::VectorXd::Random(objOP->ncols());
	double beta=start.norm();
	start/=beta;
	Vec y=start;
	Eigen::MatrixXd T(maxstep+1,maxstep+1);
	vector<double> energies;
	double alpha,diff;
	Vec last=y;
	Vec secondlast;
	int print=0;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
	
	for(int k=0;k>-1;k++){
		if (k>maxstep) break;
		if (diff<0) break;//need more stopping conditions
		 
		if (k>0) last/=beta;
		//project out converged eigenvectors
		y=(objOP->*MultOPx)(last);
		//project out again
		
		//do an iteration
		alpha=last.dot(y);
		if (k>0) y=y-beta*secondlast;
		y=y-alpha*last;
		beta=y.norm();
		T(k,k)=alpha;
		T(k,k+1)=beta;
		T(k+1,k)=beta;
		
		secondlast=last;
		last=y;
		
		//analyze t matrix
		if (k==0) energies.push_back(alpha);
		else{
			es.compute(T.topLeftCorner(k+1,k+1));
			energies.push_back(es.eigenvalues()(0));
			diff=energies[k-1]-energies[k];
			if(print){
				cout<<"---------------------------"<<endl;
				cout<<T.topLeftCorner(k+1,k+1)<<endl;
				cout<<es.eigenvalues()<<endl;
			}
		}	
	}
	//compute eigenvectors
	Vec outvec=es.eigenvectors().col(0)(0)*start;
	last=start;
	for(int k=0;k<energies.size()-1;k++){
		//project out converged eigenvectors
		y=(objOP->*MultOPx)(last);
		//project out again
		
		//do an iteration
		alpha=T(k,k);
		if (k>0) y=y-beta*secondlast;
		y=y-alpha*last;
		beta=T(k,k+1);
		y/=beta;

		outvec=outvec+es.eigenvectors().col(0)(k+1)*y;		
		secondlast=last;
		last=y;
	}
	cout<<objOP->calcVarEigen(outvec)<<endl;		
	polishtest(energies[energies.size()-1],outvec);		
	//for(int i=0;i<objOP->ncols();i++) cout<<setprecision(15)<<outvec(i)<<endl;	
	//for(int i=0;i<energies.size();i++) cout<<setprecision(15)<<energies[i]<<endl;
}

//duncan's algorithm for polishing the lowest eigenvalue of a matrix
template<class ARFLOAT, class ARFOP>
void Lanczos<ARFLOAT, ARFOP>::polish(double e, Vec start){
	int maxstep=5, print=1;
	double beta=start.norm();
	start/=beta;
	Vec y=start;
	Eigen::MatrixXd T=Eigen::MatrixXd::Zero(maxstep+1,maxstep+1);
	Eigen::MatrixXd T2(maxstep+1,maxstep+1);
	Eigen::MatrixXd TP1=Eigen::MatrixXd::Zero(maxstep+1,maxstep+1);
	//P1(0,0)=1;
	
	double alpha,diff;
	Vec last=y, secondlast;
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> gs;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
	TP1=Eigen::MatrixXd::Zero(maxstep+1,maxstep+1);
	TP1(1,1)=1;
	
	//the first part is the same as ordinary lanczos, except we subtract off the energy each time, and we started with an eigenvector
	int k;
	for(k=0;k>-1;k++){
		if (k>=maxstep) break;
		if (diff<0) break;		//need more stopping conditions
		
		if (k>0) last/=beta;
		//project out converged eigenvectors
		y=(objOP->*MultOPx)(last);
		y=y-e*last;
		//project out again
		
		//do an iteration
		alpha=last.dot(y);
		if (k>0) y=y-beta*secondlast;
		y=y-alpha*last;
		beta=y.norm();
		T(k,k)=alpha;
		T(k,k+1)=beta;
		T(k+1,k)=beta;
		
		secondlast=last;
		last=y;
		
		//analyze t matrix
		if (k==0) continue;
		else{
			T2=T*T;
			if (k>=1) TP1(1,1)=pow(T(1,1),2);
			if (k>=2){
				TP1(1,2)=T(1,2)*T(1,1);
				TP1(2,1)=T(1,2)*T(1,1);
				TP1(2,2)=pow(T(1,2),2);
			}
			gs.compute(TP1.topLeftCorner(k+1,k+1),T2.topLeftCorner(k+1,k+1));
			es.compute(T2.topLeftCorner(k+1,k+1)-TP1.topLeftCorner(k+1,k+1));
			if(gs.info()!=Eigen::Success){
				cout<<"no convergence for generalized eigenproblem"<<endl;
				exit(0);
			}
			if(print){
				cout<<"*******************************"<<endl;
				cout<<T2.topLeftCorner(k+1,k+1)<<endl;
//				cout<<TP1.topLeftCorner(k+1,k+1)<<endl;
				cout<<"---------------------------"<<endl;
				cout<<gs.eigenvalues()<<endl;
				cout<<es.eigenvalues()<<endl;
				cout<<"---------------------------"<<endl;
				cout<<gs.eigenvectors()<<endl;
				cout<<es.eigenvectors()<<endl;
			}
		}
	}
	int N=k;
	//compute eigenvectors
//	cout<<gs.eigenvectors().col(0)(0)*T(0,1)<<endl;
	Vec outvec=Vec::Zero(objOP->ncols()); 
//	Vec outvec=gs.eigenvectors().col(0)(0)*start;
	last=start;
	for(int k=0;k<N-1;k++){
		//project out converged eigenvectors
		y=(objOP->*MultOPx)(last);
		y=y-e*last;
		//project out again
		
		//do an iteration
		alpha=T(k,k);
		if (k>0) y=y-beta*secondlast;
		y=y-alpha*last;
		beta=T(k,k+1);
		y/=beta;

		outvec=outvec+es.eigenvectors().col(0)(k+1)*y;		
		secondlast=last;
		last=y;
//		cout<<objOP->calcVarEigen(start-T(0,1)*outvec)<<endl;				
	}
	//orthonormalize outvec
	double s=outvec.dot(start);
	double n=outvec.norm();
	cout<<"s="<<s<<", norm="<<n<<endl;
	outvec=outvec-s*start;
	outvec/=n;
	cout<<T(0,1)<<" "<<start.norm()<<" "<<objOP->calcVarEigen(start)<<endl;
	cout<<objOP->calcVarEigen(start-T(0,1)*outvec)<<endl;				
	
}
//tries to polish by dense solving the generalized eigenproblem
template<class ARFLOAT, class ARFOP>
void Lanczos<ARFLOAT, ARFOP>::polishtest(double e, Vec start){
	Eigen::MatrixXd H=Eigen::MatrixXd::Zero(objOP->ncols(),objOP->ncols());
	Eigen::MatrixXd HP(objOP->ncols(),objOP->ncols());
	Eigen::VectorXd v,w;
	for(int i=0;i<objOP->ncols();i++){
		v=Eigen::VectorXd::Zero(objOP->ncols());
		v(i)=1;
		w=(objOP->*MultOPx)(v);		
		for(int j=0;j<objOP->ncols();j++){
			if(i==j) H(i,j)=w(j)-e;
			else H(i,j)=w(j);	
			HP(i,j)=start(i)*start(j);
		}
	}
	cout<<"starting test"<<endl;
	cout<<"e="<<e<<endl;
	cout<<"start"<<endl<<start<<endl;
	cout<<H<<endl;
	cout<<HP<<endl;
	Eigen::MatrixXd A=H*HP*H;
	cout<<A<<endl;
	Eigen::VectorXd b=Eigen::VectorXd::Zero(objOP->ncols());
	Eigen::VectorXd x=A.ldlt().solve(b);
	cout<<x<<endl;
}
#endif
