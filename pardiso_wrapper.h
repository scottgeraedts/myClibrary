/*
 * diagonlizes with the shift-invert methd using superLU and arpack
 * it can also do regular diagonalization
*/

#ifndef PARDISOWRAPPER_H
#define PARDISOWRAPPER_H
#define SUPERLU_INC_

#include <iostream>
#include <algorithm>
#include "utils.h"
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//#define USE_COMPLEX

#include <Eigen/PardisoSupport>
using namespace std;

#include "arscomp.h"
#include "arssym.h"


extern"C"{
	//matrix multiplication
	//void zgemm_(char *transa, char *transb, int *rows, int *cols, int *k, double *alpha, complex<double> *a, int *lda, complex<double> *b, int *ldb, double *beta, complex<double> *c, int *ldc);
	//sparse matrix-vector
	void mkl_cspblas_zcsrsymv_(char *uplo, int *rows, complex<double> *vals, int *ia, int *ja, complex<double> *v, complex<double> *w);
	void mkl_cspblas_zcsrgemv_(char *uplo, int *rows, complex<double> *vals, int *ia, int *ja, complex<double> *v, complex<double> *w);


	//sparse matrix-vector
	void mkl_cspblas_dcsrsymv_(char *uplo, int *rows, double *vals, int *ia, int *ja, double *v, double *w);
	void mkl_cspblas_dcsrgemv_(char *uplo, int *rows, double *vals, int *ia, int *ja, double *v, double *w);

}

template<class ART>
class Pardiso_Wrapper {

 private:

	int n; // Number of rows and columns.
	int nonzero;

	Eigen::Matrix<ART,-1,1> sparseLU_out;	
	ART *vals,*raw_evals;
	int *ia, *ja;

	Eigen::PardisoLDLT< Eigen::SparseMatrix<ART> > sparse_solver;
 //used to 'undo' shifts on the diagonal, so we can use the same matrix for all shift-inverts
	double oldE;

	bool first_time;

 public:

	int verbose;
  int ncols() { return n; }
	void setrows(int x){ 
		n=x;
	}
	
	vector<double> eigvals;
	vector< vector<ART> > eigvecs;
//	vector<int> lowlevpos;
//	double getE(int a){return eigvals[lowlevpos[a]];} //ARPACK sometimes returns eigenvalues in the wrong order, these functions correct that

	Eigen::SparseMatrix<ART> EigenSparse;
	void CSR_from_Sparse(Eigen::SparseMatrix< ART > &sparse);
	void ShiftInvert(double E);
	void print();

  void MultMv(ART *v, ART *w); //original Matvec
  void MultInv(ART *v, ART *w); //uses precomputed LU decomposition to solve a system when doing shift-invert

	Eigen::Matrix<ART,-1,-1> EigenDense;

  int eigenvalues(int k, double E=-100); //computes eigenvalues and eigenvectors using ARPACK, E should be -100 if you want the ground state. If E is not -100, does shift-invert around E
  double single_energy(string whichp); //finds only the highest or lowest state, used for estimating energy bounds before doing shift-invert to target a section of the spectrum

  ~Pardiso_Wrapper();
  Pardiso_Wrapper(int nrows )
  // Constructor.
  {
    n = nrows;
	verbose=0;
	oldE=0.;
	first_time=true;
	if(verbose>0) cout<<"hitting the right constructor"<<endl;
  } // Constructor.
  Pardiso_Wrapper()
  // Constructor.
  {
	cout<<"hitting the wrong constructor"<<endl;
    n = 1;
  } // Constructor.
  
	
}; // MatrixWithProduct

template<class ART>
void Pardiso_Wrapper<ART>::CSR_from_Sparse(Eigen::SparseMatrix<ART > &sparse){
	EigenSparse=sparse;
	nonzero=EigenSparse.nonZeros();
	if(verbose>0) cout<<"nonzero elements "<<nonzero<<endl;
	vals=EigenSparse.valuePtr();
	ja=EigenSparse.innerIndexPtr();
	ia=EigenSparse.outerIndexPtr();
}	

template<class ART>
void Pardiso_Wrapper<ART>::ShiftInvert(double E){
//subtract a constant from the diagonal elements
	if(first_time){
		sparse_solver.analyzePattern(EigenSparse);
		first_time=false;
	}
	for(int i=0;i<ncols();i++) EigenSparse.coeffRef(i,i)-=(E-oldE);
	sparse_solver.factorize(EigenSparse);
}

template<class ART>
void Pardiso_Wrapper<ART>::MultInv(ART *v, ART *w){
	Eigen::Map <Eigen::Matrix<ART, Eigen::Dynamic, 1> > mapped_v(v,n);
	sparseLU_out=sparse_solver.solve(mapped_v);
	Eigen::Map <Eigen::Matrix<ART, -1, 1> > (w,n,1)=sparseLU_out; //using just out.data() fails for an unknown reason
}
template<class ART>
Pardiso_Wrapper<ART>::~Pardiso_Wrapper(){
	if(verbose>0) cout<<"deallocated successfully"<<endl;
}


#ifndef USE_COMPLEX
template<>
inline void Pardiso_Wrapper<double>::MultMv(double *v, double *w){
//	char uplo='u';
//	mkl_cspblas_zcsrsymv_(&uplo, &n, vals, ia, ja, v, w);
	char uplo='n';
	mkl_cspblas_dcsrgemv_(&uplo, &n, vals, ia, ja, v, w);
}



template<>
inline double Pardiso_Wrapper<double>::single_energy(string type){
	ARSymStdEig<double, Pardiso_Wrapper<double> >  dprob(ncols(), 5, this, &Pardiso_Wrapper<double>::MultMv,type,(int)0, 1e-10,1e6);
	dprob.FindEigenvalues();
	return dprob.Eigenvalue(4);

}

template<>
inline int Pardiso_Wrapper<double>::eigenvalues(int stop, double E){
	vector<double>temp(n,0);
	int Nconverged;
	if (E==-100){
		ARSymStdEig<double, Pardiso_Wrapper<double> >  dprob(ncols(), stop, this, &Pardiso_Wrapper<double>::MultMv,"SA",(int)0, 1e-14,1e6);
//		dprob.FindEigenvalues();
		dprob.FindEigenvectors();
		Nconverged=dprob.ConvergedEigenvalues();

		eigvals=vector<double>(Nconverged,0);
		eigvecs=vector<vector<double> >(dprob.ConvergedEigenvalues(),vector<double>(n,0));
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=dprob.Eigenvalue(k);
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
	}else{
		time_t walltime=time(NULL);
		clock_t CPUtime=clock();
		if(verbose>0) cout<<"making sparse"<<endl;
		ShiftInvert(E);
		walltime=time(NULL)-walltime;
		CPUtime=clock()-CPUtime;
		if(verbose>0) cout<<"the LU decomposition took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;

		walltime=time(NULL);
		CPUtime=clock();

		ARSymStdEig<double, Pardiso_Wrapper<double> > dprob(n,stop, this, &Pardiso_Wrapper<double>::MultInv,"LM");
		dprob.FindEigenvectors();
		walltime=time(NULL)-walltime;
		CPUtime=clock()-CPUtime;
		if(verbose>0) cout<<"the eigensolving took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;
		
		eigvals=vector<double>(stop,0);
		eigvecs=vector<vector<double> >(dprob.ConvergedEigenvalues(),vector<double> (n,0));
		for(int k=0;k<stop;k++){
			
			eigvals[k]=1./dprob.Eigenvalue(k)+E;
			//cout<<eigvals[k]<<endl;
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
		Nconverged=stop;
		sort(eigvals.begin(),eigvals.end());
	}
	//lowlevpos=sort_indexes(eigvals);
	return Nconverged;
//	for(int i=0;i<Nconverged;i++) cout<<dprob.Eigenvalue(i)<<endl;
	
}

#else
//does the analysis of the shape of the matrix, to speed up future solutions
template<>
inline void Pardiso_Wrapper<complex<double> >::MultMv(complex<double> *v, complex<double> *w){
//	char uplo='u';
//	mkl_cspblas_zcsrsymv_(&uplo, &n, vals, ia, ja, v, w);
	char uplo='n';
	mkl_cspblas_zcsrgemv_(&uplo, &n, vals, ia, ja, v, w);
}


template<>
inline double Pardiso_Wrapper<complex<double> >::single_energy(string type){
	ARCompStdEig<double, Pardiso_Wrapper<complex<double> > >  dprob(ncols(), 5, this, &Pardiso_Wrapper<complex<double> >::MultMv,type,(int)0, 1e-10,1e6);
	dprob.FindEigenvalues();
	return real(dprob.Eigenvalue(4));

}

template<>
inline int Pardiso_Wrapper<complex<double> >::eigenvalues(int stop, double E){
	vector<double>temp(n,0);
	int Nconverged;
	if (E==-100){
		if(verbose>0) cout<<"not doing shift invert"<<endl;
		ARCompStdEig<double, Pardiso_Wrapper<complex<double> > >  dprob(ncols(), stop, this, &Pardiso_Wrapper<complex<double> >::MultMv,"SR",(int)0, 1e-14,1e6);
//		dprob.FindEigenvalues();
		dprob.FindSchurVectors();
		Nconverged=dprob.ConvergedEigenvalues();

		eigvals=vector<double>(Nconverged,0);
		eigvecs=vector<vector< complex<double> > >(dprob.ConvergedEigenvalues(),vector<complex<double> >(n,0));
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=real(dprob.Eigenvalue(k));
			eigvecs[k]=*(dprob.StlSchurVector(k));
		}
	}else{
		time_t walltime=time(NULL);
		clock_t CPUtime=clock();
		if(verbose>0) cout<<"making sparse"<<endl;
		ShiftInvert(E);
		walltime=time(NULL)-walltime;
		CPUtime=clock()-CPUtime;
		if(verbose>0) cout<<"the LU decomposition took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;

		complex<double> shift(E,0);
		walltime=time(NULL);
		CPUtime=clock();

		ARCompStdEig<double, Pardiso_Wrapper<complex<double> > > dprob(n,stop, this, &Pardiso_Wrapper<complex<double> >::MultInv,"LM");
		dprob.FindSchurVectors();
		walltime=time(NULL)-walltime;
		CPUtime=clock()-CPUtime;
		if(verbose>0) cout<<"the eigensolving took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;
		
		eigvals=vector<double>(stop,0);
		eigvecs=vector<vector< complex<double> > >(dprob.ConvergedEigenvalues(),vector<complex<double> > (n,0));
		for(int k=0;k<stop;k++){
			
			eigvals[k]=1./real(dprob.Eigenvalue(k))+E;
			cout<<eigvals[k]<<endl;
			eigvecs[k]=*(dprob.StlSchurVector(k));
		}
		Nconverged=stop;
		sort(eigvals.begin(),eigvals.end());
	}
	//lowlevpos=sort_indexes(eigvals);
	return Nconverged;
//	for(int i=0;i<Nconverged;i++) cout<<dprob.Eigenvalue(i)<<endl;
	
}
#endif

#endif
