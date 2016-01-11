/*
this can do all kinds of things to a matrix, all it needs is a matvec that inherits from it
*/

#ifndef MATPROD_H
#define MATPROD_H

//#include "lapacke.h"
#include <iostream>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
using namespace std;

#include "arscomp.h"
#include "arssym.h"


template<class ART>
class MatrixWithProduct {

 private:

	int m, n; // Number of rows and columns.
	double E1,E2;
	ART *dense;
	Eigen::SparseMatrix<ART> sparse;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix<ART> > sparseLU_solver;
	Eigen::Matrix<ART, Eigen::Dynamic, 1> sparseLU_out; //used to store the results of solving the linear system in MultInvSparse
	
	int *ipiv;

 public:

  int nrows() { return m; }//size of matrix, ncols is kept for backwards compatibility with ARPACK but shouldnt be used
  int ncols() { return n; }
	void setrows(int x){ 
		m=x; n=x;
	}
	
	vector<double> eigvals;
	vector< vector<ART> > eigvecs;
	vector<int> lowlevpos;
	double getE(int a){return eigvals[lowlevpos[a]];} //ARPACK sometimes returns eigenvalues in the wrong order, these functions correct that
	vector<ART> getEV(int a){return eigvecs[lowlevpos[a]];}
	vector<int> sort_indexes(const vector<double> &v); //sorts the output of ARPACK so that the above functions return things in the right order
	bool compy(int,int);

  virtual void MultMv(ART* v, ART* w); //original Matvec
//  virtual Eigen::Matrix<ART, Eigen::Dynamic, 1> MultEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1>);//used so we can communicate with my lanczos, which is based on eigen

  void MultM2v(ART* v, ART* w); //squares the original matvec

  void MultInvDense(ART*v, ART* w); //deprecated
  void denseLU(); //deprecated
  void denseSolve();//deprecated

	void EigenDenseEigs(); //dense diagonalization
	Eigen::Matrix<ART,-1,-1> EigenDense;
	void makeDense(); //given a matvec, this computes a dense Eigen matrix for use with dense solvers
	void printDense();//prints the dense matrix

  void MultInvSparse(ART*v, ART* w); //uses precomputed LU decomposition to solve a system when doing shift-invert
  void makeSparse(double E); //makes a sparse matrix from a matvec and LU decomposes
  void SparseFromDense(double E); //makes a sparse matrix from a dense matrix and LU decomposes
  void sparseSolve(); //for testing purposes, dense solves sparse matrices

  int eigenvalues(int k, double E=-100); //computes eigenvalues and eigenvectors using ARPACK, E should be -100 if you want the ground state. If E is not -100, does shift-invert around E
  double find_middle();
  double single_energy(string whichp); //finds only the highest or lowest state, used for estimating energy bounds before doing shift-invert to target a section of the spectrum
  double calcVarEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1> v); //calculates the variance of an eigenvector to see if polishing is needed
    
  ~MatrixWithProduct();
  
  MatrixWithProduct(int nrows, double _E1=0, double _E2=0, int ncols = 0 )
  // Constructor.
  {
    m = nrows;
    n = (ncols?ncols:nrows);
    E1=_E1; E2=_E2;
    dense=NULL;
    ipiv=NULL;
  } // Constructor.
  MatrixWithProduct()
  // Constructor.
  {
    m = 1;
    n = 1;
    E1=0.; E2=0.;
    dense=NULL;
    ipiv=NULL;
  } // Constructor.
  
	
}; // MatrixWithProduct

template<class ART>
void MatrixWithProduct<ART>::MultM2v(ART* v, ART* w)
{
	if(E1==0 && E2==0){
		cout<<"energies were never set so can't compute the square thingy"<<endl;
		exit(0);
	}
	double *w1=new double[this->ncols()];
	double *w2=new double[this->ncols()];
	MultMv(v,w1);
	for(int i=0;i<this->ncols();i++)
		w1[i]=w1[i]-E1*v[i];
	
	MultMv(w1,w);
	for(int i=0;i<this->ncols();i++)
		w[i]=-w[i]+E2*w1[i];
	delete [] w1;
	delete [] w2;

} //  MultM2v.

//turn the matvec into a dense matrix
template<class ART>
void MatrixWithProduct<ART>::makeDense(){
    EigenDense=Eigen::Matrix<ART,-1,-1>::Zero(n,n);
//	dense=new ART[n*n];
	ART *v=new ART[n];
	ART *w=new ART[n];
	for(int i=0;i<n;i++){
		for(int j=0; j<n; j++){
			if(i==j) v[j]=1;
			else v[j]=0;
			w[j]=0;
		}
		MultMv(v,w);
		for(int j=0; j<n; j++){
//			dense[i+j*n]=w[j];
			EigenDense(j,i)=w[j];
		}
	}
	delete [] v;
	delete [] w;
//	cout<<EigenDense<<endl;
}
template<class ART>
void MatrixWithProduct<ART>::EigenDenseEigs(){
	Eigen::SelfAdjointEigenSolver< Eigen::Matrix<ART,-1,-1> > es(EigenDense);

	//store eigenvalues
	double *temp=new double[n];	
	Eigen::Map< Eigen::Matrix<double,-1,1> >(temp,n,1)=es.eigenvalues();
	eigvals=vector<double>(temp,temp+n);

	//store eigenvectors
	ART *temp2=new ART[n*n];
	lowlevpos=vector<int>(n,0);
	Eigen::Map< Eigen::Matrix<ART,-1,-1> >(temp2,n,n)=es.eigenvectors();
	eigvecs.clear();
	for(int i=0;i<n;i++){
		eigvecs.push_back(vector<ART>(temp2+i*n,temp2+(i+1)*n));
		lowlevpos[i]=i;
	}

	delete [] temp;
	delete [] temp2;
}
template<class ART>
void MatrixWithProduct<ART>::printDense(){
	cout<<EigenDense<<endl;
//	for(int i=0;i<n;i++){
//		for(int j=0;j<n;j++) cout<<dense[i+j*n]<<" ";
//		cout<<endl;
//	}
}
/*
//compute the dense LU decomposition
template<class ART>
void MatrixWithProduct<ART>::denseLU(){
	ipiv=new int[n];
	LAPACKE_dsytrf(LAPACK_COL_MAJOR,'U',n,dense,n,ipiv);
}

//compute solution to linear system (i.e. multiply by the inverse)
template<class ART>
void MatrixWithProduct<ART>::MultInvDense(ART *v, ART *w){
	for(int i=0;i<n;i++) w[i]=v[i];
	LAPACKE_dsytrs(LAPACK_COL_MAJOR,'U',n,1,dense,n,ipiv,w,n);
}
//calls lapack on the dense matrix to get the solution
//since we can't call lapack with a template this needs to be changed depending on what ART is
template<class ART>
void MatrixWithProduct<ART>::denseSolve(){
	ART *w=new double[n];
	LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,dense,n,w); 

//	for(int i=0;i<n;i++){
//		if ( abs(w[i])<1)	cout<<w[i]<<endl;
//	}

	for(int i=0;i<n;i++) cout<<w[i]<<endl;
	delete [] w;
}
*/
template<class ART>
void MatrixWithProduct<ART>::makeSparse(double E){
	sparse.resize(n,n);
	vector<Eigen::Triplet<ART> > coeff;
	Eigen::Triplet<ART> temp;
	dense=new ART[n*n];
	ART *v=new ART[n];
	ART *w=new ART[n];
	for(int i=0;i<n;i++){
		for(int j=0; j<n; j++){
			if(i==j) v[j]=1;
			else v[j]=0;
			w[j]=0;
		}
		MultMv(v,w);
		for(int j=0;j<n;j++) w[j]-=v[j]*E;
		for(int j=0; j<n; j++){
		//	if (w[j]!=0) temp=Eigen::Triplet<ART>(j,i,w[j]);
			if (abs(w[j])>1e-16) coeff.push_back( Eigen::Triplet<ART>(j,i,w[j]) );
		}
	}
	sparse.setFromTriplets(coeff.begin(), coeff.end() );
	delete [] v;
	delete [] w;
	sparseLU_solver.compute(sparse);
	if(sparseLU_solver.info()!=0) {
	  // decomposition failed
	  cout<<"decomposition failed! "<<sparseLU_solver.info()<<endl;
	}	
}

template<class ART>
void MatrixWithProduct<ART>::SparseFromDense(double E){
	sparse=EigenDense.sparseView();
	for(int i=0;i<n;i++) sparse.coeffRef(i,i)-=E;
	sparseLU_solver.compute(sparse);
	if(sparseLU_solver.info()!=0) {
	  // decomposition failed
	  cout<<"decomposition failed! "<<sparseLU_solver.info()<<endl;
	}	
}			
			
template<class ART>
void MatrixWithProduct<ART>::MultInvSparse(ART *v, ART *w){

	Eigen::Map <Eigen::Matrix<ART, Eigen::Dynamic, 1> > mapped_v(v,n);
	sparseLU_out=sparseLU_solver.solve(mapped_v);
	Eigen::Map <Eigen::Matrix<ART, -1, 1> > (w,n,1)=sparseLU_out; //using just out.data() fails for an unknown reason
}

template<class ART>
void MatrixWithProduct<ART>::MultMv(ART *v, ART *w){

	Eigen::Map <Eigen::Matrix<ART, Eigen::Dynamic, 1> > mapped_v(v,n);
	sparseLU_out=EigenDense*mapped_v;
	Eigen::Map <Eigen::Matrix<ART, -1, 1> > (w,n,1)=sparseLU_out; //using just out.data() fails for an unknown reason
}


//this doesn't do what you think it does! It converts the sparse matrix to a dense one and solves the dense one
template <class ART>
void MatrixWithProduct<ART>::sparseSolve(){
	Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> dMat=Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic>(sparse);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> > es(dMat);
	cout<<es.eigenvalues()<<endl;
}

template <class ART>
double MatrixWithProduct<ART>::calcVarEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1> evec){
	Eigen::Matrix<ART, Eigen::Dynamic, 1> w(n);
	double norm,eval;

	norm=evec.norm();
	evec/=norm;
	w=EigenDense*evec;
	eval=w.dot(evec); //should I also return this eigenvalue?
//	cout<<eval<<endl;
	w=w-eval*evec;
	return w.norm();
}

template<>
double MatrixWithProduct< complex<double> >::find_middle(){
	ARCompStdEig<double, MatrixWithProduct< complex<double> > >  dprob(ncols(), 5, this, &MatrixWithProduct< complex<double> >::MultMv,"LR",(int)0, 1e-4,1e6);
	dprob.FindEigenvalues();
	double upper=dprob.Eigenvalue(0).real();
	ARCompStdEig<double, MatrixWithProduct< complex<double> > >  dprob2(ncols(), 5, this, &MatrixWithProduct< complex<double> >::MultMv,"SR",(int)0, 1e-4,1e6);
	dprob2.FindEigenvalues();
	double lower=dprob.Eigenvalue(0).real();
	return 0.5*(upper+lower);
}

template<>
double MatrixWithProduct< double >::find_middle(){
	ARSymStdEig<double, MatrixWithProduct< double > >  dprob(ncols(), 5, this, &MatrixWithProduct<double>::MultMv,"LM",(int)0, 1e-4,1e6);
	dprob.FindEigenvalues();
	double upper=dprob.Eigenvalue(0);
	ARSymStdEig<double, MatrixWithProduct< double > >  dprob2(ncols(), 5, this, &MatrixWithProduct< double >::MultMv,"SM",(int)0, 1e-4,1e6);
	dprob2.FindEigenvalues();
	double lower=dprob.Eigenvalue(0);
	return 0.5*(upper+lower);
}
	
//template<class ART> //this generic template only serves to set the default value of E
//int MatrixWithProduct< ART >::eigenvalues(int stop, double E=-100){
//	return 0;
//}
template<>
inline int MatrixWithProduct< complex<double> >::eigenvalues(int stop, double E){
	vector< complex<double> >temp(n,0);
	int Nconverged;
	if (E==-100){
		ARCompStdEig<double, MatrixWithProduct< complex<double> > >  dprob(ncols(), stop, this, &MatrixWithProduct< complex<double> >::MultMv,"SR",(int)0, 1e-10,1e6);
		dprob.FindEigenvectors();
		Nconverged=dprob.ConvergedEigenvalues();

		eigvals=vector<double>(Nconverged,0);
		eigvecs=vector<vector< complex<double> > >(dprob.ConvergedEigenvalues(),temp);
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=dprob.Eigenvalue(k).real();
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
	}else{
		//SparseFromDense(E);
		makeSparse(E);
		ARCompStdEig<double, MatrixWithProduct< complex<double> > >  dprob(ncols(), stop, this, &MatrixWithProduct< complex<double> >::MultInvSparse,"LM");
		dprob.FindEigenvectors();
		
		eigvals=vector<double>(dprob.ConvergedEigenvalues(),0);
		eigvecs=vector<vector< complex<double> > >(dprob.ConvergedEigenvalues(),temp);
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=1./dprob.Eigenvalue(k).real()+E;
			//eigvals[k]=dprob.Eigenvalue(k).real();
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
		Nconverged=dprob.ConvergedEigenvalues();
	}
	if (Nconverged!=stop) cout<<"didn't get as many eigenvalues as expected! "<<Nconverged<<endl;
	lowlevpos=sort_indexes(eigvals);	//get the sorted indices, then put eigenvalues in the right order
	vector<double> temp_eigvals(stop,0);
	for(int i=0;i<stop;i++) temp_eigvals[i]=eigvals[lowlevpos[i]];
	eigvals=temp_eigvals;
	
	return Nconverged;
//	for(int i=0;i<Nconverged;i++) cout<<dprob.Eigenvalue(i)<<endl;
	
}
template<>
inline int MatrixWithProduct< double >::eigenvalues(int stop, double E){
	vector<double>temp(n,0);
	int Nconverged;
	if (E==-100){
		ARSymStdEig<double, MatrixWithProduct<double> >  dprob(ncols(), stop, this, &MatrixWithProduct<double>::MultMv,"SA",(int)0, 1e-10,1e6);
		dprob.FindEigenvectors();
		Nconverged=dprob.ConvergedEigenvalues();

		eigvals=vector<double>(Nconverged,0);
		eigvecs=vector<vector<double> >(dprob.ConvergedEigenvalues(),temp);
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=dprob.Eigenvalue(k);
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
	}else{
		makeSparse(E);
		//SparseFromDense(E);
		ARSymStdEig<double, MatrixWithProduct<double> >  dprob(ncols(), stop, this, &MatrixWithProduct<double>::MultInvSparse,"LM");
		dprob.FindEigenvectors();
		
		eigvals=vector<double>(dprob.ConvergedEigenvalues(),0);
		eigvecs=vector<vector<double> >(dprob.ConvergedEigenvalues(),temp);
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=1./dprob.Eigenvalue(k)+E;
//			eigvals[k]=dprob.Eigenvalue(k);
//			cout<<eigvals[k]<<endl;
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
		Nconverged=dprob.ConvergedEigenvalues();
	}
	lowlevpos=sort_indexes(eigvals);	
	return Nconverged;
//	for(int i=0;i<Nconverged;i++) cout<<dprob.Eigenvalue(i)<<endl;
	
}

template <class ART>
vector<int> MatrixWithProduct<ART>::sort_indexes(const vector<double> &v) {

  // initialize original index locations
  vector<int> idx(v.size());
  for (int i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
	//I can't use std::sort because c++ is super gay
	int temp;
	for(int j=idx.size();j>0;j--){
		for(int i=0;i<j-1;i++){
			if(v[idx[i]]>v[idx[i+1]]){
				 temp=idx[i];
				 idx[i]=idx[i+1];
				 idx[i+1]=temp;
			}
		}
	}
  return idx;
}	
//destructor, delete the dense matrices
template<class ART>
MatrixWithProduct<ART>::~MatrixWithProduct(){
	delete [] dense;
	delete [] ipiv;
}
#endif // MATPROD_H

