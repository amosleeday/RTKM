#pragma once
#include "Afxtempl.h"
using namespace std;
#include <vector>
class CMatrix
{
public:
	CMatrix(void);
	~CMatrix(void);
	
public:
	void built_mat( std::vector< std::vector<double> > &A, int row,int col);
	//construct a vector
	void built_vec(std::vector<double> &A,int row);
	bool eye(std::vector< std::vector<double> > & A,int M,int N);
	bool zeros(std::vector< std::vector<double> > & A,int M,int N);


	std::vector<int> find(std::vector<int> A,int m,int x);
	bool plus(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> >&  C,int row,int col);
	bool plus(std::vector<double> A,std::vector<double> B,std::vector<double> &C,int row);
	bool plus(std::vector< std::vector<double> > & A,std::vector< std::vector<double> >B,int row,int col);//A=A+B
	bool plusT(std::vector< std::vector<double> > & A,std::vector< std::vector<double> >B,int row,int col);//A=A+B'
	bool subtr(std::vector<double> A,std::vector<double> B,std::vector<double> &C,int row);
	bool subtr(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int row,int col);
	bool subtr(std::vector< std::vector<double> >&A,std::vector< std::vector<double> >  B,int row,int col);
	bool multip(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int M,int N,int L);
	bool multip(std::vector<double> A,std::vector<double>  B,std::vector< std::vector<double> >  &C,int M,int L);
	bool multip(std::vector< std::vector<double> >  A, std::vector<double>  B, std::vector<double>  & C,int M,int N);
	bool multip(std::vector< std::vector<double> >  A, double s, std::vector< std::vector<double> > & C,int M,int N);
	bool multip(std::vector< std::vector<double> >  &A, double s,int M,int N);

	bool multipT(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  B,std::vector< std::vector<double> > & C,int M,int N,int L); // C=A*B'

	
	double inner_prod(std::vector<double> A,int M);
	double inner_prod( std::vector< std::vector<double> > A, std::vector< std::vector<double> > B, int i, int j, int N ); //dot( A(:,i),B(:, j)
	double norm(std::vector<double> A,int M);
	bool MatrixInv(std::vector< std::vector<double> > &C, int N);
	bool MatrixInvDong( std::vector<std::vector<double> >& a,int n);
	bool MatrixInvDong( std::vector<std::vector<double> > a, std::vector<std::vector<double> >& b,int n);
	bool MatrixSovleDong( std::vector<std::vector<double> > a, std::vector<std::vector<double> > b,std::vector<std::vector<double> > &c, int n, int m);
	bool InvSymMat(std::vector< std::vector<double> > A,std::vector< std::vector<double> > &C,int N);
	bool InvGuass(std::vector< std::vector<double> > A,std::vector< std::vector<double> > &B,int M);
	bool InvGenMat(std::vector< std::vector<double> > A,std::vector< std::vector<double> > &B,int M);
	bool MatT(std::vector< std::vector<double> > A,std::vector< std::vector<double> >& B,int M,int N);
	bool form_ATPA(std::vector< std::vector<double> > A, std::vector< std::vector<double> > P,std::vector< std::vector<double> > &ATPA,int M,int N);
	bool form_APAT(std::vector< std::vector<double> > A, std::vector< std::vector<double> > P,std::vector< std::vector<double> > &ATPA,int M,int N);
	bool form_ATPA(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  &ATPA,int M,int N);
	bool form_ATPL(std::vector< std::vector<double> > A,std::vector< std::vector<double> > P,std::vector<double> L,std::vector<double> &ATPL,int M,int N);
	bool form_ATPL(std::vector< std::vector<double> > A,std::vector<double> L,std::vector<double> &ATPL,int m,int n);
	bool form_ATPAL(std::vector< std::vector<double> > A,std::vector< std::vector<double> > P,std::vector<double> L,std::vector< std::vector<double> > &ATPA,std::vector<double> & ATPL,int M,int NN);
	bool form_ATPAL(std::vector< std::vector<double> > A,std::vector<double> L,std::vector< std::vector<double> >& ATPA, std::vector<double> &ATPL,int M,int N);
	bool lPl(std::vector< std::vector<double> > P,std::vector<double> L, double &lpl);
	bool ATPAL_DD(std::vector< std::vector<double> > A,std::vector<double>L,std::vector< std::vector<double> > &ATPA,std::vector<double> &ATPL,int M,int N,int ISL1);
	bool ATPAL_2DInter(std::vector< std::vector<double> > A,std::vector<double> L,std::vector< std::vector<double> > &ATPA,std::vector<double> &ATPL,int M,int N,int ii,int ISL1);
	bool ATPAL_2DInter_Del(std::vector< std::vector<double> > A,std::vector<double> L,std::vector< std::vector<double> > &ATPA,std::vector<double> &ATPL,int M,int N,int Delete,int ii,int ISL1);
	bool Adjustment(std::vector< std::vector<double> > NN, std::vector<double > W, double lpl, int ObsNum, std::vector< double>& X, std::vector<std::vector< double> > &Dxx,double &segma0 );
	bool calc_resi(std::vector< std::vector<double> > A,std::vector<double> X,std::vector<double> L,std::vector<double> &V,int M,int N);
	double calc_vtpv( std::vector<double> V,std::vector< std::vector<double> > P,int M);
	bool del_row(std::vector< std::vector<double> > &A,int M,int N,int irow);
	bool del_elem(std::vector<double> &A,int M,int irow);
	bool del_col(std::vector< std::vector<double> > &A,int M,int N,int icol);
	bool del_row_col(std::vector< std::vector<double> > &A, int M,int icol);
	bool add_row_col(std::vector< std::vector<double> > &A, int M,int icol);
	bool add_col(std::vector< std::vector<double> > &A, int M,int N, int icol);
	bool add_row(std::vector< std::vector<double> > &A, int M,int N, int irow);
	double Mutipvec(std::vector<double> A,std::vector<double> B,int M);

	bool FprintMatrix(std::vector< std::vector<double> > A, CString output);
	bool FprintMatrix(std::vector<double>  A, CString output);
	bool FprintMatrix(std::vector< std::vector<double> > A, fstream& fout);
	bool FprintMatrix(std::vector<double>  A, fstream& fout);

	bool Parition0(std::vector< std::vector<double> >  &NN, std::vector< std::vector<double> >  &UU,std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2);
	bool Partion00(std::vector< std::vector<double> >  &NN, std::vector< std::vector<double> > &UU,std::vector< std::vector<double> > N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  U1,std::vector< std::vector<double> > U2);	
	bool Parition0(std::vector< std::vector<double> >  &Q11, std::vector< std::vector<double> >  &Q12,std::vector< std::vector<double> >  &Q22,std::vector< std::vector<double> >  &X1,std::vector< std::vector<double> >  &X2,std::vector< std::vector<double> >  N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  U1,std::vector< std::vector<double> > U2);
	bool Parition1(std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um);
	bool Parition11(std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um);
	
	bool Parition1(std::vector< std::vector<double> >  N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  U1,std::vector< std::vector<double> >  U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um,std::vector< std::vector<double> > &Gain0, std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y);
	// N=[  N11   |   N12 ]   [ x ]  = [ U1 ]         ( N22-N21*N11^(-1)*N12 )*y=N22-N21*N11^(-1)*U1 
	//      [  N12   |   N22 ]   [ y ]  = [ U2 ]			x=N11^(-1)( U1-N12*y )
	bool Parition2(std::vector< std::vector<double> >  &Q11,std::vector< std::vector<double> >  &Q12,std::vector< std::vector<double> > &Q22, std::vector< std::vector<double> >  &X1,std::vector< std::vector<double> >  &X2,std::vector< std::vector<double> >& X1bar ,std::vector< std::vector<double> > & Qhat,std::vector< std::vector<double> > & Xhat);
	bool Parition3(std::vector< std::vector<double> >  &Q11,std::vector< std::vector<double> >  &Q12,std::vector< std::vector<double> > &Q22, std::vector< std::vector<double> >  &X1,std::vector< std::vector<double> >  &X2,std::vector< std::vector<double> >& X2bar ,std::vector< std::vector<double> > & Qhat,std::vector< std::vector<double> > & Xhat);

	bool PulsNormal(std::vector< std::vector<double> >  &N11,std::vector< std::vector<double> >  &N12,std::vector< std::vector<double> > &N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um);
	bool PulsNormal1(std::vector< std::vector<double> >  N11,std::vector< std::vector<double> >  N12,std::vector< std::vector<double> > N22, std::vector< std::vector<double> >  &U1,std::vector< std::vector<double> >  &U2, std::vector< std::vector<double> > & C,std::vector< std::vector<double> > & Um);

	bool LU(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  U,int row);
	bool LDR(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  L,std::vector<double> D,std::vector< std::vector<double> >  R,int row);
	bool LDLT(std::vector< std::vector<double> >  A,std::vector< std::vector<double> >  L,std::vector<double> D,int row);
		
	//bool LDLTFactorization(std::vector< std::vector<double> > Q,std::vector< std::vector<double> > & L,std::vector< std::vector<double> > &D,int n);
	bool LTDLFactorization(std::vector< std::vector<double> > Q,std::vector< std::vector<double> > & L,std::vector< std::vector<double> >  &D,int n);
	bool free_memo(std::vector< std::vector<double> >   &A,int M,int N);


	void form_RefTMat(std::vector< std::vector<double> > &T, int AmbNum, int ip_ref);// ip_ref: the postiton where the last reference satellite located in the new prn list 
	void buble_sort(std::vector<int> list, std::vector< std::vector<double> >& T );
	bool cond(std::vector< std::vector<double> > Q, int N,double &num);
	double Regu_par(std::vector< std::vector<double> > N11, std::vector< std::vector<double> > N12,std::vector< std::vector<double> > N22);
	double norm(std::vector< std::vector<double> > Q);
	double trace(std::vector< std::vector< double > > Q);
};
