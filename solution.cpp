#pragma  once
#include "stdafx.h"//changed
#include<math.h>//calc
#include "solution.h"
#include "ExtFun.h"
#include <iomanip>
#include "Position.h"

/*********************************Normal Equation***************************************/
/*
 *Form Normal Equation 2, after the Reform Error equation
 *suppose that there are only position parameter and ambiguity parameter,here
 *
 *I:
 *	A1	the design matrix of position for code and phase	[(k*n) + (k'*n)]x 3		k=num	of phase type=typenumctrlphs  k'=num of code type=typenumctrlcode  n=num of single type obs 
 *	A2	the design matrix of ambiguity (k*n)*(k*n)
 *	Weight the weight matrix of obs here,  [(k*n) + (k'*n)] x [(k*n) + (k'*n)]
 *O:
 *	N11	the normal equ of pos
 *	N12
 *	N22 the normal equ of amb
 *	
 *	Note:
 *		consider the correlation of combination obs, Weight may not be diagonal 
 **/
void	 FormNormalEquationCoef2(math::matrix<double>	A1,math::matrix<double>A2,math::matrix<double>& N11,math::matrix<double>& N12,math::matrix<double>& N22,math::matrix<double>Weight)
{
	N11	=(~A1)*Weight*A1;
	N12	=(~A1)*Weight*A2;
	N22	=(~A2)*Weight*A2;
}

void  FormNormalEquationConst2(math::matrix<double>	A1,math::matrix<double>A2,math::matrix<double>Weight,math::matrix<double>L,math::matrix<double>&U1,math::matrix<double>&U2)
{
	U1=(~A1)*Weight*L;
	U2=(~A2)*Weight*L;
}

/*
 *Solve the equation
 *		| N11		N12 |	| x1 |		| U1 |
		|					   |		|	  |	=	|	   |
 *		| N21		N22 |	| x2 |		| U2 |
 *		ColU2   the dimension of ambiguity to fixed
 *		paraNum1	the number of parameters of type 1(position and so on)
 *	O:
 *		N1		L11
 *		N12		L12
 *		N22		Q22
 *		PartNorm		inv(N11)*N12;
 *	Note:
 *		more detail, refer to the report
 */
void SolveNormalEquationCholesky2(math::matrix<double>& N1,math::matrix<double>& N12,math::matrix<double>& N2,math::matrix<double>& U1,math::matrix<double>& U2,math::matrix<double>&PartNorm,int ColU2,int paraNum1)
{
	N1= Cholesky(N1,paraNum1);//L11 is lower matrix
	
	math::matrix<double>L1Inv;
	L1Inv=InvLowTri(N1,paraNum1);
	PartNorm=MultiplyselfUpperLower(L1Inv,paraNum1)*N12;//inv(N11)*N12;
	N12=L1Inv*N12;//L12;
	
	//math::matrix<double>L22;  
	N2=N2-(~N12)*N12;//QInv22
	N2=Cholesky(N2,ColU2);//L22

	//math::matrix<double>x1_hat;  L11* [trans(L11)*x]=U1
	U1=SolveLowerEquation(N1,U1,paraNum1);//
	//math::matrix<double>x2_hat;
	U2=SolveLowerEquation(N2,U2-(~N12)*U1,ColU2);
//	double resid_temp=-( (~U1)*U1+(~U2)*U2);
	//float solution
	U2=SolveUpperEquation2(N2,U2,ColU2);
	U1=SolveUpperEquation2(N1,U1-N12*U2,paraNum1);

	N2=InvLowTri(N2,N2.ColNo());  //inv(L22)
	N2=MultiplyselfUpperLower(N2,ColU2);  //inv(N2)=Qaa
}

void SolveNormalEquationCholesky2(math::matrix<double>& N1,math::matrix<double>& N12,math::matrix<double>& N2,math::matrix<double>& U1,math::matrix<double>& U2)
{
	int ColU2=N2.RowNo();
	int paraNum1=N1.RowNo();
	N1= Cholesky(N1,paraNum1);//L11 is lower matrix

	math::matrix<double>L1Inv;
	L1Inv=InvLowTri(N1,paraNum1);
//	PartNorm=MultiplyselfUpperLower(L1Inv,paraNum1)*N12;//inv(N11)*N12;
	N12=L1Inv*N12;//L12;

	//math::matrix<double>L22;  
	N2=N2-(~N12)*N12;//QInv22
	N2=Cholesky(N2,ColU2);//L22

	//math::matrix<double>x1_hat;  L11* [trans(L11)*x]=U1
	U1=SolveLowerEquation(N1,U1,paraNum1);//
	//math::matrix<double>x2_hat;
	U2=SolveLowerEquation(N2,U2-(~N12)*U1,ColU2);
	//	double resid_temp=-( (~U1)*U1+(~U2)*U2);
	//float solution
	U2=SolveUpperEquation2(N2,U2,ColU2);
	U1=SolveUpperEquation2(N1,U1-N12*U2,paraNum1);

	N2=InvLowTri(N2,N2.ColNo());  //inv(L22)
	N2=MultiplyselfUpperLower(N2,ColU2);  //inv(N2)=Qaa
}



/* the subscript s denotes the ambiguity  */
extern void SolveNEQCholesky2(math::matrix<double> Nk,math::matrix<double> Nks,math::matrix<double> Ns,math::matrix<double>& Qss,math::matrix<double>& Qks,math::matrix<double>& Qkk, math::matrix<double>& U1,math::matrix<double>& U2,int ColN2,int paraNum1)
{
	math::matrix<double>N1=Nk,N12=Nks,N2=Ns;
	N1= Cholesky(N1,paraNum1);//L11 is lower matrix

	math::matrix<double>L1Inv;
	L1Inv=InvLowTri(N1,paraNum1);
	math::matrix<double> invN11=MultiplyselfUpperLower(L1Inv,paraNum1);
	math::matrix<double> Part=invN11*N12;//inv(N11)*N12;
	N12=L1Inv*N12;//L12;

	//math::matrix<double>L22;  
	N2=N2-(~N12)*N12;//QInv22
	N2=Cholesky(N2,ColN2);//L22

	//math::matrix<double>x1_hat;  L11* [trans(L11)*x]=U1
	U1=SolveLowerEquation(N1,U1,paraNum1);//
	//math::matrix<double>x2_hat;
	U2=SolveLowerEquation(N2,U2-(~N12)*U1,ColN2);
	//	double resid_temp=-( (~U1)*U1+(~U2)*U2);
	//float solution
	U2=SolveUpperEquation2(N2,U2,ColN2);
	U1=SolveUpperEquation2(N1,U1-N12*U2,paraNum1);
	
	N2=InvLowTri(N2,N2.ColNo());  //inv(L22)
	Qss=MultiplyselfUpperLower(N2,ColN2);  //inv(N2)=Qaa
	Qks=-Part*Qss;
	Qkk=invN11*(EyeMat(paraNum1)-Nks*(~Qks));
}



/* subscript a denotes the amb, b is other parameter */
math::matrix<double> FixAmbQXWZ(math::matrix<double>Qaa,math::matrix<double>Qba,math::matrix<double>&Qbb,
							math::matrix<double>& a_hat,math::matrix<double>&b_hat,double threshold,double& ratio)
{
	Ambiguity AR;
	int m=2, n=Qaa.RowNo();
	double *s=new double[m];
	math::matrix<double> FixedAmb(n,m);
	AR.Lambda(n,m,a_hat,Qaa,FixedAmb,s);
	math::matrix<double>a_check=GetBlockMat(FixedAmb,1,n,1,1,2);

	if ( (ratio=s[1]/s[0])>threshold)
	{
		math::matrix<double>Qt= Qba*CholeskyInv(Qaa);
		b_hat=b_hat+Qt*(a_check-a_hat);
		Qbb=Qbb-Qt*(~Qba);
	}
	delete[] s;
	return a_check;
}
/*edit the vcm for QXWZ*/
void editVcm(math::matrix<double>&Qaa,math::matrix<double>&Qba,math::matrix<double>&Qbb,
	math::matrix<double>& a_hat,math::matrix<double>&b_hat)
{
	int dimAA=Qaa.RowNo(), dimBB=Qbb.RowNo();
	math::matrix<double>temp_Qbb(dimBB+1,dimBB+1);
	temp_Qbb(dimBB,dimBB)=Qaa(0,0);
	for (int i=0;i<dimBB;i++)
	{
		temp_Qbb(i,dimBB)=Qba(i,0);
		temp_Qbb(dimBB,i)=temp_Qbb(i,dimBB);
		for (int j=i;j<dimBB;j++)
		{
			temp_Qbb(i,j)=Qbb(i,j);
			temp_Qbb(j,i)=temp_Qbb(i,j);
		}
	}
	math::matrix<double>temp_Qba(dimBB+1,dimAA-1);
	for (int i=1;i<dimAA;i++)
	{
		temp_Qba(dimBB,i-1)=Qaa(0,i);
		for (int j=0;j<dimBB-1;j++)
		{
			temp_Qba(j,i-1)=Qba(j,i);
		}
	}
	Qaa=GetBlockMat(Qaa,2,dimAA,2,dimAA,1);
	Qba=temp_Qba;
	Qbb=temp_Qbb;
	b_hat=VecMat(1,b_hat,GetBlockMat(a_hat,1,1,1,1,2));
	a_hat=GetBlockMat(a_hat,2,dimAA,1,1,2);
}
/*
 *Fix Solution
 *I:
 *	dimAmb	dimension of ambiguity
 *	ddctrl
 *	NormalPart		inv(N11)*N12
 *	Q22					
 *	U2					amb_float
 */
int FixSolu(int dimAmb,double& ratio,DdCtrl ddctrl,math::matrix<double>ambfix,math::matrix<double>NormPart,math::matrix<double>&Q22,math::matrix<double>U2,math::matrix<double>&U1)
{
	Ambiguity AmRe;
	double s[2];s[0]=0.0;s[1]=0.0;
	math::matrix<double>fixRes(dimAmb,2);
	int mode=AmRe.Lambda(dimAmb,2,U2,Q22,fixRes,s);
	if (mode!=-1)
	{
		ratio=s[1]/s[0];
		if (ratio>=ddctrl.ddambctrl.ratiothrsd)
		{
			math::matrix<double> afix(dimAmb,1);
			for (int i=0;i<dimAmb;i++)
			{
				afix(i,0)=fixRes(i,0);
			}
			ambfix=afix;
			U1=U1-NormPart*(afix-U2);
		}
		else
		{
			return -1;
		}
	}
	return mode;
}

void GetParSet(int* pos,int numTofix,int* posRest,int restAmbNum,math::matrix<double>& N11,math::matrix<double>& N12,math::matrix<double>& U1,math::matrix<double>& N22,math::matrix<double>& U2)
{
	int i,j;
	math::matrix<double> SubN22(numTofix,numTofix);
	math::matrix<double> SubN11(restAmbNum+3,restAmbNum+3);
	math::matrix<double> SubN12(restAmbNum+3,numTofix);
	math::matrix<double> SubU1(restAmbNum+3,1);
	math::matrix<double> SubU2(numTofix,1);

	for (i=0;i<3;i++)
	{
		SubU1(i,0)=U1(i,0);
	}

	for (i=0;i<restAmbNum;i++)
	{
		SubU1(i+3,0)=U2(posRest[i],0);
		for (j=0;j<3;j++)
		{
			SubN11(j,i+3)=N12(j,posRest[i]);
		}
		SubN11(i+3,i+3)=SubN22(posRest[i],posRest[i]);
			for (j=i+1;j<restAmbNum;j++)
			{
				SubN11(i+3,j+3)=SubN22(posRest[i],posRest[j]);
				SubN11(j+3,i+3)=SubN11(i+3,j+3);
			}
	}
	PutMat(SubN11,N11,1,1,1);
	for (i=0;i<numTofix;i++)
	{
		for (j=0;j<3;j++)
		{
			SubN12(j,i)=N12(j,pos[i]);
		}
		SubU2(i,0)=U2(pos[i],0);
		SubN22(i,i)=N22(pos[i],pos[i]);
		for (j=i+1;j<numTofix;j++)
		{
			SubN22(i,j)=N22(pos[i],pos[j]);
			SubN22(j,i)=N22(pos[i],pos[j]);
		}
	}
	N11=SubN11;N12=SubN12;N22=SubN22;SubU1=U1;SubU2=U2;
}

void GetParPos(int* pos,int* posRest,int& numTofix,int& restAmbNum,int* numOfList,int* prntofix,DdCtrl ddctrl,DdObsInfo obsinfo,DdAmbInfo ambinfo)
{
	double ParMakEle=ddctrl.ddambctrl.parEleMask;
	int num=ddctrl.PhsTypeNo();
	int i,j,k,cnt=0,cnt1=0,cnt2=0,cnt3=0;
	double elesum=0.0;
	for (i=0;i<num;i++)
	{
		cnt=0;
		for (k=0;k<i;k++)
		{
			cnt+=ambinfo.NoUnfix(k);
		}
		cnt1=0;
		for (j=0;j<obsinfo.NoPhs(i);j++)
		{
		elesum=(obsinfo.eleRovRov[j]+obsinfo.eleRovBase[j])/2.0;
			if (ambinfo.fixFlag[i][j]==0)//&&elesum*R2D>ParMakEle
			{
				if (elesum*R2D>=ParMakEle)
				{
					prntofix[cnt2]=ambinfo.prnList[i][j];
					pos[cnt2++]=cnt+j;
				}
				else
				{
					posRest[cnt3++]=cnt+j;
				}

			}//end if
		}//end j
		numOfList[i]=cnt1;
	}//end i
	restAmbNum=cnt3;
	numTofix=cnt2;
}

/*
 *eliminate the parameter of x3
 *
 */
void ElimPara3(math::matrix<double>& N1,math::matrix<double>& N12,math::matrix<double>N13,math::matrix<double>& N2,math::matrix<double> N23,math::matrix<double> N3, math::matrix<double>& U1,math::matrix<double>& U2,math::matrix<double>U3)
{
	math::matrix<double> temp;
	int numMatRe=0;
	math::matrix<double> matcon;
	N3=CholeskyInv(N3,N3.RowNo());
	matcon=N13*N3;
	N1	=N1-matcon*(~N13);//````
	U1	=U1-matcon*U3;//````

	math::matrix<double> trans;
	temp=N23*N3;
	U2	=U2- temp*U3;//````
	N12	=N12-matcon*trans;//````
	N2		=N2-temp*trans;//````
	
	/*
	N1	=N1-N13*N3*(~N13);
	N12	=N12-N13*N3*(~N23);
	N2	=N2-N23*N3*(~N23);
	U1	=U1-N13*N3*U3;
	U2	=U2-N23*N3*U3;*/
}

/*
 *the observation information from DDdata
 *I:
 *	PrePrnList
 *	CurPrnList
 *	PreNum
 *	CurNum
 *O:
 *	upNum			initiating = 0
 *	downNum	initiating = 0
 *	upPrn			add those sate to the end of equations
 *	downPrn		eliminate the information of those sate
 *Note:
 *	 
 */
void SateUpDown(int* PrePrnList,int*CurPrnList,int PreNum,int CurNum,int& upNum,int& downNum,int* upPrn,int* downPrn, int* posUp,int* posDown)
{
	int flag=0;
	//for falling 
	for (int i=0;i<PreNum;i++)
	{
		flag=0;
		for (int j=0;j<CurNum;j++)
		{
			if (PrePrnList[i]==CurPrnList[j])
			{
				flag=1;
				break;
			}
		}
		if (flag==0)
		{
			posDown[downNum]=i;
			downPrn[downNum++]=PrePrnList[i];
		}
	}

	
	//for rising
	for (int i=0;i<CurNum;i++)
	{
		flag=0;
		for (int j=0;j<PreNum;j++)
		{
			if (PrePrnList[j]==CurPrnList[i])
			{
				flag=1;
				break;
			}
		}
		if (flag==0)
		{
			posUp[upNum]=i;
			upPrn[upNum++]=CurPrnList[i];
		}
	}
}

void SateUp(int& upPhsNum,int& upNum,int* upPrn,int* posUp,int* posUpPhs,DdObsInfo& preddobsinfo,DdObsInfo currddobsinfo,DdCtrl ddctrl)
{
	upNum=0;
	int PreNum, CurNum,i,j,k;
	int PrePrnList[MAXNUMSATE],CurPrnList[MAXNUMSATE];
	int flag=0,cnt=0;
	int num;
	num=(ddctrl.pseudoFlag>3)?ddctrl.pseudoFlag-3:ddctrl.pseudoFlag;
	for (k=0;k<num;k++)
	{
		cnt=0;
		for (i=0;i<k;i++)
		{
			cnt+=currddobsinfo.numCod[i];
		}
		preddobsinfo.SetList(PrePrnList,k,0);
		currddobsinfo.SetList(CurPrnList,k,0);
		PreNum=preddobsinfo.numCod[k];
		CurNum=currddobsinfo.numCod[k];
	//for rising
		for (int i=0;i<CurNum;i++)
		{
			flag=0;
			for (int j=0;j<PreNum;j++)
			{
				if (PrePrnList[j]==CurPrnList[i]&&CurPrnList[i]>0)
				{
					flag=1;
					break;
				}
			}
			if (flag==0)
			{
				posUp[upNum]=i+cnt;
				upPrn[upNum++]=CurPrnList[i];
			}
		}//end i
	}//end k

	upPhsNum=0;
	int cnt1 =preddobsinfo.SumPhs(),cnt2=0;
	num=(ddctrl.ddambctrl.flag>3)?ddctrl.ddambctrl.flag-3:ddctrl.ddambctrl.flag;
	for (k=0;k<num;k++)
	{
		cnt=cnt1;
		cnt2=0;
		for (i=0;i<k;i++)
		{
			cnt+=currddobsinfo.numPhs[i];
			cnt2+=currddobsinfo.numPhs[i];
		}
		preddobsinfo.SetList(PrePrnList,k,1);
		currddobsinfo.SetList(CurPrnList,k,1);
		PreNum=preddobsinfo.numPhs[k];
		CurNum=currddobsinfo.numPhs[k];
		//for rising
		for (int i=0;i<CurNum;i++)
		{
			flag=0;
			for (int j=0;j<PreNum;j++)
			{
				if (PrePrnList[j]==CurPrnList[i]&&CurPrnList[i]>0)
				{
					flag=1;
					break;
				}
			}
			if (flag==0)
			{
				posUpPhs[upPhsNum]=i+cnt2;
				upPhsNum++;
				posUp[upNum]=i+cnt;
				upPrn[upNum++]=CurPrnList[i];
			}
		}//end i
	}//end k
}

void ProcUp(int preNum,int currNum,int upNum,int* posUp,math::matrix<double>*preNorEqu)
{
	int preNo=preNum,currNo=currNum;
	if (upNum!=0 )
	{
		int i,j;
		for (i=0;i<upNum;i++)//for upNum>1
		{
			preNorEqu[1]=InsertZeroCol(preNorEqu[1],posUp[i],1);//ConvergeMat(3,preNorEqu[1],preNum,upNum)*
			preNorEqu[2]=InsertZeroRowCol(preNorEqu[2],posUp[i],1);//DiagMatSym(preNorEqu[2],preNum,upNum);
			preNorEqu[4]=InsertZeroRow(preNorEqu[4],posUp[i],1);//VecMat(1,preNorEqu[4],preNum,upNum);
			for (j=i+1;j<upNum;j++) posUp[j]++; 
		}

	}
}

void SateDown(int& downPhsNum,int& downNum,int* downPrn,int* posDown,int* posDownPhs,DdObsInfo& preddobsinfo,DdObsInfo currddobsinfo,DdCtrl ddctrl)
{
	DdObsInfo temp;
	downNum=0;
	int PreNum, CurNum,i,j,k;
	int PrePrnList[MAXNUMSATE],CurPrnList[MAXNUMSATE];
	int flag=0,cnt=0;
	//for falling 
	int num;
	num=(ddctrl.pseudoFlag>3)?ddctrl.pseudoFlag-3:ddctrl.pseudoFlag;
	for (k=0;k<num;k++)
	{
		cnt=0;
		for (i=0;i<k;i++)
		{
			cnt+=preddobsinfo.numCod[i];
		}
		
		preddobsinfo.SetList(PrePrnList,k,0);
		currddobsinfo.SetList(CurPrnList,k,0);
		PreNum=preddobsinfo.numCod[k];
		CurNum=currddobsinfo.numCod[k];
		for (int i=0;i<PreNum;i++)
		{
			flag=0;
			for (int j=0;j<CurNum;j++)
			{
				if (PrePrnList[i]==CurPrnList[j]&&CurPrnList[j]>0)
				{
					flag=1;
					temp.prnlistCod[k][temp.numCod[k]]=CurPrnList[j];
					temp.numCod[k]++;
					break;
				}
			}
			if (flag==0)
			{
				posDown[downNum]=i+cnt;
				downPrn[downNum++]=PrePrnList[i];
			}
		}

	}
	downPhsNum=0;
	int cnt1 =preddobsinfo.SumCod(),cnt2=0;
	num=(ddctrl.ddambctrl.flag>3)?ddctrl.ddambctrl.flag-3:ddctrl.ddambctrl.flag;
	for (k=0;k<num;k++)
	{
		cnt=cnt1;cnt2=0;
		for (i=0;i<k;i++)
		{
			cnt+=preddobsinfo.numPhs[i];
			cnt2+=preddobsinfo.numPhs[i];
		}
		preddobsinfo.SetList(PrePrnList,k,1);
		currddobsinfo.SetList(CurPrnList,k,1);
		PreNum=preddobsinfo.numPhs[k];
		CurNum=currddobsinfo.numPhs[k];
		for (int i=0;i<PreNum;i++)
		{
			flag=0;
			for (int j=0;j<CurNum;j++)
			{
				if (PrePrnList[i]==CurPrnList[j]&&CurPrnList[j]>0)
				{
					flag=1;
					temp.prnlistPhs[k][temp.numPhs[k]]=CurPrnList[j];
					temp.numPhs[k]++;
					break;
				}
			}
			if (flag==0)
			{
				posDownPhs[downPhsNum]=i+cnt2;
				downPhsNum++;
				posDown[downNum]=i+cnt;
				downPrn[downNum++]=PrePrnList[i];
			}
		}
		
	}
	preddobsinfo=temp;
}
/*
 *remove the fallen satellite in normal equation
 *just process the amb and pos, ignore the iono
 */
void ProcDown(int preNum,int currNum,int downPhsNum,int* posDownPhs,math::matrix<double>*preNorEqu)
{
	int preNo=preNum,currNo=currNum;
	if (downPhsNum!=0 )
	{
		//suppose that the falling sates are at the bottom of list in preNormEqu, and... 
		//now the only one is supported, but there is problem in the multi-freq 
		math::matrix<double>vecDown1,vecDown2,temp,vecU2;
		int i,j;
		for (i=0;i<downPhsNum;i++)
		{
			preNorEqu[1]=MoveVecColEnd(preNorEqu[1],posDownPhs[i]+1);

			preNorEqu[2]=MoveVecRowColEnd(preNorEqu[2],posDownPhs[i]+1);

			preNorEqu[4]=MoveVecRowEnd(preNorEqu[4],posDownPhs[i]+1);

			for (j=i+1;j<downPhsNum;j++) 	posDownPhs[j]--;
		}
		int numblk=preNo-downPhsNum;
		vecDown1=GetBlockMat(preNorEqu[1],1,3,numblk+1,preNo,2);

		vecDown2=GetBlockMat(preNorEqu[2],1,numblk,numblk+1,preNo,2);

		temp=GetBlockMat(preNorEqu[2],numblk+1,preNo,preNo-downPhsNum+1,preNo,1);
		temp=CholeskyInv(temp,downPhsNum);

		preNorEqu[0]=preNorEqu[0]-vecDown1*temp*(~vecDown1);

		preNorEqu[1]=preNorEqu[1]-vecDown1*temp*(~vecDown2);
		preNorEqu[2]=GetBlockMat(preNorEqu[2],1,numblk,1,numblk,1)-vecDown2*temp*(~vecDown2);
		
		vecU2=GetBlockMat(preNorEqu[4],numblk+1,preNo,1,1,2);
		preNorEqu[3]=preNorEqu[3]-vecDown1*temp*vecU2;
		preNorEqu[4]=GetBlockMat(preNorEqu[4],1,numblk,1,1,2)-vecDown2*temp*vecU2;
	}
}

/*
 *The superposition of normal equation, just in static situation and initialization
 *I:
 *  currNorEqu   currN11,currN12,currN22,currU1,currU2
 *  
 * O:
 *  preNorEqu   the order is same as currNorEqu
 */
void SuperpositionNorEqu(math::matrix<double>*currNorEqu,math::matrix<double>*preNorEqu,DdCtrl ddctrl,DdObsInfo& preddobsinfo,DdObsInfo& currddobsinfo)
{
	//first, eliminate the N11
	//currNorEqu[2]=currNorEqu[2]-(~currNorEqu[1]*CholeskyInv(currNorEqu[0])*currNorEqu[1]);
	//currNorEqu[4]=currNorEqu[4]-(~currNorEqu[1]*CholeskyInv(currNorEqu[0])*currNorEqu[3]);

	int upNum=0,downNum=0,downPhsNum=0,upPhsNum=0,i;
	int upPrn[MAXSATERISE],downPrn[MAXSATEFALL];
	int posUp[MAXSATERISE],posUpPhs[MAXSATERISE],posDownPhs[MAXSATEFALL],posDown[MAXSATEFALL];
	
	int preNum=preddobsinfo.SumCod()+preddobsinfo.SumPhs();
	int curNum=currddobsinfo.SumCod()+currddobsinfo.SumPhs();
	SateDown(downPhsNum,downNum,downPrn,posDown,posDownPhs,preddobsinfo,currddobsinfo,ddctrl);
	ProcDown(preNum,curNum,downPhsNum,posDownPhs,preNorEqu);

	SateUp(upPhsNum,upNum,upPrn,posUp,posUpPhs,preddobsinfo,currddobsinfo,ddctrl);
	ProcUp(preNum,curNum,upPhsNum,posUpPhs,preNorEqu);
	//cout<<preNorEqu[4].ColNo()<<"  "<<currNorEqu[4].ColNo()<<"  "<< preNorEqu[2].ColNo() <<"  "<< currNorEqu[2].ColNo() <<endl;
	for ( i=0;i<5;i++)	preNorEqu[i]=preNorEqu[i]+currNorEqu[i];
}

/*not recommended*/
double PartialAr(math::matrix<double>& N11,math::matrix<double>& N12,math::matrix<double>& N22,math::matrix<double>& U1,math::matrix<double>&U2,
					DdAmbInfo ambinfo,DdObsInfo obsinfo,double maskele,double& ratio)
{
	DdAmbInfo parambinfo=ambinfo;
	math::matrix<double>Ub=U1,Norm12=N12;
	/*here N11=inv(N11) */
	N22=N22-~N12*N11*N12;
	U2=U2-~N12*N11*U1;
	N22=CholeskyInv(N22);
	U2=N22*U2;
	U1=N11*(U1-N12*U2);
	N12=-N11*N12*N22;
	double ele=maskele;
	N22=SortQaa(N22,N12,U2,ambinfo,parambinfo,obsinfo);
	math::matrix<double> QbaPar,QaaPar,aPar;
	int mode;

	//while(ele<80.0)
	//{
		aPar=U2;
		QaaPar=SelectQaa(N22,N12,QbaPar,aPar,obsinfo,parambinfo,ele);
		if(QaaPar.RowNo()==1) return ele;//break;
		int dim=aPar.RowNo();
		Ambiguity ar;
		math::matrix<double>FF(dim,2);
		double  ss[2];
		ar.Lambda(dim,2,aPar,QaaPar,FF,ss);
		ratio=ss[1]/ss[0];
		if (ratio>2.0)
		{
			math::matrix<double>fixAmb=GetBlockMat(FF,1,FF.RowNo(),1,1,2);
			U1=U1-QbaPar*CholeskyInv(QaaPar)*(aPar-fixAmb);
			//cout<<~U1;
			//			U1=N11*(Ub-Norm12*fixAmb);
			int	resultmode=2;
			return ele; //break;
		}
		else
		{
			ele+=10.0;
		}
	//}
	return ele;
}

/*****************************************************************************************************************/
/*get the index of sate whose ele is bigger than elemask, for PAR
 * and the index that is less than elemask
 * elemask unit: degree
 */
void ParIndex(int typeNo,DdData curdata,int* indexToBeFixed,int& numFix,int& numLeft,int* indexLeft,double elemask)
{
	InitPtr(indexToBeFixed,MAXNUMSATE);
	InitPtr(indexLeft,MAXNUMSATE);
	int index=0;
	numFix=numLeft=0;
	int i,j;
	for (i=0;i<typeNo;i++)
	{
		for (j=0;j<curdata.pairNum;j++)
		{
			if (curdata.datarecord[j].vadFlgPhs[i]!=0)
			{
				if (curdata.ele[j]*R2D>=elemask)
				{
					indexToBeFixed[numFix++]=index;
				}
				else
				{
					indexLeft[numLeft++]=index;
				}
				index++;
			}//end if
		}//end loop j
	}
}

/*after getting the index, extract the partial vcm and a_hat
 * Qaa are seperate QFix, QLeft, QLeft_fix,and Qba is Qba_left, Qba_fix
 * a_hat --> a_left, a_fix
 */
void extractVcm(math::matrix<double>Qaa,math::matrix<double>Qba,math::matrix<double>a_hat,math::matrix<double>&QFix,
						math::matrix<double>&QLeft,math::matrix<double>&QLeft_Fix,math::matrix<double>&Qba_Left,math::matrix<double>&Qba_Fix,
						math::matrix<double>&a_fix,math::matrix<double>&a_left,int* indexToBeFixed,int numFix,int numLeft,int* indexLeft)
{
	int dimAA=Qaa.RowNo(), dimBB=Qba.RowNo();
	int i,j,Row,Col;
	for (i=0;i<numFix;i++)
	{
		Row=indexToBeFixed[i];
		a_fix(i,0)=a_hat(Row,0);
		for (j=i;j<numFix;j++)
		{
			Col=indexToBeFixed[j];
			QFix(i,j)=Qaa(Row,Col);
			QFix(j,i)=Qaa(Row,Col);
		}
		for (j=0;j<dimBB;j++)
		{
			Qba_Fix(j,i)=Qba(j,Row);
		}
		for (j=0;j<numLeft;j++)
		{
			Col=indexLeft[j];
			QLeft_Fix(j,i)=Qaa(Col,Row);
		}
	}
	for (i=0;i<numLeft;i++)
	{
		Row=indexLeft[i];
		a_left(i,0)=a_hat(Row,0);
		for (j=i;j<numLeft;j++)
		{
			Col=indexLeft[j];
			QLeft(i,j)=Qaa(Row,Col);
			QLeft(j,i)=Qaa(Row,Col);
		}
		for (j=0;j<dimBB;j++)
		{
			Qba_Left(j,i)=Qba(j,Row);
		}
	}
}


/* after extracting the partial vcm, rebuild the vcm and x_hat
 * this step is prepared for update after the LAMBDA
 * or update with block mat(simple) 
 * 
 * but this function can be used in another condition
 */
void rebuildVcm(math::matrix<double>& Qbb_re,math::matrix<double>&Qba_re,math::matrix<double>& b_re,
						math::matrix<double>Qbb,math::matrix<double>Qba_left,math::matrix<double>Qaa_left,math::matrix<double>Qba_fix,
						math::matrix<double>Qaa_Left_fix,math::matrix<double>b_hat,math::matrix<double>a_left)
{
	int row=Qbb.RowNo(),col=Qaa_left.RowNo();
	b_re=VecMat(1,b_hat,a_left);
	Qba_re=VecMat(Qba_fix.ColNo(),Qba_fix,Qaa_Left_fix);
	Qbb_re=DiagMatSym(Qbb,Qaa_left);
	for (int i=0;i<row;i++)
	{
		for (int j=0;j<col;j++)
		{
			Qbb_re(i,row+j)=Qba_left(i,j);
			Qbb_re(row+j,i)=Qbb_re(i,row+j);
		}
	}
}

/*
 after fixing the ambiguity, update the vcm with block matrix
 this function replace the rebuild

 Qbb		Qba_left	|	Qba_fix				#	b_hat
			Qaa_left	|	Qaa_left_fix		#	a_left
sym						|	Qaa_fix				#	a_fix
 */
void updateWithBlock(math::matrix<double>& Qbb,math::matrix<double>& Qba_left,math::matrix<double>&Qaa_left,
								math::matrix<double>&b_hat,math::matrix<double>&a_left,math::matrix<double> Qaa_fix,
								math::matrix<double>Qba_fix,math::matrix<double>Qaa_Left_fix,math::matrix<double>a_check,math::matrix<double>a_fix_float)
{
	math::matrix<double>invQaa=CholeskyInv(Qaa_fix);
	math::matrix<double>Qt=Qba_fix*invQaa;
	b_hat=b_hat+Qt*(a_check-a_fix_float);
	Qbb=Qbb-Qt*(~Qba_fix);
	Qba_left=Qba_left-Qt*(~Qaa_Left_fix);
	
	Qt=Qaa_Left_fix*invQaa;
	a_left=a_left+Qt*(a_check-a_fix_float);
	Qaa_left=Qaa_left-Qt*(~Qaa_Left_fix);
}

/*
 *I: 
 *	Qamb		contaions the vcm of trop and amb
 *	Qbb			the vcm of iono
 *	a_hat		the float solu of trop and amb
 *	b_hat		the float solu of iono 
 *	
 *O:
 *	a_hat_update			the updated float solu of amb
 *	a_check					the fixed solu of amb
 *	b_hat_update			the updated float solu of amb
 */
extern void QXWZPartialARwithUpdate(math::matrix<double>Qamb,math::matrix<double>Qb_amb,math::matrix<double>Qbb,
												math::matrix<double>a_hat,math::matrix<double> b_hat,math::matrix<double>&a_hat_update,
												math::matrix<double>&a_check,math::matrix<double>&b_hat_update,
												DdData curdata,double& ratio,double eleThreshold,double ratioThrs,
												DdAmbInfo& curAmb)
{
	math::matrix<double> Qamb_copy=Qamb,Qb_amb_copy=Qb_amb,Qbb_copy=Qbb,a_hat_copy=a_hat,b_hat_copy=b_hat;
	/*		*/
	editVcm(Qamb_copy,Qb_amb_copy,Qbb_copy,a_hat_copy,b_hat_copy);
	int indexToBeFixed[MAXNUMSATE],indexLeft[MAXNUMSATE];
	InitPtr(indexToBeFixed,MAXNUMSATE);	InitPtr(indexLeft,MAXNUMSATE);
	int numFix=0;
	int numLeft=0;
	int typeNo=2; /*the typeNo need to input from outsides*/
	double elemask=eleThreshold;
	ParIndex(typeNo,curdata, indexToBeFixed,numFix,numLeft, indexLeft,elemask);

	math::matrix<double> Qfix(numFix,numFix);
	math::matrix<double>QLeft(numLeft,numLeft),QLeft_Fix(numLeft,numFix),
									Qba_Left(Qb_amb_copy.RowNo(),numLeft),Qba_Fix(Qb_amb_copy.RowNo(),numFix),
									a_fix(numFix,1),a_left(numLeft,1);
	extractVcm(Qamb_copy,Qb_amb_copy,a_hat_copy,Qfix,
						QLeft,QLeft_Fix,Qba_Left,Qba_Fix,a_fix,a_left,
						indexToBeFixed,numFix, numLeft,indexLeft);

	/*Qfix and a_fix to do PAR  */
	Ambiguity AR;
	int m=2, n=Qfix.RowNo();
	double *s=new double[m];
	math::matrix<double> FixedAmb(n,m);
	AR.Lambda(n,m,a_fix,Qfix,FixedAmb,s);
	a_check=GetBlockMat(FixedAmb,1,n,1,1,2);
	
	ratio=s[1]/s[0];
	if (ratio>ratioThrs)
	{
		//cout<<a_check;
		updateWithBlock(Qbb_copy,Qba_Left,QLeft,b_hat_copy,a_left,Qfix,Qba_Fix,QLeft_Fix,a_check,a_fix);
		a_hat_update=a_left;
		b_hat_update=b_hat_copy;
		TransferFixAmb(2,curAmb,curdata,a_check,numFix,indexToBeFixed);
	}
	//return ratio;
}


/*
 * transfer the fixed amb
 * 
 * checked
 */
extern void TransferFixAmb(int typeNo,DdAmbInfo& curAmb,DdData curdata,math::matrix<double>afix,int numFix,int* indexToBeFixed)
{
	DdAmbInfo tempAmb=curAmb;
	int* group=new int[typeNo]; InitPtr(group,typeNo);
	int* cnt=new int[typeNo];InitPtr(cnt,typeNo);
	int i,j,preDim;
	for(i=0;i<typeNo;i++)	
	{
		cnt[i]=curAmb.NoSat(i);
		group[i]=i;
	}
	for (i=numFix;i>0;i--)
	{
		int ind=indexToBeFixed[i-1];
		//int typeNum=typeNo;
		int loopNum=0;
		while(ind>=curdata.pairNum)
		{
			ind-=cnt[loopNum];
			loopNum++;
		}
		tempAmb.fixFlag[loopNum][ind]=1;
		tempAmb.fixSolu[loopNum][ind]=afix(i-1,0);
	}

	delete[] group,cnt;
	curAmb=tempAmb;
}

/******************************************************check AR*********************************************************/

/*for QXWZ 
 *curdata     obs    /meter
 *Lc		the combination of IF  
 */
void IonoFreeData(DdData curdata,DdData& comData,DdObsInfo obsinfo,DdObsInfo& obsinfo_temp)
{
	int sysid=Prn2Sysid(curdata.refPrn);
	double freq1=FreqSys(1,0)*1e-9, freq2=FreqSys(1,1)*1e-9;
	double num1=SQ(freq1)/(SQ(freq1)-SQ(freq2)), num2=SQ(freq2)/(SQ(freq1)-SQ(freq2));
	int cnt1=0,cnt2=0;
	obsinfo_temp.eleRefBase=obsinfo.eleRefBase;
	obsinfo_temp.eleRefRov=obsinfo.eleRefRov;
	for (int i=0;i<curdata.pairNum;i++)
	{
		if (curdata.datarecord[i].vadFlgCod[0]*curdata.datarecord[i].vadFlgCod[1]==1)
		{
			comData.datarecord[cnt1].PsRange[0]	=num1*curdata.datarecord[i].PsRange[0]-
																			num2*curdata.datarecord[i].PsRange[1];
			comData.ele[cnt1]=curdata.ele[i];
			comData.mapWet[cnt1]=curdata.mapWet[i];
			comData.rovPrn[cnt1]=curdata.rovPrn[i];
			comData.satePosRov[cnt1]=curdata.satePosRov[i];
			comData.satePosBase[cnt1]=curdata.satePosBase[i];
			comData.tropCor[cnt1]=curdata.tropCor[i];
			comData.datarecord[cnt1].numVadCod++;
			comData.datarecord[cnt1].vadFlgCod[0]=1;
			obsinfo_temp.numCod[0]++;
			obsinfo_temp.prnlistCod[0][cnt1]=curdata.rovPrn[i];
			obsinfo_temp.eleRovRov[cnt1]=obsinfo.eleRovRov[i];
			obsinfo_temp.eleRovBase[cnt1]=obsinfo.eleRovBase[i];
			cnt1++;
		}
		if (curdata.datarecord[i].vadFlgPhs[0]*curdata.datarecord[i].vadFlgPhs[1]==1)
		{
			comData.datarecord[cnt2].Phase[0]	=num1*curdata.datarecord[i].Phase[0]-
																			num2*curdata.datarecord[i].Phase[1];
			comData.datarecord[cnt2].numVadPhs++;
			comData.datarecord[cnt2].vadFlgPhs[0]=1;
			obsinfo_temp.numPhs[0]++;
			obsinfo_temp.prnlistPhs[0][cnt1]=curdata.rovPrn[i];
			cnt2++;
		}
	}
	comData.pairNum=cnt2;
	comData.refPrn=curdata.refPrn;
	comData.sec=curdata.sec;
	comData.week=curdata.week;
	for (int i=0;i<3;i++)
	{
		comData.refRecPos[i]=curdata.refRecPos[i];
		comData.refSatPos_Base[i]=curdata.refSatPos_Base[i];
		comData.refSatPos_Rov[i]=curdata.refSatPos_Rov[i];
		comData.rovRecPos[i]=curdata.rovRecPos[i];
	}
}

/* for QXWZ
 *the unit of ambiguity in ambinfo is cycle
 *the unit of ambiguity in comAmb is meter
 */
void IonoFreeAmb(DdAmbInfo ambinfo,DdAmbInfo& comAmb,int sysid)
{
	double freq1= FreqSys(sysid,0),freq2=FreqSys(sysid,1);
	double lam1=CLIGHT/freq1,lam2=CLIGHT/freq2;
	freq1*=1e-9;
	freq2*=1e-9;
	int numsat=ambinfo.NoSat(0);
	//PtrEqual(ambinfo.prnList[0],comAmb.prnList[0],numsat);
	int ind2=0,cnt=0;
	int numsate2=ambinfo.NoSat(1);
	for (int i=0;i<numsat;i++)
	{
		ind2=FindPosInt(ambinfo.prnList[1],numsate2,ambinfo.prnList[0][i]);
		if (ind2!=-1 && ambinfo.fixFlag[0][i]*ambinfo.fixFlag[1][ind2]!=0)
		{
			comAmb.fixSolu[0][cnt]=SQ(freq1)/(SQ(freq1)-SQ(freq2))*lam1*ambinfo.fixSolu[0][i]-
												SQ(freq2)/(SQ(freq1)-SQ(freq2))*lam2*ambinfo.fixSolu[0][ind2];
			comAmb.fixFlag[0][cnt]=1;
		}
		comAmb.prnList[0][cnt]=ambinfo.prnList[0][i];
		cnt++;
	}
}


/*For QXWZ 
 *DesMatPos		from current epoch 
 *L					double difference const, modify L, 
 *here ambiguity is meter
 *return 	L+N
 */
math::matrix<double> ReFormConstWithAmb(math::matrix<double>L, DdAmbInfo ambinfo,DdCtrl ddctrl,DdData dddata)
{
	math::matrix<double>errEquConst=L;
	int num	=ddctrl.PhsTypeNo();
	int i,j,k,cnt1=0;
	int cnt=dddata.SumCod();
	double lam;
	for (i=0;i<num;i++)
	{
		//lam=CLIGHT/ddctrl.freqPhs[i];
		for (j=0;j<dddata.pairNum;j++)
		{
			if (dddata.datarecord[j].vadFlgPhs[i]==1)
			{
				int pos,numsate=ambinfo.NoSat(i);
				pos=FindPosInt(ambinfo.prnList[i],numsate,dddata.rovPrn[j]);
				if (pos!=-1 || ambinfo.fixFlag[i][pos]==1)
				{
						errEquConst(cnt,0)+=ambinfo.fixSolu[i][pos];
				}
				cnt++;
			}
		}
	}
	return errEquConst;
}


/*
 *after fix the ambiguity, check the correctness of AR with Lc 
 *	1. rebuild data;
 *	2. rebuild fixed IF-amb
 *	3. rebuild Desmat-->DesmatPos,  L, weight
 */
extern void checkAR(DdAmbInfo ambinfo,DdData curdata,DdCtrl ddctrl,DdObsInfo obsinfo)
{
	Position pt;
	DdCtrl ddctrl_temp=ddctrl;
	ddctrl_temp.codtype[1]=0; ddctrl_temp.pseudoFlag=4;
	ddctrl_temp.phstype[1]=0; ddctrl_temp.ddambctrl.flag=4;
	
	DdData data_temp;
	DdObsInfo obsinfo_temp;
	/*rebuild data */
	IonoFreeData(curdata,data_temp,obsinfo,obsinfo_temp);

	DdAmbInfo ambinfo_temp;
	/*rebuild IF-ambiguity*/
	IonoFreeAmb(ambinfo,ambinfo_temp,ddctrl.sysid);

 	int row=data_temp.SumCod()+data_temp.SumPhs();
	
	math::matrix<double> DesMatPos(row,3);
	
	math::matrix<double>	Weight(row,row);
	math::matrix<double>	L_if(row,1);

	pt.FormDesMatPos(DesMatPos,L_if,data_temp,1,1);
	pt.FormResidual(L_if,data_temp,ddctrl_temp,ambinfo_temp);
	
	int dimAmb=ambinfo_temp.NoUnfix(0);
	int numCod=data_temp.SumCod();
	int numPhs=data_temp.SumPhs();
	math::matrix<double> DesMatAmb(row,dimAmb);
	pt.FormDesMatAmb(DesMatAmb,ddctrl_temp,data_temp,obsinfo_temp,ambinfo_temp);

	/*notice the coefficient of desmatamb*/
	ddctrl_temp.weightMode=ddctrl.weightMode;
	Weight=pt.FormWeightVc(ddctrl_temp,data_temp,obsinfo_temp);
	math::matrix<double> Nt=~DesMatPos*Weight;
	math::matrix<double> N1=Nt*DesMatPos;
	math::matrix<double> N12=Nt*DesMatAmb;
	math::matrix<double> Ntt=~DesMatAmb*Weight;
	math::matrix<double> N2=Ntt*DesMatAmb;
	math::matrix<double>Lm= L_if;
	//Lm=ReFormConstWithAmb(L_if,ambinfo_temp,ddctrl_temp,data_temp);
	math::matrix<double> U1=Nt*Lm;
	math::matrix<double> U2=Ntt*Lm;	
	SolveNormalEquationCholesky2(N1,N12,N2,U1,U2);
	cout<<U1;
}

/***************************************************************************************************************

/*
 *after solution, pass the fixed ambiguity to current ambinfo
 *I:
 *	afix  the fixed ambiguity
 *	resultflag the flag of fixing 0=unfixed 1=fix all 2=partial
 */
void PassFixedAmb(DdCtrl ddctrl,DdAmbInfo& curAmb,DdObsInfo obsinfo,math::matrix<double>afix, int resultflag)
{
	if (resultflag==1)
	{			
		int i,j,k,cnt=0;
		int* prnlistPar=new int[curAmb.SumUnfix()]; //list of unfixed sate

		/* get the prnlist of all sate*/
		for (i=0;i<ddctrl.PhsTypeNo();i++)
		{
			for (j=0;j<curAmb.NoSat(i);j++)
			{
				if (curAmb.fixFlag[i][j]==0)	prnlistPar[cnt++]=curAmb.prnList[i][j];
			}
		}
	
		for (i=0;i<ddctrl.PhsTypeNo();i++)
		{
			cnt=0;
			for (k=0;k<i;k++) cnt+=curAmb.NoUnfix(i);
			int loop1=curAmb.NoUnfix(i);
			int loop2=curAmb.NoSat(i);
			int s=0;
			for (j=0;j<loop1;j++)
			{
				for (k=s;k<loop2;k++)
				{
					if (prnlistPar[j+cnt]==curAmb.prnList[i][k])
					{
						curAmb.fixFlag[i][k]=1;
						curAmb.fixSolu[i][k]=afix(j+cnt,0);
						s=k+1;
						break;
					}
				}
			}
		}
		delete [] prnlistPar;
	}
}


/*debug!!!*/
void SuperpositionNorEquAmb(DdCtrl ddctrl,math::matrix<double>& PreNorAmb,math::matrix<double> CurNorAmb,DdAmbInfo PreAmbInfo,DdAmbInfo CurAmbInfo)
{
	int i,j,k,cntpre=0,cntcur=0,cnt=0;
	int* CurPrnlist=new int[CurAmbInfo.SumUnfix()]; //list of unfixed sate
	int* PrePrnlist=new int[PreAmbInfo.SumUnfix()];
	int typeNo=ddctrl.PhsTypeNo();
	for (i=0;i<typeNo;i++)
	{
		int numcur=CurAmbInfo.NoSat(i);
		for (j=0;j<numcur;j++)
		{
			if (CurAmbInfo.fixFlag[i][j]==0)	CurPrnlist[cnt++]=CurAmbInfo.prnList[i][j];
		}
		numcur=PreAmbInfo.NoSat(i);
		for (j=0;j<numcur;j++)
		{
			if (PreAmbInfo.fixFlag[i][j]==0)	PrePrnlist[cnt++]=PreAmbInfo.prnList[i][j];
		}
	}

	int* upPos=new int[CurAmbInfo.SumNoSat()/3];
	int* downPos=new int[CurAmbInfo.SumNoSat()/3];
	int upNum=0,downNum=0;
	for (i=0;i<typeNo;i++)
	{
		cntcur=cntpre=0;
		for(j=0;j<i;j++) 
		{
			cntcur+=CurAmbInfo.NoUnfix(j);
			cntpre+=PreAmbInfo.NoUnfix(j);
		}
		int unfixc=CurAmbInfo.NoUnfix(i);
		for (j=0;j<unfixc;j++)
		{
			int upflag=0;
			for (k=0;k<unfixc;k++)
			{
				if (CurPrnlist[j+cntcur]==PrePrnlist[k+cntpre]) break;
				if(CurPrnlist[j+cntcur]<PrePrnlist[k+cntpre]) 
				{
					upflag=1;
					break;
				}
				if(CurPrnlist[j+cntcur]>PrePrnlist[k+cntpre]) continue;
			}
			if(upflag==1) upPos[upNum++]=j+cntcur;
		}
		int unfixp=PreAmbInfo.NoUnfix(i);
		for (j=0;j<unfixp;j++)
		{
			int downflag=0;
			for (k=0;k<unfixc;k++)
			{
				if (CurPrnlist[k+cntcur]==PrePrnlist[j+cntpre]) break;
				if(CurPrnlist[k+cntcur]>PrePrnlist[j+cntpre]) 
				{
					downflag=1;
					break;
				}
				if(CurPrnlist[k+cntcur]<PrePrnlist[j+cntpre]) continue;
			}
			if(downflag==1) downPos[downNum++]=j+cntpre;
		}

	}//end i

	delete[] downPos,upPos,PrePrnlist, CurPrnlist;
}

/*
 *find the rise sate and fallen sate, record the PRN and number, by comparing the current and previous AmbInfo
 */
void PosUpDown(DdCtrl ddctrl, DdAmbInfo CurAmbInfo,DdAmbInfo PreAmbInfo,int& upNum,int& downNum,int* upPos,int* downPos)
{
	int i,j,k,cntpre=0,cntcur=0,cnt1=0,cnt2=0;
	int* CurPrnlist=new int[CurAmbInfo.SumUnfix()]; //list of unfixed sate
	int* PrePrnlist=new int[PreAmbInfo.SumUnfix()];
	int typeNo=ddctrl.PhsTypeNo();
	for (i=0;i<typeNo;i++)
	{
		int numcur=CurAmbInfo.NoSat(i);
		for (j=0;j<numcur;j++)
		{
			if (CurAmbInfo.fixFlag[i][j]==0)	CurPrnlist[cnt1++]=CurAmbInfo.prnList[i][j];
		}
		numcur=PreAmbInfo.NoSat(i);
		for (j=0;j<numcur;j++)
		{
			if (PreAmbInfo.fixFlag[i][j]==0)	PrePrnlist[cnt2++]=PreAmbInfo.prnList[i][j];
		}
	}

	//int* upPos=new int[CurAmbInfo.SumNoSat()/3];
	//int* downPos=new int[CurAmbInfo.SumNoSat()/3];
	upNum=downNum=0;
	for (i=0;i<typeNo;i++)
	{
		cntcur=cntpre=0;
		for(j=0;j<i;j++) 
		{
			cntcur+=CurAmbInfo.NoUnfix(j);
			cntpre+=PreAmbInfo.NoUnfix(j);
		}
		//up
		int unfixc=CurAmbInfo.NoUnfix(i);
		for (j=0;j<unfixc;j++)
		{
			int upflag=0;
			for (k=0;k<unfixc;k++)
			{
				if (CurPrnlist[j+cntcur]==PrePrnlist[k+cntpre]) break;
				if(CurPrnlist[j+cntcur]<PrePrnlist[k+cntpre]) 
				{
					upflag=1;
					break;
				}
				if(CurPrnlist[j+cntcur]>PrePrnlist[k+cntpre]) continue;
			}
			if(upflag==1) upPos[upNum++]=j+cntcur;
		}

		//down 
		int unfixp=PreAmbInfo.NoUnfix(i);
		for (j=0;j<unfixp;j++)
		{
			int downflag=0;
			for (k=0;k<unfixc;k++)
			{
				if (CurPrnlist[k+cntcur]==PrePrnlist[j+cntpre]) break;
				if(CurPrnlist[k+cntcur]>PrePrnlist[j+cntpre]) 
				{
					downflag=1;
					break;
				}
				if(CurPrnlist[k+cntcur]<PrePrnlist[j+cntpre]) continue;
			}
			if(downflag==1) downPos[downNum++]=j+cntpre;
		}

	}//end i

	delete[] PrePrnlist, CurPrnlist;
}
/*
 * superposition the normal equ of long baseline 
 * the LSE of k th epoch need to be solved, then update the LSE with the ambiguity information of k-1 th
 * with the variance-covariance matrix (VCM) 
 * ignore the constraint of trops and ionos
 *I:
 *						the prnlist of cur
 *						the prnlist of pre
 *PreNor						the VCM of ambiguity in k-1 th, it eliminates the other parameters.
 *
 *						the float solu of k-1 th 
 *O:
 *						the updated VCM of k th epoch, and the VCM contains the Trops, Ionos and the amb,,[0,1,2]=[Atmosphere , At_Am, Amb ]
 *						the updated solu of k th epoch, and it contains the Trops, Ionos and the amb solu
 *	Note:		
 *			2015.12
 */
extern void SuperPosLonBWithFixedCrd(DdCtrl ddctrl,math::matrix<double>& PreVcmAmb,math::matrix<double>& PreSoluAmb,
							math::matrix<double>* CurVcm,math::matrix<double>* CurSolu,DdAmbInfo PreAmbInfo,DdAmbInfo CurAmbInfo,DdData curData,DdData preData)
{
	int  upNum=0, dowmNum=0;
	int upPos[6];
	int downPos[6];
	PosUpDown(ddctrl,CurAmbInfo,PreAmbInfo,upNum,dowmNum,upPos,downPos);
	int i=0;
	int paraNum=1;/*the number of parameters of trops*/
	 //if the sate falls in k th, eliminate PreVcmAmb of the certain sate in  amb part 
	for (i=dowmNum;i>0&&dowmNum>0;i--)
	{
		PreVcmAmb=ElimRowCol(PreVcmAmb,PreSoluAmb,downPos[i-1]+paraNum+1);
	}
	/*process the cycle slip*/
	math::matrix<double> Q_t;//=CholeskyInv(PreVcmAmb);

	// if the satellite rises in k th, inflate the VCM of k-1 th with 0 in ionos part and amb part 
	for (i=upNum;i>0&&upNum>0;i--)
	{
		PreVcmAmb=CholeskyInv(PreVcmAmb);// N
		PreSoluAmb=PreVcmAmb*PreSoluAmb;//U

		PreVcmAmb=InsertZeroRowCol(PreVcmAmb,upPos[i-1]+1+paraNum);
		PreSoluAmb=InsertZeroRow(PreSoluAmb,upPos[i-1]+paraNum,1);

		PreVcmAmb=CholeskyInv(PreVcmAmb);//Q
		PreSoluAmb=PreVcmAmb*PreSoluAmb;//x_hat
	}

	math::matrix<double> invAmb=CholeskyInv(PreVcmAmb+CurVcm[2]);
	CurSolu[0]=CurSolu[0]+CurVcm[1]*invAmb*(PreSoluAmb-CurSolu[1]);
	CurSolu[1]=CurSolu[1]+CurVcm[2]*invAmb*(PreSoluAmb-CurSolu[1]);

	CurVcm[0]=CurVcm[0]-CurVcm[1]*invAmb*(~CurVcm[1]);
	CurVcm[1]=CurVcm[1]-CurVcm[1]*invAmb*CurVcm[2];
	CurVcm[2]=CurVcm[2]-CurVcm[2]*invAmb*CurVcm[2];

}


extern void SuperPosLonBWithFixedCrdNEQ(DdCtrl ddctrl,math::matrix<double>& PreNEQAmb,math::matrix<double>& PreNEQUAmb,
	math::matrix<double> CurNEQAmb,math::matrix<double> CurNEQUAmb,DdAmbInfo PreAmbInfo,DdAmbInfo CurAmbInfo,
	DdData curData,DdData preData)
{
	int  upNum=0, dowmNum=0;
	int upPos[6];
	int downPos[6];
	PosUpDown(ddctrl,CurAmbInfo,PreAmbInfo,upNum,dowmNum,upPos,downPos);
	int i=0,j;
	int paraNum=1;/*the number of parameters of trops*/
	//if the sate falls in k th, eliminate PreNEQAmb of the certain sate in  amb part 
	for (i=dowmNum;i>0&&dowmNum>0;i--)
	{
		PreNEQAmb=ElimRowColNEQ(PreNEQAmb,PreNEQUAmb,downPos[i-1]+paraNum+1);
	}
	/*process the cycle slip*/
	math::matrix<double> Q_t=CholeskyInv(PreNEQAmb);
	math::matrix<double> x_t=Q_t*PreNEQUAmb;
	int typeNo=ddctrl.PhsTypeNo(),cnt=0;
	for (i=0;i<typeNo;i++)
	{
		for (j=0;j<curData.pairNum;j++)
		{
			if (curData.datarecord[j].vadFlgPhs[i]==1)
			{
				if (curData.datarecord[j].isCycleSlip[i]==1)
				{
					Q_t(cnt,cnt)=INFINUMBER;
				}
				cnt++;
			}
		}
	}

	// if the satellite rises in k th, inflate the NEQ of k-1 th with 0 in ionos part and amb part 
	for (i=upNum;i>0&&upNum>0;i--)
	{
		PreNEQAmb=InsertZeroRowCol(PreNEQAmb,upPos[i-1]+1+paraNum);
		PreNEQUAmb=InsertZeroRow(PreNEQUAmb,upPos[i-1]+paraNum,1);
	}
	PreNEQAmb+=CurNEQAmb;
	PreNEQUAmb+=CurNEQUAmb;
}


/*
 *The sequential solution of short baseline 
 *ignore the ionosphere and troposphere parameters
 *I:
 *fixMode   1=all, 2=par  
 *
 *O:
 *mode		3 =all fixed sequentially, 4=partial sequentially 0= unfixed
 *
 *Note:
 *	this is the initialization of the RTK of short baseline, too.
 */
extern math::matrix<double> SoluShortSequence(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double> Weight,math::matrix<double>L,
													DdObsInfo obsinfo,math::matrix<double>* preNorEqu,DdAmbInfo ambinfo ,DdObsInfo& preddobsinfo,DdObsInfo curddobsinfo,int& resultmode,double& ratio,int fixMode)
{
	int dimAmb=ambinfo.SumUnfix();
	int ambunfixed=dimAmb;
	dimAmb=(dimAmb==0)?1:dimAmb;
	
	math::matrix<double> N11, N12, N22, U1,U2;
	FormNormalEquationCoef2(DesMatPos,DesMatAmb,N11,N12,N22,Weight);
	FormNormalEquationConst2(DesMatPos,DesMatAmb,Weight,L,U1,U2);
	N11=CholeskyInv(N11);
	
		N22=N22-~N12*N11*N12;
		U2=U2-~N12*N11*U1;
		/* add  */
		
		N22=CholeskyInv(N22);
		U2=N22*U2;
		math::matrix<double>Ub=U1;
		U1=N11*(U1-N12*U2);
		Ambiguity ar;
		math::matrix<double>FF(N22.RowNo(),2);
		double  ss[2];
		ar.Lambda(N22.RowNo(),2,U2,N22,FF,ss);
		ratio=ss[1]/ss[0];
		if (ratio>ddctrl.ddambctrl.ratiothrsd)
		{
			math::matrix<double>fixAmb=GetBlockMat(FF,1,FF.RowNo(),1,1,2);
			U1=N11*(Ub-N12*fixAmb);
			resultmode=1;
			PassFixedAmb(ddctrl,ambinfo,obsinfo,fixAmb,resultmode);
		}
	

	resultmode=0;
	if (ambunfixed==0)
	{
		U1=CholeskyInv(N11)*U1;
		resultmode=1;
	}
	if (ambunfixed>0)
	{
		math::matrix<double>currNorEqu[5];
		currNorEqu[0]=N11;currNorEqu[1]=N12,currNorEqu[2]=N22;
		currNorEqu[3]=U1;currNorEqu[4]=U2;
		SuperpositionNorEqu(currNorEqu,preNorEqu,ddctrl,preddobsinfo,curddobsinfo);
		
		N11=preNorEqu[0];N12=preNorEqu[1];N22=preNorEqu[2];U1=preNorEqu[3];U2=preNorEqu[4];
		math::matrix<double>NormPart;
		math::matrix<double>fixAmb(dimAmb,1);
		int* prnlist=new int[dimAmb]; //list of unfixed sate
		int numOfList[3];
		numOfList[0]=0;numOfList[1]=0;numOfList[2]=0;

		if (fixMode==1)
		{// all ambiguities are fixed
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,dimAmb,3);
			//now N11=L11, N12=L12, N22=Q22=Qamb, NormPart=inv(N1)*N12 ,U1=x1_hat, U2=x2_hat
			resultmode=FixSolu(dimAmb,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
			//if(ambunfixed>0) PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
		}
		else if (fixMode==2)
		{//partial ambiguities 
			int numTofix=0,paraNum1=0,restAmbNum=0;//dimAmb+3=dimTofix+paraNum1
			int *posFix=new int[ambinfo.SumUnfix()]; //position to fix
			double*ele=new double[ambinfo.SumUnfix()];
			int *posRest=new int[ambinfo.SumUnfix()];

			GetParPos(posFix,posRest,numTofix,restAmbNum,numOfList,prnlist,ddctrl,obsinfo,ambinfo);
			GetParSet(posFix,numTofix,posRest,restAmbNum,N11,N12,U1,N22,U2);
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,numTofix,restAmbNum+3);
			resultmode=FixSolu(numTofix,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
//			PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
			delete[] posFix,posRest,ele;
		}
		delete[] prnlist;
	}
	preddobsinfo=curddobsinfo;
	return U1;
}

/*
 *consider the iono, eliminate the delay . then superposition
 */
extern math::matrix<double> SoluIonoSequence(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double>DesMatIono,math::matrix<double> Weight,math::matrix<double>L,
																		DdObsInfo obsinfo,math::matrix<double>* preNorEqu,DdAmbInfo ambinfo ,DdObsInfo& preddobsinfo,DdObsInfo curddobsinfo,int& resultmode,double& ratio,int fixMode)
{
	int dimAmb=ambinfo.SumUnfix();
	int ambunfixed=dimAmb;
	dimAmb=(dimAmb==0)?1:dimAmb;
	//int fixMode=ddctrl.ambFlag;

	math::matrix<double> N11=(~DesMatPos)*Weight*DesMatPos;
	math::matrix<double> N12=(~DesMatPos)*Weight*DesMatAmb;
	math::matrix<double> N13=(~DesMatPos)*Weight*DesMatIono;
	math::matrix<double> N23=(~DesMatAmb)*Weight*DesMatIono;
	math::matrix<double> N22=(~DesMatAmb)*Weight*DesMatAmb;
	math::matrix<double> N33=(~DesMatIono)*Weight*DesMatIono;
	math::matrix<double> U1=(~DesMatPos)*Weight*L;
	math::matrix<double> U2=(~DesMatAmb)*Weight*L;
	math::matrix<double> U3=(~DesMatIono)*Weight*L;
	//eliminate  N33
	N33=CholeskyInv(N33);
	N11	=N11-N13*N33*(~N13);
	N12	=N12-N13*N33*(~N23);
	N22	=N22-N23*N33*(~N23);
	U1	=U1-N13*N33*U3;
	U2	=U2-N23*N33*U3;

	resultmode=0;
	if (ambunfixed==0)
	{
		U1=CholeskyInv(N11)*U1;
		resultmode=1;
	}
	if (ambunfixed>0)
	{
		math::matrix<double>currNorEqu[5];
		currNorEqu[0]=N11;currNorEqu[1]=N12,currNorEqu[2]=N22;
		currNorEqu[3]=U1;currNorEqu[4]=U2;
		SuperpositionNorEqu(currNorEqu,preNorEqu,ddctrl,preddobsinfo,curddobsinfo);
		N11=preNorEqu[0];N12=preNorEqu[1];N22=preNorEqu[2];U1=preNorEqu[3];U2=preNorEqu[4];
		math::matrix<double>NormPart;
		math::matrix<double>fixAmb(dimAmb,1);
		int* prnlist=new int[dimAmb]; //list of unfixed sate
		int numOfList[3];
		numOfList[0]=0;numOfList[1]=0;numOfList[2]=0;

		if (fixMode==1)
		{// all ambiguities are fixed
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,dimAmb,3);
			//now N11=L11, N12=L12, N22=Q22=Qamb, NormPart=inv(N1)*N12 ,U1=x1_hat, U2=x2_hat

			resultmode=FixSolu(dimAmb,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
		//	PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
		}
		else if (fixMode==2)
		{//partial ambiguities 
			int numTofix=0,paraNum1=0,restAmbNum=0;//dimAmb+3=dimTofix+paraNum1
			int *posFix=new int[ambinfo.SumUnfix()]; //position to fix
			double*ele=new double[ambinfo.SumUnfix()];
			int *posRest=new int[ambinfo.SumUnfix()];

			GetParPos(posFix,posRest,numTofix,restAmbNum,numOfList,prnlist,ddctrl,obsinfo,ambinfo);
			GetParSet(posFix,numTofix,posRest,restAmbNum,N11,N12,U1,N22,U2);
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,numTofix,restAmbNum+3);
			//change the size of N11, N12, N22 U1, U2. N11 and N12 include the unfixed set 
			resultmode=FixSolu(numTofix,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
	//		PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
			delete[] posFix,posRest,ele;
		}
		delete[] prnlist;
	}
	return U1;
}
/*
 *The solution of short baseline epoch by epoch
 *ignore the ionosphere and troposphere parameters
 *I:
 *fixMode   1=all, 2=par  
 *
 *O:
 *resultmode		1 =all fixed epoch by epoch, 2=partial epoch by epoch 0= unfixed
 *
 *return x_pos
 */
extern math::matrix<double> SoluShortEpoch(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double> Weight,math::matrix<double>L,
													  DdObsInfo obsinfo,  DdAmbInfo& ambinfo, int& resultmode,double& ratio,math::matrix<double>* prenorm)
{
	int dimAmb=ambinfo.SumUnfix();
	int ambunfixed=dimAmb;
	dimAmb=(dimAmb==0)?1:dimAmb;
	int fixMode=ddctrl.ambFlag;
	
	math::matrix<double> N11(3,3);
	math::matrix<double> N12(3,dimAmb);
	math::matrix<double> N22(dimAmb,dimAmb);
	math::matrix<double> U1(3,1);
	math::matrix<double> U2(dimAmb,1);
	FormNormalEquationCoef2(DesMatPos,DesMatAmb,N11,N12,N22,Weight);
	FormNormalEquationConst2(DesMatPos,DesMatAmb,Weight,L,U1,U2);
	N11=CholeskyInv(N11);
	N22=N22-~N12*N11*N12;
	U2=U2-~N12*N11*U1;
	N22=CholeskyInv(N22);
	U2=N22*U2;
	math::matrix<double>Ub=U1;
	U1=N11*(U1-N12*U2);
	Ambiguity ar;
	math::matrix<double>FF(N22.RowNo(),2);
	double  ss[2];
	ar.Lambda(N22.RowNo(),2,U2,N22,FF,ss);
	if (ratio=ss[1]/ss[0]>ddctrl.ddambctrl.ratiothrsd)
	{
		U1=N11*(Ub-N12*GetBlockMat(FF,1,FF.RowNo(),1,1,2));
	}
	cout<<ss[1]/ss[0];

	prenorm[0]=N11;prenorm[1]=N12;prenorm[2]=N22;prenorm[3]=U1;prenorm[4]=U2;
	if (ambunfixed==0)
	{
		U1=CholeskyInv(N11)*U1;
		resultmode=1;
	}
	if (ambunfixed>0)
	{
		math::matrix<double>currN11,currN12,currN22,currU1,currU2;
		math::matrix<double>NormPart;
		math::matrix<double>fixAmb(dimAmb,1);

		int* prnlist=new int[dimAmb]; //list of unfixed sate
		int numOfList[3];
		numOfList[0]=0;numOfList[1]=0;numOfList[2]=0;
		if (fixMode==1)
		{// all ambiguities are fixed
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,dimAmb,3);
			//now N11=L11, N12=L12, N22=Q22=Qamb, NormPart=inv(N1)*N12 ,U1=x1_hat, U2=x2_hat
			resultmode=FixSolu(dimAmb,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
//			if(ambunfixed>0) PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
		}
		else if (fixMode==2)
		{//partial ambiguities 
			int numTofix=0,restAmbNum=0;//dimAmb+3=dimTofix+paraNum1
			int *posFix=new int[ambinfo.SumUnfix()]; double*ele=new double[ambinfo.SumUnfix()];
			int *posRest=new int[ambinfo.SumUnfix()];

			GetParPos(posFix,posRest,numTofix,restAmbNum,numOfList,prnlist,ddctrl,obsinfo,ambinfo);
			GetParSet(posFix,numTofix,posRest,restAmbNum,N11,N12,U1,N22,U2);
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,numTofix,restAmbNum+3);
			resultmode=FixSolu(numTofix,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
//			PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
			
			delete[] posFix,posRest,ele;
		}
		delete[] prnlist;
	}
	
	return U1;//x_hat 
}

extern math::matrix<double> SoluShortEpoch(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double> Weight,math::matrix<double>L,
	DdObsInfo obsinfo,  DdAmbInfo& ambinfo, int& resultmode,double& ratio)
{
	int dimAmb=ambinfo.SumUnfix();
	int ambunfixed=dimAmb;
	dimAmb=(dimAmb==0)?1:dimAmb;
	int fixMode=ddctrl.ambFlag;

	math::matrix<double> N11, N12, N22, U1,U2;
	FormNormalEquationCoef2(DesMatPos,DesMatAmb,N11,N12,N22,Weight);
	FormNormalEquationConst2(DesMatPos,DesMatAmb,Weight,L,U1,U2);
	N11=CholeskyInv(N11);
	int s=1;
	if (ambunfixed==0)
	{
		resultmode=10;
		U1=(N11)*U1;
	}
	else if(s==1)
	{
		N22=N22-~N12*N11*N12;
		U2=U2-~N12*N11*U1;
		N22=CholeskyInv(N22);
		U2=N22*U2;
		math::matrix<double>Ub=U1;
		U1=N11*(U1-N12*U2);
		Ambiguity ar;
		math::matrix<double>FF(N22.RowNo(),2);
		double  ss[2];
		ar.Lambda(N22.RowNo(),2,U2,N22,FF,ss);
		ratio=ss[1]/ss[0];
		if (ratio>ddctrl.ddambctrl.ratiothrsd)
		{
			math::matrix<double>fixAmb=GetBlockMat(FF,1,FF.RowNo(),1,1,2);
			for (int i=0;i<6;i++)
			{
				//cout<<setiosflags(ios::fixed)<<setprecision(7)<<setw(20)<<fixAmb(i,0)*CLIGHT*FREQ1/(SQ(FREQ1)-SQ(FREQ2))-fixAmb(i+6,0)*CLIGHT*FREQ2/(SQ(FREQ1)-SQ(FREQ2))<<endl;
				cout<<setiosflags(ios::fixed)<<setprecision(7)<<setw(20)<<fixAmb(i,0)-fixAmb(i+6,0)<<endl;
			}
			U1=N11*(Ub-N12*fixAmb);
			resultmode=1;
			//PassFixedAmb(ddctrl,ambinfo,obsinfo,fixAmb,resultmode);
		}
	}
	else if(s==0)
	{
		double ele;
		ele=PartialAr(N11,N12,N22,U1,U2,ambinfo,obsinfo,15.0,ratio);
	}
	
	return U1;
}

extern math::matrix<double> SoluShortEpoch(DdCtrl ddctrl,math::matrix<double>* NormEqu,DdAmbInfo& ambinfo, int& resultmode,double& ratio)
{
	int flagmode=ddctrl.ambFlag;		double ratiothrsd=ddctrl.ddambctrl.ratiothrsd;
	math::matrix<double> N11=NormEqu[0], N12=NormEqu[1], N22=NormEqu[2], U1=NormEqu[3],U2=NormEqu[4];
	int dimAmb=N22.RowNo();
	int ambunfixed=dimAmb;
	dimAmb=(dimAmb==0)?1:dimAmb;
	int fixMode=flagmode;
	N11=CholeskyInv(N11);
	if (ambunfixed==0)
	{
		resultmode=10;
		U1=(N11)*U1;
	}
	else if(fixMode==1)
	{
		N22=N22-~N12*N11*N12;
		U2=U2-~N12*N11*U1;
		N22=CholeskyInv(N22);
		U2=N22*U2;
		math::matrix<double>Ub=U1;
		U1=N11*(U1-N12*U2);
		Ambiguity ar;
		math::matrix<double>FF(N22.RowNo(),2);
		double  ss[2];
		ar.Lambda(N22.RowNo(),2,U2,N22,FF,ss);
		ratio=ss[1]/ss[0];
		if (ratio>ratiothrsd)
		{
			math::matrix<double>fixAmb=GetBlockMat(FF,1,FF.RowNo(),1,1,2);
			U1=N11*(Ub-N12*fixAmb);
			resultmode=1;
			//PassFixedAmb(ddctrl,ambinfo,obsinfo,fixAmb,resultmode);
		}
	}
	else if(fixMode==2)
	{
		double ele;
	//	ele=PartialAr(N11,N12,N22,U1,U2,ambinfo,obsinfo,15.0,ratio);
	}

	return U1;
}

/*
 *take the ionosphere delay into account 
 *eliminate the iono 
 */
extern math::matrix<double> SoluEpochIono(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double>DesMatIono,math::matrix<double> Weight,math::matrix<double>L,
																		DdObsInfo obsinfo,  DdAmbInfo& ambinfo, int& resultmode,double& ratio,math::matrix<double>* prenorm)
{
	int dimAmb=ambinfo.SumUnfix();
	int ambunfixed=dimAmb;
	dimAmb=(dimAmb==0)?1:dimAmb;
	int fixMode=ddctrl.ambFlag;

	math::matrix<double> N11=(~DesMatPos)*Weight*DesMatPos;
	math::matrix<double> N12=(~DesMatPos)*Weight*DesMatAmb;
	math::matrix<double> N13=(~DesMatPos)*Weight*DesMatIono;
	math::matrix<double> N23=(~DesMatAmb)*Weight*DesMatIono;
	math::matrix<double> N22=(~DesMatAmb)*Weight*DesMatAmb;
	math::matrix<double> N33=(~DesMatIono)*Weight*DesMatIono;
	math::matrix<double> U1=(~DesMatPos)*Weight*L;
	math::matrix<double> U2=(~DesMatAmb)*Weight*L;
	math::matrix<double> U3=(~DesMatIono)*Weight*L;
	//eliminate  N33
	N33=CholeskyInv(N33);
	N11	=N11-N13*N33*(~N13);
	N12	=N12-N13*N33*(~N23);
	N22	=N22-N23*N33*(~N23);
	U1	=U1-N13*N33*U3;
	U2	=U2-N23*N33*U3;

	prenorm[0]=N11;prenorm[1]=N12;prenorm[2]=N22;prenorm[3]=U1;prenorm[4]=U2;
	if (ambunfixed==0)
	{
		U1=CholeskyInv(N11)*U1;
		resultmode=1;
	}
	if (ambunfixed>0)
	{
		math::matrix<double>currN11,currN12,currN22,currU1,currU2;
		math::matrix<double>NormPart;
		math::matrix<double>fixAmb(dimAmb,1);

		int* prnlist=new int[dimAmb]; //list of unfixed sate
		int numOfList[3];
		numOfList[0]=0;numOfList[1]=0;numOfList[2]=0;
		if (fixMode==1)
		{// all ambiguities are fixed
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,dimAmb,3);
			//now N11=L11, N12=L12, N22=Q22=Qamb, NormPart=inv(N1)*N12 ,U1=x1_hat, U2=x2_hat
			resultmode=FixSolu(dimAmb,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
//			if(ambunfixed>0) PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);
		}
		else if (fixMode==2)
		{//partial ambiguities 
			int numTofix=0,restAmbNum=0;//dimAmb+3=dimTofix+paraNum1
			int *posFix=new int[ambinfo.SumUnfix()]; double*ele=new double[ambinfo.SumUnfix()];
			int *posRest=new int[ambinfo.SumUnfix()];

			GetParPos(posFix,posRest,numTofix,restAmbNum,numOfList,prnlist,ddctrl,obsinfo,ambinfo);
			GetParSet(posFix,numTofix,posRest,restAmbNum,N11,N12,U1,N22,U2);
			SolveNormalEquationCholesky2(N11,N12,N22,U1,U2,NormPart,numTofix,restAmbNum+3);
			resultmode=FixSolu(numTofix,ratio,ddctrl,fixAmb,NormPart,N22,U2,U1);
			resultmode=(resultmode==0)?fixMode:0;
//			PassFixedAmb(ambinfo,obsinfo,fixAmb,prnlist,numOfList,resultmode);

			delete[] posFix,posRest,ele;
		}
		delete[] prnlist;
	}

	return U1;//x_hat 
}


extern void SoluCrossCode( int refpos,SdData lastSdData,int nEpoch,math::matrix<double> DesMatPos,math::matrix<double>L,math::matrix<double>*Ne,math::matrix<double>*Le,math::matrix<double>*Nr,math::matrix<double>*Lr)
{
	double c0=0.330,c1=0.380;
	math::matrix<double>Qzd(lastSdData.satnum,lastSdData.satnum);
	for (int ik=0;ik<lastSdData.satnum;ik++) Qzd(ik,ik)=SQ( c0/(sin(lastSdData.ele[ik])+c1));
	math::matrix<double>Df;Df=EyeMat(lastSdData.satnum);
	Df=RemoveRow(Df,refpos+1);
	SetMatClo(Df,refpos+1,-1.0);
	math::matrix<double>Qdd0;
	math::matrix<double>I;
	I=EyeMat(2);
	Qdd0=Kronecker( I,Df*Qzd*(~Df)*2.0, 1);
	Qdd0=CholeskyInv(Qdd0);

	I(0,1)=0.53;I(1,0)=0.53;
	math::matrix<double>Qsd1;
	Qsd1=2.0*Kronecker(I,Qzd,2);
	I(0,1)=0.0;I(1,0)=0.0;
	math::matrix<double>Qdd1;
	Qdd1=Kronecker(I,Df,2)*Qsd1*Kronecker(I,Df,2);
	Qdd1=CholeskyInv(Qdd1);
	Ne[nEpoch]=DesMatPos*Qdd0*(~DesMatPos);
	Le[nEpoch]=DesMatPos*Qdd0*L;

	Ne[nEpoch]=DesMatPos*Qdd1*(~DesMatPos);
	Le[nEpoch]=DesMatPos*Qdd1*L;
}

extern void AddCode(int num,int interval,math::matrix<double>*Ne,math::matrix<double>*Le,math::matrix<double>*Nr,math::matrix<double>*Lr)
{
	int i,j;
	math::matrix<double>Ne1(3,3),Nte(3,3),Le1(3,1),Lte(3,1);
	math::matrix<double>Nr1(3,3),Ntr(3,3),Lr1(3,1),Ltr(3,1);
	for (i=0;i<num;i++)
	{
		if(i==num-interval)	break;
		//ZeroMat(Ne1);ZeroMat(Le1);
		//ZeroMat(Nr1);ZeroMat(Lr1);
		if (i==0)
		{
			for (j=i;j<i+interval;j++)
			{
				Ne1=Ne1+Ne[j];
				Le1=Le1+Le[j];

				Nr1=Nr1+Nr[j];
				Lr1=Lr1+Lr[j];
			}
		}
		else
		{
			Ne1=Ne1+Ne[i]-Nte;
			Le1=Le1+Le[i]-Lte;

			Nr1=Nr1+Nr[i]-Ntr;
			Lr1=Lr1+Lr[i]-Ltr;
		}
		Nte=Ne[i];
		Ne[i]=CholeskyInv(Ne1);//N11^-1
		Lte=Le[i];
		Le[i]=Ne1*Le1;//store the xhat

		Ntr=Nr[i];
		Nr[i]=CholeskyInv(Nr1);//N11^-1
		Ltr=Lr[i];
		Lr[i]=Nr1*Lr1;//store the xhat
	}
}

extern void IonoFreeData(DdData curdata,DdData& comData,DdObsInfo* obsinfo,DdObsInfo* obsinfo_temp)
{
	int sysid=Prn2Sysid(curdata.refPrn);
	double freq1=FreqSys(1,0)*1e-9, freq2=FreqSys(1,1)*1e-9;
	double num1=SQ(freq1)/(SQ(freq1)-SQ(freq2)), num2=SQ(freq2)/(SQ(freq1)-SQ(freq2));
	int cnt1=0,cnt2=0;
	if (obsinfo!=NULL)
	{
		obsinfo_temp->eleRefBase=obsinfo->eleRefBase;
		obsinfo_temp->eleRefRov=obsinfo->eleRefRov;
	}

	for (int i=0;i<curdata.pairNum;i++)
	{
		if (curdata.datarecord[i].vadFlgCod[0]*curdata.datarecord[i].vadFlgCod[1]==1)
		{
			comData.datarecord[cnt1].PsRange[0]	=num1*curdata.datarecord[i].PsRange[0]-
				num2*curdata.datarecord[i].PsRange[1];
			comData.ele[cnt1]=curdata.ele[i];
			comData.mapWet[cnt1]=curdata.mapWet[i];
			comData.rovPrn[cnt1]=curdata.rovPrn[i];
			comData.satePosRov[cnt1]=curdata.satePosRov[i];
			comData.satePosBase[cnt1]=curdata.satePosBase[i];
			comData.tropCor[cnt1]=curdata.tropCor[i];
			comData.datarecord[cnt1].numVadCod++;
			comData.datarecord[cnt1].vadFlgCod[0]=1;
			if (obsinfo!=NULL)
			{
				obsinfo_temp->numCod[0]++;
				obsinfo_temp->prnlistCod[0][cnt1]=curdata.rovPrn[i];
				obsinfo_temp->eleRovRov[cnt1]=obsinfo->eleRovRov[i];
				obsinfo_temp->eleRovBase[cnt1]=obsinfo->eleRovBase[i];
			}
			cnt1++;
		}
		if (curdata.datarecord[i].vadFlgPhs[0]*curdata.datarecord[i].vadFlgPhs[1]==1)
		{
			comData.datarecord[cnt2].Phase[0]	=num1*curdata.datarecord[i].Phase[0]-
				num2*curdata.datarecord[i].Phase[1];
			comData.datarecord[cnt2].numVadPhs++;
			comData.datarecord[cnt2].vadFlgPhs[0]=1;
			if (obsinfo!=NULL)
			{
				obsinfo_temp->numPhs[0]++;
				obsinfo_temp->prnlistPhs[0][cnt1]=curdata.rovPrn[i];
			}
			cnt2++;
		}
	}
	comData.pairNum=cnt2;
	comData.refPrn=curdata.refPrn;
	comData.sec=curdata.sec;
	comData.week=curdata.week;
	for (int i=0;i<3;i++)
	{
		comData.refRecPos[i]=curdata.refRecPos[i];
		comData.refSatPos_Base[i]=curdata.refSatPos_Base[i];
		comData.refSatPos_Rov[i]=curdata.refSatPos_Rov[i];
		comData.rovRecPos[i]=curdata.rovRecPos[i];
	}
}
