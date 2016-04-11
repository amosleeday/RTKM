#include "stdafx.h"//changed
#include "VariableDef.h"
#include "MStruct.h"
#include<math.h>//calc
#include "ExtFun.h"
#include "Ambiguity.h"
#include "matrix.h"

/*		this file includes the solution method 	*/
/*		after the class of Position						*/


/*
 *FUNCTION:
 *						Form the normal equations
 *						Solve the normal equations
 *
 *MODE:				
 *						static						kinematic
 *		
 *METHOD:
 *						least-square				kalman filter
 *
 *PARAMETERIZATION:
 *						Ionosphere  (fix, float, weighted)
 *						Troposphere(fix, float[const, random walk], decorrelation)
 *
 */

extern math::matrix<double> SoluShortSequence(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double> Weight,math::matrix<double>L, DdObsInfo obsinfo,math::matrix<double>* preNorEqu,DdAmbInfo ambinfo ,DdObsInfo& preddobsinfo,DdObsInfo curddobsinfo,int& resultmode,double& ratio,int fixMode);

extern math::matrix<double> SoluShortEpoch(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double> Weight,math::matrix<double>L, DdObsInfo obsinfo, DdAmbInfo& ambinfo, int& resultmode,double& ratio,math::matrix<double>* prenorm);

extern math::matrix<double> SoluShortEpoch(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double> Weight,math::matrix<double>L, DdObsInfo obsinfo, DdAmbInfo& ambinfo, int& resultmode,double& ratio);

extern math::matrix<double> SoluShortEpoch(DdCtrl ddctrl,math::matrix<double>* NormEqu,DdAmbInfo& ambinfo, int& resultmode,double& ratio);

extern math::matrix<double> SoluEpochIono(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double>DesMatIono,math::matrix<double> Weight,math::matrix<double>L, DdObsInfo obsinfo, DdAmbInfo& ambinfo, int& resultmode,double& ratio,math::matrix<double>* prenorm);

extern math::matrix<double> SoluIonoSequence(DdCtrl ddctrl,math::matrix<double> DesMatPos,math::matrix<double> DesMatAmb,math::matrix<double>DesMatIono,math::matrix<double> Weight,math::matrix<double>L, DdObsInfo obsinfo,math::matrix<double>* preNorEqu,DdAmbInfo ambinfo ,DdObsInfo& preddobsinfo,DdObsInfo curddobsinfo,int& resultmode,double& ratio,int fixMode);

extern void SoluCrossCode(int nEpoch,math::matrix<double> DesMAtPos,math::matrix<double>L,math::matrix<double>*Ne,math::matrix<double>*Le,math::matrix<double>*Nr,math::matrix<double>*Lr);

extern void PosUpDown(DdCtrl ddctrl, DdAmbInfo CurAmbInfo,DdAmbInfo PreAmbInfo,int& upNum,int& downNum,int* upPos,int* downPos);
extern void SuperPosLonBWithFixedCrd(DdCtrl ddctrl,math::matrix<double>& PreVcmAmb,math::matrix<double>& PreSoluAmb, math::matrix<double>* CurVcm,math::matrix<double>* CurSolu,DdAmbInfo PreAmbInfo,DdAmbInfo CurAmbInfo,DdData curData,DdData preData);
extern void SuperPosLonBWithFixedCrdNEQ(DdCtrl ddctrl,math::matrix<double>& PreNEQAmb,math::matrix<double>& PreNEQUAmb, math::matrix<double> CurNEQ,math::matrix<double> CurNEQU,DdAmbInfo PreAmbInfo,DdAmbInfo CurAmbInfo,DdData curData,DdData preData);
//extern void SolveNEQCholesky2(math::matrix<double> Nk,math::matrix<double> Nks,math::matrix<double> Ns,math::matrix<double>& Qss,math::matrix<double>& U1,math::matrix<double>& U2,int ColN2,int paraNum1);
extern void SolveNEQCholesky2(math::matrix<double> Nk,math::matrix<double> Nks,math::matrix<double> Ns,math::matrix<double>& Qss,math::matrix<double>& Qks,math::matrix<double>& Qkk, math::matrix<double>& U1,math::matrix<double>& U2,int ColN2,int paraNum1);
extern math::matrix<double> FixAmbQXWZ(math::matrix<double>Qaa,math::matrix<double>Qba,math::matrix<double>&Qbb, math::matrix<double>& a_hat,math::matrix<double>&b_hat,double threshold,double& ratio);
extern void editVcm(math::matrix<double>&Qaa,math::matrix<double>&Qba,math::matrix<double>&Qbb, math::matrix<double>& a_hat,math::matrix<double>&b_hat);
//extern void QXWZPartialARwithUpdate(math::matrix<double>Qamb,math::matrix<double>Qb_amb,math::matrix<double>Qbb, math::matrix<double>a_hat,math::matrix<double> b_hat,math::matrix<double>&a_hat_update, math::matrix<double>&a_check,math::matrix<double>&b_hat_update, DdData curdata,double& ratio,double eleThreshold,double ratioThrs);
extern void QXWZPartialARwithUpdate(math::matrix<double>Qamb,math::matrix<double>Qb_amb,math::matrix<double>Qbb, math::matrix<double>a_hat,math::matrix<double> b_hat,math::matrix<double>&a_hat_update, math::matrix<double>&a_check,math::matrix<double>&b_hat_update, DdData curdata,double& ratio,double eleThreshold,double ratioThrs, DdAmbInfo& curAmb);
extern void TransferFixAmb(int typeNo,DdAmbInfo& curAmb,DdData curdata,math::matrix<double>afix,int numFix,int* indexToBeFixed);
extern void checkAR(DdAmbInfo ambinfo,DdData curdata,DdCtrl ddctrl,DdObsInfo obsinfo);
extern void SolveNormalEquationCholesky2(math::matrix<double>& N1,math::matrix<double>& N12,math::matrix<double>& N2,math::matrix<double>& U1);
extern void IonoFreeData(DdData curdata,DdData& comData,DdObsInfo* obsinfo,DdObsInfo* obsinfo_temp);