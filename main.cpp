#include "stdafx.h"//changed
#include "atlstr.h"
#include <stdlib.h>
#include <stdio.h>
#include "RTCM_NMEA.h"
#include <fstream>//file
#include <iomanip>
#include <iostream>
#include <windows.h>
#include<math.h>//calc
#include "matrix.h"
#include <winsock2.h>//net
#pragma comment(lib,"WS2_32.lib")//net
#include "Rinex.h"
#include "Position.h"
#include "solution.h"
#include <time.h>

void testLambda(int dimAmb)
{

	ifstream infile;
	infile.open("C:\\Users\\503\\Desktop\\te.txt",ios::in);
	if (!infile)
	{
		cerr<<"open file err!"<<endl;
		exit(1);
	}
	math::matrix<double>Ms(dimAmb+1,dimAmb),xt,Q,F(dimAmb,2);
	for (int m=0;m<dimAmb+1;m++)
	{
		for (int n=0;n<dimAmb;n++)
		{
			infile>>Ms(m,n);
		}
	}
	xt=~GetBlockMat(Ms,1,1,1,dimAmb,2);
	Q=GetBlockMat(Ms,2,dimAmb+1,1,dimAmb,2);
	Ambiguity amb;
	double s[2];
	amb.Lambda(dimAmb,2,xt,Q,F,s);
	cout<<F;
	cout<<s[1]/s[0]<<endl;
}

void testDll()
{
	//
	//-----------------DLL GNSSTB------------
	/*xyz2neu neu;
	HINSTANCE hInst=LoadLibrary(_T("GNSSTB.dll"));
	xyz2neu func_neu;
	func_neu=(xyz2neu)GetProcAddress(hInst, "XYZ2NEU");
	if(func_neu==NULL)
	{
		FreeLibrary(hInst);
	}
	double XYZ[3]={-2870673.6948,4623083.5405,3315621.9756};
	double B0=31;
	double L0=121.31;
	double NEU[3];
	func_neu(XYZ,B0,L0,NEU);
	FreeLibrary(hInst);
*/

}



void testMat()
{
	int k=FreqSys(1,1)/FreqSys(1,2);
	//cout<<k<<endl;
	math::matrix<double> yu(4,4);
	for (int i=0;i<4;i++)
	{
		yu(0,i)=i+1;
		yu(i,0)=i+1;
		if (i==1)
		{
			yu(1,1)=5;yu(1,i+1)=9;yu(1,i+2)=10;
			yu(1+i,1)=9;yu(i+2,1)=10;
		}
		yu(2,2)=22;yu(2,3)=20;
		yu(3,2)=20;
		yu(3,3)=37;

	}
	//cout<<yu<<endl;

	math::matrix<double> yk(4,4);
	math::matrix<double> ol(4,1);
	math::matrix<double>* yyy=new math::matrix<double>[2];

	for (int i=0;i<4;i++)
	{
		for (int j=0;j<1;j++)
		{
			ol(i,j)=(i+1.2)+j*0.2;
		}
		//math::matrix<double> ols(4,1);
	}
	yyy[0]=ol;
	math::matrix<double>I(2,2);
	//cout<<yu<<endl;
	//cout<<InsertZeroCol(yu,0,1);
	for (int i=0;i<2;i++)
	{
		I(i,i)=12.0;
	}
	//cout<<Kronecker(yu,I,2);
}

void testneu()
{
	double coor[3],co[3];
	coor[0]=-2500946.203;
	coor[1]=-4670472.824;
	coor[2]=3539500.578;
	co[0]=-20507778.8799;
	co[1]=-9658597.467;
	co[2]=13900102.049;
	double ele,azi;
	XYZ2RAH(coor,1,co,ele,azi);
	cout<<ele<<"  "<<azi<<endl;
}

void SingleSys(Position spp,DdData& dddataPre,DdData& dddataCurr,SdData& lastSdData,ObsEpochData& roverData1,ObsEpochData& baseData1,
						SppCtrl sppctrl,DdCtrl ddctrl, DdAmbInfo& ambinfo,DdAmbInfo& preambinfo,SppInfo sppinfo,SppInfo baseInfo,DdObsInfo ddobsinfo1,
						BroadEphHeader brdephheader,BroadEphData* broadephdata,BroadEphDataGlo* gloeph,int nSate,int nSateGlo,
						ObsHeader obsHeader1,ObsEpochData	epochData1,ObsHeader obsHeader2,ObsEpochData	epochData2,
						double* baseCrd,double* roverCrd,int& nEpoch, int sysidUse,math::matrix<double>* NormEqu)
{
	lastSdData.ZeroElem();roverData1.ZeroElem();	baseData1.ZeroElem();  sppinfo.ZeroElem(); baseInfo.ZeroElem();
	if (sysidUse==1 || sysidUse==3 ||sysidUse==5)
	{
		spp.StandPosition(brdephheader,broadephdata,epochData1,nSate,sppctrl,sppinfo,roverData1);
		spp.baseStn(baseInfo,broadephdata,epochData2,baseData1,nSate,sppctrl);
	}
	else if(sysidUse==2)
	{
		//spp.StandPosition(brdephheader,gloeph,epochData1,nSateGlo,sppctrl,sppinfoGlo,roverData);//Glo
		//spp.baseStn(baseInfoGlo,gloeph,epochData2,baseData,nSateGlo,sppctrl);//glo

	}


	if(Norm(roverCrd,3)>0.0) PtrEqual(roverCrd,sppinfo.recPos,3);
	int refPrn=0;
	ddobsinfo1.ZeroEle();
	int tprn=0;
	refPrn=spp.SelectRefSate(baseInfo,sppinfo,5.0,baseData1,roverData1,lastSdData,ddobsinfo1,0,tprn);
	int refpos;
	refpos=GetPos(lastSdData.prn,refPrn,lastSdData.satnum);

	dddataCurr.ZeroElem(); ambinfo.ZeroElem();
	spp.DoubleDiff(refPrn,refpos,lastSdData,dddataCurr);
	if(dddataCurr.pairNum<4) exit(1);
	//DdData ts=dddataCurr;
	/*	-----------------------------------------end data preparation part	---------------------------------------*/
	dddataCurr=spp.ComObsPhsCod(ddctrl,ddobsinfo1,ambinfo,dddataCurr);
	//if(nEpoch>1) spp.PassPreAmb(preambinfo,ambinfo,ddctrl.PhsTypeNo());
	/* -----------------------------------------end cycle slip---------------------------------------- */

	int obsNum	=ddobsinfo1.SumCod()+ddobsinfo1.SumPhs();//the number of all obs in one system at current epoch
	int phsNum	=ddobsinfo1.SumPhs();
	int ionoNum=dddataCurr.pairNum;
	int ambNum=ambinfo.SumUnfix();
	//if(nEpoch>1)ambNum=preambinfo.SumUnfix();
	if(ambNum==0) ambNum=1;
	/*	equation part	*/
	math::matrix<double> DesMatPos(obsNum,3);
	math::matrix<double> DesMatAmb(obsNum,ambNum);
	math::matrix<double> Weight(obsNum,obsNum);
	math::matrix<double> L(obsNum,1);
	math::matrix<double> DesMatTrop(obsNum,1);
	math::matrix<double> DesMatIono(obsNum,obsNum);

	spp.FormDdErrEq(DesMatPos,DesMatTrop,DesMatIono,DesMatAmb,L,dddataCurr,ddctrl,ambinfo,ddobsinfo1);
	Weight=spp.FormWeightVc(ddctrl,dddataCurr,ddobsinfo1);
	NormEqu[0]=~DesMatPos*Weight*(DesMatPos);		NormEqu[1]=~DesMatPos*Weight*(DesMatAmb);
	NormEqu[2]=~DesMatAmb*Weight*(DesMatAmb);	NormEqu[3]=~DesMatPos*Weight*L;
	NormEqu[4]=~DesMatAmb*Weight*L;
}

void RtkGlo(Position spp,ObsEpochData& baseData1,SppCtrl sppctrl,DdCtrl ddctrl, BroadEphHeader brdephheader,BroadEphData* broadephdata,
					BroadEphDataGlo* gloeph,int nSateGlo, ObsHeader obsHeader1,ObsEpochData	epochData1,ObsHeader obsHeader2,ObsEpochData	epochData2,
					double* baseCrd,double* roverCrd,math::matrix<double>& N11)
{
	DdData dddataCurr;SdData lastSdData; ObsEpochData roverData,baseData;
	DdObsInfo ddobsinfo;SppInfoGlo sppinfo, baseInfoGlo;DdAmbInfo ambinfo;
	spp.StandPosition(brdephheader,gloeph,epochData1,nSateGlo,sppctrl,sppinfo,roverData);
	spp.baseStn(baseInfoGlo,gloeph,epochData2,baseData,nSateGlo,sppctrl);//glo
	roverData.Cycle2MeterGlo(sppinfo.freqNum);
	baseData.Cycle2MeterGlo(baseInfoGlo.freqNum);
	
	if(Norm(roverCrd,3)>0.0) PtrEqual(roverCrd,sppinfo.recPos,3);
	
	int* dNum=new int[max(baseInfoGlo.validnum, sppinfo.validnum)];
	int tprn=0;
	int refPrn= spp.SelectRefSate(baseInfoGlo,sppinfo,5.0,baseData,roverData,lastSdData,ddobsinfo,dNum);
	int refpos=GetPos(lastSdData.prn,refPrn,lastSdData.satnum);
	
	dddataCurr.ZeroElem(); ambinfo.ZeroElem();
	spp.DoubleDiff(refPrn,refpos,lastSdData,dddataCurr);
	if(dddataCurr.pairNum<4) exit(1);
	
	/*	-----------------------------------------end data preparation part	---------------------------------------*/
//	dddataCurr=spp.ComObsPhsCod(ddctrl,ddobsinfo,ambinfo,dddataCurr,sppinfo);// units of obs : m
	/* -----------------------------------------end cycle slip---------------------------------------- */


	int obsNum	=ddobsinfo.SumCod()+ddobsinfo.SumPhs();//the number of all obs in one system at current epoch
	int phsNum	=ddobsinfo.SumPhs();
	int ionoNum=dddataCurr.pairNum;
	int ambNum=ambinfo.SumUnfix();
	//if(nEpoch>1)ambNum=preambinfo.SumUnfix();
	if(ambNum==0) ambNum=1;
	/*	equation part	*/
	math::matrix<double> DesMatPos(obsNum,3);
	math::matrix<double> DesMatAmb(obsNum,ambNum);
	math::matrix<double> Weight(obsNum,obsNum);
	math::matrix<double> L(obsNum,1);
	math::matrix<double> DesMatTrop(obsNum,1);
	math::matrix<double> DesMatIono(obsNum,obsNum);
	spp.FormDdErrEq(DesMatPos,DesMatTrop,DesMatIono,DesMatAmb,L,dddataCurr,ddctrl,ambinfo,ddobsinfo);

}

void ProcessSingle(CString commandFile,int sysidUse,const char* outFile)
{
	ProcessRinex    readdata;
	CString CommadFile	=_T("CommandFile.txt");
	//control block
	DdData dddataPre;
	DdData dddataCurr;
	SdData lastSdData;
	ObsEpochData roverData,baseData;	//single system
	DdObsInfo ddobsinfo,preddobsinfo;
	DdAmbInfo ambinfo,preambinfo;
	int numCod=0,numPhs=0;
	math::matrix<double> prenorm[5];
	double dt;

	SppCtrl sppctrl;
	Position spp;
	SppInfo sppinfo,baseInfo;
	SppInfoGlo sppinfoGlo,baseInfoGlo;

	DdCtrl ddctrl;
	InPutFileSet  inputfile;
	double baseCrd[3],roverCrd[3];
	DdCtrl ddCtrlPtr[MAXSYS];
	ReadCommFile(ddCtrlPtr,inputfile,_T("CommandFile.txt"),baseCrd,roverCrd);
	inputfile.CheckPath();
	ddctrl=ddCtrlPtr[sysidUse-1];
	sppctrl.sysid=ddctrl.sysid;
	sppctrl.maskEle=ddctrl.maskele;
	//read broadcast eph
	CString NFileName=inputfile.fileEph[0];
	BroadEphHeader brdephheader;
	BroadEphData*		broadephdata=new BroadEphData[MAXEPHNUM];
	BroadEphDataGlo* gloeph=new BroadEphDataGlo[MAXEPHNUM_GLO];
	int nSate	=0;
	int nSateGlo	=0;
	readdata.ReadBroadEph(NFileName, brdephheader,broadephdata,gloeph, nSate,nSateGlo);
	

	//Rover station
	ObsHeader			obsHeader1;
	ObsEpochData		epochData1;
	CString OFileName1	=inputfile.fileRover[0];//E:\\数据\\SHCORS\\066(0307)\\SHCH\\0001066J00.14O""0002066C00.14O"E:\\CUTB2014120124.obs
	CStdioFile	Ofile1;
	Ofile1.Open(OFileName1,CFile::modeRead);
	CString	line1;
	readdata.ReadObsHeader(Ofile1,line1,obsHeader1);

	//Base station
	for (int i=0;i<3;i++) baseInfo.recPos[i]=baseCrd[i];

	ObsHeader			obsHeader2;
	ObsEpochData		epochData2;

	CString OFileName2	=inputfile.fileBase[0];
	CStdioFile	Ofile2;
	Ofile2.Open(OFileName2,CFile::modeRead);
	CString	line2;
	readdata.ReadObsHeader(Ofile2,line2,obsHeader2);
	int	nEpoch1=0,nEpoch2=0;

	int	eventFlag1=0,eventFlag2=0;
	int nEpoch=0;
	//output to file
	fstream fout;
	fout.open(outFile,ios::out);
	while (eventFlag1!=10 && eventFlag2!=10)
	{
		epochData1.ZeroElem();	epochData2.ZeroElem();
		epochData1=readdata.GetEpochData(Ofile1,line1,obsHeader1,nEpoch1,eventFlag1);
		epochData1.Sort();
		epochData2=readdata.GetEpochData(Ofile2,line2,obsHeader2,nEpoch2,eventFlag2);
		epochData2.Sort();
		if(epochData1.sateNum>=5&&ddctrl.CodTypeNo()>1)	epochData1=epochData1.AutoShrink(2);
		if(epochData2.sateNum>=5&&ddctrl.CodTypeNo()>1)	epochData2=epochData2.AutoShrink(2);
		dt=(epochData1.week-epochData2.week)*7*86400.0+(epochData1.sec-epochData2.sec);
		if (1)//mark  fabs(dt)<0.5
		{
			nEpoch++;
			lastSdData.ZeroElem();
			roverData.ZeroElem();	baseData.ZeroElem();  sppinfo.ZeroElem(); baseInfo.ZeroElem();
			if (sysidUse==1 || sysidUse==3 ||sysidUse==5)
			{
				spp.StandPosition(brdephheader,broadephdata,epochData1,nSate,sppctrl,sppinfo,roverData);
				spp.baseStn(baseInfo,broadephdata,epochData2,baseData,nSate,sppctrl);
			}
			else if(sysidUse==2)
			{
				spp.StandPosition(brdephheader,gloeph,epochData1,nSateGlo,sppctrl,sppinfoGlo,roverData);//Glo
				spp.baseStn(baseInfoGlo,gloeph,epochData2,baseData,nSateGlo,sppctrl);//glo

			}
			int t0,t1,t2,t3;
			t0=((int)roverData.sec%86400);
			t1=t0/3600; //hour
			t2=(t0-t1*3600)/60;  //hour
			t3=(t0-t1*3600-t2*60)%60;
#if FileOutSpp
			//fout<<"epoch:   "<<t1<<"  "<<t2<<"  "<<t3<<endl;
			//SppFileOut(fout,sppinfo);
			
#endif
			if (nEpoch>288)
			{
				nEpoch=nEpoch;
			}
			//continue;
			if(Norm(roverCrd,3)>0.0) PtrEqual(roverCrd,sppinfo.recPos,3);
			int refPrn=0;
			ddobsinfo.ZeroEle();
			int tprn=0;
			refPrn=spp.SelectRefSate(baseInfo,sppinfo,5.0,baseData,roverData,lastSdData,ddobsinfo,0,tprn);
			int refpos;
			refpos=GetPos(lastSdData.prn,refPrn,lastSdData.satnum);

			dddataCurr.ZeroElem(); ambinfo.ZeroElem();
			spp.DoubleDiff(refPrn,refpos,lastSdData,dddataCurr);
			if(dddataCurr.pairNum<4) exit(1);
			DdData ts=dddataCurr;
			/*	-----------------------------------------end data preparation part	---------------------------------------*/
			dddataCurr=spp.ComObsPhsCod(ddctrl,ddobsinfo,ambinfo,dddataCurr);
			if(nEpoch>1) spp.PassPreAmb(preambinfo,ambinfo,ddctrl.PhsTypeNo());
			/* -----------------------------------------end cycle slip---------------------------------------- */

			int obsNum	=ddobsinfo.SumCod()+ddobsinfo.SumPhs();//the number of all obs in one system at current epoch
			int phsNum	=ddobsinfo.SumPhs();
			int ionoNum=dddataCurr.pairNum;
			int ambNum=ambinfo.SumUnfix();
			//if(nEpoch>1)ambNum=preambinfo.SumUnfix();
			if(ambNum==0) ambNum=1;
			/*	equation part	*/
			math::matrix<double> DesMatPos(obsNum,3);
			math::matrix<double> DesMatAmb(obsNum,ambNum);
			math::matrix<double> Weight(obsNum,obsNum);
			math::matrix<double> L(obsNum,1);
			math::matrix<double> DesMatTrop(obsNum,1);
			math::matrix<double> DesMatIono(obsNum,obsNum);

			spp.FormDdErrEq(DesMatPos,DesMatTrop,DesMatIono,DesMatAmb,L,dddataCurr,ddctrl,ambinfo,ddobsinfo);
			Weight=spp.FormWeightVc(ddctrl,dddataCurr,ddobsinfo);

			/*---------------------	end	processing method-------------------	*/
			int resultflag=0; double ratio=0.0;
			if(nEpoch%1000==0) cout<<nEpoch<<endl;
			if (nEpoch>0)
			{
				double xxb[3];
				math::matrix<double>xb=SoluShortEpoch(ddctrl,DesMatPos,DesMatAmb,Weight,L,ddobsinfo,ambinfo,resultflag,ratio);
				preambinfo=ambinfo;
				for(int i=0;i<3;i++) xxb[i]=xb(i,0)+sppinfo.recPos[i];  
				fout<<resultflag<<"   "<<dddataCurr.pairNum<<"   "<<~XYZ2NEU(sppinfo.recPos,xxb);
				//cout<<~XYZ2NEU(sppinfo.recPos,xxb);
				preddobsinfo=ddobsinfo;
			}
		//	int t0,t1,t2,t3;
			//t0=((int)dddataCurr.sec%86400);
			//t1=t0/3600; //hour
			//t2=(t0-t1*3600)/60;  //hour
			//t3=(t0-t1*3600-t2*60)%60;
			//cout<<"ratio    "<<ratio<<"   "<<t1<<" "<<t2<<" "<<t3<<" "<<dddataCurr.pairNum<<"  "<<nEpoch<<endl;
			//cout<<endl;
			dddataPre		=dddataCurr;
			/*	end equation part	*/
		}
		else
		{
			if (dt>0.0) epochData2=readdata.GetEpochData(Ofile2,line2,obsHeader2,nEpoch2,eventFlag2);
			else epochData1=readdata.GetEpochData(Ofile1,line1,obsHeader1,nEpoch1,eventFlag1);
		}
		if (eventFlag1==10 || eventFlag2==10) nEpoch1++;
	}//end while

	fout.close();
	Ofile1.Close();
	Ofile2.Close();
	delete[] broadephdata,gloeph;
}

void ProcessDouble(CString commandFile,int* sysidUse)
{
	ProcessRinex    readdata;
	Position spp;
	CString CommadFile	=_T("CommandFile.txt");
	
	/*                            data                                  */
	DdCtrl ddCtrlPtr[MAXSYS]; DdData dddataPre,dddataCurr;  SdData lastSdData;
	ObsEpochData roverData1,baseData1;DdObsInfo ddobsinfo1,preddobsinfo[2];
	DdAmbInfo ambinfo[2],preambinfo[2];
	int numCod=0,numPhs=0;
	double dt;// the time diff of obs time from two files

	SppCtrl sppctrl[2];

	SppInfo sppinfo[2],baseInfo[2];
	SppInfoGlo sppinfoGlo,baseInfoGlo;

	double baseCrd[3],roverCrd[3];
	InPutFileSet  inputfile;
	ReadCommFile(ddCtrlPtr,inputfile,_T("CommandFile.txt"),baseCrd,roverCrd);
	inputfile.CheckPath();
	for(int i=0;i<2;i++) sppctrl[i].sysid=ddCtrlPtr[ sysidUse[i]-1 ].sysid;

	//read broadcast eph
	
	BroadEphHeader brdephheader;
	BroadEphData*		broadephdata=new BroadEphData[MAXEPHNUM];
	BroadEphDataGlo* gloeph=new BroadEphDataGlo[MAXEPHNUM_GLO];
	int nSate	=0,nSateGlo	=0;
	CString NFileName=inputfile.fileEph[0];
	readdata.ReadBroadEph(NFileName, brdephheader,broadephdata,gloeph, nSate,nSateGlo);
	

	//Rover station
	ObsHeader			obsHeader1,obsHeader2;
	ObsEpochData		epochData1,epochData2;
	CString OFileName1	=inputfile.fileRover[0];//E:\\数据\\SHCORS\\066(0307)\\SHCH\\0001066J00.14O""0002066C00.14O"E:\\CUTB2014120124.obs
	CString OFileName2	=inputfile.fileBase[0];
	CStdioFile	Ofile1;CStdioFile	Ofile2;CString	line1;CString	line2;

	Ofile1.Open(OFileName1,CFile::modeRead);
	readdata.ReadObsHeader(Ofile1,line1,obsHeader1);
	Ofile2.Open(OFileName2,CFile::modeRead);
	readdata.ReadObsHeader(Ofile2,line2,obsHeader2);
	//Base station
	for (int i=0;i<3;i++) baseInfo[0].recPos[i]=baseInfo[1].recPos[i]=baseCrd[i];
	for (int i=0;i<3;i++) baseInfo[0].recPos[i]=baseInfo[1].recPos[i]=baseCrd[i];
	int	nEpoch1=0,nEpoch2=0;
	int	eventFlag1=0,eventFlag2=0;
	int nEpoch=0;
	//output to file
	CString outs;
	GetOutPath(CommadFile,outs);
	char* outchar=new char[outs.GetLength()+1];
	fstream fout;
	fout.open(outchar,ios::out);
	while (eventFlag1!=10 && eventFlag2!=10)
	{
		epochData1.ZeroElem();	epochData2.ZeroElem();
		epochData1=readdata.GetEpochData(Ofile1,line1,obsHeader1,nEpoch1,eventFlag1);
		epochData1.Sort();
		epochData2=readdata.GetEpochData(Ofile2,line2,obsHeader2,nEpoch2,eventFlag2);
		epochData2.Sort();

		dt=(epochData1.week-epochData2.week)*7*86400.0+(epochData1.sec-epochData2.sec);
		if (fabs(dt)<0.5)
		{
			nEpoch++;
			math::matrix<double> NormEqu[5];
			math::matrix<double> N11,N12,N22,U1,U2;
			for (int i=0;i<2;i++)
			{
				SingleSys(spp,dddataPre,dddataCurr,lastSdData,roverData1,baseData1,sppctrl[i],ddCtrlPtr[sysidUse[i]-1],
					ambinfo[i],preambinfo[i],sppinfo[i],baseInfo[i],ddobsinfo1,brdephheader,broadephdata,gloeph,nSate,nSateGlo,
					obsHeader1,epochData1,obsHeader2,epochData2,baseCrd,roverCrd,nEpoch,sysidUse[i],NormEqu);
				
				if (i==0)
				{
					N11=NormEqu[0];N12=NormEqu[1];N22=NormEqu[2];
					U1=NormEqu[3];	U2=NormEqu[4];
				}
				if (i==1)
				{
					N11 +=NormEqu[0];N12=ConvergeMat(N12.RowNo(),N12,NormEqu[1]);
					N22=DiagMatSym(N22,NormEqu[2]);
					U1 +=NormEqu[3];	U2=VecMat(1,U2, NormEqu[4]);
				}
			}
			NormEqu[0]=N11;NormEqu[1]=N12;NormEqu[2]=N22;
			NormEqu[3]=U1;	NormEqu[4]=U2;
			/*---------------------	end	processing method-------------------	*/
			int resultflag=0; double ratio=0.0;
			if(nEpoch%1000==0) cout<<nEpoch<<endl;

				double xxb[3];
				math::matrix<double>xb=SoluShortEpoch(ddCtrlPtr[0],NormEqu,ambinfo[0],resultflag,ratio);
				preambinfo[0]=ambinfo[0];
				for(int i=0;i<3;i++) sppinfo[0].recPos[i]=roverCrd[i];  
				for(int i=0;i<3;i++) xxb[i]=xb(i,0)+sppinfo[0].recPos[i];  

				fout<<~XYZ2NEU(sppinfo[0].recPos,xxb);
				//cout<<~XYZ2NEU(sppinfo[0].recPos,xxb);

			int t0,t1,t2,t3;
			t0=((int)dddataCurr.sec%86400);
			t1=t0/3600; //hour
			t2=(t0-t1*3600)/60;  //hour
			t3=(t0-t1*3600-t2*60)%60;
			 if(fabs(xb(0,0))>1.0 || fabs(xb(1,0))>1.0) cout<<t1<<" "<<t2<<" "<<t3<<" "<<dddataCurr.pairNum<<"  "<<nEpoch<<endl;
			//cout<<endl;
			dddataPre		=dddataCurr;
			/*	end equation part	*/
		}
		else
		{
			if (dt>0.0) epochData2=readdata.GetEpochData(Ofile2,line2,obsHeader2,nEpoch2,eventFlag2);
			else epochData1=readdata.GetEpochData(Ofile1,line1,obsHeader1,nEpoch1,eventFlag1);
		}
		if (eventFlag1==10 || eventFlag2==10) nEpoch1++;
	}//end while

	fout.close();
	Ofile1.Close();
	Ofile2.Close();
	delete[] broadephdata,gloeph,outchar;
}


/*for QXWZ,   2015.12*/
void ProcessSingleFixPos(int sysidUse)
{
	ProcessRinex    readdata;
	CString CommadFile	=_T("CommandFile.txt");
	//control block
	DdData dddataPre;
	DdData dddataCurr;
	SdData lastSdData;
	ObsEpochData roverData,baseData;	//single system
	DdObsInfo ddobsinfo,preddobsinfo;
	DdAmbInfo ambinfo,preambinfo;
	int numCod=0,numPhs=0;
	math::matrix<double> prenorm[5];
	double dt;

	SppCtrl sppctrl;
	Position spp;
	SppInfo sppinfo,baseInfo;
	SppInfoGlo sppinfoGlo,baseInfoGlo;

	DdCtrl ddctrl;
	InPutFileSet  inputfile;
	double baseCrd[3],roverCrd[3];
	DdCtrl ddCtrlPtr[MAXSYS];
	ReadCommFile(ddCtrlPtr,inputfile,_T("CommandFile.txt"),baseCrd,roverCrd);
	inputfile.CheckPath();
	ddctrl=ddCtrlPtr[sysidUse-1];
	sppctrl.sysid=ddctrl.sysid;

	//read broadcast eph
	CString NFileName=inputfile.fileEph[0];
	BroadEphHeader brdephheader;
	BroadEphData*		broadephdata=new BroadEphData[MAXEPHNUM];
	BroadEphDataGlo* gloeph=new BroadEphDataGlo[MAXEPHNUM_GLO];
	int nSate	=0;
	int nSateGlo	=0;
	readdata.ReadBroadEph(NFileName, brdephheader,broadephdata,gloeph, nSate,nSateGlo);
	
	

	//Rover station
	ObsHeader			obsHeader1;
	ObsEpochData		epochData1;
	CString OFileName1	=inputfile.fileRover[0];
	CStdioFile	Ofile1;
	Ofile1.Open(OFileName1,CFile::modeRead);
	CString	line1;
	readdata.ReadObsHeader(Ofile1,line1,obsHeader1);

	//Base station
	for (int i=0;i<3;i++) baseInfo.recPos[i]=baseCrd[i];

	ObsHeader			obsHeader2;
	ObsEpochData		epochData2;

	CString OFileName2	=inputfile.fileBase[0];
	CStdioFile	Ofile2;
	Ofile2.Open(OFileName2,CFile::modeRead);
	CString	line2;
	readdata.ReadObsHeader(Ofile2,line2,obsHeader2);
	int	nEpoch1=0,nEpoch2=0;

	int	eventFlag1=0,eventFlag2=0;
	int nEpoch=0;
	//output to file
	fstream fout;
	fout.open("QXWZ.txt",ios::out);
	math::matrix<double>Npre_amb,Upre_amb;

	int cycleslip_index_L1[MAXNUMSATE];
	int cycleslip_index_L2[MAXNUMSATE];

	double ptrMapPre[MAXNUMSATE];
	double ptrMapCur[MAXNUMSATE];

	int tprn=0;
	math::matrix<double> PreVcmAmb;
	math::matrix<double> PreSoluAmb;
	math::matrix<double>* CurVcm=new math::matrix<double>[3];
	math::matrix<double>* CurSolu=new math::matrix<double>[2];

	math::matrix<double>Nneq;
	math::matrix<double>Uneq;
	while (eventFlag1!=10 && eventFlag2!=10)
	{
		epochData1.ZeroElem();	epochData2.ZeroElem();
		epochData1=readdata.GetEpochData(Ofile1,line1,obsHeader1,nEpoch1,eventFlag1);
		epochData1.Sort();
		epochData2=readdata.GetEpochData(Ofile2,line2,obsHeader2,nEpoch2,eventFlag2);
		epochData2.Sort();

		dt=(epochData1.week-epochData2.week)*7*86400.0+(epochData1.sec-epochData2.sec);
		if (epochData1.sec<260112.0)
		{
			//continue;
		}
		if (fabs(dt)<0.5)
		{
			clock_t start,finish;
			double dura=0.0;
			nEpoch++;
			start=clock();
			lastSdData.ZeroElem();
			roverData.ZeroElem();	baseData.ZeroElem();  sppinfo.ZeroElem(); baseInfo.ZeroElem();
			if (sysidUse==1 || sysidUse==3 ||sysidUse==5)
			{
				spp.StandPosition(brdephheader,broadephdata,epochData1,nSate,sppctrl,sppinfo,roverData);
				spp.baseStn(baseInfo,broadephdata,epochData2,baseData,nSate,sppctrl);
			}
			else if(sysidUse==2)
			{
				spp.StandPosition(brdephheader,gloeph,epochData1,nSateGlo,sppctrl,sppinfoGlo,roverData);//Glo
				spp.baseStn(baseInfoGlo,gloeph,epochData2,baseData,nSateGlo,sppctrl);//glo

			}

			

			if(Norm(roverCrd,3)>0.0) PtrEqual(roverCrd,sppinfo.recPos,3);
			int refPrn=0;
			ddobsinfo.ZeroEle();
			
			refPrn=spp.SelectRefSate(baseInfo,sppinfo,5.0,baseData,roverData,lastSdData,ddobsinfo,1,tprn);
			int refpos;
			refpos=GetPos(lastSdData.prn,refPrn,lastSdData.satnum);

			dddataCurr.ZeroElem(); ambinfo.ZeroElem();
			spp.DoubleDiff(refPrn,refpos,lastSdData,dddataCurr);
			if(dddataCurr.pairNum<4) continue;//exit(1)
			DdData ts=dddataCurr;
			/*	-----------------------------------------end data preparation part	---------------------------------------*/
			dddataCurr=spp.ComObsPhsCod(ddctrl,ddobsinfo,ambinfo,dddataCurr);
			//if(nEpoch>1) spp.PassPreAmb(preambinfo,ambinfo,ddctrl.PhsTypeNo());
			/* -----------------------------------------end cycle slip---------------------------------------- */
			for (int s=0;s<dddataCurr.pairNum;s++) dddataCurr.tropCor[s]=ts.tropCor[s];
			int obsNum	=ddobsinfo.SumCod()+ddobsinfo.SumPhs();//the number of all obs in one system at current epoch
			int phsNum	=ddobsinfo.SumPhs();
			int ionoNum=dddataCurr.pairNum;
			int ambNum=ambinfo.SumUnfix();
			//if(nEpoch>1)ambNum=preambinfo.SumUnfix();
			if(ambNum==0) ambNum=1;
			/*	equation part	*/
			math::matrix<double> DesMatPos(obsNum,3);
			math::matrix<double> DesMatAmb(obsNum,ambNum);
			math::matrix<double> Weight(obsNum,obsNum);
			math::matrix<double> L(obsNum,1);
			math::matrix<double> DesMatTrop(obsNum,1);
			math::matrix<double> DesMatIono(obsNum,dddataCurr.pairNum);

			spp.FormDdErrEq(DesMatPos,DesMatTrop,DesMatIono,DesMatAmb,L,dddataCurr,ddctrl,ambinfo,ddobsinfo);
			math::matrix<double> DesMatTrop_Hour(obsNum,1);
			double firstTime=0.0,curTime=0.0;
			spp.FormDesMatTrop_ReNew(ddctrl, nEpoch, firstTime, curTime,  roverData,dddataCurr, dddataPre,ptrMapCur, ptrMapPre,DesMatTrop, DesMatTrop_Hour,3600.0);
			Weight=spp.FormWeightVc(ddctrl,dddataCurr,ddobsinfo);
			math::matrix<double>DesAtmos=DesMatIono;
			DesMatAmb=ConvergeMat(DesMatTrop_Hour.RowNo(),DesMatTrop_Hour,DesMatAmb);
			math::matrix<double>Natmos=~DesMatIono*Weight*DesMatIono;
			math::matrix<double>Namb=~DesMatAmb*Weight*DesMatAmb;
			math::matrix<double>Natsam=~DesMatIono*Weight*DesMatAmb;
			math::matrix<double>Uatmos=~DesMatIono*Weight*L;
			math::matrix<double>Ua=~DesMatAmb*Weight*L;
			
			math::matrix<double> Uaaa,Uatat,Qaaa;

			math::matrix<double>Qamb,Qbb,Qb_amb;
			if (nEpoch==1)
			{
				Nneq=Namb;
				Uneq=Ua;
				Uaaa=Ua;
				Uatat=Uatmos;
				SolveNEQCholesky2(Natmos,Natsam,Nneq,Qamb,Qb_amb,Qbb,Uatat,Uaaa,Nneq.RowNo(),Natmos.RowNo());

				/*reduction of atmos parameters*/
				Nneq=Namb-~Natsam*CholeskyInv(Natmos)*Natsam;
				Uneq=Ua-~Natsam*CholeskyInv(Natmos)*Uatmos;
				
			}
			
			/****************************method 1**********************************↓*/
			/*the LSE of current epoch*/
			math::matrix<double>invNatmos=CholeskyInv(Natmos);
			math::matrix<double>Qaa=CholeskyInv( Namb-~Natsam*invNatmos*Natsam );
			math::matrix<double>Qatsam=-invNatmos*Natsam*Qaa;
			math::matrix<double>Qats=CholeskyInv(Natmos-Natsam*CholeskyInv(Namb)*~Natsam);
			CurSolu[0]=Qats*Uatmos+Qatsam*Ua; CurSolu[1]=~Qatsam*Uatmos+Qaa*Ua;
			CurVcm[0]=Qats;  CurVcm[1]=Qatsam;	CurVcm[2]=Qaa;
			/****************************method 1******************* ***************↑*/

			if (nEpoch>1)
			{
				/*method 1: update the LSE solution with the previous ambiguity constraint*/
				SuperPosLonBWithFixedCrd(ddctrl,PreVcmAmb,PreSoluAmb,CurVcm,CurSolu,preambinfo,ambinfo,dddataCurr,dddataPre);

				/*method 2: Normal equ*/
				SuperPosLonBWithFixedCrdNEQ(ddctrl,Nneq,Uneq,Namb,Ua,preambinfo,ambinfo,dddataCurr,dddataPre);
				Uaaa= Uneq;
				Uatat=Uatmos;
				SolveNEQCholesky2(Natmos,Natsam,Nneq,Qamb,Qb_amb,Qbb,Uatat,Uaaa,Nneq.RowNo(),Natmos.RowNo());

				/*now, Uatat is the solu of ionosphere, Uaaa is the solu of trop and amb*/
				/*reduction of atmos parameters*/
				Nneq-=~Natsam*CholeskyInv(Natmos)*Natsam;
				Uneq-=~Natsam*CholeskyInv(Natmos)*Uatmos;
			}


			/*****************  try to fix the ambiguity******************/
			math::matrix<double>a_hat=Uaaa,b_hat=Uatat;
			double ratio=0.0;
			//editVcm(Qamb,Qb_amb,Qbb,a_hat,b_hat);


			math::matrix<double>a_check;
				//a_check=FixAmbQXWZ(Qamb,Qb_amb,Qbb,a_hat,b_hat,2.0,ratio);/*fix all*/
			
			math::matrix<double>a_hat_update,b_hat_update;
			double ratiothrs=1.5;
			double elemaskPar=35.0;
			DdAmbInfo ambinfo_temp=ambinfo;
			QXWZPartialARwithUpdate(Qamb,Qb_amb,Qbb,a_hat,b_hat,a_hat_update,a_check,b_hat_update,
													dddataCurr,ratio,elemaskPar,ratiothrs,ambinfo_temp);
			

			if (ratio>ratiothrs)
			{
				checkAR(ambinfo_temp,dddataCurr,ddctrl,ddobsinfo);
				//cout<<b_hat-GetBlockMat(b_hat_update,1,b_hat.RowNo(),1,1,2);
				ratio=ratio;
			}
			/*pass VCM and Solu*/
			PreSoluAmb=CurSolu[1];
			PreVcmAmb=CurVcm[2];

			/*output the solu*/
			finish=clock();

			int t0,t1,t2,t3;
			t0=((int)roverData.sec%86400);
			t1=t0/3600; //hour
			t2=(t0-t1*3600)/60;  //min
			t3=(t0-t1*3600-t2*60)%60;//sec
			cout<<t1<<"  "<<t2<<"  "<<t3<<"  "<<(double)(finish-start)<<" ms  "<<ratio<<endl;
			if (ratio<ratiothrs)
			{
				FileOutTropsIonoFloat(fout,Uatat,Uaaa,dddataCurr,ptrMapPre,1);//CurSolu[0],CurSolu[1]
			} 
			else
			{
				FileOutTropsIonoFix(fout,b_hat_update,dddataCurr,ptrMapPre,1);
			}
			
			/*---------------------	end	processing method-------------------	*/
			/*pass the AmbInfo*/
			preambinfo=ambinfo;
			dddataPre=dddataCurr;
		}
		else
		{
			if (dt>0.0) epochData2=readdata.GetEpochData(Ofile2,line2,obsHeader2,nEpoch2,eventFlag2);
			else epochData1=readdata.GetEpochData(Ofile1,line1,obsHeader1,nEpoch1,eventFlag1);
		}
		if (eventFlag1==10 || eventFlag2==10) nEpoch1++;
	}//end while

	fout.close();
	Ofile1.Close();
	Ofile2.Close();
	delete[] broadephdata,gloeph;
}



void main()
{
	CString CommadFile	=_T("CommandFile.txt");
	char* outfile="BSP.txt";
	ProcessSingle( CommadFile, 1,outfile);
	//ProcessSingleFixPos(1);
	

}




