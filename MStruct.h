//#pragma once
//#include "stdafx.h"//changed
//#include "VariableDef.h"
//
//
//
//
///*-----------------------------------------------------------
//This h-file includes the Structs/Classes of rinex 
//	 1.header(obs and ephemeris files) 
//		1.1 ephemeris (broadcast and precise)
//			
//	    1.2 observable
//	 2. obs
//
//Also, it contains the reading class of all files above
//Maybe it needs to seperate this file.
//
//
//Watch out that the reading Functions/Classes are based on
//	MFC-Class and CPP-Class
//
//
//
//	all "valid" ,1 is OK. but eph sathealth 0 is ok;
//
//
//Reference:
//	1. "E:\北斗\Rinex_V2.12版本数据说明.txt"
//	2. "E:\北斗\rinex301.pdf"
//
//Define:
//	
//
//					GPS			Glonass				Galileo			SBAS				BDS				Mixed
//
//  prn(int)	1-50			51-100				101-150		151-200		201-250		---
//
//sysid(int)		1				2						3					4					5				6
//
//	 -------------------------------------------------------------- */
//
//
//
//	//-------------broadcast ephemeris------------------
//			// without GLO
//
//class tests
//{
//public:
//	tests();
//	tests(int a);
//	tests & operator=(const tests& m);
//	~tests();
//public:
//	int* alpha;
//	double* beta;
//	int num;
//};
//
//class BroadEphIonoCorr{
//	public:
//		BroadEphIonoCorr();
//		BroadEphIonoCorr& operator=(const BroadEphIonoCorr& m);
//public:
//			int				validA;
//			int				validB;
//			double			alpha1;
//			double			alpha2;
//			double			alpha3;
//			double			alpha4;
//			double			beta1;
//			double			beta2;
//			double			beta3;
//			double			beta4;
//};
//
//class BroadTimeCorr{	
//		// GPUT =GPS to UTC a0,a1
//		// GLUT =GLO to UTC a0,a0=TauC
//	public:
//		BroadTimeCorr();
//		public:
//			int				valid;
//			int				transTimeType;// 0 : no time corr in the eph-file 
//			double			a0;  //	sec
//			double			a1; //	sec/sec
//			double			refTime;
//			int				refWeek;
//	};
//
//class BroadEphHeader{
//	public:
//		BroadEphHeader();
//		public:
//			double							version;
//			int								sysid;
//			BroadEphIonoCorr		ionoCorrGPS;
//			BroadEphIonoCorr		ionoCorrBDS;
//			BroadEphIonoCorr		ionoCorrGAL;
//			BroadTimeCorr				GAUT;
//			BroadTimeCorr				GPUT;
//			BroadTimeCorr				SBUT;
//			BroadTimeCorr				GLUT;
//			BroadTimeCorr				GPGA;
//			BroadTimeCorr				GLGP;
//			//add     BDS
//			double							leapSecond;
//	};
//
//
//class BroadEphData
//{		//2.x	and 3.x      32 parameters,but 18 is necessary;
//		public:
//			BroadEphData();
//		public:
//			int				prn;
//			double			toc;						
//			double			satClkBs;
//			double			satClkDrt;
//			double			satClkDrtRe;
//			double			iode;								//1
//			double			crs;
//			double			deltan;
//			double			m0;
//			double			cuc;								//2
//			double			e;	//Eccentricity
//			double			cus;
//			double			sqA;
//			double			toe; //sec of GPS week //3
//			double			cic;
//			double			Omega;
//			double			cis;
//			double			i0;									//4
//			double			crc;
//			double			omega;
//			double			OmegaDot;
//			double			idot;								//5
//			double			codeL2;
//			double			gpsw;//with toe
//			double			flagL2;
//			double			satAccuracy;					//6
//			double			satHealth;//0 is ok 
//			double			tgd;
//			double			iodClk;
//			double			transTime;						//7
//			double			fitInterval;
//			double			spare1;
//			double			spare2;
//
//	};
//
//class BroadEphDataGlo
//{
//public:
//	BroadEphDataGlo();
//	~BroadEphDataGlo();
//public:
//	int		prn;
//	double	toc;
//	double toe;
//	double	NegativeTauN;
//	double PositiveGammaN;
//	double	UTCSec;
//	double*	Pos;
//	double*	Vel;
//	double*	Acc;
//	int		health;
//	int		freqNum;
//	int		Age;
//
//
//
//};
//
//	//-----------Observable      header---------------------------
//
//class SingleObsType
//{	//2.x and 3.x
//	public:
//			SingleObsType();
//	public:
//			int				sysid;		//for 3.x. single system id  //useless
//			CString			obsCode;
//			int				freCode;
//			CString			channel;//for 3.x
//	
//};
//class ObsType2
//{
//		public:
//			ObsType2();
//		public:
//			int										sysid;
//			int										sysNum;
//			int										typeNum;
//			SingleObsType						allType[MAXTYPE2];
//};
//
//class ObsType3
//{
//		public:
//			ObsType3();
//		public:
//			int						sysid;
//			int						typeNum;
//			SingleObsType  	allType[MAXTYPE3];
//};
//
//
//class ObsHeader
//{
//		public:
//			ObsHeader();
//		public:
//			double			version;
//			int				sysid;	//for file
//			double			appPos[3];
//			double			HEN[3];
//			double			interval;
//			ObsType2		obsType2;
//			ObsType3		obsType3[5];
//};
//
//
//
//	
//
////---------------time---------------------------
//class SecWeek
//{//system time, Second and week
//public:
//	SecWeek();
//public:
//	double sec;
//	double week;
//};
//class YMDHMS
//{//system time year mon day hour min sec 
//public:
//	YMDHMS();
//public:
//	int year,mon,day,hour,min;
//	double sec;
//};
//
//class Gtime			//GNSS time class
//{
//public:
//		void weeksec(int& Week,double& Second,YMDHMS gnsstime,int sysid);
//		void ymdhms(int Week,double Sec,YMDHMS& gnsstime);
//		void doy(int Year,int Month,int Day,int& YearofDOY,int& DayofYear);
//};
//	//-----------Observable      data---------------------------
//
//
//class ObsDataRecord
//{//observation data record  /  single satellite for 2.x and 3.x in the same struct  
//public:
//		ObsDataRecord();
//		ObsDataRecord& operator=(const ObsDataRecord& m);
//public:
//		int		sysid;			
//		int		PRN;
//		int		LLI[NFREQ+NEXOBS];
//		int		SNR[NFREQ+NEXOBS];
//		double	Phase[NFREQ+NEXOBS];
//		double	PsRange[NFREQ+NEXOBS];
//		double	Doppler[NFREQ+NEXOBS];
//		/*Define:   to be modified
//					pos		1				2				3
//					freq		L1				L2				L5
//								G1			G2			G3
//								E1=L1		E5b=B2		E5a=L5		ext:/E5a+b	E6			
//								B1				B2				B3				B2<B3
//		
//		In rinex 2.x		
//							the data from Chinese manufacturer is special :(BDS L1=B1 ; L2=B2)
//		
//		In rinex 3.x				
//								1   L1/G1/E2-L1-E1/B1
//								2	 L2/G2
//								5	 L5	 /E5a
//								6			 /E6			/B3
//								7			 /E5b		/B2
//								8			 /E5a+b           
//
//*/
//};
//
//class ObsEpochData{
//public:
//	ObsEpochData();
//	//ObsEpochData(int num);
//	ObsEpochData& operator=(const ObsEpochData& m);
//	//~ObsEpochData();
//public:
//	int						sateNum;
//	int						flag;
//	//YMDHMS				obsTime;//receiver sampling time
//	int						week;
//	double					sec;
//	ObsDataRecord 	obsdatarecord[MAXNUMSATE];		
//
//};
//
///*
// 
//include the flag of ionosphere dalay  	ionoflag	0=IF(defalut)  1=Klobuchar  2=fixed 3=float
//	the mask elevation (default:8.0 deg)
//	sysnum  default:1
//	sysid		 default:1
//	freqindex		=0,1,2  the index of frequency  , this is meaningful if ionoflag !=0
//	IFIndex  the index of dual-freq , meaningful if ionoflag=0;
//	*/
//class SppCtrl
//{
//public:
//	SppCtrl();
//public:
//	int		ionoFlag;//
//	int		freqIndex;
//	int		IFIndex[2];
//	int		sysNum;
//	int		sysid;
//	double	maskEle;
//};
//
//class SatePos
//{
//public:
//	SatePos();
//	SatePos& operator =(const SatePos & m);
//
//public:
//	double sateXYZ[3];
//};
//
//class SppInfo
//{
//public:
//	//SppInfo();
//	SppInfo(int sateNum);
//	SppInfo& operator=(const SppInfo& m);
//	~SppInfo();
//public:
//	double*		recPos;
//	SatePos*		satePos;
//	SatePos*		sateVel;
//	double*		sateclkVel;
//	double*		sateclkerr;
////	double*		emiTime;
//	int*				prnList;
//	double*		ele;
//	double*		azi;
//	double*		mapWet;
//	double*		tropCorr;
//	double*		codeCorr;
//	int				validnum;
//	double			dtr;
//	double*		residual;
//};
//
//class DataRecord
//{
//	public:
//	DataRecord();
//	double Phase[NFREQ+NEXOBS];
//	double PsRange[NFREQ+NEXOBS];
//	DataRecord& operator=(const DataRecord& m);
//};
//
//class SdDataRecord
//{
//public:
//	SdDataRecord();
//public:
//	double Phase[NFREQ+NEXOBS];
//	double PsRange[NFREQ+NEXOBS];
//};
//
//class SdData
//{
//public:
//	SdData();
//	SdData(int num);
//	SdData& operator=(const SdData & m);
//	void rememo(int num);
//	~SdData();
//public:
//	int			satnum;
//	int				week;
//	double			sec;
//	int*			prn;
//	double*		refRecPos;
//	double*		rovRecPos;
//	SatePos*		satePosBase;// the satepos for base station
//	SatePos*		satePosRov;
//	double*		mapWet;
//	DataRecord* sddatarecord;
//};
//
//class DdData
//{
//public:
//	DdData(int num);
//	DdData& operator=(const DdData& m);
//	void rememo(int num);
//	~DdData();
//	
//public:
//	int				pairNum;
//	int*			rovPrn;
//	int				refPrn;
//	int				week;
//	double			sec;
//	double*		refRecPos;// the base station's pos
//	double*		rovRecPos;//the rover station's pos  unknown parameter
//	double*		refSatPos_Rov;//the ref satepos to rover station
//	double*		refSatPos_Base;//the ref satepos to base station
//	SatePos*		satePosBase;// the list of satpos to base station
//	SatePos*		satePosRov;//the list of satpos to rover station
//	DataRecord* datarecord;
//	double*		mapWet;
//};
//
//
///*
//	flag	the num of freq, 0=combination,	1=single,	2=double,	3=triple
//	coef  the coef of frequencies ,[a,b,c] 
//		flag=0:	coef stands the comb coef
//		flag=1:	first one stands the index of freq for different system eg.[2 0 0]
//		flag=2: first two stand the index of freq for ...eg[1 3 0]
//		flag=3: coef stands nothing
//	sysid
//	revised:
//		on 2015.04.13
//		add coef2 and coef3 for combination , the max combination number is 3
//		coef is changed from 1*3 to 3*3, every row stands for different combination
//		and if flag>3,the number of types(including combinations) is (flag-3)
//		coef1 and coef2 stand for the 2nd and 3rd combination coef 
//		
//		flag=1:	coef[0][0] stands the index of freq for different system eg.[2 0 0]
//		flag=2: first two(coef[0][0] and coef[1][0]) stand the index of freq for ... eg[1 3 0]
//		flag=3: coef[0][0] , coef[1][0] ,coef[2][0] are meaningful ;  IE  if flag=1~3  there is no combination
//
//		flag=4: coef[0][0:2] is meaningful
//		flag=5: coef[0:1][0:2] is meaningful
//		flag=6: 
//*/
//class DdAmbCtrl
//{
//public:
//	int flag;
//	int sysid;
//	int coef[3][3];
//	int numAmb[3];
//	
//};
//
///*
//ddctrl		
//	tropFlag	0=fixed;	1=est;	2=reserved;
//	ionoFlag	0=comb;	1=IF;		2=float;	3=fixed; 
//	ambFlag	0=float;	1=fix;	2=par;
//	pseudoFlag the number of pseudoranges used in resolution 0=combine, 1=single, 2=......... max=3	(2015.4.10)  
//	
//	pseduoCoef(3*1) is synchronized to pseduoFlag, the first one indicates the pseudorang uesd,eg [2,0,0] only L2    
//	pseduoCoef is used to determine the ionoDesMat,which is also determined by phase uesd
//	phaseFalg is removed, because the ambctrl included the phase ctrl
//	otherFlag	0=spare for the rest parameters , more to be defined
//*/
//class DdCtrl
//{
//public:
//	DdCtrl();
//public:
//	int sysid;
//	int tropFlag;
//	int tropStep;
//	int ionoFlag;
//	int ambFlag;	
//	int otherFlag;
//	int pseudoFlag;
//	int pseudoCoef[3];//for combination
//	DdAmbCtrl ddambctrl;
//	//int* phaseFlag;
//	//int	combcoef[3];
//};
//
//
//
//
//
//
//
///*	
//	contain the ambiguity information of last epoch
//	num is the pair number of satellites
//	prnList list of rover sates	num*1
//	fixFlag	0=unfixed;	1=fixed;	num*k	k=1,2,3 equals to (ddambctrl.flag-3)?:;
//	fixSolu	the value of ambiguity resolution		num*k
//	here, the list of fixFlag[0:num-1] and fixSolu[0:num-1] indicates the first frequency, fixFlag[num:num*2-1]......
//*/
//class DdAmbInfo
//{
//public:
//	DdAmbInfo();
//	DdAmbInfo(int num,int nfreq);
//	~DdAmbInfo();
//	DdAmbInfo& operator=(const DdAmbInfo& m);
//public:
//	int			pairNum;
//	int			freqNum;
//	int			refSate;
//	int*			prnList;
//	int*			fixFlag;
//	double*	fixSolu;
//};
//
///*
// *Epoch information
// *include the information of data and ambiguity
// *
// */
//class EpochInfo
//{
//public:
//	EpochInfo();
//	~EpochInfo();
//	EpochInfo& operator=(const EpochInfo& m);
//public:
//	DdAmbInfo	AmbInfo;
//	DdData			ObsDataInfo;
//};
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////// 


#pragma once
#include "stdafx.h"//changed
#include "VariableDef.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>


/*-----------------------------------------------------------
This h-file includes the Structs/Classes of rinex 
	 1.header(obs and ephemeris files) 
		1.1 ephemeris (broadcast and precise)
			
	    1.2 observable
	 2. obs

Also, it contains the reading class of all files above
Maybe it needs to seperate this file.


Watch out that the reading Functions/Classes are based on
	MFC-Class and CPP-Class



	all "valid" ,1 is OK. but eph sathealth 0 is ok;



Define:
			   如切如磋，如琢如磨
			 
				GPS			Glonass				Galileo			SBAS				BDS				Mixed

  prn(int)	1-50			51-100				101-150		151-200		201-250		---

sysid(int)		1				2						3					4					5				6

	 -------------------------------------------------------------- */



	//-------------broadcast ephemeris------------------


class tests
{
public:
	tests();
	tests(int a);
	tests & operator=(const tests& m);
	~tests();
public:
	int* alpha;
	double* beta;
	int num;
};

class BroadEphIonoCorr{
	public:
		BroadEphIonoCorr();
		BroadEphIonoCorr& operator=(const BroadEphIonoCorr& m);
public:
			int				validA;
			int				validB;
			double			alpha1;
			double			alpha2;
			double			alpha3;
			double			alpha4;
			double			beta1;
			double			beta2;
			double			beta3;
			double			beta4;
};

class BroadTimeCorr{	
		// GPUT =GPS to UTC a0,a1
		// GLUT =GLO to UTC a0,a0=TauC
	public:
		BroadTimeCorr();
		public:
			int				valid;
			int				transTimeType;// 0 : no time corr in the eph-file 
			double			a0;  //	sec
			double			a1; //	sec/sec
			double			refTime;
			int				refWeek;
	};

class BroadEphHeader{
	public:
		BroadEphHeader();
		public:
			double							version;
			int								sysid;
			BroadEphIonoCorr		ionoCorrGPS;
			BroadEphIonoCorr		ionoCorrBDS;
			BroadEphIonoCorr		ionoCorrGAL;
			BroadTimeCorr				GAUT;
			BroadTimeCorr				GPUT;
			BroadTimeCorr				SBUT;
			BroadTimeCorr				GLUT;
			BroadTimeCorr				GPGA;
			BroadTimeCorr				GLGP;
			//add     BDS
			double							leapSecond;
	};


class BroadEphData
{		//2.x	and 3.x      32 parameters,but 18 is necessary;
		public:
			BroadEphData();
		public:
			int				prn;
			double			toc;						
			double			satClkBs;
			double			satClkDrt;
			double			satClkDrtRe;
			double			iode;								//1
			double			crs;
			double			deltan;
			double			m0;
			double			cuc;								//2
			double			e;	//Eccentricity
			double			cus;
			double			sqA;
			double			toe; //sec of GPS week //3
			double			cic;
			double			Omega;
			double			cis;
			double			i0;									//4
			double			crc;
			double			omega;
			double			OmegaDot;
			double			idot;								//5
			double			codeL2;
			double			gpsw;//with toe
			double			flagL2;
			double			satAccuracy;					//6
			double			satHealth;//0 is ok 
			double			tgd;
			double			iodClk;
			double			transTime;						//7
			double			fitInterval;
			double			spare1;
			double			spare2;

};



class BroadEphDataGlo
{
public:
	BroadEphDataGlo();
public:
	int		prn;
	double	toc;
	double toe;
	int week;
	int idoe;
	double	NegativeTauN;
	double PositiveGammaN;
	double	Pos[3];
	double	Vel[3];
	double	Acc[3];
	int		health;
	int		freqNum;
	int		Age;
};

class BroadEphDataSBAS
{
public:
	BroadEphDataSBAS();

public:
	int prn;            /* satellite number */
	double t0;         /* reference epoch time (GPST)   sec*/
	double week;
	double tof;        /* time of message frame (GPST) */
	int sva;            /* SV accuracy (URA index) */
	int svh;            /* SV health (0:ok) */
	double pos[3];      /* satellite position (m) (ecef) */
	double vel[3];      /* satellite velocity (m/s) (ecef) */
	double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
	double a0,a1;     /* satellite clock-offset/drift (s,s/s) */
	double iodn;
};

	//-----------Observable      header---------------------------

class SingleObsType
{	//2.x and 3.x
	public:
			SingleObsType();
	public:
			int				sysid;		//for 3.x. single system id  //useless
			CString			obsCode;
			int				freCode;
			CString			channel;//for 3.x
	
};
class ObsType2
{
		public:
			ObsType2();
		public:
			int										sysid;
			int										sysNum;
			int										typeNum;
			SingleObsType						allType[MAXTYPE2];
};

class ObsType3
{
		public:
			ObsType3();
		public:
			int						sysid;
			int						typeNum;
			SingleObsType  	allType[MAXTYPE3];
};


class ObsHeader
{
		public:
			ObsHeader();
		public:
			double			version;
			int				sysid;	//for file
			double			appPos[3];
			double			HEN[3];
			double			interval;
			ObsType2		obsType2;
			ObsType3		obsType3[5];
};



	

//---------------time---------------------------
class SecWeek
{//system time, Second and week
public:
	SecWeek();
public:
	double sec;
	double week;
};
class YMDHMS
{//system time year mon day hour min sec 
public:
	YMDHMS();
public:
	int year,mon,day,hour,min;
	double sec;
};

class Gtime			//GNSS time class
{
public:
		void weeksec(int& Week,double& Second,YMDHMS gnsstime,int sysid);
		void ymdhms(int Week,double Sec,YMDHMS& gnsstime);
		void doy(int Year,int Month,int Day,int& YearofDOY,int& DayofYear);
};
	//-----------Observable      data---------------------------

/*
 *rec
 */
class ObsDataRecord
{//observation data record  /  single satellite for 2.x and 3.x in the same struct  
public:
		ObsDataRecord();
		ObsDataRecord& operator=(const ObsDataRecord& m);
public:
		void		ZeroElem();
		int		sysid;			
		int		PRN;
		int		LLI[NFREQ+NEXOBS];
		int		SNR[NFREQ+NEXOBS];
		int		numVadPhs;
		int		vadFlgPhs[NFREQ+NEXOBS];
		int		numVadCod;
		int		vadFlgCod[NFREQ+NEXOBS];
		double	Phase[NFREQ+NEXOBS];
		int		isCycleSlip[NFREQ+NEXOBS];
		double	PsRange[NFREQ+NEXOBS];
		double	Doppler[NFREQ+NEXOBS];
		/*Define:   to be modified
					pos		1				2				3
					freq		L1				L2				L5
								G1			G2			G3
								E1=L1		E5b=B2		E5a=L5		ext:/E5a+b	E6			
								B1				B2				B3				B2<B3
		
		In rinex 2.x		
							the data from Chinese manufacturer is special :(BDS L1=B1 ; L2=B2)
		
		In rinex 3.x				
								1   L1/G1/E2-L1-E1/B1
								2	 L2/G2
								5	 L5	 /E5a
								6			 /E6			/B3
								7			 /E5b		/B2
								8			 /E5a+b           

*/
};

class ObsEpochData{
public:
	ObsEpochData();
	//ObsEpochData(int num);
	ObsEpochData& operator=(const ObsEpochData& m);
	//~ObsEpochData();
public:
	void						ZeroElem();
	void						ZeroElemSingle(int index);
	void						Sort();
	ObsEpochData&	Shrink();
	ObsEpochData&	AutoShrink(int typeNo);
	void						Cycle2MeterGlo(int* dNum);
	int						sateNum;
	int						flag;
	//YMDHMS				obsTime;//receiver sampling time
	int						week;
	double					sec;
	ObsDataRecord 	obsdatarecord[MAXNUMSATE];		

};

/*
 
include the flag of ionosphere dalay  	ionoflag	0=IF(defalut)  1=Klobuchar  2=fixed 3=float
	the mask elevation (default:8.0 deg)
	sysnum  default:1
	sysid		 default:1
	freqindex		=0,1,2  the index of frequency  , this is meaningful if ionoflag !=0
	IFIndex  the index of dual-freq , meaningful if ionoflag=0;
	*/
class SppCtrl
{
public:
	SppCtrl();
public:
	int		ionoFlag;//
	int		freqIndex;
	int		IFIndex[2];
	int		sysNum;
	int		sysid;
	double	maskEle;
};

class SatePos
{
public:
	SatePos();
	SatePos& operator =(const SatePos & m);
	
public:
	void		ZeroElem();
	double sateXYZ[3];
};
class PrecEphData
{
public:
	int week;
	double sec;
	SatePos satePos[MAXNUMSATE];
};
class SppInfo
{
public:
	//SppInfo();
	SppInfo();
	SppInfo& operator=(const SppInfo& m);
	//~SppInfo();
public:
	void				ZeroElem();
	double			gdop;
	double			recPos[3];
	SatePos		satePos[MAXNUMSATE];
	SatePos		sateVel[MAXNUMSATE];
	double			sateclkVel[MAXNUMSATE];
	double		sateclkerr[MAXNUMSATE];
//	double*		emiTime;
	int				prnList[MAXNUMSATE];
	double		ele[MAXNUMSATE];
	double		azi[MAXNUMSATE];
	double		mapWet[MAXNUMSATE];
	double		tropCorr[MAXNUMSATE];
	double		codeCorr[MAXNUMSATE];
	int				validnum;
	double			dtr;
	double		residual[MAXNUMSATE];
};

class SppInfoGlo
{
public:
	//SppInfo();
	SppInfoGlo();
	SppInfoGlo& operator=(const SppInfoGlo& m);
	//~SppInfo();
public:
	void				ZeroElem();
	double			gdop;
	double			recPos[3];
	SatePos		satePos[MAXNUMSATE];
	SatePos		sateVel[MAXNUMSATE];
	double			sateclkVel[MAXNUMSATE];
	double		sateclkerr[MAXNUMSATE];
	//	double*		emiTime;
	int				prnList[MAXNUMSATE];
	int				freqNum[MAXNUMSATE];
	double		ele[MAXNUMSATE];
	double		azi[MAXNUMSATE];
	double		mapWet[MAXNUMSATE];
	double		tropCorr[MAXNUMSATE];
	double		codeCorr[MAXNUMSATE];
	int				validnum;
	double			dtr;
	double		residual[MAXNUMSATE];
};



class DataRecord
{
	public:
	DataRecord();
	void		ZeroElem();
	int		numVadPhs;
	int		vadFlgPhs[NFREQ+NEXOBS];
	int		numVadCod;
	int		vadFlgCod[NFREQ+NEXOBS];
	double Phase[NFREQ+NEXOBS];
	int		isCycleSlip[NFREQ+NEXOBS];
	double PsRange[NFREQ+NEXOBS];
	DataRecord& operator=(const DataRecord& m);
};

class SdDataRecord
{
public:
	SdDataRecord();
	void		ZeroElem();
public:
	double Phase[NFREQ+NEXOBS];
	double PsRange[NFREQ+NEXOBS];
};

class SdData
{
public:
	SdData();
	SdData& operator=(const SdData & m);
	void	ZeroElem();
	//void rememo(int num);
	//~SdData();
public:
	int			satnum;
	int				week;
	double			sec;
	int			prn[MAXNUMSATE];
	double		refRecPos[3];
	double		rovRecPos[3];
	SatePos		satePosBase[MAXNUMSATE];// the satepos for base station
	SatePos 		satePosRov[MAXNUMSATE];
	double 		mapWet[MAXNUMSATE];
	double		tropCor[MAXNUMSATE];
	DataRecord  sddatarecord[MAXNUMSATE];
	double			ele[MAXNUMSATE];
};

class DdData
{
public:
	DdData();
	DdData& operator=(const DdData& m);
	DdData&		Sort();
	DdData&		Check(int numFreq);
	DdData&		Shrink();
	void		ClearTail(int index);
	void		ZeroElem();
	int   SumCod();
	int   SumPhs();
	int NoCod(int ind);
	int NoPhs(int ind);
	//void rememo(int num);
	//~DdData();
	double distBaseRef();
	double distRoverRef();
	double distBaseRover(int index);
	double distRoverRover(int index);
	double	distRecSate_DD(int index);
public:
	int				pairNum;
	int			rovPrn[MAXNUMSATE];
	int				refPrn;
	int				week;
	double			sec;
	double		refRecPos[3];// the base station's pos
	double		rovRecPos[3];//the rover station's pos  unknown parameter
	double		refSatPos_Rov[3];//the ref satepos to rover station
	double		refSatPos_Base[3];//the ref satepos to base station
	SatePos		satePosBase[MAXNUMSATE];// the list of satpos to base station
	SatePos		satePosRov[MAXNUMSATE];//the list of satpos to rover station
	DataRecord datarecord[MAXNUMSATE];
	double		mapWet[MAXNUMSATE];
	double		tropCor[MAXNUMSATE];
	double		ele[MAXNUMSATE];
};


/*
	flag	the num of freq, 0=combination,	1=single,	2=double,	3=triple
	sysid
		add coef2 and coef3 for combination , the max combination number is 3
		coef is changed from 1*3 to 3*3, every row stands for different combination. note that 
		and if flag>3,the number of types(including combinations) is (flag-3)
		coef1 and coef2 stand for the 2nd and 3rd combination coef 
		
		flag=1:	coef[0][0] stands the index of freq for different system eg.[2 0 0]
		flag=2: first two(coef[0][0] and coef[1][0]) stand the index of freq for ... eg[1 3 0]
		flag=3: coef[0][0] , coef[1][0] ,coef[2][0] are meaningful ;  IE  if flag=1~3  there is no combination

		flag=4: coef[0][0:2] is meaningful
		flag=5: coef[0:1][0:2] is meaningful
		flag=6: 
	ratiothrsd  ratio threshold ,default 2.5
	parEleMask	mask elevation of par set
*/
class DdAmbCtrl
{
public:
	DdAmbCtrl();
public:
	int flag;
	int sysid;
	int coef[3][3];
	int numAmb[3];
	double ratiothrsd;
	double parEleMask;
};

/*
ddctrl		
	tropFlag	0=fixed;	1=est;	2=reserved;
	ionoFlag	0=comb;	1=IF;		2=float;	3=fixed; 
	ambFlag	0=float;	1=fix;	2=par;
	pseudoFlag the number of pseudoranges used in resolution 0=combine, 1=single, 2=......... max=3	(2015.4.10)  
	
	pseduoCoef(3*1) is synchronized to pseduoFlag, the first one indicates the pseudorang uesd,eg [2,0,0] only L2    
	pseduoCoef is used to determine the ionoDesMat,which is also determined by phase uesd
	phaseFalg is removed, because the ambctrl included the phase ctrl
	otherFlag	0=spare for the rest parameters , more to be defined

	mode	the 
*/
class DdCtrl
{
public:
	DdCtrl();
public:
	int CodTypeNo();
	int PhsTypeNo();
	int sysid;
	int mode;
	double maskele;
	int weightMode;
	double sigmaCod;
	double sigmaPhs;
	int tropFlag;
	int tropStep;
	int ionoFlag;
	int ambFlag;	
	int otherFlag;
	double thresGF;
	double thresMW;
	int pseudoFlag;
	int codtype[3];
	int phstype[3];
	int pseudoCoef[3][3];//for combination same as ambcoef
	double freqCod[3];
	double freqPhs[3];
	DdAmbCtrl ddambctrl;
	
	//int* phaseFlag;
	//int	combcoef[3];
};







/*	
	contain the ambiguity information of last epoch
	num is the pair number of satellites
	prnList list of rover sates	num*1
	fixFlag	0=unfixed;	1=fixed;	num*k	k=1,2,3 equals to (ddambctrl.flag-3)?:;
	fixSolu	the value of ambiguity resolution		num*k
	here, the list of fixFlag[0:num-1] and fixSolu[0:num-1] indicates the first frequency, fixFlag[num:num*2-1]......
*/
class DdAmbInfo
{
public:
	DdAmbInfo();
	DdAmbInfo& operator=(const DdAmbInfo& m);
	DdAmbInfo& Update(DdData curdata);
public:
	int			pairNum;
	int			freqNum;
	int			refSate;
	double		freq[NFREQ];
	int			prnList[NFREQ+NEXOBS][MAXNUMSATE];
	int			fixFlag[NFREQ+NEXOBS][MAXNUMSATE];
	double		fixSolu[NFREQ+NEXOBS][MAXNUMSATE];
	int			FindSat(int ind,int prn);
	int			NoSat(int ind);
	void			ZeroElem();
	void			GetUnfixPrnlist(int* prnptr,int flag);
	int			NoFix(int ind);
	int			NoUnfix(int ind);
	int			SumUnfix();
	int			SumNoSat();
	void			InterChange(int index,int pos1,int pos2);
};

class DdObsInfo//:class DdData
{
public:
	DdObsInfo();
	DdObsInfo& operator=(const DdObsInfo& m);
	void ZeroEle();
	void SetList(int* prnlist,int ind,int flag);//flag=0,cod; flag=1,phs
	void	GetUnfixElelist(double* eleptr,int flag,DdAmbInfo ambinfo);
	double GetEle(int prn);
	int   SumCod();
	int   SumPhs();

	int NoCod(int ind);
	int NoPhs(int ind);
	int numCod[3];
	int numPhs[3];
	double eleRefBase;
	double eleRefRov;
	double eleRovRov[MAXNUMSATE];
	double eleRovBase[MAXNUMSATE];
	int prnlistCod[3][MAXNUMSATE];
	int prnlistPhs[3][MAXNUMSATE];
};

/*
 *Epoch information
 *include the information of data and ambiguity
 *
 */
class EpochInfo
{
public:
	EpochInfo();
	//~EpochInfo();
	EpochInfo& operator=(const EpochInfo& m);
public:
	DdAmbInfo	AmbInfo;
	DdData			ObsDataInfo;
};

class AmbData
{
public:
	AmbData();
	AmbData& operator=(const AmbData& m);
	void		ZeroElem();
	int		Prn;
	double FirstObs;
	double LastObs;
	int		CurrentIndex;
	double Cycle[MAXOBSEPOCH];
	int		VadFlag[MAXOBSEPOCH];
	int		isWLInteger[MAXOBSEPOCH];
	int    	NumData();
	int		NumVadData();
	double MeanAmb();
	double CheckAmb();
	double BiasCycle();
	double StdofAver();
};

class EwlData
{
public:
	EwlData();
	EwlData& operator=(const EwlData& m);
	void		ZeroElem();
	int		Prn;
	double FirstObs;
	double LastObs;
	int		CurrentIndex;
	double EwlCycle[2][MAXOBSEPOCH];
	int		VadFlag[2][MAXOBSEPOCH];
	int    	NumData(int index);
	int		NumVadData(int index);
	double EwlMean(int index);
	double EwlCheck(int index);
	double BiasCycle(int index);
	double EwlStdofAverWl(int index);
};

class fixinfo
{
public:
	fixinfo();
	int prn;
	int valid;
	int fixepoch1;
	int fixepoch2;
	double bias1;
	double check1;
	double check2;
	double bias2;
};

class UnionEph
{
public:
	UnionEph();
	BroadEphData ephdata[500];
	int ephNum;
	BroadEphDataGlo ephdataGlo[200];
	int ephNumGlo;
};

class InPutFileSet
{
public:
	InPutFileSet();
	CString fileEph[MAXFILEEPH];
	CString fileSp3[MAXFILESP3];
	CString fileClk[MAXFILECLK];
	CString fileBase[MAXFILEBASE];
	CString fileRover[MAXFILEROVER];
	int numEph,numSp3,numClk,numBase,numRover;
	void CheckPath();
private:
   void	ModifyPath(CString& filepath);
};