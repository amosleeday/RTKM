// /*
// 
// 	this file defines all the structs of class used
// */
// 
// 
// 
// #include "stdafx.h"//changed
// #include "MStruct.h"
// 
// 
// const static double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
// const static double gst0 []={1999,8,22,0,0,0}; /* galileo system time reference */
// const static double bdt0 []={2006,1, 1,0,0,0}; /* beidou time reference */
// */
// tests::tests()
// {
// 	alpha=NULL;
// 	
// 	beta= NULL;
// }
// 
// tests::tests(int a)
// {
// 	num=a;
// 	alpha=new int[a];
// 	alpha[0]=9;
// 	beta= new double[a];
// }
// 
// tests::~tests()
// {
// 	delete[] alpha;
// 	delete[] beta;
// }
// 
// tests & tests::operator=( const tests& m ) 
// {
// 	if (m.num!=num)
// 	{
// 		delete[] alpha;
// 		delete[] beta;
// 	}
// 	num=m.num;
// 	int* alpha=new int[m.num] ;
// 	double* beta=new double	[m.num];
// 	for (int i=0;i<m.num;i++)
// 	{
// 		alpha[i]		=m.alpha[i];
// 		beta[i]		=m.beta[i];
// 	}
// 	return *this;
// }
// 
// 
// BroadEphHeader::BroadEphHeader()
// {
// 	version=0.0;
// 	sysid=0;
// 	leapSecond=0.0;
// }
// 
// /*
// 	 the parameters of ionosphere correction
// 	 default:
// 			validA=0;
// 			validB=0;
// 			1=valid; 0=invalid
// */
// BroadEphIonoCorr::BroadEphIonoCorr()
// {
// 	validA=0;
// 	validB=0;
// 	alpha1=0.0;
// 	alpha2=0.0;
// 	alpha3=0.0;
// 	alpha4=0.0;
// 	beta1=0.0;
// 	beta2=0.0;
// 	beta3=0.0;
// 	beta4=0.0;
// }
// BroadEphIonoCorr& BroadEphIonoCorr::operator=(const BroadEphIonoCorr& m)
// {
// 	validA=m.validA;
// 	validB=m.validB;
// 	alpha1=m.alpha1;
// 	alpha2=m.alpha2;
// 	alpha3=m.alpha3;
// 	alpha4=m.alpha4;
// 	beta1=m.beta1;
// 	beta2=m.beta2;
// 	beta3=m.beta3;
// 	beta4=m.beta4;
// 	return *this;
// }
// 
// 
// BroadTimeCorr::BroadTimeCorr()
// {	// GPUT =GPS to UTC a0,a1
// 			valid=0;
// 			transTimeType=0; // 0 : no time corr in the eph-file 
// 			a0=0.0;
// 			a1=0.0;
// 			refTime=0;
// 			refWeek=0;
// 	}
// 
// BroadEphData::BroadEphData()
// {
// 	prn=0;
// 	toc=0.0;						
// 	satClkBs=0.0;
// 	satClkDrt=0.0;
// 	satClkDrtRe=0.0;
// 	iode=0.0;
// 	crs=0.0;
// 	deltan=0.0;
// 	m0=0.0;
// 	cuc=0.0;
// 	e=0.0;	//Eccentricity
// 	cus=0.0;
// 	sqA=0.0;
// 	toe=0.0; //sec of GPS week
// 	cic=0.0;
// 	Omega=0.0;
// 	cis=0.0;
// 	i0=0.0;
// 	crc=0.0;
// 	omega=0.0;
// 	OmegaDot=0.0;
// 	idot=0.0;
// 	codeL2=0.0;
// 	gpsw=0.0;//with toe
// 	flagL2=0.0;
// 	satAccuracy=0.0;
// 	satHealth=0.0;
// 	tgd=0.0;
// 	iodClk=0.0;
// 	transTime=0.0;
// 	fitInterval=0.0;
// 	spare1=0.0;
// 	spare2=0.0;
// }
// 
// 
// 
// ObsHeader::ObsHeader()
// {
// 
// 			version=0.0;
// 			sysid=0;
// 			interval=1.0;
// 			for (int i=0;i<3;i++)
// 			{
// 				appPos[i]=0.0;
// 				HEN[i]=0.0;
// 			}
// }
// 
// ObsType2::ObsType2()
// {
// 
// 		sysid=0;
// 		sysNum=0;
// 		typeNum=0;
// 	}
// ObsType3::ObsType3()
// 	{
// 
// 		sysid=0;
// 		typeNum=0;
// 	}
// SingleObsType::SingleObsType()
// {	//2.x and 3.x
// 
// 	sysid=0;//for 3.x. single system id
// 	obsCode="F";// F is failure
// 	freCode=0;
// 	channel="F";//for 3.x
// 	}
// 
// SecWeek::SecWeek()
// {
// 	sec=0.0;
// 	week=0.0;
// }
// YMDHMS::YMDHMS()
// {
// 	year=0;mon=0;day=0;hour=0;min=0;
// 	sec=0.0;
// }
// 
// ObsDataRecord::ObsDataRecord()
// {
// 	sysid=0;
// 	PRN=0;
// 	for(int i=0;i<NFREQ+NEXOBS;i++)
// 	{
// 		LLI[i]=0;
// 		SNR[i]=0;
// 		Phase[i]=0.0;
// 		PsRange[i]=0.0;
// 		Doppler[i]=0.0;
// 	}
// }
// ObsDataRecord& ObsDataRecord::operator=(const ObsDataRecord& m)
// {
// 	sysid=	m.sysid;
// 	PRN=	m.PRN;
// 	for(int i=0;i<NFREQ+NEXOBS;i++)
// 	{
// 		LLI[i]			=m.LLI[i];
// 		SNR[i]		=m.SNR[i];
// 		Phase[i]		=m.Phase[i];
// 		PsRange[i]=m.PsRange[i];
// 		Doppler[i]	=m.Doppler[i];
// 	}
// 	return *this;
// }
// 
// ObsEpochData::ObsEpochData()
// {
// 	sateNum=0;
// 	flag=0;
// 	week=0;
// 	sec=0.0;
// //	obsdatarecord=NULL;
// }
// 
// 
// 
// ObsEpochData& ObsEpochData::operator=(const ObsEpochData& m)
// {
// 	//if (sateNum!=m.sateNum)
// 	//{
// 	//	delete[] obsdatarecord;
// 	//	ObsDataRecord* obsdatarecord=new ObsDataRecord[m.sateNum];
// 	//}
// 	sateNum			=m.sateNum;
// 	flag					=m.flag;
// 	week				=m.week;
// 	sec					=m.sec;
// 	for(int i=0;i<MAXNUMSATE;i++)
// 		obsdatarecord[i]	=m.obsdatarecord[i];
// 	return *this;
// }
// 
// //ObsEpochData::~ObsEpochData()
// //{
// //	delete[] obsdatarecord;
// //}
// 
// 
// 
// void Gtime::weeksec(int& Week,double& Second,YMDHMS gnsstime,int sysid)
// {
// 	//返回值为GPS周秒
// 	
// 	int y;
// 	
//    int m;
// 	if(gnsstime.mon<=2)
// 	{
// 		y=gnsstime.year-1;
// 		m=gnsstime.mon+12;
// 	}
// 	else
// 	{
// 		y=gnsstime.year;
// 		m=gnsstime.mon;
// 	}
// 	double JD;
// 	JD=((int)(365.25*y))+((int)(30.6001*(m+1.0)))+gnsstime.day+gnsstime.hour/24.0+1720981.5;
// 	int weekday=(((int)(JD+0.5))%7);//计算结果中0代表星期一，1代表星期二
// 	 weekday++;
// 	 if (weekday==7)
// 	 {
// 		 weekday=0;
// 	 }
// 	 Second=(weekday*24.0+(double)gnsstime.hour)*3600.0+(double)gnsstime.min*60.0+gnsstime.sec;
// 	 Week=(int)((JD-2444244.5)/7.0);
// }
// void Gtime::ymdhms(int Week,double Sec,YMDHMS& gnsstime)
// {
// 	 double MJD = Week*7.0+44244.0+Sec/3600.0/24.0;  //modified julian day
// 	 int a   =(int)(MJD +0.5+0.5+1.0e-10)+2400000;
// 	 double FRAC=MJD +0.5+2400000.5-a;
// 	 int b = a + 1537;
// 	 int Cc = (int)((b - 122.1)/365.25 + 1.0e-10);
// 	 int d = (int)(365.25*Cc+1.0e-16);
// 	 int Ee = (int)((b-d)/30.6001 + 1.0e-10);
// 	 gnsstime.day=b-d-(int)(30.6001*Ee);
// 	 gnsstime.mon=Ee-1-12*((int)(Ee/14.0+1.0e-10));
// 	 gnsstime.year=Cc - 4715 -(int)((7.0+gnsstime.mon)/10.0 + 1.0e-10);
// 	 double TmHour = FRAC*24.0;
// 	 gnsstime.hour= (int)(TmHour + 1.0e-10);
// 	 double TmMin=(TmHour-gnsstime.hour)*60;
// 	 gnsstime.min=(int)(TmMin+1.0e-10);
// 	 gnsstime.sec= (TmMin-gnsstime.min)*60.0;	
// }
// void Gtime::doy(int Year,int Month,int Day,int& YearofDOY,int& DayofYear)
// 	{
// 		int DayNums[12];
// 	 DayNums[0]=31;
// 	 DayNums[2]=31;
// 	 DayNums[4]=31;
// 	 DayNums[6]=31;
// 	 DayNums[7]=31;
// 	 DayNums[9]=31;
// 	 DayNums[11]=31;
// 	 DayNums[3]=30;
// 	 DayNums[5]=30;
// 	 DayNums[8]=30;
// 	 DayNums[10]=30;
// 
// 	 if ((Year%400==0) || (Year%4==0 && Year%100!=0))
// 		 DayNums[1]=29;
// 	 else
// 		 DayNums[1]=28;
// 	 DayofYear=0;
// 	 for (int i=0;i<Month-1;i++)
// 		  DayofYear+=DayNums[i];
// 	 DayofYear+=Day;
// 
// 	 YearofDOY=Year;
// 	}
// 
// 
// 
// SppCtrl::SppCtrl()
// {
// 	freqIndex	=0;
// 	ionoFlag	=0;
// 	sysNum	=1;
// 	sysid			=1;
// 	maskEle	=5.0;
// 	IFIndex[0]	=0;
// 	IFIndex[1]	=1;
// }
// SatePos::SatePos()
// {
// 	for(int i=0;i<3;i++)
// 		sateXYZ[i]	=0.0;
// }
// SatePos& SatePos::operator=(const SatePos& m)
// {
// 	for(int i=0;i<3;i++)
// 		sateXYZ[i]	=m.sateXYZ[i];
// 	return *this;
// }
// 
// 
// 
// //SppInfo::SppInfo()
// //{
// //			/*recPos=NULL;
// //			satePos=NULL;
// //			sateVel=NULL;
// //			sateclkVel=NULL;
// //			sateclkerr=NULL;
// //			emiTime=NULL;
// //			prnList=NULL;
// //			ele=NULL;
// //			azi=NULL;
// //			mapWet=NULL;
// //			tropCorr=NULL;
// //			codeCorr=NULL;*/
// //}
// SppInfo::SppInfo(int sateNum)
// {
// 	recPos		=new double[3];
// 	satePos		=new SatePos[sateNum];
// 	sateVel		=new SatePos[sateNum];
// //	emiTime	=new double[sateNum];
// 	sateclkVel =new double[sateNum];
// 	prnList		=new int[sateNum];
// 	ele			=new double[sateNum];
// 	azi			=new double[sateNum];
// 	tropCorr	=new double[sateNum];//add to h file 03.18
// 	mapWet	=new double[sateNum];
// 	codeCorr	=new double[sateNum];
// 	sateclkerr	=new double[sateNum];
// 	residual		=new double[sateNum];
// 	validnum	=sateNum;
// 	dtr			=0.0;
// 	for(int i=0;i<3;i++)
// 	{
// 		recPos[i]=	0.0;
// 	}
// 	for(int i=0;i<sateNum;i++)
// 	{
// 		sateclkVel[i]	=0.0;
// //		emiTime[i]		=	0.0;
// 		prnList[i]		=	0;
// 		ele[i]				=0.0;
// 		azi[i]				=0.0;
// 		sateclkerr[i]	=0.0;
// 		tropCorr[i]		=0.0;
// 		mapWet[i]		=0.0;
// 		codeCorr[i]	=0.0;
// 		residual[i]		=0.0;
// 
// 		for(int j=0;j<3;j++)
// 		{
// 			satePos[i].sateXYZ[j]=0.0;	
// 			sateVel[i].sateXYZ[j]=0.0;	
// 		}
// 	}
// }
// 
// SppInfo::~SppInfo()
// {
// 	delete[]		recPos;
// 	delete[]		satePos;
// 	delete[]		sateVel;
// 	delete[]		sateclkVel;
// 	delete[]		sateclkerr;
// //	delete[]		emiTime;
// 	delete[]		prnList;
// 	delete[]		ele;
// 	delete[]		azi;
// 	delete[]		mapWet;
// 	delete[]		tropCorr;
// 	delete[]		codeCorr;
// 	delete[]		residual;
// }
// 
// SppInfo& SppInfo::operator=( const SppInfo& m ) 
// {
// 	validnum	=m.validnum;
// 	dtr			=m.dtr;
// 	for(int i=0;i<3;i++)
// 	{
// 		recPos[i]=	m.recPos[i];
// 	}
// 	int sateNum=m.validnum;
// 	for(int i=0;i<sateNum;i++)
// 	{
// 		sateclkVel[i]	=m.sateclkVel[i];
// //		emiTime[i]		=	m.emiTime[i];
// 		prnList[i]		=	m.prnList[i];
// 		ele[i]				=m.ele[i];
// 		azi[i]				=m.azi[i];
// 		sateclkerr[i]	=m.sateclkerr[i];
// 		tropCorr[i]		=m.tropCorr[i]	;
// 		mapWet[i]		=m.mapWet[i];
// 		codeCorr[i]	=m.codeCorr[i];
// 		for(int j=0;j<3;j++)
// 		{
// 			satePos[sateNum].sateXYZ[j]=m.satePos[sateNum].sateXYZ[j];	
// 			sateVel[sateNum].sateXYZ[j]=m.sateVel[sateNum].sateXYZ[j];	
// 		}
// 	}
// 
// 	return *this;
// }
// 
// 
// SdDataRecord::SdDataRecord()
// {
// 	for(int i=0;i<NFREQ+NEXOBS;i++)
// 	{
// 		Phase[i]		=0.0;
// 		PsRange[i]=0.0;
// 	}
// }
// SdData::SdData(int num)
// {
// 	week	=0;
// 	sec		=0.0;
// 	satnum			=num;
// 	prn		=new int[num];
// 	satePosBase	=new SatePos[num];
// 	satePosRov	=new SatePos[num];;
// 	sddatarecord=new DataRecord[num];
// 	mapWet		=new double[num];
// 	refRecPos	=new double[3];
// 	rovRecPos=new double[3];
// 	for (int i=0;i<3;i++)
// 	{
// 		refRecPos[i]		=0.0;
// 		rovRecPos[i]		=0.0;
// 	}
// 	for (int i=0;i<num;i++)
// 	{
// 		prn[i]			=0;
// 		mapWet[i]		=0.0;
// 	}
// }
// 
// SdData::SdData()
// {
// 	week=0;
// 	sec=0.0;
// 	satnum=0;
// 	prn=NULL;
// 	satePosBase=NULL;
// 	satePosRov=NULL;
// 	sddatarecord=NULL;
// 	mapWet=NULL;
// 	refRecPos=NULL;
// 	rovRecPos=NULL;
// }
// 
// SdData& SdData::operator=(const SdData& m)
// {
// 	if (m.satnum!=satnum)
// 	{
// 		delete[] prn;
// 		delete[] satePosBase;
// 		delete[] satePosRov;
// 		delete[] sddatarecord;
// 		delete[]	mapWet;
// 		prn		=new int[m.satnum];
// 		satePosBase	=new SatePos[m.satnum];
// 		satePosRov	=new SatePos[m.satnum];
// 		sddatarecord=new DataRecord[m.satnum];
// 		mapWet		=new double[m.satnum];
// 	}
// 	
// 	week	=m.week;
// 	sec		=m.sec;
// 	satnum			=m.satnum;
// 	for (int i=0;i<3;i++)
// 	{
// 		refRecPos[i]		=m.refRecPos[i];
// 		rovRecPos[i]		=m.rovRecPos[i];
// 	}
// 	for (int i=0;i<m.satnum;i++)
// 	{
// 		satePosBase[i]	=m.satePosBase[i];
// 		satePosRov[i]	=m.satePosRov[i];
// 		sddatarecord[i]=m.sddatarecord[i];
// 		prn[i]		=m.prn[i];
// 		mapWet[i]		=m.mapWet[i];
// 	}
// 	return *this;
// }
// SdData::~SdData()
// {
// 	delete[] prn;					//prn=NULL;
// 	delete[]	mapWet;			//mapWet=NULL;
// 	delete[] satePosBase;		//satePosBase=NULL;
// 	delete[] satePosRov;		//satePosRov=NULL;
// 	delete[] refRecPos;		//refRecPos=NULL;
// 	delete[] rovRecPos;		//rovRecPos=NULL;
// 	delete[] sddatarecord;	//sddatarecord=NULL;
// 	
// }
// 
// void SdData::rememo( int num )
// {
// 	week	=0;
// 	sec		=0.0;
// 	satnum			=num;
// 
// 	delete[] prn;
// 	delete[] satePosBase;
// 	delete[] satePosRov;
// 	delete[] sddatarecord;
// 	delete[]	mapWet;
// 	prn		=new int[num];
// 	satePosBase	=new SatePos[num];
// 	satePosRov	=new SatePos[num];
// 	sddatarecord=new DataRecord[num];
// 	mapWet		=new double[num];
// 
// 	for (int i=0;i<3;i++)
// 	{
// 		refRecPos[i]		=0.0;
// 		rovRecPos[i]		=0.0;
// 	}
// 	for (int i=0;i<num;i++)
// 	{
// 		prn[i]			=0;
// 		mapWet[i]		=0.0;
// 	}
// 
// }
// 
// 
// DataRecord:: DataRecord()
// {
// 	for (int i=0;i<NFREQ+NEXOBS;i++)
// 	{
// 		Phase[i]		=0.0;
// 		PsRange[i]=0.0;
// 	}
// }
// DataRecord& DataRecord::operator=(const DataRecord& m)
// {
// 	for (int i=0;i<NFREQ+NEXOBS;i++)
// 	{
// 		Phase[i]		=m.Phase[i];
// 		PsRange[i]=m.PsRange[i];
// 	}
// 	return *this;
// }
// 
// DdData::DdData(int num)
// {
// 	pairNum			=num;
// 	rovPrn				=new int[num];
// 	refPrn				=0;
// 	week				=0;
// 	sec					=0.0;
// 	refRecPos			=new double[3];
// 	rovRecPos			=new double[3];
// 	refSatPos_Rov	=new double[3];
// 	refSatPos_Base	=new double[3];
// 	satePosBase		=new SatePos[num];
// 	satePosRov		=new SatePos[num];
// 	datarecord		=new DataRecord[num];
// 	mapWet			=new double[num];
// 	for (int i=0;i<3;i++)
// 	{
// 		refRecPos[i]		=0.0;
// 		rovRecPos[i]		=0.0;
// 		refSatPos_Rov[i]	=0.0;
// 		refSatPos_Base[i]	=0.0;
// 	}
// 		for (int i=0;i<num;i++)
// 		{
// 			mapWet[i]		=0.0;	
// 			rovPrn[i]=0;
// 		}
// 	
// }
// 
// DdData& DdData::operator=(const DdData& m)
// {
// 	
// 	 refPrn	=m.refPrn;
// 	 week	=m.week;
// 	 sec		=m.sec;
// 	for (int i=0;i<3;i++)
// 	{
// 		 refRecPos[i]		=m.refRecPos[i];
// 		 rovRecPos[i]		=m.rovRecPos[i];
// 		 refSatPos_Rov[i]	=m.refSatPos_Rov[i];
// 		 refSatPos_Base[i]	=m.refSatPos_Base[i];
// 	}
// 	if ( pairNum==m.pairNum)//initiate new dddata variable with m.pairNum 
// 	{
// 		for (int i=0;i<pairNum;i++)
// 		{
// 			 satePosBase[i]	=m.satePosBase[i];
// 			 satePosRov[i]	=m.satePosRov[i];
// 			 datarecord[i]	=m.datarecord[i];
// 			 mapWet[i]		=m.mapWet[i];	
// 			 rovPrn[i]			=m.rovPrn[i];
// 		}
// 	}
// 	if (pairNum!=m.pairNum)
// 	{
// 		pairNum=m.pairNum;
// 		delete[] rovPrn,satePosBase,satePosRov,datarecord,mapWet;
// 
// 		rovPrn				=new int[pairNum];
// 		satePosBase		=new SatePos[pairNum];
// 		satePosRov		=new SatePos[pairNum];
// 		datarecord		=new DataRecord[pairNum];
// 		mapWet			=new double[pairNum];
// 		for (int i=0;i<pairNum;i++)
// 		{
// 			satePosBase[i]	=m.satePosBase[i];
// 			satePosRov[i]	=m.satePosRov[i];
// 			datarecord[i]	=m.datarecord[i];
// 			mapWet[i]		=m.mapWet[i];	
// 			rovPrn[i]			=m.rovPrn[i];
// 		}
// 	}
// 	
// 	return *this;
// }
// DdData::~DdData()
// {
// 	delete[] rovPrn;
// 	delete[] refRecPos;
// 	delete[] rovRecPos;
// 	delete[]	refSatPos_Rov;
// 	delete[]	refSatPos_Base;
// 	delete[]	satePosBase;
// 	delete[]	satePosRov;
// 	delete[] datarecord;
// }
// 
// void DdData::rememo( int num )
// {
// 
// }
// 
// 
// 
// 
// DdCtrl::DdCtrl()
// {
// 	tropFlag		=0;
// 	tropStep		=3600;
// 	ionoFlag		=0;
// 	ambFlag		=0;
// 	otherFlag		=0;
// 	//combcoef		=new int[3];
// 	pseudoFlag	=0;
// 
// 	//phaseFlag		=new int[3];
// }
// 
// 
// 
// 
// 
// 
// DdAmbInfo::DdAmbInfo()
// {
// 	prnList	=NULL;
// 	fixFlag	=NULL;
// 	fixSolu	=NULL;
// }
// 
// DdAmbInfo::DdAmbInfo( int num,int nfreq )
// {
// 	pairNum=num;
// 	freqNum=nfreq;
// 	prnList	=new int[num];
// 	fixSolu	=new double[num*freqNum];
// 	fixFlag	=new int[num*freqNum];
// 	for (int i=0;i<num;i++)
// 	{
// 		prnList[i]=0;
// 		for (int j=0;j<freqNum;j++)
// 		{
// 			fixFlag[i*freqNum+j]=0;
// 			fixSolu[i*freqNum+j]=0;
// 		}
// 	}
// }
// 
// DdAmbInfo& DdAmbInfo::operator=(const DdAmbInfo& m)
// {
// 	pairNum		=m.pairNum;
// 	refSate			=m.refSate;
// 	for (int i=0;i<pairNum;i++)
// 	{
// 		prnList[i]	=m.prnList[i];
// 		for (int j=0;j<freqNum;j++)
// 		{
// 			fixFlag[j*pairNum+i]	=m.fixFlag[j*pairNum+i];
// 			fixSolu[j*pairNum+i]	=m.fixSolu[j*pairNum+i];
// 		}
// 	}
// 	return *this;
// }
// DdAmbInfo::~DdAmbInfo()
// {
// 	delete[] prnList;
// 	delete[] fixFlag;
// 	delete[] fixSolu;
// }
// 
// BroadEphDataGlo::BroadEphDataGlo()
// {
// 	Pos=new double[3];
// 	Vel=new double[3];
// 	Acc=new double[3];
// 	prn			=0;
// 	for (int i=0;i<3;i++)
// 	{
// 		Pos[i]=0.0;
// 		Vel[i]=0.0;
// 		Acc[i]=0.0;
// 	}
// 	toc			=0.0;
// 	NegativeTauN=0.0;
// 	PositiveGammaN=0.0;
// 	UTCSec=0.0;
// 	health=0;
// 	Age=0;
// 	freqNum=	0;
// }
// BroadEphDataGlo::~BroadEphDataGlo()
// {
// 	delete[] Pos;
// 	delete[] Vel;
// 	delete[] Acc;
// }
//////////////////////////////////////////////
////////////////////////////////////////////// 
/*

	this file defines all the structs of class used
*/


#pragma once
#include "stdafx.h"//changed
#include "MStruct.h"
#include "ExtFun.h"
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <iostream>


const static double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
const static double gst0 []={1999,8,22,0,0,0}; /* galileo system time reference */
const static double bdt0 []={2006,1, 1,0,0,0}; /* beidou time reference */

tests::tests()
{
	alpha=NULL;
	
	beta= NULL;
}

tests::tests(int a)
{
	num=a;
	alpha=new int[a];
	alpha[0]=9;
	beta= new double[a];
}

tests::~tests()
{
	delete[] alpha;
	delete[] beta;
}

tests & tests::operator=( const tests& m ) 
{
	if (m.num!=num)
	{
		delete[] alpha;
		delete[] beta;
	}
	num=m.num;
	int* alpha=new int[m.num] ;
	double* beta=new double	[m.num];
	for (int i=0;i<m.num;i++)
	{
		alpha[i]		=m.alpha[i];
		beta[i]		=m.beta[i];
	}
	return *this;
}

template<typename Ts>
extern void exchange(Ts* x, Ts* y)
{
	Ts temp;
	temp=*x;
	*x=*y;
	*y=temp;
}

BroadEphHeader::BroadEphHeader()
{
	version=0.0;
	sysid=0;
	leapSecond=0.0;
}

/*
	 the parameters of ionosphere correction
	 default:
			validA=0;
			validB=0;
			1=valid; 0=invalid
*/
BroadEphIonoCorr::BroadEphIonoCorr()
{
	validA=0;
	validB=0;
	alpha1=0.0;
	alpha2=0.0;
	alpha3=0.0;
	alpha4=0.0;
	beta1=0.0;
	beta2=0.0;
	beta3=0.0;
	beta4=0.0;
}
BroadEphIonoCorr& BroadEphIonoCorr::operator=(const BroadEphIonoCorr& m)
{
	validA=m.validA;
	validB=m.validB;
	alpha1=m.alpha1;
	alpha2=m.alpha2;
	alpha3=m.alpha3;
	alpha4=m.alpha4;
	beta1=m.beta1;
	beta2=m.beta2;
	beta3=m.beta3;
	beta4=m.beta4;
	return *this;
}


BroadTimeCorr::BroadTimeCorr()
{	// GPUT =GPS to UTC a0,a1
			valid=0;
			transTimeType=0; // 0 : no time corr in the eph-file 
			a0=0.0;
			a1=0.0;
			refTime=0;
			refWeek=0;
	}

BroadEphData::BroadEphData()
{
	prn=0;
	toc=0.0;						
	satClkBs=0.0;
	satClkDrt=0.0;
	satClkDrtRe=0.0;
	iode=0.0;
	crs=0.0;
	deltan=0.0;
	m0=0.0;
	cuc=0.0;
	e=0.0;	//Eccentricity
	cus=0.0;
	sqA=0.0;
	toe=0.0; //sec of GPS week
	cic=0.0;
	Omega=0.0;
	cis=0.0;
	i0=0.0;
	crc=0.0;
	omega=0.0;
	OmegaDot=0.0;
	idot=0.0;
	codeL2=0.0;
	gpsw=0.0;//with toe
	flagL2=0.0;
	satAccuracy=0.0;
	satHealth=0.0;
	tgd=0.0;
	iodClk=0.0;
	transTime=0.0;
	fitInterval=0.0;
	spare1=0.0;
	spare2=0.0;
}



ObsHeader::ObsHeader()
{

			version=0.0;
			sysid=0;
			interval=1.0;
			for (int i=0;i<3;i++)
			{
				appPos[i]=0.0;
				HEN[i]=0.0;
			}
}

ObsType2::ObsType2()
{

		sysid=0;
		sysNum=0;
		typeNum=0;
	}
ObsType3::ObsType3()
	{

		sysid=0;
		typeNum=0;
	}
SingleObsType::SingleObsType()
{	//2.x and 3.x

	sysid=0;//for 3.x. single system id
	obsCode="F";// F is failure
	freCode=0;
	channel="F";//for 3.x
	}

SecWeek::SecWeek()
{
	sec=0.0;
	week=0.0;
}
YMDHMS::YMDHMS()
{
	year=0;mon=0;day=0;hour=0;min=0;
	sec=0.0;
}

ObsDataRecord::ObsDataRecord()
{
	sysid=0;
	PRN=0;numVadPhs=0;numVadCod=0;
	for(int i=0;i<NFREQ+NEXOBS;i++)
	{
		LLI[i]=0;
		SNR[i]=0;
		Phase[i]=0.0;
		vadFlgPhs[i]=vadFlgCod[i]=0;
		PsRange[i]=0.0;
		Doppler[i]=0.0;
		isCycleSlip[i]=0;
	}
}
ObsDataRecord& ObsDataRecord::operator=(const ObsDataRecord& m)
{
	sysid=	m.sysid;
	PRN=	m.PRN;
	numVadPhs=m.numVadPhs;
	numVadCod=m.numVadCod;
	for(int i=0;i<NFREQ+NEXOBS;i++)
	{
		LLI[i]			=m.LLI[i];
		SNR[i]		=m.SNR[i];
		Phase[i]		=m.Phase[i];
		vadFlgPhs[i]=m.vadFlgPhs[i];
		vadFlgCod[i]=m.vadFlgCod[i];
		PsRange[i]=m.PsRange[i];
		Doppler[i]	=m.Doppler[i];
		isCycleSlip[i]=m.isCycleSlip[i];
	}
	return *this;
}

void ObsDataRecord::ZeroElem()
{
	sysid=0;
	PRN=0;numVadPhs=0;numVadCod=0;
	for(int i=0;i<NFREQ+NEXOBS;i++)
	{
		LLI[i]=0;
		SNR[i]=0;
		Phase[i]=0.0;
		vadFlgPhs[i]=vadFlgCod[i]=0;
		PsRange[i]=0.0;
		Doppler[i]=0.0;
		isCycleSlip[i]=0;
	}
}

ObsEpochData::ObsEpochData()
{
	sateNum=0;
	flag=0;
	week=0;
	sec=0.0;
//	obsdatarecord=NULL;
}



ObsEpochData& ObsEpochData::operator=(const ObsEpochData& m)
{
	//if (sateNum!=m.sateNum)
	//{
	//	delete[] obsdatarecord;
	//	ObsDataRecord* obsdatarecord=new ObsDataRecord[m.sateNum];
	//}
	sateNum			=m.sateNum;
	flag					=m.flag;
	week				=m.week;
	sec					=m.sec;
	for(int i=0;i<MAXNUMSATE;i++)
		obsdatarecord[i]	=m.obsdatarecord[i];
	return *this;
}

void ObsEpochData::ZeroElem()
{
	sateNum=0;
	flag=0;
	week=0;
	sec=0.0;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		obsdatarecord[i].ZeroElem();
	}
}

void ObsEpochData::Sort()
{
	ObsDataRecord temp;
	int i,j,pos,minelem;
	for (i=0;i<MAXNUMSATE;i++)
	{
		pos=i;
		minelem=obsdatarecord[i].PRN;
		for (j=i;j<MAXNUMSATE;j++)
		{
			if (obsdatarecord[j].PRN<minelem && obsdatarecord[j].PRN!=0)
			{
				minelem=obsdatarecord[j].PRN;
				pos=j;
			}
		}
		if (i!=pos)
		{
			temp=obsdatarecord[pos];
			obsdatarecord[pos]=obsdatarecord[i];
			obsdatarecord[i]=temp;
		}
	}
}

void ObsEpochData::ZeroElemSingle( int index )
{
	obsdatarecord[index].ZeroElem();
}

ObsEpochData& ObsEpochData::Shrink()
{
	ObsEpochData temp;
	int cnt=0;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		if (obsdatarecord[i].PRN!=0)
		{
			temp.obsdatarecord[cnt++]=obsdatarecord[i];
		}
	}
	temp.sateNum=cnt;
	temp.sec=sec;
	temp.week=week;
	temp.flag=flag;
	return temp;
}

/*this function is for QXWZ company,
 *
 *extract the obs that lacks the required types 
 */
ObsEpochData& ObsEpochData::AutoShrink(int typeNo )
{
	int num=typeNo;
	int* flag=new int[num];
	InitPtr(flag,num);
	for (int i=0;i<sateNum;i++)
	{
		for (int j=0;j<num;j++)   
		{
			if (obsdatarecord[i].vadFlgPhs[j]==0 || obsdatarecord[i].vadFlgCod[j]==0)
				obsdatarecord[i].PRN=0;
		}
	}
	delete[] flag;
	return Shrink();
}

void ObsEpochData::Cycle2MeterGlo(int* dNum)
{
	for(int i=0;i<sateNum&&Prn2Sysid(obsdatarecord[i].PRN)==2;i++)
	{
		for (int j=0;j<NFREQ;j++)
		{
			if (obsdatarecord[i].vadFlgPhs[j]==1)
			{
				obsdatarecord[i].Phase[j] *=CLIGHT/FreqSysGlo(j,dNum[i]);
			}
		}
	}
}



//ObsEpochData::~ObsEpochData()
//{
//	delete[] obsdatarecord;
//}



void Gtime::weeksec(int& Week,double& Second,YMDHMS gnsstime,int sysid)
{
	//返回值为GPS周秒
	
	int y;
	
    int m;
	if(gnsstime.mon<=2)
	{
		y=gnsstime.year-1;
		m=gnsstime.mon+12;
	}
	else
	{
		y=gnsstime.year;
		m=gnsstime.mon;
	}
	double JD;
	JD=((int)(365.25*y))+((int)(30.6001*(m+1.0)))+gnsstime.day+gnsstime.hour/24.0+1720981.5;
	int weekday=(((int)(JD+0.5))%7);//计算结果中0代表星期一，1代表星期二
	 weekday++;
	 if (weekday==7)
	 {
		 weekday=0;
	 }
	 Second=(weekday*24.0+(double)gnsstime.hour)*3600.0+(double)gnsstime.min*60.0+gnsstime.sec;
	 Week=(int)((JD-2444244.5)/7.0);
}
void Gtime::ymdhms(int Week,double Sec,YMDHMS& gnsstime)
{
	 double MJD = Week*7.0+44244.0+Sec/3600.0/24.0;  //modified julian day
	 int a   =(int)(MJD +0.5+0.5+1.0e-10)+2400000;
	 double FRAC=MJD +0.5+2400000.5-a;
	 int b = a + 1537;
	 int Cc = (int)((b - 122.1)/365.25 + 1.0e-10);
	 int d = (int)(365.25*Cc+1.0e-16);
	 int Ee = (int)((b-d)/30.6001 + 1.0e-10);
	 gnsstime.day=b-d-(int)(30.6001*Ee);
	 gnsstime.mon=Ee-1-12*((int)(Ee/14.0+1.0e-10));
	 gnsstime.year=Cc - 4715 -(int)((7.0+gnsstime.mon)/10.0 + 1.0e-10);
	 double TmHour = FRAC*24.0;
	 gnsstime.hour= (int)(TmHour + 1.0e-10);
	 double TmMin=(TmHour-gnsstime.hour)*60;
	 gnsstime.min=(int)(TmMin+1.0e-10);
	 gnsstime.sec= (TmMin-gnsstime.min)*60.0;	
}
void Gtime::doy(int Year,int Month,int Day,int& YearofDOY,int& DayofYear)
	{
		int DayNums[12];
	 DayNums[0]=31;
	 DayNums[2]=31;
	 DayNums[4]=31;
	 DayNums[6]=31;
	 DayNums[7]=31;
	 DayNums[9]=31;
	 DayNums[11]=31;
	 DayNums[3]=30;
	 DayNums[5]=30;
	 DayNums[8]=30;
	 DayNums[10]=30;

	 if ((Year%400==0) || (Year%4==0 && Year%100!=0))
		 DayNums[1]=29;
	 else
		 DayNums[1]=28;
	 DayofYear=0;
	 for (int i=0;i<Month-1;i++)
		  DayofYear+=DayNums[i];
	 DayofYear+=Day;

	 YearofDOY=Year;
	}



SppCtrl::SppCtrl()
{
	freqIndex	=0;
	ionoFlag	=0;
	sysNum	=1;
	sysid			=1;
	maskEle	=5.0;
	IFIndex[0]	=0;
	IFIndex[1]	=1;
}
SatePos::SatePos()
{
	for(int i=0;i<3;i++)
		sateXYZ[i]	=0.0;
}
SatePos& SatePos::operator=(const SatePos& m)
{
	for(int i=0;i<3;i++)
		sateXYZ[i]	=m.sateXYZ[i];
	return *this;
}

void SatePos::ZeroElem()
{
	for(int i=0;i<3;i++)
		sateXYZ[i]	=0.0;
}



SppInfo::SppInfo()
{
	validnum	=0;
	dtr			=0.0;
	for(int i=0;i<3;i++)
	{
		recPos[i]=	0.0;
	}
	for(int i=0;i<MAXNUMSATE;i++)
	{
		sateclkVel[i]	=0.0;
		//		emiTime[i]		=	0.0;
		prnList[i]		=	0;
		ele[i]				=0.0;
		azi[i]				=0.0;
		sateclkerr[i]	=0.0;
		tropCorr[i]		=0.0;
		mapWet[i]		=0.0;
		codeCorr[i]	=0.0;
		residual[i]		=0.0;

		for(int j=0;j<3;j++)
		{
			satePos[i].sateXYZ[j]=0.0;	
			sateVel[i].sateXYZ[j]=0.0;	
		}
	}
}

SppInfo& SppInfo::operator=( const SppInfo& m ) 
{
	validnum	=m.validnum;
	dtr			=m.dtr;
	for(int i=0;i<3;i++)
	{
		recPos[i]=	m.recPos[i];
	}
	int sateNum=m.validnum;
	for(int i=0;i<sateNum;i++)
	{
		sateclkVel[i]	=m.sateclkVel[i];
//		emiTime[i]		=	m.emiTime[i];
		prnList[i]		=	m.prnList[i];
		ele[i]				=m.ele[i];
		azi[i]				=m.azi[i];
		sateclkerr[i]	=m.sateclkerr[i];
		tropCorr[i]		=m.tropCorr[i]	;
		mapWet[i]		=m.mapWet[i];
		codeCorr[i]	=m.codeCorr[i];
		for(int j=0;j<3;j++)
		{
			satePos[sateNum].sateXYZ[j]=m.satePos[sateNum].sateXYZ[j];	
			sateVel[sateNum].sateXYZ[j]=m.sateVel[sateNum].sateXYZ[j];	
		}
	}

	return *this;
}

void SppInfo::ZeroElem()
{
	validnum	=0;
	dtr			=0.0;
	
	for(int i=0;i<MAXNUMSATE;i++)
	{
		satePos[i].ZeroElem();
		sateVel[i].ZeroElem();
		sateclkVel[i]	=0.0;
		//		emiTime[i]		=	0.0;
		prnList[i]		=	0;
		ele[i]				=0.0;
		azi[i]				=0.0;
		sateclkerr[i]	=0.0;
		tropCorr[i]		=0.0;
		mapWet[i]		=0.0;
		codeCorr[i]	=0.0;
		residual[i]		=0.0;
	}
}


SppInfoGlo::SppInfoGlo()
{
	validnum	=0;
	dtr			=0.0;
	for(int i=0;i<3;i++)
	{
		recPos[i]=	0.0;
	}
	for(int i=0;i<MAXNUMSATE;i++)
	{
		sateclkVel[i]	=0.0;
		//		emiTime[i]		=	0.0;
		prnList[i]		=	0;
		freqNum[i]=0;
		ele[i]				=0.0;
		azi[i]				=0.0;
		sateclkerr[i]	=0.0;
		tropCorr[i]		=0.0;
		mapWet[i]		=0.0;
		codeCorr[i]	=0.0;
		residual[i]		=0.0;

		for(int j=0;j<3;j++)
		{
			satePos[i].sateXYZ[j]=0.0;	
			sateVel[i].sateXYZ[j]=0.0;	
		}
	}
}

SppInfoGlo& SppInfoGlo::operator=( const SppInfoGlo& m ) 
{
	validnum	=m.validnum;
	dtr			=m.dtr;
	for(int i=0;i<3;i++)
	{
		recPos[i]=	m.recPos[i];
	}
	int sateNum=m.validnum;
	for(int i=0;i<sateNum;i++)
	{
		sateclkVel[i]	=m.sateclkVel[i];
		//		emiTime[i]		=	m.emiTime[i];
		prnList[i]		=	m.prnList[i];
		freqNum[i]	=m.freqNum[i];
		ele[i]				=m.ele[i];
		azi[i]				=m.azi[i];
		sateclkerr[i]	=m.sateclkerr[i];
		tropCorr[i]		=m.tropCorr[i]	;
		mapWet[i]		=m.mapWet[i];
		codeCorr[i]	=m.codeCorr[i];
		for(int j=0;j<3;j++)
		{
			satePos[sateNum].sateXYZ[j]=m.satePos[sateNum].sateXYZ[j];	
			sateVel[sateNum].sateXYZ[j]=m.sateVel[sateNum].sateXYZ[j];	
		}
	}

	return *this;
}

void SppInfoGlo::ZeroElem()
{
	validnum	=0;
	dtr			=0.0;

	for(int i=0;i<MAXNUMSATE;i++)
	{
		satePos[i].ZeroElem();
		sateVel[i].ZeroElem();
		sateclkVel[i]	=0.0;
		//		emiTime[i]		=	0.0;
		prnList[i]		=	0;
		freqNum[i]=0;
		ele[i]				=0.0;
		azi[i]				=0.0;
		sateclkerr[i]	=0.0;
		tropCorr[i]		=0.0;
		mapWet[i]		=0.0;
		codeCorr[i]	=0.0;
		residual[i]		=0.0;
	}
}


SdDataRecord::SdDataRecord()
{
	for(int i=0;i<NFREQ+NEXOBS;i++)
	{
		Phase[i]		=0.0;
		PsRange[i]=0.0;
	}
}

void SdDataRecord::ZeroElem()
{
	for(int i=0;i<NFREQ+NEXOBS;i++)
	{
		Phase[i]		=0.0;
		PsRange[i]=0.0;
	}
}


SdData::SdData()
{
	week=0;
	sec=0.0;
	satnum=0;
	for (int i=0;i<3;i++)
	{
		refRecPos[i]	=0.0;
		rovRecPos[i]	=0.0;
	}
	for (int i=0;i<MAXNUMSATE;i++)
	{
		prn[i]=0;
		mapWet[i]=0.0;
		ele[i]=0.0;
		tropCor[i]=0.0;
	}
	
}

SdData& SdData::operator=(const SdData& m)
{
	week	=m.week;
	sec		=m.sec;
	satnum			=m.satnum;
	for (int i=0;i<3;i++)
	{
		refRecPos[i]		=m.refRecPos[i];
		rovRecPos[i]		=m.rovRecPos[i];
	}
	for (int i=0;i<m.satnum;i++)
	{
		satePosBase[i]	=m.satePosBase[i];
		satePosRov[i]	=m.satePosRov[i];
		sddatarecord[i]=m.sddatarecord[i];
		prn[i]		=m.prn[i];
		mapWet[i]		=m.mapWet[i];
		tropCor[i]=m.tropCor[i];
		ele[i]=m.ele[i];
	}
	return *this;
}

void SdData::ZeroElem()
{
	week=0;
	sec=0.0;
	satnum=0;
	for (int i=0;i<3;i++)
	{
		refRecPos[i]	=0.0;
		rovRecPos[i]	=0.0;
	}
	for (int i=0;i<MAXNUMSATE;i++)
	{
		prn[i]=0;
		mapWet[i]=0.0;
		ele[i]=0.0;
		tropCor[i]=0.0;
		satePosBase[i].ZeroElem();
		satePosRov[i].ZeroElem();
		sddatarecord[i].ZeroElem();
	}

}





DataRecord:: DataRecord()
{
	numVadPhs=0;
	numVadCod=0;
	
	for (int i=0;i<NFREQ+NEXOBS;i++)
	{
		Phase[i]		=0.0;
		PsRange[i]=0.0;
		vadFlgPhs[i]=vadFlgCod[i]=0.0;
		isCycleSlip[i]=0;
	}
}
DataRecord& DataRecord::operator=(const DataRecord& m)
{
	numVadPhs=m.numVadPhs;
	numVadCod=m.numVadCod;
	for (int i=0;i<NFREQ+NEXOBS;i++)
	{
		Phase[i]		=m.Phase[i];
		PsRange[i]=m.PsRange[i];
		vadFlgPhs[i]=m.vadFlgPhs[i];
		vadFlgCod[i]=m.vadFlgCod[i];
		isCycleSlip[i]=m.isCycleSlip[i];
	}
	return *this;
}

void DataRecord::ZeroElem()
{
	numVadPhs=0;
	numVadCod=0;

	for (int i=0;i<NFREQ+NEXOBS;i++)
	{
		Phase[i]		=0.0;
		PsRange[i]=0.0;
		vadFlgPhs[i]=vadFlgCod[i]=0.0;
		isCycleSlip[i]=0;
	}
}



DdData& DdData::operator=(const DdData& m)
{
	
	 refPrn	=m.refPrn;
	 week	=m.week;
	 sec		=m.sec;
	 pairNum=m.pairNum;
	for (int i=0;i<3;i++)
	{
		 refRecPos[i]		=m.refRecPos[i];
		 rovRecPos[i]		=m.rovRecPos[i];
		 refSatPos_Rov[i]	=m.refSatPos_Rov[i];
		 refSatPos_Base[i]	=m.refSatPos_Base[i];
	}
		for (int i=0;i<MAXNUMSATE;i++)
		{
			 satePosBase[i]	=m.satePosBase[i];
			 satePosRov[i]	=m.satePosRov[i];
			 datarecord[i]	=m.datarecord[i];
			 mapWet[i]		=m.mapWet[i];	
			 rovPrn[i]			=m.rovPrn[i];
			 tropCor[i]			=m.tropCor[i];
			 ele[i]				=m.ele[i];
		}
	return *this;
}

DdData::DdData()
{
	refPrn	=0;
	week	=0;
	sec		=0.0;
	pairNum=0;
	for (int i=0;i<3;i++)
	{
		refRecPos[i]		=0.0;
		rovRecPos[i]		=0.0;
		refSatPos_Rov[i]	=0.0;
		refSatPos_Base[i]	=0.0;
	}
	for (int i=0;i<MAXNUMSATE;i++)
	{
		mapWet[i]		=0.0;
		rovPrn[i]			=0;
		tropCor[i]=0.0;
		ele[i]=0.0;
	}
}

void DdData::ZeroElem()
{
	refPrn	=0;
	week	=0;
	sec		=0.0;
	pairNum=0;
	for (int i=0;i<3;i++)
	{
		refRecPos[i]		=0.0;
		rovRecPos[i]		=0.0;
		refSatPos_Rov[i]	=0.0;
		refSatPos_Base[i]	=0.0;
	}
	for (int i=0;i<MAXNUMSATE;i++)
	{
		mapWet[i]		=0.0;
		rovPrn[i]			=0;
		datarecord[i].ZeroElem();
		satePosBase[i].ZeroElem();
		satePosRov[i].ZeroElem();
		tropCor[i]=0.0;
		ele[i]=0.0;
	}
}

int DdData::SumCod()
{
	int cnt=0;
	for (int i=0;i<NFREQ;i++)
	{
		cnt+=NoCod(i);
	}
	return cnt;
}

int DdData::NoCod( int ind )
{
	int cnt=0;
	for (int i=0;i<pairNum;i++)
	{
		if (datarecord[i].vadFlgCod[ind]==1)
		{
			cnt++;
		}
	}
	return cnt;
}

int DdData::SumPhs()
{
	int cnt=0;
	for (int i=0;i<NFREQ;i++)
	{
		cnt+=NoPhs(i);
	}
	return cnt;
}

int DdData::NoPhs( int ind )
{
	int cnt=0;
	for (int i=0;i<pairNum;i++)
	{
		if (datarecord[i].vadFlgPhs[ind]==1)
		{
			cnt++;
		}
	}
	return cnt;
}
DdData& DdData::Sort()
{
	int i,j,pos,minelem;
	for (i=0;i<MAXNUMSATE;i++)
	{
		pos=i;
		minelem=rovPrn[i];
		for (j=i;j<MAXNUMSATE;j++)
		{
			if (rovPrn[j]<minelem && rovPrn[j]!=0)
			{
				minelem=rovPrn[j];
				pos=j;
			}
		}
		if (i!=pos)
		{
			exchange(&datarecord[pos],&datarecord[i]);

			exchange(&satePosRov[pos],&satePosRov[i]);

			exchange(&satePosBase[pos],&satePosBase[i]);
		    
			exchange(&mapWet[pos],&mapWet[i]);

			exchange(&tropCor[pos],&tropCor[i]);

			exchange(&ele[pos],&ele[i]);
		}
	}

	return *this;
}

double DdData::distBaseRef()
{
	return DistofVector(refRecPos,refSatPos_Base,3);
}

double DdData::distRoverRef()
{
	return DistofVector(rovRecPos,refSatPos_Rov,3);
}

double DdData::distBaseRover( int index )
{
	double s[3];
	for (int i=0;i<3;i++)
	{
		s[i]=refRecPos[i]-satePosBase[index].sateXYZ[i];
	}
	return Norm(s,3);
}

double DdData::distRoverRover( int index )
{
	double s[3];
	for (int i=0;i<3;i++)
	{
		s[i]=rovRecPos[i]-satePosRov[index].sateXYZ[i];
	}
	return Norm(s,3);
}


/*
 *num indicates the num of freq uesd
 */
DdData& DdData::Check( int numFreq )
{
	for (int i=0;i<pairNum;i++)
	{
		if (datarecord[i].numVadCod<numFreq || datarecord[i].numVadPhs<numFreq)
		{
			rovPrn[i]=0;
		}
	}
	return Shrink();
}

DdData& DdData::Shrink()
{
	DdData temp;
	int cnt=0;
	int curnum=pairNum;
	for (int i=0;i<curnum;i++)
	{
		if (rovPrn[i]>0)
		{
			  datarecord[cnt]=datarecord[i];
			  rovPrn[cnt]=rovPrn[i];
			  mapWet[cnt]=mapWet[i];
			  tropCor[cnt]=tropCor[i];
			  ele[cnt]=ele[i];
			  satePosBase[cnt]=satePosBase[i];
			  satePosRov[cnt++]=satePosRov[i];
		}
	}
	  pairNum=cnt;
	  refPrn=refPrn;
	  week=week;
	  sec=sec;
	  ClearTail(cnt);
	return *this/*temp*/;
}

void DdData::ClearTail( int index )
{
	for (int i=index;i<MAXNUMSATE;i++)
	{
		datarecord[i].ZeroElem();
		rovPrn[i]=0;
		mapWet[i]=tropCor[i]=ele[i]=0.0;
		satePosBase[i].ZeroElem();
		satePosRov[i].ZeroElem();
	}
}

double DdData::distRecSate_DD( int index )
{
	return distRoverRover(index)-distBaseRover(index)-distRoverRef()+distBaseRef();
}

DdAmbCtrl::DdAmbCtrl()
{
	flag=1;
	sysid=1;
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			coef[i][j]=0;
		}
		numAmb[i]=0;
	}
	ratiothrsd=2.5;
	parEleMask=20.0;
}



DdCtrl::DdCtrl()
{
	weightMode=0;
	sigmaCod=0.0;
	sigmaPhs=0.0;
	tropFlag		=0;
	tropStep		=0;
	maskele		=0.0;
	thresGF=0.0;
	thresMW=0.0;
	ionoFlag		=0;
	ambFlag		=0;
	otherFlag		=0;
	for (int i=0;i<3;i++)
	{
		freqCod[i]=freqPhs[i]=0;
		for (int j=0;j<3;j++)
		{
			pseudoCoef[i][j]=0;
		}
	}
	pseudoFlag	=0;
	//phaseFlag		=new int[3];
}

int DdCtrl::CodTypeNo()
{
	return (pseudoFlag>3)?pseudoFlag-3:pseudoFlag;
}

int DdCtrl::PhsTypeNo()
{
	return (ddambctrl.flag>3)?ddambctrl.flag-3:ddambctrl.flag;
}




DdAmbInfo::DdAmbInfo()
{
	pairNum=0;
	freqNum=0;
	refSate=0;
	for (int i=0;i<NEXOBS;i++)	freq[i]=0.0;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		for (int j=0;j<3;j++)
		{
			prnList[j][i]=fixFlag[j][i]=0;
			fixSolu[j][i]=0.0;
		}
	}
}



DdAmbInfo& DdAmbInfo::operator=(const DdAmbInfo& m)
{
	ZeroElem();
	pairNum		=m.pairNum;
	refSate			=m.refSate;
	freqNum		=m.freqNum;
	for (int i=0;i<NEXOBS;i++)	freq[i]=m.freq[i];
	for (int i=0;i<pairNum;i++)
	{
		for (int j=0;j<3;j++)
		{
			prnList[j][i]	=m.prnList[j][i];
			fixFlag[j][i]	=m.fixFlag[j][i];
			fixSolu[j][i]	=m.fixSolu[j][i];
		}
	}
	return *this;
}

int DdAmbInfo::NoSat( int ind )
{
	int i;
	int cc=0;
	for (i=0;i<MAXNUMSATE;i++)
	{
		if (prnList[ind][i]!=0)
		{
			cc++;
		}
	}
	return cc;
}

int DdAmbInfo::FindSat( int ind,int prn )
{
	int pos=-1;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		if (prnList[ind][i]=prn)
		{
			pos=i;
			break;
		}
		
	}
	return pos;
}

void DdAmbInfo::ZeroElem()
{
	pairNum=0;
	freqNum=0;
	refSate=0;
	for (int i=0;i<NEXOBS;i++)	freq[i]=0.0;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		for (int j=0;j<3;j++)
		{
			prnList[j][i]=fixFlag[j][i]=0;
			fixSolu[j][i]=0.0;
		}
	}
}

int DdAmbInfo::NoUnfix( int ind )
{
	return NoSat(ind)-NoFix(ind);
}

int DdAmbInfo::NoFix( int ind )
{
	int cnt=0;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		if (fixFlag[ind][i]==1)
		{
			cnt++;
		}
	}
	return cnt;
}

int DdAmbInfo::SumUnfix()
{
	int cnt=0;
	for (int i=0;i<3;i++)
	{
		cnt+=NoUnfix(i);
	}
	return cnt;
}

int DdAmbInfo::SumNoSat()
{
	int cnt=0;
	for (int i=0;i<3;i++)
	{
		cnt+=NoSat(i);
	}
	return cnt;
}

/* flag =10,get all unfix sate*/
void DdAmbInfo::GetUnfixPrnlist( int* prnptr,int flag )
{
	if (flag!=10)
	{
		int cnt=0;
		for (int i=0;i<NoSat(flag);i++)
		{
			if(fixFlag[flag][i]==0)  prnptr[cnt++]=prnList[flag][i];
		}
	}
}
/*change the data of ptr[pos1] and ptr[pos2] */
void DdAmbInfo::InterChange(int index,int pos1,int pos2)
{
	int dat;
	dat=prnList[index][pos1];
	prnList[index][pos1]=prnList[index][pos2];
	prnList[index][pos2]=dat;

	dat=fixFlag[index][pos1];
	fixFlag[index][pos1]=fixFlag[index][pos2];
	fixFlag[index][pos2]=dat;

	double solu;
	solu=fixSolu[index][pos1];
	fixSolu[index][pos1]=fixSolu[index][pos2];
	fixSolu[index][pos2]=solu;
}

DdAmbInfo& DdAmbInfo::Update( DdData curdata )
{
	DdAmbInfo minfo;
	minfo.pairNum=curdata.pairNum;
	minfo.refSate=curdata.refPrn;
	int i,j;
	for (i=0;i<NFREQ+NEXOBS;i++)
	{
		int cnt=0;
		int num=curdata.pairNum;
		for (j=0;j<num;j++)
		{
			if(curdata.datarecord[j].vadFlgPhs[i]==1)	
			{
				minfo.prnList[i][cnt]=curdata.rovPrn[j];
				// watch out that the index is 0
				if (curdata.datarecord[j].isCycleSlip[0]==1)
				{
					minfo.fixFlag[i][cnt]=0;
					minfo.fixSolu[i][cnt]=0.0;
				}
				cnt++;
			}
		}

	}
	
	return minfo;
}




BroadEphDataGlo::BroadEphDataGlo()
{

	prn			=0;
	for (int i=0;i<3;i++)
	{
		Pos[i]=0.0;
		Vel[i]=0.0;
		Acc[i]=0.0;
	}
	toc			=0.0;
	toe=0.0;
	week=0;
	idoe=0;
	NegativeTauN=0.0;
	PositiveGammaN=0.0;
	health=0;
	Age=0;
	freqNum=	0;
}



BroadEphDataSBAS::BroadEphDataSBAS()
{
	week=prn=sva=svh=0;
	t0=tof=a0=a1=0.0;
	for (int i=0;i<3;i++)
	{
		pos[i]=vel[i]=acc[i]=0.0;
	}
}

DdObsInfo::DdObsInfo()
{
	eleRefBase=0.0;
	eleRefRov=0.0;

	for (int i=0;i<MAXNUMSATE;i++)
	{
		eleRovBase[i]=eleRovRov[i]=0.0;
	}
	for (int i=0;i<3;i++)
	{
		numCod[i]=numPhs[i]=0;
		for (int j=0;j<MAXNUMSATE;j++)
		{
			prnlistPhs[i][j]=prnlistCod[i][j]=0;
		}
	}
}

void DdObsInfo::ZeroEle()
{
	eleRefBase=0.0;
	eleRefRov=0.0;

	for (int i=0;i<MAXNUMSATE;i++)
	{
		eleRovBase[i]=eleRovRov[i]=0.0;
	}
	for (int i=0;i<3;i++)
	{
		numCod[i]=numPhs[i]=0;
		for (int j=0;j<MAXNUMSATE;j++)
		{
			prnlistPhs[i][j]=prnlistCod[i][j]=0;
		}
	}
}
int DdObsInfo:: NoCod(int ind)
{
	return numCod[ind];
}
int DdObsInfo::SumCod()
{
	int temp=0;
	for (int i=0;i<3;i++)
	{
		temp+=numCod[i];
	}
	return temp;
}

int DdObsInfo::SumPhs()
{
	int temp=0;
	for (int i=0;i<3;i++)
	{
		temp+=numPhs[i];
	}
	return temp;
}

int DdObsInfo::NoPhs( int ind )
{
	return numPhs[ind];
}
void DdObsInfo:: SetList(int* prnlist,int ind,int flag)
{
	if (flag==0)
	{
		for (int i=0;i<numCod[ind];i++)
		{
			prnlist[i]=prnlistCod[ind][i];
		}
	}
	if (flag==1)
	{
		for (int i=0;i<numPhs[ind];i++)
		{
			prnlist[i]=prnlistPhs[ind][i];
		}
	}
}

DdObsInfo& DdObsInfo::operator=( const DdObsInfo& m )
{
	ZeroEle();
	eleRefBase=m.eleRefBase;
	eleRefRov=m.eleRefRov;

	for (int i=0;i<MAXNUMSATE;i++)
	{
		eleRovBase[i]=m.eleRovBase[i];
		eleRovRov[i]=m.eleRovRov[i];
	}
	for (int i=0;i<3;i++)
	{
		numCod[i]=m.numCod[i];
		numPhs[i]=m.numPhs[i];
		for (int j=0;j<MAXNUMSATE;j++)
		{
			prnlistCod[i][j]=m.prnlistCod[i][j];
			prnlistPhs[i][j]=m.prnlistPhs[i][j];
		}
	}
	return *this;
}

/* flag =10,get all unfix sate*/
void DdObsInfo::GetUnfixElelist( double* eleptr,int flag ,DdAmbInfo ambinfo)
{
	if (flag!=10)
	{
		int cnt=0;
		for (int i=0;i<ambinfo.NoSat(flag);i++)
		{
			if(ambinfo.fixFlag[flag][i]==0)  eleptr[cnt++]=(eleRovRov[i]+eleRovBase[i])/2.0;
		}
	}
}

/* return average ele of rover sat to 2 station*/
double DdObsInfo::GetEle( int prn )
{
	double ele=-1.0;
	for (int i=0;i<MAXNUMSATE;i++)
	{
		if(prn==prnlistCod[0][i] || prn==prnlistPhs[0][i]) 
		{
			ele=(eleRovRov[i]+eleRovBase[i])/2.0;
			break;
		}
	}
	return ele;
}

AmbData::AmbData()
{
	Prn=0;
	CurrentIndex=-1;
	FirstObs=LastObs=0.0;
	InitPtr(Cycle,MAXOBSEPOCH);
	InitPtr(VadFlag,MAXOBSEPOCH);
	InitPtr(isWLInteger,MAXOBSEPOCH);
}

AmbData& AmbData::operator=( const AmbData& m )
{
	Prn=m.Prn;
	FirstObs=m.FirstObs;
	LastObs=m.LastObs;
	CurrentIndex=m.CurrentIndex;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		Cycle[i]=m.Cycle[i];
		VadFlag[i]=m.VadFlag[i];
	}
	return *this;
}

int AmbData::NumData()
{
	int num=0;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if (Cycle[i]!=0.0)
		{
			num++;
		}
	}
	return num;
}

int AmbData::NumVadData()
{
	int num=0;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if (VadFlag[i]!=0)
		{
			num++;
		}
	}
	return num;
}

double AmbData::MeanAmb()
{
	double data_t[MAXOBSEPOCH];
	int cnt=0;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if (VadFlag[i]!=0)
			data_t[cnt++]=Cycle[i];
	}

	return  Mean(data_t,cnt);
}

double AmbData::CheckAmb()
{
	double mean=MeanAmb();
	return ROUND( mean );
}

double AmbData::StdofAver()
{
	double r=0.0;
	double mwCheck=CheckAmb();
	int      vadnum=NumVadData();
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if(VadFlag[i]!=0)
			r+= SQ(Cycle[i]-mwCheck);
	}
	return sqrt(r)/vadnum;
}

double AmbData::BiasCycle()
{
	double mean=MeanAmb(), check=ROUND(mean);
	return fabs(mean-check);
}

void AmbData::ZeroElem()
{
	Prn=0;
	FirstObs=LastObs=0.0;
	CurrentIndex=-1;
	InitPtr(Cycle,MAXOBSEPOCH);
	InitPtr(VadFlag,MAXOBSEPOCH);
}



UnionEph::UnionEph()
{
	ephNum=ephNumGlo=0;
}

InPutFileSet::InPutFileSet()
{
	numEph=numSp3=numClk=numBase=numRover=0;
	int i;
	for (i=0;i<MAXFILEEPH;i++)		fileEph[i]=_T("FileEph");
	for (i=0;i<MAXFILESP3;i++)		fileSp3[i]=_T("FileSP3");
	for (i=0;i<MAXFILECLK;i++)		fileSp3[i]=_T("FileCLK");
	for (i=0;i<MAXFILEBASE;i++)	fileSp3[i]=_T("FileBase");
	for (i=0;i<MAXFILEROVER;i++)	fileSp3[i]=_T("FileRover");
}

void InPutFileSet::ModifyPath(CString& filepath)
{
	int pos1=0, pos2=0;
	pos1=filepath.Find(_T("\\"),pos2+1);
	pos2=filepath.Find(_T("\\"),pos1+1);
	if (pos2-pos1!=1)
	{
		filepath.Replace(_T("\\"),_T("/"));
	}
}


void InPutFileSet::CheckPath()
{
	int i;
	for (i=0;i<numEph;i++) ModifyPath(fileEph[i]);
	for (i=0;i<numSp3;i++) ModifyPath(fileSp3[i]);
	for (i=0;i<numClk;i++) ModifyPath(fileClk[i]);
	for (i=0;i<numBase;i++) ModifyPath(fileBase[i]);
	for (i=0;i<numRover;i++) ModifyPath(fileRover[i]);
}


EwlData::EwlData()
{
	Prn=0;
	CurrentIndex=-1;
	FirstObs=LastObs=0.0;
	InitPtr(EwlCycle[0],MAXOBSEPOCH);
	InitPtr(VadFlag[0],MAXOBSEPOCH);
	InitPtr(EwlCycle[1],MAXOBSEPOCH);
	InitPtr(VadFlag[1],MAXOBSEPOCH);
}

EwlData& EwlData::operator=( const EwlData& m )
{
	Prn=m.Prn;
	FirstObs=m.FirstObs;
	LastObs=m.LastObs;
	CurrentIndex=m.CurrentIndex;
	for (int j=0;j<2;j++)
	{
		for (int i=0;i<MAXOBSEPOCH;i++)
		{
			EwlCycle[j][i]=m.EwlCycle[j][i];
			VadFlag[j][i]=m.VadFlag[j][i];
		}
	}
	return *this;
}

int EwlData::NumData(int index)
{
	int num=0;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if (EwlCycle[index][i]!=0.0)
		{
			num++;
		}
	}
	return num;
}

int EwlData::NumVadData(int index)
{
	int num=0;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if (VadFlag[index][i]!=0)
		{
			num++;
		}
	}
	return num;
}

double EwlData::EwlMean(int index)
{
	double data_t[MAXOBSEPOCH];
	int cnt=0;
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if (VadFlag[index][i]!=0)
			data_t[cnt++]=EwlCycle[index][i];
	}

	return  Mean(data_t,cnt);
}

double EwlData::EwlCheck(int index)
{
	double mean=EwlMean(index);
	return ROUND( mean );
}

double EwlData::EwlStdofAverWl(int index)
{
	double r=0.0;
	double mwCheck=EwlCheck(index);
	int      vadnum=NumVadData( index);
	for (int i=0;i<MAXOBSEPOCH;i++)
	{
		if(VadFlag[index][i]!=0)
			r+= SQ(EwlCycle[index][i]-mwCheck);
	}
	return sqrt(r)/vadnum;
}

double EwlData::BiasCycle(int index)
{
	double mean=EwlMean(index), check=ROUND(mean);
	return fabs(mean-check);
}

void EwlData::ZeroElem()
{
	Prn=0;
	FirstObs=LastObs=0.0;
	CurrentIndex=-1;
	InitPtr(EwlCycle[0],MAXOBSEPOCH);
	InitPtr(VadFlag[0],MAXOBSEPOCH);
	InitPtr(EwlCycle[1],MAXOBSEPOCH);
	InitPtr(VadFlag[1],MAXOBSEPOCH);
}


fixinfo::fixinfo()
{
	prn=fixepoch2=fixepoch1=valid=0;
	check1=check2=bias1=bias2=0.0;
}
