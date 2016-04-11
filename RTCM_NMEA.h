/* 
Reference :
	1.Protocol of NMEA 0183
	2.U-blox  receiver description
	3.Trimble OME receiver help
Frame
	|	$	|	<Address> |		{,<value>}		|	*<checksum>	|	<CR><LF>
	Address =talker identifier + sentence formatter
	"|"is the self-defining delimiter,and it doesn't exist in the sentence


	the talker identifiers are slightly different according to the different chips' manufacturers
	for example
	GP is GPS;		GL=GLONASS;	BD(Samsung Galaxy)&GB(ublox)=BDS

	GNSS capable receivers will always output the message with the GN talker ID
	GNSS capable receivers will also output the message with the GP and/or GL talker ID when using more than one constellation for the position fix
	GN may appear repeatedly
	see the definition of struct
Satellite numbering in NMEA
	GPS		SBAS				QZSS		GLONASS		BDS
	1-32		33-64				-			65-96			200+		STRICT
	...		 &152-158	193-197		...				 ...			EXTENDED
	"..."is the same 
	unknown satellite numbers are always 255 in ublox!
	more, see manufacturer's protocol

Written by:	Li Zhen Collage of Surveying and Geo-informatics of Tongji University
Date:			2014.10.15
Revised		struct to class 
by:	Li Zhen
Date:			2014.11.05
Email:			14lizhen_csgi@tongji.edu.cn
*/

//#pragma once
using namespace std;
#include<fstream>
#include<iomanip>
#include "stdafx.h"
#include "atlstr.h"


#define CLIGHT				299792458.0					//speed of light
#define PI							3.1415926535897932		
#define D2R						PI/180						//deg 2 rad
#define R2D						180/PI						//rad 2 deg

#define RE_WGS84			6378137.0						// earth semimajor axis (WGS84) (m) 
#define FE_WGS84			1.0/298.257223563	   // earth flattening (WGS84)
#define ee						0.00669437999013
#define RE_CGCS2000		6378137.0
#define FE_CGCS2000		1.0/298.257222101


//DTAI is changeable. See IERS bulletin
#define DTAI					35	.0								// since 2013.05			
#define BDT_GPST			-14.0								//BDT-GPST=-14.0
#define GPST_UTC			16.0						//GPST-UTC	=DTAI-19.0		
#define BDT_UTC				2.0						//BDT-UTC=DTAI-33.0


#define ParaNum				4									//parameter number in the equation


//----------data------------------
#define MAXNUMGSV							20				//max number in GPGSV ,maxnum is 12 in GSA  
#define MAXNUMGSA							12
#define MAXCORSNUMSAT					25
#define MAXEPOCH								50000

/*	Global Positioning System Fix Data	*/
class gpgga{				//  GPGGA struct

public:
	gpgga();

public:
	CString		GP_ID;
	double		utc;				//	hhnnss.ssss
	double		lat;				//	Latitude:ddmm.mmmm
	CString		ns_ind;			//	N/S indicator
	double		lon;				//	Longtitude:dddmm.mmmm
	CString		ew_ind;		//	E/W indicator
	int			pos_fix;			//	position fix 0=invalid;1=SPS;2=DGPS;3=PPS fix;4=rtk;5=float rtk;6=unknown
	int			num_sate;		//	number of satellites used
	double		hdop;			
	double		alt;				//altitude
	CString		alt_unit;		//altitude units
	double		hei;		//geoid separation, same as height of geoid
	CString		sep_unit;	//seperation units
	double		time_frac;		//time fracrion since DGPS:seconds
	CString		diff_id;		//DGPS station ID:0000~1023	
};

class gngns{		//GNGNS	same as GPGGA

public:
	CString		GP_ID;
	double		utc;				//1	hhnnss.ssss
	double		lat;				//	Latitude:ddmm.mmmm
	CString		ns_ind;			//	N/S indicator
	double		lon;				//	Longtitude:dddmm.mmmm
	CString		ew_ind;		//5	E/W indicator
	int			mode_ind;	//one character stands for one constellation,1st is GPS 2nd is GLO ,more in the future
	//N;A;D;P;R;F;E;M;S.see the procotol
	int			num_sate;		//	number of satellites used
	double		hdop;			
	double		orth_height;	//orthometric height in meters(MSL reference)
	double		geo_sep;		//10	geoid separation ,ellipsoid surface-geoid surface
	double		data_age;	//the age of differential data ,NULL if tID is GN
	CString		ref_sta;		//DGPS station ID:0000~1023	
};


/*	GPS DOP and Active Satellites	*/
class gpgsa{				//  GPGSA  struct

public:
		gpgsa();

public:
	CString		GP_ID;
	CString		mode;		//	 A=auto;M=manual
	int			pos_dim;		//	postion dimension 1=invalid;2=2D;3=3D
	int			PRN[MAXNUMGSA];
	double		pdop;			
	double		hdop;
	double		vdop;
	int			validnum;
} ;//gngsa,gbgsn;;

/*	GPS Satellites in View	*/
class sate_info{		// 1 satellite information
public:
	sate_info();
public:
	int		PRN;
	double	ele;		//elevation in degrees
	double azi	;		//azimuth in degrees
	double	snr;		//signal noise ratio dBHZ 0-99
};
class gpgsv{

public:
	gpgsv();
public:
	CString		GP_ID;
	int			numsat;
	sate_info sat[MAXNUMGSV];	//
};


/*	Recommended Minimum Specific GPS/TRANSIT Data 	*/
class gprmc{		//GPRMC

public:
	gprmc();

public:
	CString		GP_ID;
	double utc;				//	1.hhmmss.ssss
	CString		sta;			// 2.status A=valid;V=invalid
	double lat;				//	3.Latitude:ddmm.mmmm
	CString		ns_ind;		//	4.N/S indicator
	double lon;				//	5.Longtitude:dddmm.mmmm
	CString		ew_ind;		//	6.E/W indicator
	double	vel;				//	7.velocity over ground:knots
	double cou;				//	8.course over ground:degrees
	int	utc_date;//9.UTC Date:DDMMYY		
	double mag_var;		// 10.magnetic varitation
	CString		mv_ind;		//	11.dirction of magnetic variation E/W
} ;


/*	Recommended Minimum Specific GPS/TRANSIT Data for 3.0	*/
class gprmc3{		//GPRMC for 3.0

public:
	CString		GP_ID;
	double utc;				//	hhmmss.ssss
	CString		sta;			//status A=valid;V=invalid
	double lat;				//	Latitude:ddmm.mmmm
	CString		ns_ind;		//	N/S indicator
	double lon;				//	Longtitude:dddmm.mmmm
	CString		ew_ind;		//	E/W indicator
	double	vel;				//	velocity over ground:knots
	double cou;				//	course over ground:degrees
	int	utc_date;			//UTC Date:DDMMYY		
	double mag_var;		// magnetic varitation
	CString		mv_ind;		//	dirction of magnetic variation E/W
	CString		mode_ind;		//A=auto;D=difference;E=estimation;N=null
};

/*	Track Made Good and Ground Speed	*/
class gpvtg{		//GPVTG
public:
	CString		GP_ID;
	double	cou1;		//1	course in degrees
	CString		ref1;		//2	reference T=true heading;M=magnetic heading;
	double	cou2;		//3
	CString		ref2;		//4
	double	vel1;			//5horizontal speed
	CString		unit1;	//6N=knots;K=km/h
	double	vel2;			//7
	CString		unit2;	//8
} ;

/*	Track Made Good and Ground Speed	for 3.0*/
class gpvtg3{		//GPVTG for 3.0
public:
	CString		GP_ID;
	double	cou1;	//	course in degrees
	CString		ref1;		//reference T=true heading;M=magnetic heading;
	double	cou2;
	CString		ref2;
	double	vel1;		//horizontal speed
	CString		unit1;//N=knots;K=km/h
	double	vel2;		
	CString		unit2;
	CString		mode_ind;		//A=auto;D=difference;E=estimation;N=null
};

/*	Geographic position-latitude/longitude	*/
class gpgll{		//GPGLL
public:
	CString		GP_ID;
	double lat;				//	Latitude:ddmm.mmmm
	CString		ns_ind;		//	N/S indicator
	double lon;				//	Longtitude:dddmm.mmmm
	CString		ew_ind;		//	E/W indicator
	double utc;				//	hhnnss.ssss
	CString		sta;				// Status A=valid;V=invalid;
};

/*	UTC and local date/time data	*/
class gpzda{		//GPZDA
public:
	CString		GP_ID;
	double		time;		//hhmmss.ss
	int			day;
	int			month;
	int			year;
	int			loc_h; //offset to local time zone in hours
	int			loc_m;// offset to local in minutes
};

/*	GPS pseudorange noise statistics	*/
class gpgst{				//GPGST
public:
	CString		GP_ID;
	double		utc;			//hhmmss.ss
	double		rms;
	double		std_maj;	//STD (meters) of semi-major axis of error ellipse
	double		std_min;	//STD (meters) of semi-minor axis of error ellipse
	double		maj_orien;//Orientation of semi-major axis of error ellipse (true north:degrees)
	double		std_lat;		//STD (meters) of latitude error
	double		std_lon;		//STD (meters) of longitude error
	double		std_alt;		//STD (meters) of altitude error
};

//invalid    data
//		Type						double/int								CString	
//		Value						-1												' '
//						(-100 in gpzda.loc_h/loc_m)



class NMEA_process
{
public:
	bool read_data(CString filename);

	bool check_sum(CString line,int &flag);
	bool Get_gpgga(CString line,gpgga &gga);
	bool Get_gpgsa(CString line,gpgsa &gsa);
	bool Get_gpgsv(CString line,gpgsv &gsv);
	bool Get_gprmc(CString line,gprmc &rmc);
	bool Get_gprmc3(CString line,gprmc3 &rmc3);
	bool Get_gpvtg(CString line,gpvtg &rmc);
	bool Get_gpvtg3(CString line,gpvtg3 &rmc);
	bool Get_gpgll(CString line,gpgll &gll);
	bool Get_gpzda(CString line,gpzda &zda);
	bool Get_gpgst(CString line,gpgst &gst);
	bool Get_gngns(CString line,gngns &gns);
	//bool Calc_para(gpgsv gsv,rtcm diff);
private:
	gpgga		gga_pub;				//	1
	gpgsa		gsa_pub;				//	2
	gpgsv		gsv_pub;				//	3
	gprmc		rmc_pub;				//	4
	gprmc3		rmc3_pub;				//	5
	gpvtg		vtg_pub;				//	6
	gpvtg3		vtg3_pub;				//	7
	gpgll		gll_pub;				//	8
	gpzda		zda_pub;				//	9
	gpgst		gst_pub;				//	10
	gngns		gns_pub;				//	11
};
	//	more in the future




//-----------------------------------------------------------------------------------



class SingleCorrInfo{

public:
	SingleCorrInfo();

public:
	int			sysid;	//	1=GPS;5=BDS;
	int			prn;
	double		corrCom;
	double		corrIono;//com=iono+trop
	double		corrTrop;//
	double		corrDot;// not support in FAN_corr
};
class EpochCorrInfo{

public:
	EpochCorrInfo();

public:
	int						GPSW;
	double					GPSS;
	int						NumSat;
	SingleCorrInfo		CorrSatInfo[MAXCORSNUMSAT];
};
	
class process_CORS
{
public:
	bool read_rtcm(CString filename,EpochCorrInfo *SatCorr);

};


//--------------------------------------------------------------------------------------


class GNSSTBext
{
public:
	//--------TimeTrans--------
	int			DayofYear(int GpsWeek,double GpsSec);
    void			Gps_Week_Sec(int* Gps_Week,double* Gps_Second,int Y,int M,int D,int hour,int minute,double second);
	bool			GPStoYMDHMS(int GPSWeek,double GPSSecond,int& Year,int& Month,int& Day,int& Hour,int& Minute,double& Second,
								double& JD,double& MJD,int& Weekday);
				
	bool			 BLH2XYZ(double *BLH,double *XYZ);
	void			 GPFormat2time1(double gp_utc,int& utc_hour,int& utc_min,double &utc_sec);
	void			 GPFormat2time2(int gp_utc,int& utc_day,int& utc_mon,int& utc_year);
	double		 GPFormat2Degree(double gpl);
};






