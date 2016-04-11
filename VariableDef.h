/*This file is used to define constant variables

	const is unchangeable£¬so ¡°#define¡±is recommended


---------------------------------------------------------------------------*/
 

#define ROUND(x)			(floor(x+0.5))
#define SQ(x)					(pow(x,2))
#define SWAP(x,y)				do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define SGN(x)					((x)<=0.0?-1.0:1.0)
#define LOOPMAX			10000		/* max loop of lambda */

#define NFREQ					3
#define NEXOBS				0
#define MAXTYPE2			45// max number of obstypes in 2.x
#define MAXTYPE3			26// max number of obstypes in 3.x  single system 
#define MAXNUMSATE		35     /* max sate number of gps   32*/
#define MAXNUMSATE_BDS	  35
#define MAXNUMSATE_GLO	 30
#define MAXNUMSATE_GAL 33
#define MAXNUMSATE_SBAS 10
#define MAXOBSEPOCH	30

#define CLIGHT				299792458.0					//speed of light
#define PI							3.1415926535897932		
#define D2R						(PI/180)						//deg 2 rad
#define R2D						(180/PI)						//rad 2 deg

#define RE_WGS84			6378137.0						// earth semimajor axis (WGS84) (m) 
#define FE_WGS84			(1.0/298.257223563)	   // earth flattening (WGS84)
#define OMEGAE				7.2921151467E-5			//earth angular velocity GPS  rad/s
#define GM						3.986005E14					//earth gravitational constant
#define ee						0.00669437999013			//eccentricity ^2

#define RE_CGCS2000		6378137.0						//	the parameters of BDS and GAL are same,see reference
#define FE_CGCS2000		(1.0/298.257222101)
#define OMEGAE_BDS		7.2921151467E-5			//earth angular velocity BDS
#define GM_BDS				3.986004418E14			//
#define ee_BDS					0.006694380022901	

#define RE_GLO				6378136.0        /* radius of earth (m)             */
#define MU_GLO				3.9860044E14     /* gravitational constant          */
#define FE_Glo					(1.0/298.25784)
#define J2_GLO					1.0826257E-3     /* 2nd zonal harmonic of geopot   ] */
#define OMGE_GLO			7.292115E-5      /* earth angular velocity (rad/s) */
#define ee_GLO				0.006694366177482  
#define TSTEP					60.0             /* integration step glonass ephemeris (s) */

#define SIN_5					-0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5					 0.9961946980917456 /* cos(-5.0 deg) */

#define RE_GAL				6378137.0
#define FE_GAL					(1.0/298.257223563)
#define MU_GAL				3.986004418E14
#define OMEGA_GAL		7.2921151467E-5
#define ee_GAL					0.00669437999013

#define ERREPH_GLO		5.0            /* error of glonass ephemeris (m) */
#define RTOL_KEPLER		1E-14         /* relative tolerance for Kepler equation */

#define  INFINUMBER		1.0E10

#define RE_BJ54				6378245.0
#define FE_BJ54				(1.0/298.3)
#define ee_BJ54				0.006693421622966

#define RE_XA80				6378140.0
#define FE_XA80				(1.0/298.257)
#define ee_XA80				0.006694384999588

#define HIONOS				350000.0						// model height of ionosphere



#define FREQ1					1.57542E9           /* L1/E1  frequency (Hz) */
#define FREQ2					1.22760E9           /* L2     frequency (Hz) */
#define FREQ5					1.17645E9           /* L5/E5a frequency (Hz) */
#define FREQ6					1.27875E9           /* E6/LEX frequency (Hz) */
#define FREQ7					1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8					1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ1_GLO			1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO			0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO			1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO			0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO			1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ2_BDS			1.561098E9          /* BeiDou B1 frequency (Hz) */
#define FREQ7_BDS			1.20714E9           /* BeiDou B2 frequency (Hz) */
#define FREQ6_BDS			1.26852E9           /* BeiDou B3 frequency (Hz) */

#define FREQ1RATIO			1.57542           /* L1/E1  frequency (Hz) */
#define FREQ2RATIO			1.22760           /* L2     frequency (Hz) */
#define FREQ5RATIO			1.17645           /* L5/E5a frequency (Hz) */
#define FREQ6RATIO			1.27875           /* E6/LEX frequency (Hz) */
#define FREQ7RATIO			1.20714           /* E5b    frequency (Hz) */
#define FREQ8RATIO			1.191795          /* E5a+b  frequency (Hz) */
#define FREQ1_GLORATIO	1.60200           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLORATIO   0.56250E-3           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLORATIO   1.24600          /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLORATIO   0.43750E-3           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLORATIO   1.202025         /* GLONASS G3 frequency (Hz) */
#define FREQ2_BDSRATIO   1.561098          /* BeiDou B1 frequency (Hz) */
#define FREQ7_BDSRATIO   1.20714           /* BeiDou B2 frequency (Hz) */
#define FREQ6_BDSRATIO   1.26852           /* BeiDou B3 frequency (Hz) */

#define  MAXSATERISE			10					/* the max of rising number*/
#define  MAXSATEFALL			10					/* the max of falling number*/

#define  INTERPOLYORDER			10		    /* order of interpolation */

#define MAXSYS						6				/*	max number of system   */
#define MAXEPHNUM				5000			/* max number of GPS+GAL+BDS+QZSS*/
#define MAXEPHNUM_GLO		1500
#define MAXEPHNUM_SBAS		200		

#define MAXFILEEPH				10				/*max number of ephemeris files  */
#define MAXFILESP3					10
#define MAXFILECLK				10
#define MAXFILEBASE				10
#define MAXFILEROVER			10



/*  for DEBUG  */

#define  DEBUGTEST		1

#define OUTSPP			0

#define FileOutSpp		1