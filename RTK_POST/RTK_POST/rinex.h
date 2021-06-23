#ifndef RINEX_H
#define RINEX_H

#include "RTK.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*------------------------------------------------------------------------------
* rtcm3.c :rtkcmn functions
*
*          Copyright (C) 2021 by linjun,HuiTian Xpeng, All rights reserved.
* references:rtklib
* history : 2021/05/19 1.0  new
*-----------------------------------------------------------------------------*/
#define _POSIX_C_SOURCE 199506
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>
//#include "rtklib.h"

/* constants -----------------------------------------------------------------*/
#define MAXFREQ     7                   /* max NFREQ */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define INT_SWAP_TRAC 86400.0           /* swap interval of trace file (s) */
#define INT_SWAP_STAT 86400.0           /* swap interval of solution status file (s) */

#define MAXEXFILE   1024                /* max number of expanded files */
#define MAXSBSAGEF  30.0                /* max age of SBAS fast correction (s) */
#define MAXSBSAGEL  1800.0              /* max age of SBAS long term corr (s) */
#define MAXSBSURA   8                   /* max URA of SBAS satellite */
#define MAXBAND     10                  /* max SBAS band of IGP */
#define MAXNIGP     201                 /* max number of IGP in SBAS band */
#define MAXNGEO     4                   /* max number of GEO satellites */
#define MAXCOMMENT  100                 /* max number of RINEX comments */
#define MAXSTRPATH  1024                /* max length of stream path */
#define MAXSTRMSG   1024                /* max length of stream message */
#define MAXSTRRTK   8                   /* max number of stream in RTK server */
#define MAXSBSMSG   32                  /* max number of SBAS msg in RTK server */
#define MAXSOLMSG   8191                /* max length of solution message */
#define MAXRAWLEN   16384               /* max length of receiver raw message */
#define MAXERRMSG   4096                /* max length of error/warning message */
#define MAXANT      64                  /* max length of station name/antenna type */
#define MAXSOLBUF   256                 /* max number of solution buffer */
#define MAXOBSBUF   128                 /* max number of observation data buffer */
#define MAXNRPOS    16                  /* max number of reference positions */
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define MAXGISLAYER 32                  /* max number of GIS data layers */
#define MAXRCVCMD   4096                /* max length of receiver commands */

#define RNX2VER     2.10                /* RINEX ver.2 default output version */
#define RNX3VER     3.00                /* RINEX ver.3 default output version */

#define OBSTYPE_PR  0x01                /* observation type: pseudorange */
#define OBSTYPE_CP  0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP 0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR 0x08                /* observation type: SNR */
#define OBSTYPE_ALL 0xFF                /* observation type: all */

#define FREQTYPE_L1 0x01                /* frequency type: L1/E1/B1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/E5b/B2 */
#define FREQTYPE_L3 0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_L4 0x08                /* frequency type: L6/E6/B3 */
#define FREQTYPE_L5 0x10                /* frequency type: E5ab */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
#define TSYS_IRN    6                   /* time system: IRNSS time */

#define POLYCRC32   0xEDB88320u /* CRC32 polynomial */
#define POLYCRC24Q  0x1864CFBu  /* CRC24Q polynomial */

#define MAX_VAR_EPH SQR(300.0)  /* max variance eph to reject satellite (m^2) */

#ifdef OBS_100HZ
#define DTTOL       0.005               /* tolerance of time difference (s) */
#else
#define DTTOL       0.025               /* tolerance of time difference (s) */
#endif

#define LLI_SLIP    0x01                /* LLI: cycle-slip */
#define LLI_HALFC   0x02                /* LLI: half-cycle not resovled */
#define LLI_BOCTRK  0x04                /* LLI: boc tracking of mboc signal */
#define LLI_HALFA   0x40                /* LLI: half-cycle added */
#define LLI_HALFS   0x80                /* LLI: half-cycle subtracted */
/* type definitions ----------------------------------------------------------*/

typedef struct {        /* time struct */
	time_t time;        /* time (s) expressed by standard time_t */
	double sec;         /* fraction of second under 1 s */
} gtime_t;

typedef struct {        /* observation data record */
	gtime_t time;       /* receiver sampling time (GPST) */
	uint8_t sat, rcv;    /* satellite/receiver number */
	uint16_t SNR[NFREQ + NEXOBS]; /* signal strength (0.001 dBHz) */
	uint8_t  LLI[NFREQ + NEXOBS]; /* loss of lock indicator */
	uint8_t code[NFREQ + NEXOBS]; /* code indicator (CODE_???) */
	double L[NFREQ + NEXOBS]; /* observation data carrier-phase (cycle) */
	double P[NFREQ + NEXOBS]; /* observation data pseudorange (m) */
	float  D[NFREQ + NEXOBS]; /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {        /* observation data */
	int n, nmax;         /* number of obervation data/allocated */
	obsd_t *data;       /* observation data records */
} obs_t;

typedef struct {        /* earth rotation parameter data type */
	double mjd;         /* mjd (days) */
	double xp, yp;       /* pole offset (rad) */
	double xpr, ypr;     /* pole offset rate (rad/day) */
	double ut1_utc;     /* ut1-utc (s) */
	double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
	int n, nmax;         /* number and max number of data */
	erpd_t *data;       /* earth rotation parameter data */
} erp_t;

typedef struct {        /* antenna parameter type */
	int sat;            /* satellite number (0:receiver) */
	char type[MAXANT];  /* antenna type */
	char code[MAXANT];  /* serial number or satellite code */
	gtime_t ts, te;      /* valid time start and end */
	double off[NFREQ][3]; /* phase center offset e/n/u or x/y/z (m) */
	double var[NFREQ][19]; /* phase center variation (m) */
	/* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} pcv_t;

typedef struct {        /* antenna parameters type */
	int n, nmax;         /* number of data/allocated */
	pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

typedef struct {        /* almanac type */
	int sat;            /* satellite number */
	int svh;            /* sv health (0:ok) */
	int svconf;         /* as and sv config */
	int week;           /* GPS/QZS: gps week, GAL: galileo week */
	gtime_t toa;        /* Toa */
	/* SV orbit parameters */
	double A, e, i0, OMG0, omg, M0, OMGd;
	double toas;        /* Toa (s) in week */
	double f0, f1;       /* SV clock parameters (af0,af1) */
} alm_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
	int sat;            /* satellite number */
	int iode, iodc;      /* IODE,IODC */
	int sva;            /* SV accuracy (URA index) */
	int svh;            /* SV health (0:ok) */
	int week;           /* GPS/QZS: gps week, GAL: galileo week */
	int code;           /* GPS/QZS: code on L2 */
	/* GAL: data source defined as rinex 3.03 */
	/* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
	int flag;           /* GPS/QZS: L2 P data flag */
	/* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
	gtime_t toe, toc, ttr; /* Toe,Toc,T_trans */
	/* SV orbit parameters */
	double A, e, i0, OMG0, omg, M0, deln, OMGd, idot;
	double crc, crs, cuc, cus, cic, cis;
	double toes;        /* Toe (s) in week */
	double fit;         /* fit interval (h) */
	double f0, f1, f2;    /* SV clock parameters (af0,af1,af2) */
	double tgd[6];      /* group delay parameters */
	/* GPS/QZS:tgd[0]=TGD */
	/* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
	/* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
	/*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
	double Adot, ndot;   /* Adot,ndot for CNAV */
} eph_t;

typedef struct {        /* GLONASS broadcast ephemeris type */
	int sat;            /* satellite number */
	int iode;           /* IODE (0-6 bit of tb field) */
	int frq;            /* satellite frequency number */
	int svh, sva, age;    /* satellite health, accuracy, age of operation */
	gtime_t toe;        /* epoch of epherides (gpst) */
	gtime_t tof;        /* message frame time (gpst) */
	double pos[3];      /* satellite position (ecef) (m) */
	double vel[3];      /* satellite velocity (ecef) (m/s) */
	double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
	double taun, gamn;   /* SV clock bias (s)/relative freq bias */
	double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

typedef struct {        /* precise ephemeris type */
	gtime_t time;       /* time (GPST) */
	int index;          /* ephemeris index for multiple files */
	double pos[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
	float  std[MAXSAT][4]; /* satellite position/clock std (m|s) */
	double vel[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
	float  vst[MAXSAT][4]; /* satellite velocity/clk-rate std (m/s|s/s) */
	float  cov[MAXSAT][3]; /* satellite position covariance (m^2) */
	float  vco[MAXSAT][3]; /* satellite velocity covariance (m^2) */
} peph_t;

typedef struct {        /* precise clock type */
	gtime_t time;       /* time (GPST) */
	int index;          /* clock index for multiple files */
	double clk[MAXSAT][1]; /* satellite clock (s) */
	float  std[MAXSAT][1]; /* satellite clock std (s) */
} pclk_t;

typedef struct {        /* SBAS ephemeris type */
	int sat;            /* satellite number */
	gtime_t t0;         /* reference epoch time (GPST) */
	gtime_t tof;        /* time of message frame (GPST) */
	int sva;            /* SV accuracy (URA index) */
	int svh;            /* SV health (0:ok) */
	double pos[3];      /* satellite position (m) (ecef) */
	double vel[3];      /* satellite velocity (m/s) (ecef) */
	double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
	double af0, af1;     /* satellite clock-offset/drift (s,s/s) */
} seph_t;

typedef struct {        /* NORAL TLE data type */
	char name[32];     /* common name */
	char alias[32];     /* alias name */
	char satno[16];     /* satellilte catalog number */
	char satclass;      /* classification */
	char desig[16];     /* international designator */
	gtime_t epoch;      /* element set epoch (UTC) */
	double ndot;        /* 1st derivative of mean motion */
	double nddot;       /* 2st derivative of mean motion */
	double bstar;       /* B* drag term */
	int etype;          /* element set type */
	int eleno;          /* element number */
	double inc;         /* orbit inclination (deg) */
	double OMG;         /* right ascension of ascending node (deg) */
	double ecc;         /* eccentricity */
	double omg;         /* argument of perigee (deg) */
	double M;           /* mean anomaly (deg) */
	double n;           /* mean motion (rev/day) */
	int rev;            /* revolution number at epoch */
} tled_t;

typedef struct {        /* NORAD TLE (two line element) type */
	int n, nmax;         /* number/max number of two line element data */
	tled_t *data;       /* NORAD TLE data */
} tle_t;

typedef struct {        /* TEC grid type */
	gtime_t time;       /* epoch time (GPST) */
	int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
	double rb;          /* earth radius (km) */
	double lats[3];     /* latitude start/interval (deg) */
	double lons[3];     /* longitude start/interval (deg) */
	double hgts[3];     /* heights start/interval (km) */
	double *data;       /* TEC grid data (tecu) */
	float *rms;         /* RMS values (tecu) */
} tec_t;

typedef struct {        /* SBAS message type */
	int week, tow;       /* receiption time */
	uint8_t prn, rcv;    /* SBAS satellite PRN,receiver number */
	uint8_t msg[29];    /* SBAS message (226bit) padded by 0 */
} sbsmsg_t;

typedef struct {        /* SBAS messages type */
	int n, nmax;         /* number of SBAS messages/allocated */
	sbsmsg_t *msgs;     /* SBAS messages */
} sbs_t;

typedef struct {        /* SBAS fast correction type */
	gtime_t t0;         /* time of applicability (TOF) */
	double prc;         /* pseudorange correction (PRC) (m) */
	double rrc;         /* range-rate correction (RRC) (m/s) */
	double dt;          /* range-rate correction delta-time (s) */
	int iodf;           /* IODF (issue of date fast corr) */
	int16_t udre;       /* UDRE+1 */
	int16_t ai;         /* degradation factor indicator */
} sbsfcorr_t;

typedef struct {        /* SBAS long term satellite error correction type */
	gtime_t t0;         /* correction time */
	int iode;           /* IODE (issue of date ephemeris) */
	double dpos[3];     /* delta position (m) (ecef) */
	double dvel[3];     /* delta velocity (m/s) (ecef) */
	double daf0, daf1;   /* delta clock-offset/drift (s,s/s) */
} sbslcorr_t;

typedef struct {        /* SBAS satellite correction type */
	int sat;            /* satellite number */
	sbsfcorr_t fcorr;   /* fast correction */
	sbslcorr_t lcorr;   /* long term correction */
} sbssatp_t;

typedef struct {        /* SBAS satellite corrections type */
	int iodp;           /* IODP (issue of date mask) */
	int nsat;           /* number of satellites */
	int tlat;           /* system latency (s) */
	sbssatp_t sat[MAXSAT]; /* satellite correction */
} sbssat_t;

typedef struct {        /* SBAS ionospheric correction type */
	gtime_t t0;         /* correction time */
	int16_t lat, lon;    /* latitude/longitude (deg) */
	int16_t give;       /* GIVI+1 */
	float delay;        /* vertical delay estimate (m) */
} sbsigp_t;

typedef struct {        /* IGP band type */
	int16_t x;          /* longitude/latitude (deg) */
	const int16_t *y;   /* latitudes/longitudes (deg) */
	uint8_t bits;       /* IGP mask start bit */
	uint8_t bite;       /* IGP mask end bit */
} sbsigpband_t;

typedef struct {        /* SBAS ionospheric corrections type */
	int iodi;           /* IODI (issue of date ionos corr) */
	int nigp;           /* number of igps */
	sbsigp_t igp[MAXNIGP]; /* ionospheric correction */
} sbsion_t;

typedef struct {        /* DGPS/GNSS correction type */
	gtime_t t0;         /* correction time */
	double prc;         /* pseudorange correction (PRC) (m) */
	double rrc;         /* range rate correction (RRC) (m/s) */
	int iod;            /* issue of data (IOD) */
	double udre;        /* UDRE */
} dgps_t;

typedef struct {        /* SSR correction type */
	gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
	double udi[6];      /* SSR update interval (s) */
	int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
	int iode;           /* issue of data */
	int iodcrc;         /* issue of data crc for beidou/sbas */
	int ura;            /* URA indicator */
	int refd;           /* sat ref datum (0:ITRF,1:regional) */
	double deph[3];    /* delta orbit {radial,along,cross} (m) */
	double ddeph[3];    /* dot delta orbit {radial,along,cross} (m/s) */
	double dclk[3];    /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
	double hrclk;       /* high-rate clock corection (m) */
	float  cbias[MAXCODE]; /* code biases (m) */
	double pbias[MAXCODE]; /* phase biases (m) */
	float  stdpb[MAXCODE]; /* std-dev of phase biases (m) */
	double yaw_ang, yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
	uint8_t update;     /* update flag (0:no update,1:update) */
} ssr_t;

typedef struct {        /* navigation data type */
	int n, nmax;         /* number of broadcast ephemeris */
	int ng, ngmax;       /* number of glonass ephemeris */
	int ns, nsmax;       /* number of sbas ephemeris */
	int ne, nemax;       /* number of precise ephemeris */
	int nc, ncmax;       /* number of precise clock */
	int na, namax;       /* number of almanac data */
	int nt, ntmax;       /* number of tec grid data */
	eph_t *eph;         /* GPS/QZS/GAL/BDS/IRN ephemeris */
	geph_t *geph;       /* GLONASS ephemeris */
	seph_t *seph;       /* SBAS ephemeris */
	peph_t *peph;       /* precise ephemeris */
	pclk_t *pclk;       /* precise clock */
	alm_t *alm;         /* almanac data */
	tec_t *tec;         /* tec grid data */
	erp_t  erp;         /* earth rotation parameters */
	double utc_gps[8];  /* GPS delta-UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
	double utc_glo[8];  /* GLONASS UTC time parameters {tau_C,tau_GPS} */
	double utc_gal[8];  /* Galileo UTC parameters */
	double utc_qzs[8];  /* QZS UTC parameters */
	double utc_cmp[8];  /* BeiDou UTC parameters */
	double utc_irn[9];  /* IRNSS UTC parameters {A0,A1,Tot,...,dt_LSF,A2} */
	double utc_sbs[4];  /* SBAS UTC parameters */
	double ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	double ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
	double ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	double ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	int glo_fcn[32];    /* GLONASS FCN + 8 */
	double cbias[MAXSAT][3]; /* satellite DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
	double rbias[MAXRCV][2][3]; /* receiver DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
	pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
	sbssat_t sbssat;    /* SBAS satellite corrections */
	sbsion_t sbsion[MAXBAND + 1]; /* SBAS ionosphere corrections */
	dgps_t dgps[MAXSAT]; /* DGPS corrections */
	ssr_t ssr[MAXSAT];  /* SSR corrections */
} nav_t;

typedef struct {        /* station parameter type */
	char name[MAXANT]; /* marker name */
	char marker[MAXANT]; /* marker number */
	char antdes[MAXANT]; /* antenna descriptor */
	char antsno[MAXANT]; /* antenna serial number */
	char rectype[MAXANT]; /* receiver type descriptor */
	char recver[MAXANT]; /* receiver firmware version */
	char recsno[MAXANT]; /* receiver serial number */
	int antsetup;       /* antenna setup id */
	int itrf;           /* ITRF realization year */
	int deltype;        /* antenna delta type (0:enu,1:xyz) */
	double pos[3];      /* station position (ecef) (m) */
	double del[3];      /* antenna position delta (e/n/u or x/y/z) (m) */
	double hgt;         /* antenna height (m) */
	int glo_cp_align;   /* GLONASS code-phase alignment (0:no,1:yes) */
	double glo_cp_bias[4]; /* GLONASS code-phase biases {1C,1P,2C,2P} (m) */
} sta_t;

typedef struct {        /* solution type */
	gtime_t time;       /* time (GPST) */
	double rr[6];       /* position/velocity (m|m/s) */
	/* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
	float  qr[6];       /* position variance/covariance (m^2) */
	/* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
	/* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
	float  qv[6];       /* velocity variance/covariance (m^2/s^2) */
	double dtr[6];      /* receiver clock bias to time systems (s) */
	uint8_t type;       /* type (0:xyz-ecef,1:enu-baseline) */
	uint8_t stat;       /* solution status (SOLQ_???) */
	uint8_t ns;         /* number of valid satellites */
	float age;          /* age of differential (s) */
	float ratio;        /* AR ratio factor for valiation */
	float thres;        /* AR ratio threshold for valiation */
} sol_t;

typedef struct {        /* solution buffer type */
	int n, nmax;         /* number of solution/max number of buffer */
	int cyclic;         /* cyclic buffer flag */
	int start, end;      /* start/end index */
	gtime_t time;       /* current solution time */
	sol_t *data;        /* solution data */
	double rb[3];       /* reference position {x,y,z} (ecef) (m) */
	uint8_t buff[MAXSOLMSG + 1]; /* message buffer */
	int nb;             /* number of byte in message buffer */
} solbuf_t;

typedef struct {        /* solution status type */
	gtime_t time;       /* time (GPST) */
	uint8_t sat;        /* satellite number */
	uint8_t frq;        /* frequency (1:L1,2:L2,...) */
	float az, el;        /* azimuth/elevation angle (rad) */
	float resp;         /* pseudorange residual (m) */
	float resc;         /* carrier-phase residual (m) */
	uint8_t flag;       /* flags: (vsat<<5)+(slip<<3)+fix */
	uint16_t snr;       /* signal strength (*SNR_UNIT dBHz) */
	uint16_t lock;      /* lock counter */
	uint16_t outc;      /* outage counter */
	uint16_t slipc;     /* slip counter */
	uint16_t rejc;      /* reject counter */
} solstat_t;

typedef struct {        /* solution status buffer type */
	int n, nmax;         /* number of solution/max number of buffer */
	solstat_t *data;    /* solution status data */
} solstatbuf_t;

typedef struct {        /* RTCM control struct type */
	int staid;          /* station id */
	int stah;           /* station health */
	int seqno;          /* sequence number for rtcm 2 or iods msm */
	int outtype;        /* output message type */
	gtime_t time;       /* message time */
	gtime_t time_s;     /* message start time */
	obs_t obs;          /* observation data (uncorrected) */
	nav_t nav;          /* satellite ephemerides */
	sta_t sta;          /* station parameters */
	dgps_t *dgps;       /* output of dgps corrections */
	ssr_t ssr[MAXSAT];  /* output of ssr corrections */
	char msg[128];      /* special message */
	char msgtype[256];  /* last message type */
	char msmtype[7][128]; /* msm signal types */
	int obsflag;        /* obs data complete flag (1:ok,0:not complete) */
	int ephsat;         /* input ephemeris satellite number */
	int ephset;         /* input ephemeris set (0-1) */
	double cp[MAXSAT][NFREQ + NEXOBS]; /* carrier-phase measurement */
	uint16_t lock[MAXSAT][NFREQ + NEXOBS]; /* lock time */
	uint16_t loss[MAXSAT][NFREQ + NEXOBS]; /* loss of lock count */
	gtime_t lltime[MAXSAT][NFREQ + NEXOBS]; /* last lock time */
	int nbyte;          /* number of bytes in message buffer */
	int nbit;           /* number of bits in word buffer */
	int len;            /* message length (bytes) */
	uint8_t buff[1200]; /* message buffer */
	uint32_t word;      /* word buffer for rtcm 2 */
	uint32_t nmsg2[100]; /* message count of RTCM 2 (1-99:1-99,0:other) */
	uint32_t nmsg3[400]; /* message count of RTCM 3 (1-299:1001-1299,300-329:4070-4099,0:ohter) */
	char opt[256];      /* RTCM dependent options */
} rtcm_t;

typedef struct {        /* RINEX control struct type */
	gtime_t time;       /* message time */
	double ver;         /* RINEX version */
	char   type;        /* RINEX file type ('O','N',...) */
	int    sys;         /* navigation system */
	int    tsys;        /* time system */
	char   tobs[8][MAXOBSTYPE][4]; /* rinex obs types */
	obs_t  obs;         /* observation data */
	nav_t  nav;         /* navigation data */
	sta_t  sta;         /* station info */
	int    ephsat;      /* input ephemeris satellite number */
	int    ephset;      /* input ephemeris set (0-1) */
	char   opt[256];    /* rinex dependent options */
} rnxctr_t;

typedef struct {        /* download URL type */
	char type[32];      /* data type */
	char path[1024];    /* URL path */
	char dir[1024];    /* local directory */
	double tint;        /* time interval (s) */
} url_t;

typedef struct {        /* option type */
	const char *name;   /* option name */
	int format;         /* option format (0:int,1:double,2:string,3:enum) */
	void *var;          /* pointer to option variable */
	const char *comment; /* option comment/enum labels/unit */
} opt_t;

typedef struct {        /* SNR mask type */
	int ena[2];         /* enable flag {rover,base} */
	double mask[NFREQ][9]; /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;

typedef struct {        /* processing options type */
	int mode;           /* positioning mode (PMODE_???) */
	int soltype;        /* solution type (0:forward,1:backward,2:combined) */
	int nf;             /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
	int navsys;         /* navigation system */
	double elmin;       /* elevation mask angle (rad) */
	snrmask_t snrmask;  /* SNR mask */
	int sateph;         /* satellite ephemeris/clock (EPHOPT_???) */
	int modear;         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
	int glomodear;      /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
	int bdsmodear;      /* BeiDou AR mode (0:off,1:on) */
	int maxout;         /* obs outage count to reset bias */
	int minlock;        /* min lock count to fix ambiguity */
	int minfix;         /* min fix count to hold ambiguity */
	int armaxiter;      /* max iteration to resolve ambiguity */
	int ionoopt;        /* ionosphere option (IONOOPT_???) */
	int tropopt;        /* troposphere option (TROPOPT_???) */
	int dynamics;       /* dynamics model (0:none,1:velociy,2:accel) */
	int tidecorr;       /* earth tide correction (0:off,1:solid,2:solid+otl+pole) */
	int niter;          /* number of filter iteration */
	int codesmooth;     /* code smoothing window size (0:none) */
	int intpref;        /* interpolate reference obs (for post mission) */
	int sbascorr;       /* SBAS correction options */
	int sbassatsel;     /* SBAS satellite selection (0:all) */
	int rovpos;         /* rover position for fixed mode */
	int refpos;         /* base position for relative mode */
	/* (0:pos in prcopt,  1:average of single pos, */
	/*  2:read from file, 3:rinex header, 4:rtcm pos) */
	double eratio[NFREQ]; /* code/phase error ratio */
	double err[5];      /* measurement error factor */
	/* [0]:reserved */
	/* [1-3]:error factor a/b/c of phase (m) */
	/* [4]:doppler frequency (hz) */
	double std[3];      /* initial-state std [0]bias,[1]iono [2]trop */
	double prn[6];      /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv [5] pos */
	double sclkstab;    /* satellite clock stability (sec/sec) */
	double thresar[8];  /* AR validation threshold */
	double elmaskar;    /* elevation mask of AR for rising satellite (deg) */
	double elmaskhold;  /* elevation mask to hold ambiguity (deg) */
	double thresslip;   /* slip threshold of geometry-free phase (m) */
	double maxtdiff;    /* max difference of time (sec) */
	double maxinno;     /* reject threshold of innovation (m) */
	double maxgdop;     /* reject threshold of gdop */
	double baseline[2]; /* baseline length constraint {const,sigma} (m) */
	double ru[3];       /* rover position for fixed mode {x,y,z} (ecef) (m) */
	double rb[3];       /* base position for relative mode {x,y,z} (ecef) (m) */
	char anttype[2][MAXANT]; /* antenna types {rover,base} */
	double antdel[2][3]; /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
	pcv_t pcvr[2];      /* receiver antenna parameters {rov,base} */
	uint8_t exsats[MAXSAT]; /* excluded satellites (1:excluded,2:included) */
	int  maxaveep;      /* max averaging epoches */
	int  initrst;       /* initialize by restart */
	int  outsingle;     /* output single by dgps/float/fix/ppp outage */
	char rnxopt[2][256]; /* rinex options {rover,base} */
	int  posopt[6];     /* positioning options */
	int  syncsol;       /* solution sync mode (0:off,1:on) */
	double odisp[2][6 * 11]; /* ocean tide loading parameters {rov,base} */
	int  freqopt;       /* disable L2-AR */
	char pppopt[256];   /* ppp option */
} prcopt_t;

typedef struct {        /* solution options type */
	int posf;           /* solution format (SOLF_???) */
	int times;          /* time system (TIMES_???) */
	int timef;          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
	int timeu;          /* time digits under decimal point */
	int degf;           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
	int outhead;        /* output header (0:no,1:yes) */
	int outopt;         /* output processing options (0:no,1:yes) */
	int outvel;         /* output velocity options (0:no,1:yes) */
	int datum;          /* datum (0:WGS84,1:Tokyo) */
	int height;         /* height (0:ellipsoidal,1:geodetic) */
	int geoid;          /* geoid model (0:EGM96,1:JGD2000) */
	int solstatic;      /* solution of static mode (0:all,1:single) */
	int sstat;          /* solution statistics level (0:off,1:states,2:residuals) */
	int trace;          /* debug trace level (0:off,1-5:debug) */
	double nmeaintv[2]; /* nmea output interval (s) (<0:no,0:all) */
	/* nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
	char sep[64];       /* field separator */
	char prog[64];      /* program name */
	double maxsolstd;   /* max std-dev for solution output (m) (0:all) */
} solopt_t;

typedef struct {        /* file options type */
	char satantp[MAXSTRPATH]; /* satellite antenna parameters file */
	char rcvantp[MAXSTRPATH]; /* receiver antenna parameters file */
	char stapos[MAXSTRPATH]; /* station positions file */
	char geoid[MAXSTRPATH]; /* external geoid data file */
	char iono[MAXSTRPATH]; /* ionosphere data file */
	char dcb[MAXSTRPATH]; /* dcb data file */
	char eop[MAXSTRPATH]; /* eop data file */
	char blq[MAXSTRPATH]; /* ocean tide loading blq file */
	char tempdir[MAXSTRPATH]; /* ftp/http temporaly directory */
	char geexe[MAXSTRPATH]; /* google earth exec file */
	char solstat[MAXSTRPATH]; /* solution statistics file */
	char trace[MAXSTRPATH]; /* debug trace file */
} filopt_t;

typedef struct {        /* RINEX options type */
	gtime_t ts, te;      /* time start/end */
	double tint;        /* time interval (s) */
	double ttol;        /* time tolerance (s) */
	double tunit;       /* time unit for multiple-session (s) */
	int rnxver;         /* RINEX version (x100) */
	int navsys;         /* navigation system */
	int obstype;        /* observation type */
	int freqtype;       /* frequency type */
	char mask[7][64];   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
	char staid[32];    /* station id for rinex file name */
	char prog[32];    /* program */
	char runby[32];    /* run-by */
	char marker[64];    /* marker name */
	char markerno[32];  /* marker number */
	char markertype[32]; /* marker type (ver.3) */
	char name[2][32];   /* observer/agency */
	char rec[3][32];   /* receiver #/type/vers */
	char ant[3][32];   /* antenna #/type */
	double apppos[3];   /* approx position x/y/z */
	double antdel[3];   /* antenna delta h/e/n */
	double glo_cp_bias[4]; /* GLONASS code-phase biases (m) */
	char comment[MAXCOMMENT][64]; /* comments */
	char rcvopt[256];   /* receiver dependent options */
	uint8_t exsats[MAXSAT]; /* excluded satellites */
	int glofcn[32];     /* glonass fcn+8 */
	int outiono;        /* output iono correction */
	int outtime;        /* output time system correction */
	int outleaps;       /* output leap seconds */
	int autopos;        /* auto approx position */
	int phshift;        /* phase shift correction */
	int halfcyc;        /* half cycle correction */
	int sep_nav;        /* separated nav files */
	gtime_t tstart;     /* first obs time */
	gtime_t tend;       /* last obs time */
	gtime_t trtcm;      /* approx log start time for rtcm */
	char tobs[7][MAXOBSTYPE][4]; /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
	double shift[7][MAXOBSTYPE]; /* phase shift (cyc) {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
	int nobs[7];        /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
} rnxopt_t;

typedef struct {        /* satellite status type */
	uint8_t sys;        /* navigation system */
	uint8_t vs;         /* valid satellite flag single */
	double azel[2];     /* azimuth/elevation angles {az,el} (rad) */
	double resp[NFREQ]; /* residuals of pseudorange (m) */
	double resc[NFREQ]; /* residuals of carrier-phase (m) */
	uint8_t vsat[NFREQ]; /* valid satellite flag */
	uint16_t snr[NFREQ]; /* signal strength (*SNR_UNIT dBHz) */
	uint8_t fix[NFREQ]; /* ambiguity fix flag (1:fix,2:float,3:hold) */
	uint8_t slip[NFREQ]; /* cycle-slip flag */
	uint8_t half[NFREQ]; /* half-cycle valid flag */
	int lock[NFREQ];   /* lock counter of phase */
	uint32_t outc[NFREQ]; /* obs outage counter of phase */
	uint32_t slipc[NFREQ]; /* cycle-slip counter */
	uint32_t rejc[NFREQ]; /* reject counter */
	double gf[NFREQ - 1]; /* geometry-free phase (m) */
	double mw[NFREQ - 1]; /* MW-LC (m) */
	double phw;         /* phase windup (cycle) */
	gtime_t pt[2][NFREQ]; /* previous carrier-phase time */
	double ph[2][NFREQ]; /* previous carrier-phase observable (cycle) */
} ssat_t;

typedef struct {        /* ambiguity control type */
	gtime_t epoch[4];   /* last epoch */
	int n[4];           /* number of epochs */
	double LC[4];      /* linear combination average */
	double LCv[4];      /* linear combination variance */
	int fixcnt;         /* fix count */
	char flags[MAXSAT]; /* fix flags */
} ambc_t;

typedef struct {        /* RTK control/result type */
	sol_t  sol;         /* RTK solution */
	double rb[6];       /* base position/velocity (ecef) (m|m/s) */
	int nx, na;          /* number of float states/fixed states */
	double tt;          /* time difference between current and previous (s) */
	double *x, *P;      /* float states and their covariance */
	double *xa, *Pa;     /* fixed states and their covariance */
	int nfix;           /* number of continuous fixes of ambiguity */
	ambc_t ambc[MAXSAT]; /* ambibuity control */
	ssat_t ssat[MAXSAT]; /* satellite status */
	int neb;            /* bytes in error message buffer */
	char errbuf[MAXERRMSG]; /* error message buffer */
	prcopt_t opt;       /* processing options */
} rtk_t;

typedef struct {        /* receiver raw data control type */
	gtime_t time;       /* message time */
	gtime_t tobs[MAXSAT][NFREQ + NEXOBS]; /* observation data time */
	obs_t obs;          /* observation data */
	obs_t obuf;         /* observation data buffer */
	nav_t nav;          /* satellite ephemerides */
	sta_t sta;          /* station parameters */
	int ephsat;         /* update satelle of ephemeris (0:no satellite) */
	int ephset;         /* update set of ephemeris (0-1) */
	sbsmsg_t sbsmsg;    /* SBAS message */
	char msgtype[256];  /* last message type */
	uint8_t subfrm[MAXSAT][380]; /* subframe buffer */
	double lockt[MAXSAT][NFREQ + NEXOBS]; /* lock time (s) */
	double icpp[MAXSAT], off[MAXSAT], icpc; /* carrier params for ss2 */
	double prCA[MAXSAT], dpCA[MAXSAT]; /* L1/CA pseudrange/doppler for javad */
	uint8_t halfc[MAXSAT][NFREQ + NEXOBS]; /* half-cycle add flag */
	char freqn[MAXOBS]; /* frequency number for javad */
	int nbyte;          /* number of bytes in message buffer */
	int len;            /* message length (bytes) */
	int iod;            /* issue of data */
	int tod;            /* time of day (ms) */
	int tbase;          /* time base (0:gpst,1:utc(usno),2:glonass,3:utc(su) */
	int flag;           /* general purpose flag */
	int outtype;        /* output message type */
	uint8_t buff[MAXRAWLEN]; /* message buffer */
	char opt[256];      /* receiver dependent options */
	int format;         /* receiver stream format */
	void *rcv_data;     /* receiver dependent data */
} raw_t;


typedef struct {        /* stream converter type */
	int itype, otype;    /* input and output stream type */
	int nmsg;           /* number of output messages */
	int msgs[32];       /* output message types */
	double tint[32];    /* output message intervals (s) */
	uint32_t tick[32];  /* cycle tick of output message */
	int ephsat[32];     /* satellites of output ephemeris */
	int stasel;         /* station info selection (0:remote,1:local) */
	rtcm_t rtcm;        /* rtcm input data buffer */
	raw_t raw;          /* raw  input data buffer */
	rtcm_t out;         /* rtcm output data buffer */
} strconv_t;


typedef struct {        /* GIS data point type */
	double pos[3];      /* point data {lat,lon,height} (rad,m) */
} gis_pnt_t;

typedef struct {        /* GIS data polyline type */
	int npnt;           /* number of points */
	double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
	double *pos;        /* position data (3 x npnt) */
} gis_poly_t;

typedef struct {        /* GIS data polygon type */
	int npnt;           /* number of points */
	double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
	double *pos;        /* position data (3 x npnt) */
} gis_polygon_t;

typedef struct gisd_tag { /* GIS data list type */
	int type;           /* data type (1:point,2:polyline,3:polygon) */
	void *data;         /* data body */
	struct gisd_tag *next; /* pointer to next */
} gisd_t;

typedef struct {        /* GIS type */
	char name[MAXGISLAYER][256]; /* name */
	int flag[MAXGISLAYER];     /* flag */
	gisd_t *data[MAXGISLAYER]; /* gis data list */
	double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
} gis_t;

typedef void fatalfunc_t(const char *); /* fatal callback function type */

/* global variables ----------------------------------------------------------*/
//extern const double chisqr[];        /* chi-sqr(n) table (alpha=0.001) */
//extern const prcopt_t prcopt_default; /* default positioning options */
//extern const solopt_t solopt_default; /* default solution output options */

/* satellites, systems, codes functions --------------------------------------*/
int  satno(int sys, int prn);
int  satsys(int sat, int *prn);
int  satid2no(const char *id);
void satno2id(int sat, char *id);
uint8_t obs2code(const char *obs);
char *code2obs(uint8_t code);
double code2freq(int sys, uint8_t code, int fcn);
double sat2freq(int sat, uint8_t code, const nav_t *nav);
int  code2idx(int sys, uint8_t code);
int  satexclude(int sat, double var, int svh, const prcopt_t *opt);
int  testsnr(int base, int freq, double el, double snr,
	const snrmask_t *mask);
void setcodepri(int sys, int idx, const char *pri);
int  getcodepri(int sys, uint8_t code, const char *opt);

/* matrix and vector functions -----------------------------------------------*/
double *mat(int n, int m);
int    *imat(int n, int m);
double *zeros(int n, int m);
double *eye(int n);
double dot(const double *a, const double *b, int n);
double norm(const double *a, int n);
void cross3(const double *a, const double *b, double *c);
int  normv3(const double *a, double *b);
void matcpy(double *A, const double *B, int n, int m);
void matmul(const char *tr, int n, int k, int m, double alpha,
	const double *A, const double *B, double beta, double *C);
int  matinv(double *A, int n);
int  solve(const char *tr, const double *A, const double *Y, int n,
	int m, double *X);
int  lsq(const double *A, const double *y, int n, int m, double *x,
	double *Q);
int  filter(double *x, double *P, const double *H, const double *v,
	const double *R, int n, int m);
int  smoother(const double *xf, const double *Qf, const double *xb,
	const double *Qb, int n, double *xs, double *Qs);
void matprint(const double *A, int n, int m, int p, int q);
void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);

void add_fatal(fatalfunc_t *func);

/* time and string functions -------------------------------------------------*/
double  str2num(const char *s, int i, int n);
int     str2time(const char *s, int i, int n, gtime_t *t);
void    time2str(gtime_t t, char *str, int n);
gtime_t epoch2time(const double *ep);
void    time2epoch(gtime_t t, double *ep);
gtime_t gpst2time(int week, double sec);
double  time2gpst(gtime_t t, int *week);
gtime_t gst2time(int week, double sec);
double  time2gst(gtime_t t, int *week);
gtime_t bdt2time(int week, double sec);
double  time2bdt(gtime_t t, int *week);
char    *time_str(gtime_t t, int n);

gtime_t timeadd(gtime_t t, double sec);
double  timediff(gtime_t t1, gtime_t t2);
gtime_t gpst2utc(gtime_t t);
gtime_t utc2gpst(gtime_t t);
gtime_t gpst2bdt(gtime_t t);
gtime_t bdt2gpst(gtime_t t);
gtime_t timeget(void);
void    timeset(gtime_t t);
void    timereset(void);
double  time2doy(gtime_t t);
double  utc2gmst(gtime_t t, double ut1_utc);
int read_leaps(const char *file);

int adjgpsweek(int week);
uint32_t tickget(void);
void sleepms(int ms);

int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
	const char *base);
int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts,
	gtime_t te, const char *rov, const char *base);

/* coordinates transformation ------------------------------------------------*/
void ecef2pos(const double *r, double *pos);
void pos2ecef(const double *pos, double *r);
void ecef2enu(const double *pos, const double *r, double *e);
void enu2ecef(const double *pos, const double *e, double *r);
void covenu(const double *pos, const double *P, double *Q);
void covecef(const double *pos, const double *Q, double *P);
void xyz2enu(const double *pos, double *E);
void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
void deg2dms(double deg, double *dms, int ndec);
double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
void readpos(const char *file, const char *rcv, double *pos);
int  sortobs(obs_t *obs);
void uniqnav(nav_t *nav);
int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
int  readnav(const char *file, nav_t *nav);
int  savenav(const char *file, const nav_t *nav);
void freeobs(obs_t *obs);
void freenav(nav_t *nav, int opt);
int  readblq(const char *file, const char *sta, double *odisp);
int  readerp(const char *file, erp_t *erp);
int  geterp(const erp_t *erp, gtime_t time, double *val);

/* debug trace functions -----------------------------------------------------*/
void traceopen(const char *file);
void traceclose(void);
void tracelevel(int level);
void trace(int level, const char *format, ...);
void tracet(int level, const char *format, ...);
void tracemat(int level, const double *A, int n, int m, int p, int q);
void traceobs(int level, const obsd_t *obs, int n);
void tracenav(int level, const nav_t *nav);
void tracegnav(int level, const nav_t *nav);
void tracehnav(int level, const nav_t *nav);
void tracepeph(int level, const nav_t *nav);
void tracepclk(int level, const nav_t *nav);
void traceb(int level, const uint8_t *p, int n);

/* platform dependent functions ----------------------------------------------*/
int execcmd(const char *cmd);
int expath(const char *path, char *paths[], int nmax);
void createdir(const char *path);

/* positioning models --------------------------------------------------------*/
double satazel(const double *pos, const double *e, double *azel);
double geodist(const double *rs, const double *rr, double *e);
void dops(int ns, const double *azel, double elmin, double *dop);

/* atmosphere models ---------------------------------------------------------*/
double ionmodel(gtime_t t, const double *ion, const double *pos,
	const double *azel);
double ionmapf(const double *pos, const double *azel);
double ionppp(const double *pos, const double *azel, double re,
	double hion, double *pppos);
double tropmodel(gtime_t time, const double *pos, const double *azel,
	double humi);
double tropmapf(gtime_t time, const double *pos, const double *azel,
	double *mapfw);
int iontec(gtime_t time, const nav_t *nav, const double *pos,
	const double *azel, int opt, double *delay, double *var);
void readtec(const char *file, nav_t *nav, int opt);
int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
	const double *azel, int ionoopt, double *ion, double *var);
int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
	const double *azel, int tropopt, double *trp, double *var);

/* antenna models ------------------------------------------------------------*/
int  readpcv(const char *file, pcvs_t *pcvs);
pcv_t *searchpcv(int sat, const char *type, gtime_t time,
	const pcvs_t *pcvs);
void antmodel(const pcv_t *pcv, const double *del, const double *azel,
	int opt, double *dant);
void antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
	double *rmoon, double *gmst);
void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
	const double *odisp, double *dr);

/* geiod models --------------------------------------------------------------*/
int opengeoid(int model, const char *file);
void closegeoid(void);
double geoidh(const double *pos);

/* datum transformation ------------------------------------------------------*/
int loaddatump(const char *file);
int tokyo2jgd(double *pos);
int jgd2tokyo(double *pos);

/* rinex functions -----------------------------------------------------------*/
int readrnx(const char *file, int rcv, const char *opt, obs_t *obs,
	nav_t *nav, sta_t *sta);
int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
	double tint, const char *opt, obs_t *obs, nav_t *nav,
	sta_t *sta);
int readrnxc(const char *file, nav_t *nav);
int outrnxobsh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxobsb(FILE *fp, const rnxopt_t *opt, const obsd_t *obs, int n,
	int epflag);
int outrnxnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxgnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxhnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxlnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxqnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxcnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxinavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
int outrnxnavb(FILE *fp, const rnxopt_t *opt, const eph_t *eph);
int outrnxgnavb(FILE *fp, const rnxopt_t *opt, const geph_t *geph);
int outrnxhnavb(FILE *fp, const rnxopt_t *opt, const seph_t *seph);
int rtk_uncompress(const char *file, char *uncfile);
int convrnx(int format, rnxopt_t *opt, const char *file, char **ofile);
int  init_rnxctr(rnxctr_t *rnx);
void free_rnxctr(rnxctr_t *rnx);
int  open_rnxctr(rnxctr_t *rnx, FILE *fp);
int  input_rnxctr(rnxctr_t *rnx, FILE *fp);
/* application defined functions ---------------------------------------------*/
extern int satno(int sys, int prn);

#endif