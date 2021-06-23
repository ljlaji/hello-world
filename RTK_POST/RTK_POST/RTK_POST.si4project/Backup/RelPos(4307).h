#ifndef RELPOS_H
#define RELPOS_H

#include "rtklib.h "

//typedef unsigned char             BOOL;      /* 布尔类型 */
typedef          char             char_t;      /* 字符类型 */
typedef unsigned char             uint8_t;     /* 无符号8位整形数 */
typedef signed   char             sint8_t;     /* 有符号8位整形数 */
typedef unsigned short            uint16_t;    /* 无符号16位整形数 */
typedef signed   short            sint16_t;    /* 有符号16位整形数 */
typedef unsigned int              uint32_t;    /* 无符号32位整形数 */
typedef signed   int              sint32_t;    /* 有符号32位整形数 */
typedef          float            float32_t;      /* 单精度浮点型数 */
typedef          double           float64_t;      /* 双精度浮点型数 */
typedef signed   long long        sint64_t;    /* 有符号64位长整形数 */
typedef unsigned long long        uint64_t;    /* 无符号64位长整形数 */


#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define NX          (4+4)       /* # of estimated parameters for spp*/
#define EPSILON             (1e-13)                    /*最小的正数*/
#define MAXDIM 50
#define MAXITR      10          /* max number of iteration for point pos */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
#define NUM_SYS		2	/* number of  navigation system */
#define NUM_FRE     2



#define COORDIM 3
typedef struct {        /* processing options type */
    int mode;           /* positioning mode (PMODE_???) */
    int soltype;        /* solution type (0:forward,1:backward,2:combined) */
    int nf;             /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
    int navsys;         /* navigation system */
    double elmin;       /* elevation mask angle (rad) */
    double cnrmask;  /* SNR mask */
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
    int niter;          /* number of filter iteration */
    int codesmooth;     /* code smoothing window size (0:none) */


                        /*  2:read from file, 3:rinex header, 4:rtcm pos) */
    double eratio[NFREQ]; /* code/phase error ratio */
    double err[5];      /* measurement error factor */
                        /* [0]:reserved */
                        /* [1-3]:error factor a/b/c of phase (m) */
                        /* [4]:doppler frequency (hz) */
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

} T_PRCOPT;


typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
	uint8_t sat;            /* satellite number */
	sint8_t iode, iodc;      /* IODE,IODC */
	sint8_t sva;            /* SV accuracy (URA index) */
	uint8_t svh;            /* SV health (0:ok) */
	sint32_t week;           /* GPS/QZS: gps week, GAL: galileo week */
	sint16_t code;			/* GPS/QZS: code on L2 */

	float64_t toe, toc, ttr; /* Toe,Toc,T_trans */
	/* SV orbit parameters */
	float64_t A, e, i0, OMG0, omg, M0, deln, OMGd, idot;
	float64_t crc, crs, cuc, cus, cic, cis;
	float64_t toes;        /* Toe (s) in week */
	float64_t fit;         /* fit interval (h) */
	float64_t f0, f1, f2;    /* SV clock parameters (af0,af1,af2) */
	float64_t tgd[6];      /* group delay parameters */
	/* GPS/QZS:tgd[0]=TGD */
	/* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
	/* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
	/*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
} T_EPH;

typedef struct {        /* GLONASS broadcast ephemeris type */
	sint32_t sat;            /* satellite number */
	sint32_t iode;           /* IODE (0-6 bit of tb field) */
	sint32_t frq;            /* satellite frequency number */
	sint32_t svh, sva, age;    /* satellite health, accuracy, age of operation */
	float64_t toe;        /* epoch of epherides (gpst) */
	float64_t tof;        /* message frame time (gpst) */
	float64_t pos[3];      /* satellite position (ecef) (m) */
	float64_t vel[3];      /* satellite velocity (ecef) (m/s) */
	float64_t acc[3];      /* satellite acceleration (ecef) (m/s^2) */
	float64_t taun, gamn;   /* SV clock bias (s)/relative freq bias */
	float64_t dtaun;       /* delay between L1 and L2 (s) */
} T_GEPH;

typedef struct {        /* observation data record */
	float64_t time;       /* receiver sampling time (GPST),second from 1980*/
	uint8_t sat, rcv;    /* satellite/receiver number */
	uint16_t SNR[NFREQ]; /* signal strength (0.001 dBHz) */
	uint8_t  LLI[NFREQ]; /* loss of lock indicator */
	uint8_t code[NFREQ]; /* code indicator (CODE_???) */
	float64_t L[NFREQ]; /* observation data carrier-phase (cycle) */
	float64_t P[NFREQ]; /* observation data pseudorange (m) */
	float32_t  D[NFREQ]; /* observation data doppler frequency (Hz) */
} T_OBS;

typedef struct {        /*epoch observation data record */
	uint8_t  n,nr,nb;     /*# of observation satellite*/
	T_OBS data[MAXOBS]; 
} T_OBSS;


typedef struct {        /* navigation data type */
	sint32_t n, nmax;         /* number of broadcast ephemeris */
	sint32_t ng, ngmax;       /* number of glonass ephemeris */
	T_EPH eph[MAXSAT - NSATGLO];         /* GPS/QZS/GAL/BDS/IRN ephemeris */
	T_GEPH geph[NSATGLO + 1];       /* GLONASS ephemeris,+1 for no use glonass*/
	float64_t utc_gps[8];  /* GPS delta-UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
	float64_t utc_glo[8];  /* GLONASS UTC time parameters {tau_C,tau_GPS} */
	float64_t utc_gal[8];  /* Galileo UTC parameters */
	float64_t utc_qzs[8];  /* QZS UTC parameters */
	float64_t utc_cmp[8];  /* BeiDou UTC parameters */
	float64_t utc_irn[9];  /* IRNSS UTC parameters {A0,A1,Tot,...,dt_LSF,A2} */
	float64_t ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
	float64_t ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	sint32_t glo_fcn[32];    /* GLONASS FCN + 8 */
	pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
} T_NAV;
typedef struct {        /* solution type */
	float64_t time;       /* time (GPST) */
	float64_t rr[6];       /* position/velocity (m|m/s) */
	/* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
	float64_t  qr[6];       /* position variance/covariance (m^2) */
	/* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
	/* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
	float64_t  qv[6];       /* velocity variance/covariance (m^2/s^2) */
	float64_t dtr[6];      /* receiver clock bias to time systems (s) */
	uint8_t type;       /* type (0:xyz-ecef,1:enu-baseline) */
	uint8_t stat;       /* solution status (SOLQ_???) */
	uint8_t ns;         /* number of valid satellites */
	float age;          /* age of differential (s) */
	float ratio;        /* AR ratio factor for valiation */
	float thres;        /* AR ratio threshold for valiation */
} T_SOL;

typedef struct {
	uint8_t n;
	uint8_t sat[MAXSAT];
	float64_t pos[MAXSAT][COORDIM];
	float64_t vel[MAXSAT][COORDIM];
	float64_t ele[MAXSAT];
	float64_t azi[MAXSAT];
	float64_t dis[MAXSAT];		 
	float64_t los[MAXSAT][COORDIM];//line of sight 
} T_SSAT;

typedef struct {        /* double-different type */
	float64_t time;

//reference satellite of different system and frequency, [GPS_L1 GPS_L2 GPS_L5
//														BDS_B1 BDS_B2 BDS_B3]

	uint8_t refsat[NUM_SYS][NUM_FRE];
	uint8_t ddn;
	uint8_t ddnum[NUM_SYS][NUM_FRE];
	uint8_t ddsat[MAXOBS*NUM_FRE];
	float64_t wavelen[MAXOBS*NUM_FRE];
	float64_t ddL[MAXOBS*NUM_FRE];
	float64_t ddP[MAXOBS*NUM_FRE];
	float64_t triddL[MAXOBS*NUM_FRE];
	float64_t triddP[MAXOBS*NUM_FRE];
	float64_t ddmat[MAXOBS*NUM_FRE*COORDIM];
	

} T_DDOBS;


extern uint8_t PntPos(const T_PRCOPT *opt, const T_OBS *obs, const uint8_t n, const T_NAV *nav, T_SOL *sol, T_SSAT  *ssat);
extern void DDObs(const T_PRCOPT opt, const T_OBSS *obss, const T_SOL *solR, const T_SSAT  *ssatR, const T_SOL *solB, const T_SSAT  *ssatB, T_DDOBS *ddobs);
extern void RtkPos(const T_PRCOPT opt,const T_DDOBS ddobs);

#endif

