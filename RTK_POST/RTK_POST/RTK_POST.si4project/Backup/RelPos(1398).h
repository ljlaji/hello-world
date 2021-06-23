#ifndef RELPOS_H
#define RELPOS_H

#include "rtklib.h "

//typedef unsigned char             BOOL;      /* 布尔类型 */
typedef          char             CHAR;      /* 字符类型 */
typedef unsigned char             UINT8;     /* 无符号8位整形数 */
typedef signed   char             SINT8;     /* 有符号8位整形数 */
typedef unsigned short            UINT16;    /* 无符号16位整形数 */
typedef signed   short            SINT16;    /* 有符号16位整形数 */
typedef unsigned int              UINT32;    /* 无符号32位整形数 */
typedef signed   int              SINT32;    /* 有符号32位整形数 */
typedef          float            FP32;      /* 单精度浮点型数 */
typedef          double           FP64;      /* 双精度浮点型数 */
typedef signed   long long        SINT64;    /* 有符号64位长整形数 */
typedef unsigned long long        UINT64;    /* 无符号64位长整形数 */

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define NX          (4+4)       /* # of estimated parameters for spp*/
#define EPSILON             (1e-13)                    /*最小的正数*/
#define MAXDIM MAXSAT
#define MAXITR      10          /* max number of iteration for point pos */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

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
	UINT8 sat;            /* satellite number */
	SINT8 iode, iodc;      /* IODE,IODC */
	SINT8 sva;            /* SV accuracy (URA index) */
	UINT8 svh;            /* SV health (0:ok) */
	SINT32 week;           /* GPS/QZS: gps week, GAL: galileo week */
	SINT16 code;			/* GPS/QZS: code on L2 */

	FP64 toe, toc, ttr; /* Toe,Toc,T_trans */
	/* SV orbit parameters */
	FP64 A, e, i0, OMG0, omg, M0, deln, OMGd, idot;
	FP64 crc, crs, cuc, cus, cic, cis;
	FP64 toes;        /* Toe (s) in week */
	FP64 fit;         /* fit interval (h) */
	FP64 f0, f1, f2;    /* SV clock parameters (af0,af1,af2) */
	FP64 tgd[6];      /* group delay parameters */
	/* GPS/QZS:tgd[0]=TGD */
	/* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
	/* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
	/*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
} T_EPH;

typedef struct {        /* GLONASS broadcast ephemeris type */
	SINT32 sat;            /* satellite number */
	SINT32 iode;           /* IODE (0-6 bit of tb field) */
	SINT32 frq;            /* satellite frequency number */
	SINT32 svh, sva, age;    /* satellite health, accuracy, age of operation */
	FP64 toe;        /* epoch of epherides (gpst) */
	FP64 tof;        /* message frame time (gpst) */
	FP64 pos[3];      /* satellite position (ecef) (m) */
	FP64 vel[3];      /* satellite velocity (ecef) (m/s) */
	FP64 acc[3];      /* satellite acceleration (ecef) (m/s^2) */
	FP64 taun, gamn;   /* SV clock bias (s)/relative freq bias */
	FP64 dtaun;       /* delay between L1 and L2 (s) */
} T_GEPH;

typedef struct {        /* observation data record */
	FP64 time;       /* receiver sampling time (GPST),second from 1980*/
	UINT8 sat, rcv;    /* satellite/receiver number */
	UINT16 SNR[NFREQ]; /* signal strength (0.001 dBHz) */
	UINT8  LLI[NFREQ]; /* loss of lock indicator */
	UINT8 code[NFREQ]; /* code indicator (CODE_???) */
	FP64 L[NFREQ]; /* observation data carrier-phase (cycle) */
	FP64 P[NFREQ]; /* observation data pseudorange (m) */
	FP32  D[NFREQ]; /* observation data doppler frequency (Hz) */
} T_OBS;

typedef struct {        /* navigation data type */
	SINT32 n, nmax;         /* number of broadcast ephemeris */
	SINT32 ng, ngmax;       /* number of glonass ephemeris */
	T_EPH eph[MAXSAT - NSATGLO];         /* GPS/QZS/GAL/BDS/IRN ephemeris */
	T_GEPH geph[NSATGLO + 1];       /* GLONASS ephemeris,+1 for no use glonass*/
	FP64 utc_gps[8];  /* GPS delta-UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
	FP64 utc_glo[8];  /* GLONASS UTC time parameters {tau_C,tau_GPS} */
	FP64 utc_gal[8];  /* Galileo UTC parameters */
	FP64 utc_qzs[8];  /* QZS UTC parameters */
	FP64 utc_cmp[8];  /* BeiDou UTC parameters */
	FP64 utc_irn[9];  /* IRNSS UTC parameters {A0,A1,Tot,...,dt_LSF,A2} */
	FP64 ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	FP64 ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
	FP64 ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	FP64 ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	FP64 ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
	SINT32 glo_fcn[32];    /* GLONASS FCN + 8 */
	pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
} T_NAV;
typedef struct {        /* solution type */
	FP64 time;       /* time (GPST) */
	FP64 rr[6];       /* position/velocity (m|m/s) */
	/* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
	FP64  qr[6];       /* position variance/covariance (m^2) */
	/* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
	/* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
	FP64  qv[6];       /* velocity variance/covariance (m^2/s^2) */
	FP64 dtr[6];      /* receiver clock bias to time systems (s) */
	UINT8 type;       /* type (0:xyz-ecef,1:enu-baseline) */
	UINT8 stat;       /* solution status (SOLQ_???) */
	UINT8 ns;         /* number of valid satellites */
	float age;          /* age of differential (s) */
	float ratio;        /* AR ratio factor for valiation */
	float thres;        /* AR ratio threshold for valiation */
} T_SOL;

typedef struct {
	UINT8 sat[MAXOBS];
	FP64 pos[3 * MAXOBS];
	FP64 vel[3 * MAXOBS];
	FP64 ele[MAXOBS];
	FP64 azi[MAXOBS];
	FP64 e[3 * MAXOBS];
} T_SSAT;


extern UINT8 PntPos(const T_PRCOPT *opt, const T_OBS *obs, const UINT8 n, const T_NAV *nav, T_SOL *sol);


#endif

