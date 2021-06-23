
#ifndef RTK_H
#define RTK_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>



/* constants -----------------------------------------------------------------*/

#define ENACMP


#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */




#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2         frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a/B2a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/L6  frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9          /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ1a_GLO  1.600995E9          /* GLONASS G1a frequency (Hz) */
#define FREQ2a_GLO  1.248060E9          /* GLONASS G2a frequency (Hz) */
#define FREQ1_CMP   1.561098E9          /* BDS B1I     frequency (Hz) */
#define FREQ2_CMP   1.20714E9           /* BDS B2I/B2b frequency (Hz) */
#define FREQ3_CMP   1.26852E9           /* BDS B3      frequency (Hz) */

#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   1.5                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_QZS   1.0                 /* error factor: QZSS */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_IRN   1.5                 /* error factor: IRNSS */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_CMP     0x20                /* navigation system: BeiDou */
#define SYS_IRN     0x40                /* navigation system: IRNS */
#define SYS_LEO     0x80                /* navigation system: LEO */
#define SYS_ALL     0xFF                /* navigation system: all */


#ifndef NFREQ
#define NFREQ       3                   /* number of carrier frequencies */
#endif
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */

#ifndef NEXOBS
#define NEXOBS      0                   /* number of extended obs codes */
#endif

#define SNR_UNIT    0.001               /* SNR unit (dBHz) */

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#ifdef ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif
#ifdef ENAGAL
#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   36                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1
#else
#define MINPRNGAL   0
#define MAXPRNGAL   0
#define NSATGAL     0
#define NSYSGAL     0
#endif
#ifdef ENAQZS
#define MINPRNQZS   193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS   202                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                 /* min satellite PRN number of QZSS L1S */
#define MAXPRNQZS_S 191                 /* max satellite PRN number of QZSS L1S */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif
#ifdef ENACMP
#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   63                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1
#else
#define MINPRNCMP   0
#define MAXPRNCMP   0
#define NSATCMP     0
#define NSYSCMP     0
#endif
#ifdef ENAIRN
#define MINPRNIRN   1                   /* min satellite sat number of IRNSS */
#define MAXPRNIRN   14                  /* max satellite sat number of IRNSS */
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1) /* number of IRNSS satellites */
#define NSYSIRN     1
#else
#define MINPRNIRN   0
#define MAXPRNIRN   0
#define NSATIRN     0
#define NSYSIRN     0
#endif
#ifdef ENALEO
#define MINPRNLEO   1                   /* min satellite sat number of LEO */
#define MAXPRNLEO   10                  /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif
#define NSYS        (NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP+NSYSIRN+NSYSLEO) /* number of systems */

#ifdef ENASBS
#define MINPRNSBS   120                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   158                 /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */
#else
#define MINPRNSBS   0                /* min satellite PRN number of SBAS */
#define MAXPRNSBS   0                 /* max satellite PRN number of SBAS */
#define NSATSBS     0 /* number of SBAS satellites */
#endif


#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
                                        /* max satellite number (1 to MAXSAT) */
#define MAXSTA      255

#ifndef MAXOBS
#define MAXOBS      96                  /* max number of obs in an epoch */
#endif
#define MAXRCV      64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */

#define MAXDTOE     7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_IRN 7200.0              /* max time difference to IRNSS Toe (s) */
#define MAXDTOE_SBS 360.0               /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S   86400.0             /* max time difference to ephem toe (s) for other */
#define MAXGDOP     300.0               /* max GDOP */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless,B1codeless (GPS,BDS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A,B1A    (GAL,BDS) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1S (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2bI  (GAL,BDS) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2bQ  (GAL,BDS) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
#define CODE_L6A    30                  /* obs code: E6A,B3A    (GAL,BDS) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
#define CODE_L6S    35                  /* obs code: L6S        (QZS) */
#define CODE_L6L    36                  /* obs code: L6L        (QZS) */
#define CODE_L8I    37                  /* obs code: E5abI      (GAL) */
#define CODE_L8Q    38                  /* obs code: E5abQ      (GAL) */
#define CODE_L8X    39                  /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
#define CODE_L2I    40                  /* obs code: B1_2I      (BDS) */
#define CODE_L2Q    41                  /* obs code: B1_2Q      (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) (obsolute) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) (obsolute) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define CODE_L1D    56                  /* obs code: B1D        (BDS) */
#define CODE_L5D    57                  /* obs code: L5D(L5S),B2aD (QZS,BDS) */
#define CODE_L5P    58                  /* obs code: L5P(L5S),B2aP (QZS,BDS) */
#define CODE_L5Z    59                  /* obs code: L5D+P(L5S) (QZS) */
#define CODE_L6E    60                  /* obs code: L6E        (QZS) */
#define CODE_L7D    61                  /* obs code: B2bD       (BDS) */
#define CODE_L7P    62                  /* obs code: B2bP       (BDS) */
#define CODE_L7Z    63                  /* obs code: B2bD+P     (BDS) */
#define CODE_L8D    64                  /* obs code: B2abD      (BDS) */
#define CODE_L8P    65                  /* obs code: B2abP      (BDS) */
#define CODE_L4A    66                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4B    67                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4X    68                  /* obs code: G1al1OCd+p (GLO) */
#define MAXCODE     68                  /* max number of obs code */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */
#define SOLF_ENU    2                   /* solution format: e/n/u-baseline */

#define SOLQ_NONE   0                   /* solution status: no solution */
#define SOLQ_GIVEN   1                   /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 2                   /* solution status: single */
#define SOLQ_DGPS   3                   /* solution status: DGPS/DGNSS */
#define SOLQ_FIX    4                   /* solution status: fix */
#define SOLQ_FLOAT  5                   /* solution status: float */


#define IONOOPT_OFF 0                   /* ionosphere option: correction off */
#define IONOOPT_BRDC 1                  /* ionosphere option: broadcast model */
#define IONOOPT_IFLC 3                  /* ionosphere option: L1/L2 iono-free LC */

#define TROPOPT_OFF 0                   /* troposphere option: correction off */
#define TROPOPT_SAAS 1                  /* troposphere option: Saastamoinen model */

#define EPHOPT_BRDC 0                   /* ephemeris option: broadcast ephemeris */

#define ROUND(x)    (int)floor((x)+0.5)

/* constants and macros  on ephemeris------------------------------------------------------*/
#define SQR(x)   ((x)*(x))
#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

/*常用常量定义*/

#define NX          (4+4)       /* 单点定位估计的状态量*/
#define EPSILON     (1e-13) /*最小的正数*/
#define MAXITR      10          /* 单点定位最大迭代数 */
#define ERR_ION     5.0         /* 电离层延时误差 (m) */
#define ERR_TROP    3.0         /* 对流层延时误差 (m) */
#define ERR_SAAS    0.3         /* Saastamoinen 模型误差 (m) */
#define ERR_BRDCI   0.5         /* 广播星历电离层模型误差因子 */
#define ERR_CBIAS   0.3         /* 伪码偏移误差 (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* 测量误差估计的最小高度角 (rad) */
#define NUM_SYS		2	/* RTK中导航系统数量 */
#define NUM_FRE     2	/* RTK中导航频率数量 */
#define AMBSOLNUM           (2)    /*模糊度组合个数*/
#define MAXDIM MAXOBS*NUM_FRE
#define COORDIM 3	/*维度数量*/
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x)) /*开方运算*/
#define MOD_SF 0		/*最小二乘单频组合*/
#define MOD_WL 1		/*最小二乘宽巷组合*/	
#define PROB_FAIL 0.01	/*模糊度固定失败先验概率*/

/*数据类型定义*/
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

typedef struct {        /* 解算配置结构体 */
	int mode;           /* 定位模式 */
	int nf;             /* 频率数量 */
	int navsys;         /* 卫星系统 */
	double elmin;       /* elevation mask angle (rad) */
	double cnrmask;  /* SNR mask */
	int modear;         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
	int ionoopt;        /* ionosphere option (IONOOPT_???) */
	int tropopt;        /* troposphere option (TROPOPT_???) */
	int niter;          /* number of filter iteration */
	double thresslip;   /* slip threshold of geometry-free phase (m) */
	double maxtdiff;    /* max difference of time (sec) */
	double maxinno;     /* reject threshold of innovation (m) */
	double maxgdop;     /* reject threshold of gdop */
	double baseline[2]; /* baseline length constraint {const,sigma} (m) */
	double ru[3];       /* rover position for fixed mode {x,y,z} (ecef) (m) */
	double rb[3];       /* base position for relative mode {x,y,z} (ecef) (m) */
	double err[5];      /* measurement error factor */
	/* [0]:reserved */
	/* [1-3]:error factor a/b/c of phase (m) */
	/* [4]:doppler frequency (hz) */

	float64_t sigL[NUM_SYS][NUM_FRE];
	float64_t sigP[NUM_SYS][NUM_FRE];
	float64_t lambdaffrtmask;

} T_PRCOPT;


typedef struct {        /* GPS/QZS/GAL广播星历结构体 */
	uint8_t sat;            /* 卫星号 */
	sint8_t iode, iodc;      /* IODE,IODC */
	sint8_t sva;            /* 卫星精度（URA索引） */
	uint8_t svh;            /* 健康因子 (0:ok) */
	sint32_t week;           /* GPS/QZS: GPS周, GAL: galileo周 */
	sint16_t code;			/* galileo选星历时有用（SelEph） */

	float64_t toe, toc, ttr; /* Toe,Toc,T_trans */
	float64_t A, e, i0, OMG0, omg, M0, deln, OMGd, idot;/* 卫星轨道参数 */
	float64_t crc, crs, cuc, cus, cic, cis;
	float64_t toes;        /* Toe (s) 周内秒*/
	float64_t f0, f1, f2;    /* 卫星钟差参数 (af0,af1,af2) */
	float64_t tgd[6];      /* 群延时参数 */
	/* GPS/QZS:tgd[0]=TGD */
	/* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
	/* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
	/*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
} T_EPH;

typedef struct {        /* GLONASS 广播星历结构体 */
	sint32_t sat;            /* 卫星数 */
	sint32_t iode;           /* IODE  */
	sint32_t frq;            /* 频率数量 */
	sint32_t svh, sva, age;    /* 卫星健康, 精度, 龄期参数 */
	float64_t toe;        /* epoch of epherides (gpst) */
	float64_t tof;        /* message frame time (gpst) */
	float64_t pos[3];      /* 卫星位置 (ecef) (m) */
	float64_t vel[3];      /* 卫星速度 (ecef) (m/s) */
	float64_t acc[3];      /* 卫星加速度 (ecef) (m/s^2) */
	float64_t taun, gamn;   /* SV clock bias (s)/relative freq bias */
	float64_t dtaun;       /* delay between L1 and L2 (s) */
} T_GEPH;

typedef struct {        /* 单颗观测量结构体 */
	float64_t time;       /* 接收时刻 (GPST),1980以来的累计秒*/
	uint8_t sat, rcv;    /* 卫星/接收机编号 */
	uint16_t SNR[NFREQ]; /* 载噪比 (0.001 dBHz) */
	uint8_t code[NFREQ]; /* 观测量频段标识符，全局变量obscodes索引值 */
	float64_t L[NFREQ]; /* 载波相位 (cycle) */
	float64_t P[NFREQ]; /*伪距 (m) */
	float32_t  D[NFREQ]; /*多普勒 (Hz) */
} T_OBS;

typedef struct {        /*多颗卫星观测量结构体 */
	uint8_t  n, nr, nb;     /*总卫星，移动站卫星、基准站卫星数量*/
	T_OBS data[MAXSAT];  /*卫星观测数据*/
} T_OBSS;


typedef struct {        /* 导航数据 */
	sint32_t n, nmax;         /* 广播星历总数 */
	sint32_t ng, ngmax;       /* glonass的星历数量 */
	T_EPH eph[(MAXSAT - NSATGLO)*20];         /* GPS/QZS/GAL/BDS/IRN 星历 */
	T_GEPH geph[(NSATGLO + 1)*20];       /* GLONASS 星历,无glonass星历时，+1防止出错*/
	sint32_t glo_fcn[32];    /* GLONASS FCN + 8 */
	float64_t ion_gps[8];  /* GPS 电离层模型参数 {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_gal[4];  /* Galileo 电离层模型参数 {ai0,ai1,ai2,0} */
	float64_t ion_qzs[8];  /* QZSS 电离层模型参数 {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_cmp[8];  /* BeiDou 电离层模型参数 {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_irn[8];  /* IRNSS 电离层模型参数 {a0,a1,a2,a3,b0,b1,b2,b3} */
} T_NAV;
typedef struct {        /* solution type */
	float64_t time;       /* time (GPST) */
	uint8_t ns;         /* number of valid satellites */
	float64_t rr[6];       /* position/velocity (m|m/s) */
	/* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
	float64_t  qr[6];       /* position variance/covariance (m^2) */
	/* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
	/* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
	float64_t  qv[6];       /* velocity variance/covariance (m^2/s^2) */
	float64_t dtr[6];      /* receiver clock bias to time systems (s) */
	uint8_t type;       /* type (0:xyz-ecef,1:enu-baseline) */
	uint8_t stat;       /* solution status (SOLQ_???) */
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

typedef struct {        /* 双差观测量结构体 */
	float64_t time;
	uint8_t refsat[NUM_SYS][NUM_FRE];/*每个系统、频点的参考星, [GPS_L1，GPS_L2，GPS_L5；BDS_B1，BDS_B2，BDS_B3]*/
	uint8_t num;		//双差观测量总数量
	uint8_t ddnum[NUM_SYS][NUM_FRE];//每个系统、频点双差数量
	uint8_t ddsat[MAXOBS*NUM_FRE];//双差观测量对应卫星
	float64_t wavelen[MAXOBS*NUM_FRE];//双差观测量对应的波长
	float64_t ddL[MAXOBS*NUM_FRE];//双差载波相位(原始值-估计值)
	float64_t ddP[MAXOBS*NUM_FRE];//双差伪距(原始值-估计值)
	float64_t triddL[MAXOBS*NUM_FRE];//原始双差载波相位
	float64_t triddP[MAXOBS*NUM_FRE];//原始双差伪距
	float64_t ddmat[MAXOBS*NUM_FRE*COORDIM];//单差单位矢量（主星-副星）
	float64_t ddAmb[MAXOBS*NUM_FRE];//双差模糊度

} T_DDOBS;

typedef struct {        /* 双差模糊度结构体 */
	float64_t time;
	uint8_t refsat[NUM_SYS][NUM_FRE];/*每个系统、频点的参考星, [GPS_L1，GPS_L2，GPS_L5；BDS_B1，BDS_B2，BDS_B3]*/
	uint8_t num;		//双差模糊度总数量
	uint8_t ddnum[NUM_SYS][NUM_FRE];//每个系统、频点双差数量
	uint8_t ddsat[MAXOBS*NUM_FRE];//双差模糊度对应卫星
	float64_t wavelen[MAXOBS*NUM_FRE];//双差模糊度对应的波长
	float64_t ddamb[MAXOBS*NUM_FRE];//双差模糊度
	float64_t ratio;//ratio值
	float64_t ps;//固定成功率
	float64_t ratioThreshold;//固定成功门限
	uint8_t stat;//固定状态，SOLQ_FLOAT,SOLQ_FIX	
} T_DDAmb;

typedef struct
{
	float64_t time;
	float64_t baseline[COORDIM];
	float64_t baselinelen;
	T_SOL solrover, solbase;
	T_DDAmb ambfix;

}T_RTK;


/*----------------------------------外部函数声明-------------------------------------------------*/
/*矩阵相关函数*/
extern double Dot(const float64_t *a, const float64_t *b, sint16_t n);
extern float64_t Norm(const float64_t *a, uint16_t n);
extern void MulMat(const float64_t A[], const float64_t B[], uint16_t iM, uint16_t iN, uint16_t iK, float64_t C[]);
extern sint32_t InvMat(const float64_t A[], uint32_t m, float64_t B[]);
extern void Eye(float64_t *A, uint32_t iRow);
extern void Diag(const float64_t A[], uint32_t iRow, float64_t B[]);
extern void TransMat(const float64_t A[], uint32_t iRow, uint32_t iCol, float64_t B[]);
extern void Zeros(float64_t *A, uint16_t n, uint16_t m);
extern void AddMat(const float64_t *A, const float64_t *B, uint16_t iRow, uint16_t iCol, float64_t *C);
extern void SubMat(const float64_t *A, const float64_t *B, uint16_t iRow, uint16_t iCol, float64_t *C);
extern void PrintMat(double *A, int n, int m);
extern void CpyMat(double *A, const float64_t *B, uint32_t n, uint32_t m);
extern void CombineMat(uint8_t flag, const float64_t *A, const uint16_t n, const uint16_t m, const float64_t *B, const uint16_t p, const uint16_t q, float64_t *C);
extern void BlkDiag(const float64_t *A, uint16_t n, uint16_t m, const float64_t *B, uint16_t p, uint16_t q, float64_t *C);
extern void SelectMat(const float64_t *A, uint16_t n, uint16_t m, uint16_t row0, uint16_t col0, uint16_t p, uint16_t q, float64_t *B);
extern uint8_t Lsq(const float64_t *A, const float64_t *y, const float64_t *Q, uint16_t n, uint16_t m, float64_t *x, float64_t *Qx);

extern void Ecef2Pos(const float64_t *r, float64_t *pos);
extern uint8_t PntPos(const T_PRCOPT *opt, const T_OBS *obs, const uint8_t n, const T_NAV *nav, T_SOL *sol, T_SSAT  *ssat);
sint8_t LambdaFFRT(uint16_t n, uint16_t m, const float64_t a[], const float64_t Q[], float64_t Pf, float64_t F[], float64_t s[], float64_t mu[]);
extern float64_t Epoch2time(const double *ep);
extern void Time2epoch(float64_t t, float64_t *ep);
extern uint8_t RtkPos(const T_PRCOPT *opt, const T_OBSS *obst, const T_NAV *navt, T_RTK *rtk);

FILE *ifp[3], *ofp[2];

#endif /* RTKLIB_H */
