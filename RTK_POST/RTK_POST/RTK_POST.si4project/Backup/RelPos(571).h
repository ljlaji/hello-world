#ifndef RELPOS_H
#define RELPOS_H

#include "rtklib.h "

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
#define SFREL 0
#define WLREL 1
#define PROB_FAIL 0.01


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

	float64_t sigL[NUM_SYS][NUM_FRE];
	float64_t sigP[NUM_SYS][NUM_FRE];

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
	uint8_t  n,nr,nb;     /*总卫星，移动站卫星、基准站卫星数量*/
	T_OBS data[MAXSAT];  /*卫星观测数据*/
} T_OBSS;


typedef struct {        /* 导航数据 */
	sint32_t n, nmax;         /* 广播星历总数 */
	sint32_t ng, ngmax;       /* glonass的星历数量 */
	T_EPH eph[MAXSAT - NSATGLO];         /* GPS/QZS/GAL/BDS/IRN 星历 */
	T_GEPH geph[NSATGLO + 1];       /* GLONASS 星历,无glonass星历时，+1防止出错*/
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

typedef struct {        /* double-different type */
	float64_t time;

//reference satellite of different system and frequency, [GPS_L1 GPS_L2 GPS_L5
//														BDS_B1 BDS_B2 BDS_B3]

	uint8_t refsat[NUM_SYS][NUM_FRE];
	uint8_t num;
	uint8_t ddnum[NUM_SYS][NUM_FRE];
	uint8_t ddsat[MAXOBS*NUM_FRE];
	float64_t wavelen[MAXOBS*NUM_FRE];
	float64_t ddL[MAXOBS*NUM_FRE];
	float64_t ddP[MAXOBS*NUM_FRE];
	float64_t triddL[MAXOBS*NUM_FRE];
	float64_t triddP[MAXOBS*NUM_FRE];
	float64_t ddmat[MAXOBS*NUM_FRE*COORDIM];
	float64_t ddWL[MAXOBS*NUM_FRE];
	float64_t wavelenWL[MAXOBS*NUM_FRE*COORDIM];
	

} T_DDOBS;

typedef struct
{
	float64_t time;
	float64_t baseline[COORDIM];
	float64_t baselinelen;
	float64_t rr[COORDIM],rb[COORDIM];
	uint8_t stat;
	
	
}T_RTK;


/*----------------------------------外部函数声明-------------------------------------------------*/
/*矩阵相关函数*/
extern double Dot(const float64_t *a, const float64_t *b, sint16_t n);
extern float64_t Norm(const float64_t *a, uint16_t n);
extern void MulMat(const float64_t A[],const float64_t B[],uint16_t iM,uint16_t iN,uint16_t iK,float64_t C[]);
extern sint32_t InvMat(const float64_t A[],uint32_t m,float64_t B[]);
extern void Eye(float64_t *A,uint32_t iRow);
extern void Diag(const float64_t A[],uint32_t iRow,float64_t B[]);
extern void TransMat(const float64_t A[],uint32_t iRow,uint32_t iCol,float64_t B[]);
extern void Zeros(float64_t *A,uint16_t n,uint16_t m);
extern void AddMat(const float64_t *A,const float64_t *B,uint16_t iRow,uint16_t iCol,float64_t *C);
extern void SubMat(const float64_t *A,const float64_t *B,uint16_t iRow,uint16_t iCol,float64_t *C);
extern void PrintMat(double *A, int n, int m);
extern void CpyMat(double *A, const float64_t *B, uint32_t n, uint32_t m);
extern void CombineMat(uint8_t flag,const float64_t *A,const uint16_t n,const uint16_t m,const float64_t *B,const uint16_t p,const uint16_t q,float64_t *C);
extern void BlkDiag(const float64_t *A,uint16_t n,uint16_t m,const float64_t *B,uint16_t p, uint16_t q ,float64_t *C);
extern void SelectMat(const float64_t *A,uint16_t n,uint16_t m, uint16_t row0,uint16_t col0,uint16_t p,uint16_t q,float64_t *B);
extern uint8_t Lsq(const float64_t *A, const float64_t *y, const float64_t *Q,uint16_t n, uint16_t m, float64_t *x,float64_t *Qx);
/*星历相关函数*/







extern uint8_t PntPos(const T_PRCOPT *opt, const T_OBS *obs, const uint8_t n, const T_NAV *nav, T_SOL *sol, T_SSAT  *ssat);
sint8_t LambdaFFRT(uint16_t n, uint16_t m, const float64_t a[], const float64_t Q[], float64_t Pf, float64_t F[], float64_t s[], float64_t mu[]);
extern float64_t Epoch2time(const double *ep);
extern void Time2epoch(float64_t t, float64_t *ep);
extern uint8_t RtkPos(const T_PRCOPT *opt, const T_OBSS *obst, const T_NAV *navt, T_RTK *rtk);



#endif

