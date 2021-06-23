#ifndef RELPOS_H
#define RELPOS_H

#include "rtklib.h "

/*���ó�������*/

#define NX          (4+4)       /* ���㶨λ���Ƶ�״̬��*/
#define EPSILON     (1e-13) /*��С������*/
#define MAXITR      10          /* ���㶨λ�������� */
#define ERR_ION     5.0         /* �������ʱ��� (m) */
#define ERR_TROP    3.0         /* ��������ʱ��� (m) */
#define ERR_SAAS    0.3         /* Saastamoinen ģ����� (m) */
#define ERR_BRDCI   0.5         /* �㲥���������ģ��������� */
#define ERR_CBIAS   0.3         /* α��ƫ����� (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* ���������Ƶ���С�߶Ƚ� (rad) */
#define NUM_SYS		2	/* RTK�е���ϵͳ���� */
#define NUM_FRE     2	/* RTK�е���Ƶ������ */
#define AMBSOLNUM           (2)    /*ģ������ϸ���*/
#define MAXDIM MAXOBS*NUM_FRE
#define COORDIM 3	/*ά������*/
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x)) /*��������*/
#define SFREL 0
#define WLREL 1
#define PROB_FAIL 0.01


/*�������Ͷ���*/
//typedef unsigned char             BOOL;      /* �������� */
typedef          char             char_t;      /* �ַ����� */
typedef unsigned char             uint8_t;     /* �޷���8λ������ */
typedef signed   char             sint8_t;     /* �з���8λ������ */
typedef unsigned short            uint16_t;    /* �޷���16λ������ */
typedef signed   short            sint16_t;    /* �з���16λ������ */
typedef unsigned int              uint32_t;    /* �޷���32λ������ */
typedef signed   int              sint32_t;    /* �з���32λ������ */
typedef          float            float32_t;      /* �����ȸ������� */
typedef          double           float64_t;      /* ˫���ȸ������� */
typedef signed   long long        sint64_t;    /* �з���64λ�������� */
typedef unsigned long long        uint64_t;    /* �޷���64λ�������� */


typedef struct {        /* �������ýṹ�� */
    int mode;           /* ��λģʽ */
    int nf;             /* Ƶ������ */
    int navsys;         /* ����ϵͳ */
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


typedef struct {        /* GPS/QZS/GAL�㲥�����ṹ�� */
	uint8_t sat;            /* ���Ǻ� */
	sint8_t iode, iodc;      /* IODE,IODC */
	sint8_t sva;            /* ���Ǿ��ȣ�URA������ */
	uint8_t svh;            /* �������� (0:ok) */
	sint32_t week;           /* GPS/QZS: GPS��, GAL: galileo�� */
	sint16_t code;			/* galileoѡ����ʱ���ã�SelEph�� */

	float64_t toe, toc, ttr; /* Toe,Toc,T_trans */	
	float64_t A, e, i0, OMG0, omg, M0, deln, OMGd, idot;/* ���ǹ������ */
	float64_t crc, crs, cuc, cus, cic, cis;
	float64_t toes;        /* Toe (s) ������*/
	float64_t f0, f1, f2;    /* �����Ӳ���� (af0,af1,af2) */
	float64_t tgd[6];      /* Ⱥ��ʱ���� */
	/* GPS/QZS:tgd[0]=TGD */
	/* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
	/* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
	/*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
} T_EPH;

typedef struct {        /* GLONASS �㲥�����ṹ�� */
	sint32_t sat;            /* ������ */
	sint32_t iode;           /* IODE  */
	sint32_t frq;            /* Ƶ������ */
	sint32_t svh, sva, age;    /* ���ǽ���, ����, ���ڲ��� */
	float64_t toe;        /* epoch of epherides (gpst) */
	float64_t tof;        /* message frame time (gpst) */
	float64_t pos[3];      /* ����λ�� (ecef) (m) */
	float64_t vel[3];      /* �����ٶ� (ecef) (m/s) */
	float64_t acc[3];      /* ���Ǽ��ٶ� (ecef) (m/s^2) */
	float64_t taun, gamn;   /* SV clock bias (s)/relative freq bias */
	float64_t dtaun;       /* delay between L1 and L2 (s) */
} T_GEPH;

typedef struct {        /* ���Ź۲����ṹ�� */
	float64_t time;       /* ����ʱ�� (GPST),1980�������ۼ���*/
	uint8_t sat, rcv;    /* ����/���ջ���� */
	uint16_t SNR[NFREQ]; /* ����� (0.001 dBHz) */
	uint8_t code[NFREQ]; /* �۲���Ƶ�α�ʶ����ȫ�ֱ���obscodes����ֵ */
	float64_t L[NFREQ]; /* �ز���λ (cycle) */
	float64_t P[NFREQ]; /*α�� (m) */
	float32_t  D[NFREQ]; /*������ (Hz) */
} T_OBS;

typedef struct {        /*������ǹ۲����ṹ�� */
	uint8_t  n,nr,nb;     /*�����ǣ��ƶ�վ���ǡ���׼վ��������*/
	T_OBS data[MAXSAT];  /*���ǹ۲�����*/
} T_OBSS;


typedef struct {        /* �������� */
	sint32_t n, nmax;         /* �㲥�������� */
	sint32_t ng, ngmax;       /* glonass���������� */
	T_EPH eph[MAXSAT - NSATGLO];         /* GPS/QZS/GAL/BDS/IRN ���� */
	T_GEPH geph[NSATGLO + 1];       /* GLONASS ����,��glonass����ʱ��+1��ֹ����*/
	sint32_t glo_fcn[32];    /* GLONASS FCN + 8 */
	float64_t ion_gps[8];  /* GPS �����ģ�Ͳ��� {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_gal[4];  /* Galileo �����ģ�Ͳ��� {ai0,ai1,ai2,0} */
	float64_t ion_qzs[8];  /* QZSS �����ģ�Ͳ��� {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_cmp[8];  /* BeiDou �����ģ�Ͳ��� {a0,a1,a2,a3,b0,b1,b2,b3} */
	float64_t ion_irn[8];  /* IRNSS �����ģ�Ͳ��� {a0,a1,a2,a3,b0,b1,b2,b3} */
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


/*----------------------------------�ⲿ��������-------------------------------------------------*/
/*������غ���*/
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
/*������غ���*/







extern uint8_t PntPos(const T_PRCOPT *opt, const T_OBS *obs, const uint8_t n, const T_NAV *nav, T_SOL *sol, T_SSAT  *ssat);
sint8_t LambdaFFRT(uint16_t n, uint16_t m, const float64_t a[], const float64_t Q[], float64_t Pf, float64_t F[], float64_t s[], float64_t mu[]);
extern float64_t Epoch2time(const double *ep);
extern void Time2epoch(float64_t t, float64_t *ep);
extern uint8_t RtkPos(const T_PRCOPT *opt, const T_OBSS *obst, const T_NAV *navt, T_RTK *rtk);



#endif

