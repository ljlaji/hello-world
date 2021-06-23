#include "rtklib.h "

typedef unsigned char             BOOL;      /* 布尔类型 */
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
    double odisp[2][6*11]; /* ocean tide loading parameters {rov,base} */
    int  freqopt;       /* disable L2-AR */
    char pppopt[256];   /* ppp option */
} T_PRCOPT;


typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    UINT8 sat;            /* satellite number */
    SINT8 iode,iodc;      /* IODE,IODC */
    SINT8 sva;            /* SV accuracy (URA index) */
    UINT8 svh;            /* SV health (0:ok) */
    SINT32 week;           /* GPS/QZS: gps week, GAL: galileo week */
	SINT16 code;			/* GPS/QZS: code on L2 */

    FP64 toe,toc,ttr; /* Toe,Toc,T_trans */
                        /* SV orbit parameters */
    FP64 A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    FP64 crc,crs,cuc,cus,cic,cis;
    FP64 toes;        /* Toe (s) in week */
    FP64 fit;         /* fit interval (h) */
    FP64 f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
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
    SINT32 svh,sva,age;    /* satellite health, accuracy, age of operation */
    FP64 toe;        /* epoch of epherides (gpst) */
    FP64 tof;        /* message frame time (gpst) */
    FP64 pos[3];      /* satellite position (ecef) (m) */
    FP64 vel[3];      /* satellite velocity (ecef) (m/s) */
    FP64 acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    FP64 taun,gamn;   /* SV clock bias (s)/relative freq bias */
    FP64 dtaun;       /* delay between L1 and L2 (s) */
} T_GEPH;

typedef struct {        /* observation data record */
    FP64 time;       /* receiver sampling time (GPST),second from 1980*/
    UINT8 sat,rcv;    /* satellite/receiver number */
    UINT16 SNR[NFREQ]; /* signal strength (0.001 dBHz) */
    UINT8  LLI[NFREQ]; /* loss of lock indicator */
    UINT8 code[NFREQ]; /* code indicator (CODE_???) */
    FP64 L[NFREQ]; /* observation data carrier-phase (cycle) */
    FP64 P[NFREQ]; /* observation data pseudorange (m) */
    FP32  D[NFREQ]; /* observation data doppler frequency (Hz) */
} T_OBS;

typedef struct {        /* navigation data type */
    SINT32 n,nmax;         /* number of broadcast ephemeris */
    SINT32 ng,ngmax;       /* number of glonass ephemeris */
    T_EPH eph[MAXSAT-NSATGLO];         /* GPS/QZS/GAL/BDS/IRN ephemeris */
    T_GEPH geph[NSATGLO+1];       /* GLONASS ephemeris,+1 for no use glonass*/
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
	FP64 pos[3*MAXOBS];
	FP64 vel[3*MAXOBS];
	FP64 ele[MAXOBS];
	FP64 azi[MAXOBS];
	FP64 e[3*MAXOBS];
} T_SSAT;
/* ephemeris selections ------------------------------------------------------*/
static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double Dot(const FP64 *a, const FP64 *b, UINT16 n)
{
    FP64 c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}


/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system
* args   : int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern SINT16 SatSys(SINT16 sat, SINT16 *prn)
{
    SINT16 sys=SYS_NONE;
    if (sat<=0||MAXSAT<sat) sat=0;
    else if (sat<=NSATGPS) {
        sys=SYS_GPS; sat+=MINPRNGPS-1;
    }
    else if ((sat-=NSATGPS)<=NSATGLO) {
        sys=SYS_GLO; sat+=MINPRNGLO-1;
    }
    else if ((sat-=NSATGLO)<=NSATGAL) {
        sys=SYS_GAL; sat+=MINPRNGAL-1;
    }
    else if ((sat-=NSATGAL)<=NSATQZS) {
        sys=SYS_QZS; sat+=MINPRNQZS-1; 
    }
    else if ((sat-=NSATQZS)<=NSATCMP) {
        sys=SYS_CMP; sat+=MINPRNCMP-1; 
    }
    else if ((sat-=NSATCMP)<=NSATIRN) {
        sys=SYS_IRN; sat+=MINPRNIRN-1; 
    }
    else if ((sat-=NSATIRN)<=NSATLEO) {
        sys=SYS_LEO; sat+=MINPRNLEO-1; 
    }
    else if ((sat-=NSATLEO)<=NSATSBS) {
        sys=SYS_SBS; sat+=MINPRNSBS-1; 
    }
    else sat=0;
    if (prn) *prn=sat;
    return sys;
}
/* variance by ura ephemeris -------------------------------------------------*/
static FP64 Var_Uraeph(SINT16 sys, SINT16 ura)
{
    const FP64 ura_value[]={   
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys==SYS_GAL) { /* galileo sisa (ref [7] 5.1.11) */
        if (ura<= 49) return SQR(ura*0.01);
        if (ura<= 74) return SQR(0.5+(ura- 50)*0.02);
        if (ura<= 99) return SQR(1.0+(ura- 75)*0.04);
        if (ura<=125) return SQR(2.0+(ura-100)*0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
        return ura<0||14<ura?SQR(6144.0):SQR(ura_value[ura]);
    }
}
/* get selected satellite ephemeris -------------------------------------------
* Get the selected satellite ephemeris.
* args   : int    sys       I   satellite system (SYS_???)
* return : selected ephemeris
*            refer setseleph()
*-----------------------------------------------------------------------------*/
extern SINT16 GetSelEph(SINT16 sys)
{
    switch (sys) {
        case SYS_GPS: return eph_sel[0];
        case SYS_GLO: return eph_sel[1];
        case SYS_GAL: return eph_sel[2];
        case SYS_QZS: return eph_sel[3];
        case SYS_CMP: return eph_sel[4];
        case SYS_IRN: return eph_sel[5];
        case SYS_SBS: return eph_sel[6];
    }
    return 0;
}

/* select ephememeris --------------------------------------------------------*/
static T_EPH *SelEph(FP64 time, SINT16 sat, SINT16 iode, T_NAV *nav)
{
    FP64 t,tmax,tmin;
    SINT16 i,j=-1,sys,sel;    
    
    sys=SatSys(sat,NULL);
    switch (sys) {
        case SYS_GPS: tmax=MAXDTOE+1.0    ; sel=eph_sel[0]; break;
        case SYS_GAL: tmax=MAXDTOE_GAL    ; sel=eph_sel[2]; break;
        case SYS_QZS: tmax=MAXDTOE_QZS+1.0; sel=eph_sel[3]; break;
        case SYS_CMP: tmax=MAXDTOE_CMP+1.0; sel=eph_sel[4]; break;
        case SYS_IRN: tmax=MAXDTOE_IRN+1.0; sel=eph_sel[5]; break;
        default: tmax=MAXDTOE+1.0; break;
    }
    tmin=tmax+1.0;
    
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        if (iode>=0&&nav->eph[i].iode!=iode) continue;
        if (sys==SYS_GAL) {
            sel=GetSelEph(SYS_GAL);
            if (sel==0&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
            if (sel==1&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
            if ((nav->eph[i].toe-time)>=0.0) continue; /* AOD<=0 */
        }
        if ((t=fabs((nav->eph[i].toe-time)))>tmax) continue;
        if (iode>=0) return nav->eph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        return NULL;
    }
    return nav->eph+j;
}

/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : FP64 time     I   time (gpst)
*          T_EPH *eph       I   broadcast ephemeris
*          FP64 *rs       O   satellite position (ecef) {x,y,z} (m)
*          FP64 *dts      O   satellite clock bias (s)
*          FP64 *var      O   satellite position and clock variance (m^2)
* return : none
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void Eph2Pos(const FP64 time, const T_EPH *eph,  FP64 *rs,  FP64 *dts,
                     FP64 *var)
{
     FP64 tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
     FP64 xg,yg,zg,sino,coso;
	 
     SINT16 n,sys,prn;
    
    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=time-eph->toe;
    
    switch ((sys=SatSys(eph->sat,&prn))) {
        case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
        default:      mu=MU_GPS; omge=OMGE;     break;
    }
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;
    
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        return;
    }
    sinE=sin(E); cosE=cos(E);
    
    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=eph->A*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);
    
    /* beidou geo satellite */
    if (sys==SYS_CMP&&(prn<=5||prn>=59)) { /* ref [9] table 4-1 */
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=time-eph->toc;
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;
    
    /* relativity correction */
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);
    
    /* position and clock error variance */
    *var=Var_Uraeph(sys,eph->sva);
}
/* glonass orbit differential equations --------------------------------------*/
static void Deq(const FP64 *x, FP64 *xdot, const FP64 *acc)
{
    FP64 a,b,c,r2=Dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);
    
    if (r2<=0.0) {
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
    xdot[5]=(c-2.0*a)*x[2]+acc[2];
}

/* glonass position and velocity by numerical integration --------------------*/
static void GlOrbit(FP64 t, FP64 *x, const FP64 *acc)
{
    FP64 k1[6],k2[6],k3[6],k4[6],w[6];
    UINT8 i;
    
    Deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    Deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    Deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    Deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}					 
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : FP64 time     I   time (gpst)
*          T_GEPH *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void Geph2Pos(FP64 time, const T_GEPH *geph, FP64 *rs, FP64 *dts,
                     FP64 *var)
{
    FP64 t,tt,x[6];
    UINT8  i;
    
    t=time-geph->toe;
    
    *dts=-geph->taun+geph->gamn*t;
    
    for (i=0;i<3;i++) {
        x[i  ]=geph->pos[i];
        x[i+3]=geph->vel[i];
    }
    for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        GlOrbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];
    
    *var=SQR(ERREPH_GLO);
}

/* select glonass ephememeris ------------------------------------------------*/
static T_GEPH *SelGeph(FP64 time, SINT16 sat, SINT16 iode, const T_NAV *nav)
{
    FP64 t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    SINT16 i,j=-1;

    
    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs((nav->geph[i].toe-time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        return NULL;
    }
    return nav->geph+j;
}					 
/* satellite position and clock by broadcast ephemeris -----------------------*/
static SINT8 EphPos(FP64 time, FP64 teph, int sat, const T_NAV *nav,
                  SINT16 iode, FP64 *rs, FP64 *dts, FP64 *var, SINT8 *svh)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    double rst[3],dtst[1],tt=1E-3;
    int i,sys;
    
    sys=SatSys(sat,NULL);
    
    *svh=-1;
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=SelEph(teph,sat,iode,nav))) return 0;
        Eph2Pos(time,eph,rs,dts,var);
        time=(time+tt);
        Eph2Pos(time,eph,rst,dtst,var);
        *svh=eph->svh;
    }
    else if (sys==SYS_GLO) {
        if (!(geph=SelGeph(teph,sat,iode,nav))) return 0;
        Geph2Pos(time,geph,rs,dts,var);
        time=time+tt;
        Geph2Pos(time,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else return 0;
    
    /* satellite velocity and clock drift by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
    dts[1]=(dtst[0]-dts[0])/tt;
    
    return 1;
}



/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock
* args   : gtime_t time     I   time (gpst)
*          gtime_t teph     I   time to select ephemeris (gpst)
*          int    sat       I   satellite number
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   sat position and velocity (ecef)
*                               {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts      O   sat clock {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variance (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern UINT8 SatPos(FP64 time, FP64 teph, SINT16 sat, UINT8 ephopt,
                  const T_NAV *nav, FP64 *rs, FP64 *dts, FP64 *var,
                  FP64 *svh)
{
    
    *svh=0;
    
    switch (ephopt) {
        case EPHOPT_BRDC  : return EphPos     (time,teph,sat,nav,-1,rs,dts,var,svh);
      //  case EPHOPT_SBAS  : return satpos_sbas(time,teph,sat,nav,   rs,dts,var,svh);
     //   case EPHOPT_SSRAPC: return satpos_ssr (time,teph,sat,nav, 0,rs,dts,var,svh);
     //   case EPHOPT_SSRCOM: return satpos_ssr (time,teph,sat,nav, 1,rs,dts,var,svh);
       /* case EPHOPT_PREC  :
            if (!peph2pos(time,sat,nav,1,rs,dts,var)) break; else return 1;*/
    }
    *svh=-1;
    return 0;
}
/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern FP64 Eph2Clk(FP64 time, T_EPH *eph)
{
    FP64 t,ts;
    UINT8 i;
    
    t=ts=time-eph->toc;
    
    for (i=0;i<2;i++) {
        t=ts-(eph->f0+eph->f1*t+eph->f2*t*t);
    }
    return eph->f0+eph->f1*t+eph->f2*t*t;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern FP64 Geph2Clk(FP64 time, T_GEPH *geph)
{
    FP64 t,ts;
    UINT8 i;    
    
    t=ts=time-geph->toe;
    
    for (i=0;i<2;i++) {
        t=ts-(-geph->taun+geph->gamn*t);
    }
    return -geph->taun+geph->gamn*t;
}

/* satellite clock with broadcast ephemeris ----------------------------------*/
static UINT8 EphClk(FP64 time, FP64 teph, SINT16 sat, T_NAV *nav,
                  FP64 *dts)
{
    T_EPH  *eph;
    T_GEPH *geph;
    SINT16 sys;
    
    sys=SatSys(sat,NULL);
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=SelEph(teph,sat,-1,nav))) return 0;
        *dts=Eph2Clk(time,eph);
    }
    else if (sys==SYS_GLO) {
        if (!(geph=SelGeph(teph,sat,-1,nav))) return 0;
        *dts=Geph2Clk(time,geph);
    }

    else return 0;
    
    return 1;
}

/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : FP64 teph     I   time to select ephemeris (gpst)
*          obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*           FP64 *rs       O   satellite positions and velocities (ecef)
*           FP64 *dts      O   satellite clocks
*           FP64 *var      O   sat position and clock error variances (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          svh[i]        = obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern void SatPoss(const FP64 teph, const T_OBS *obs, const SINT16 n, const T_NAV *nav, FP64 *rs,  FP64 *dts,  FP64 *var, SINT8 *svh)
{
    FP64 time[MAXOBS]={0};
     FP64 dt,pr;
    SINT32 i,j;
    
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        
        if (j>=NFREQ) continue;
        /* transmission time by satellite clock */
        time[i]=obs[i].time-pr/CLIGHT;
        
        /* satellite clock bias by broadcast ephemeris */
        if (!EphClk(time[i],teph,obs[i].sat,nav,&dt)) {
            continue;
        }
        time[i]=time[i]-dt;
        
        /* satellite position and clock at transmission time */
		if (!SatPos(time[i], teph, obs[i].sat, EPHOPT_BRDC, nav, rs + i * 6, dts + i * 2, var + i,
                    svh+i)) {
            continue;
        }
        /* if no precise clock available, use broadcast clock instead */
        if (dts[i*2]==0.0) {
            if (!EphClk(time[i],teph,obs[i].sat,nav,dts+i*2)) continue;
            dts[1+i*2]=0.0;
            *var=SQR(STD_BRDCCLK);
        }
    }
}


/* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void Ecef2Pos(const FP64 *r, FP64 *pos)
{
    FP64 e2=FE_WGS84*(2.0-FE_WGS84),r2=Dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}

/* pseudorange residuals -----------------------------------------------------*/
static UINT8 rescode(UINT8 iter, const T_OBS*obs, UINT8 n, const FP64 *rs,
                   const FP64 *dts, const FP64 *vare, const SINT8 *svh,
                   const T_NAV *nav, const FP64 *x, const T_PRCOPT *opt,
                   FP64 *v, FP64 *H, FP64 *var, FP64 *azel, UINT8 *vsat,
                   FP64 *resp, UINT8 *ns)
{
    FP64 time=0;
    FP64 r=0,freq=0,dion=0,dtrp=0,vmeas=0,vion=0,vtrp=0,rr[3]={0},pos[3],dtr={0},e[3]={0},P=0;
    SINT16 i=0,j=0,nv=0,sat=0,sys=0,mask[NX-3]={0};
    
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            i++;
            continue;
        }
        /* excluded satellite? */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* geometric distance */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        
        if (iter>0) {
            /* test elevation mask */
            if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
            
            /* test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt)) continue;
            
            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                continue;
            }
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            dion*=SQR(FREQ1/freq);
            vion*=SQR(FREQ1/freq);
            
            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }
        }
        /* psendorange with code bias correction */
        if ((P=prange(obs+i,nav,opt,&vmeas))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        /* design matrix */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }
        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
#if 0 /* enable QZS-GPS time offset estimation */
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* variance of pseudorange error */
        var[nv++]=varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
        
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}
				   

/* estimate receiver position ------------------------------------------------*/
static UINT8 EstPos(const T_OBS *obs, const UINT8 n, const FP664 *rs, const FP64 *dts,
                  const FP64 *vare, const SINT8 *svh, const T_NAV *nav,
                  T_SOL *sol, FP64 *azel, UINT8 *vsat,
                  FP64 *resp)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
    
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    
    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals (m) */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns);
#if 0
		PrintMat(v, nv, 1);
		printf("\n");
		PrintMat(H, NX, nv);


#endif
        
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }

        /* weighted by Std */
       for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }

        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }

/*		printf("\n");
		PrintMat(v, 1, nv);
		printf("\n");
		PrintMat(H, NX, nv);
		printf("\n");
		PrintMat(dx, 1,NX);
		printf("\n");
		PrintMat(Q, NX, NX);*/
        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }
        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            
            /* validate solution */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    return 0;
}

#if 0

/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : T_OBS *obs      I   observation data
*          UINT8    n         I   number of observation data
*          T_NAV  *nav      I   navigation data
*          T_SOL  *sol      IO  solution
*          T_SSAT *ssat     IO  satellite data              (NULL: no output)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern UINT8 PntPos(const T_OBS *obs, const UINT8 n, const T_NAV *nav,T_SOL *sol,  T_SSAT ssat)
{
    FP64 rs[MAXOBS*6]={0},dts[MAXOBS*2]={0},var[MAXOBS]={0},resp[MAXOBS]={0};
	UINT8 i,stat;
	SINT8 svh[MAXOBS];
	SINT32 vsat[MAXOBS]={0},
    	
    sol->stat=SOLQ_NONE;
    
    if (n<=0) 	return 0;
    sol->time=obs[0].time;    

    /* satellite positons, velocities and clocks */
    SatPoss(sol->time,obs,n,nav,rs,dts,var,svh);
	//for (i = 0; i < n; i++)
	//{
	//	printf("\n%d,%.2f %.2f  %.2f",obs[i].sat,rs[0 + i * 6], rs[1 + i * 6], rs[2 + i * 6]);
	//}
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    /* RAIM FDE */
    if (!stat&&n>=6) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    }
    /* estimate receiver velocity with Doppler */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    return stat;
}
#endif
