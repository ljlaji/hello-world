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


typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    UINT8 sat;            /* satellite number */
    SINT8 iode,iodc;      /* IODE,IODC */
    SINT8 sva;            /* SV accuracy (URA index) */
    UINT8 svh;            /* SV health (0:ok) */
    SINT32 week;           /* GPS/QZS: gps week, GAL: galileo week */

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
    FP64 L[NFREQ+NEXOBS]; /* observation data carrier-phase (cycle) */
    FP64 P[NFREQ+NEXOBS]; /* observation data pseudorange (m) */
    FP32  D[NFREQ+NEXOBS]; /* observation data doppler frequency (Hz) */
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

/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const FP64 *a, const FP64 *b, UINT16 n)
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
extern SINT16 satsys(SINT16 sat, SINT16 *prn)
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
static FP64 var_uraeph(SINT16 sys, SINT16 ura)
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
/* select ephememeris --------------------------------------------------------*/
static T_EPH *seleph(FP64 time, SINT16 sat, SINT16 iode, T_NAV *nav)
{
    FP64 t,tmax,tmin;
    SINT16 i,j=-1,sys,sel;    
    
    sys=satsys(sat,NULL);
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
            sel=getseleph(SYS_GAL);
            if (sel==0&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
            if (sel==1&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
            if (timediff(nav->eph[i].toe,time)>=0.0) continue; /* AOD<=0 */
        }
        if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
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
extern void eph2pos(const FP64 time, const T_EPH *eph,  FP64 *rs,  FP64 *dts,
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
    
    switch ((sys=satsys(eph->sat,&prn))) {
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
    *var=var_uraeph(sys,eph->sva);
}
/* glonass orbit differential equations --------------------------------------*/
static void deq(const FP64 *x, FP64 *xdot, const FP64 *acc)
{
    FP64 a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);
    
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
static void glorbit(FP64 t, FP64 *x, const FP64 *acc)
{
    FP64 k1[6],k2[6],k3[6],k4[6],w[6];
    UINT8 i;
    
    deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    deq(w,k4,acc);
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
extern void geph2pos(FP64 time, const T_GEPH *geph, FP64 *rs, FP64 *dts,
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
        glorbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];
    
    *var=SQR(ERREPH_GLO);
}

					 
/* satellite position and clock by broadcast ephemeris -----------------------*/
static SINT8 ephpos(FP64 time, FP64 teph, int sat, const T_NAV *nav,
                  SINT16 iode, FP64 *rs, FP64 *dts, FP64 *var, SINT8 *svh)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    double rst[3],dtst[1],tt=1E-3;
    int i,sys;
    
    sys=satsys(sat,NULL);
    
    *svh=-1;
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=seleph(teph,sat,iode,nav))) return 0;
        eph2pos(time,eph,rs,dts,var);
        time=timeadd(time,tt);
        eph2pos(time,eph,rst,dtst,var);
        *svh=eph->svh;
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,iode,nav))) return 0;
        geph2pos(time,geph,rs,dts,var);
        time=timeadd(time,tt);
        geph2pos(time,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else if (sys==SYS_SBS) {
        if (!(seph=selseph(teph,sat,nav))) return 0;
        seph2pos(time,seph,rs,dts,var);
        time=timeadd(time,tt);
        seph2pos(time,seph,rst,dtst,var);
        *svh=seph->svh;
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
extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                  const nav_t *nav, double *rs, double *dts, double *var,
                  int *svh)
{
    
    *svh=0;
    
    switch (ephopt) {
        case EPHOPT_BRDC  : return ephpos     (time,teph,sat,nav,-1,rs,dts,var,svh);
      //  case EPHOPT_SBAS  : return satpos_sbas(time,teph,sat,nav,   rs,dts,var,svh);
        case EPHOPT_SSRAPC: return satpos_ssr (time,teph,sat,nav, 0,rs,dts,var,svh);
        case EPHOPT_SSRCOM: return satpos_ssr (time,teph,sat,nav, 1,rs,dts,var,svh);
       /* case EPHOPT_PREC  :
            if (!peph2pos(time,sat,nav,1,rs,dts,var)) break; else return 1;*/
    }
    *svh=-1;
    return 0;
}


/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : gtime_t teph     I   time to select ephemeris (gpst)
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
    FP64 time[2*MAXOBS]={{0}};
     FP64 dt,pr;
    SINT32 i,j;
    
    
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        
        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        
        if (j>=NFREQ) continue;
        /* transmission time by satellite clock */
        time[i]=obs[i].time-pr/CLIGHT;
        
        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i],teph,obs[i].sat,nav,&dt)) {
            continue;
        }
        time[i]=timeadd(time[i],-dt);
        
        /* satellite position and clock at transmission time */
        if (!satpos(time[i],teph,obs[i].sat,ephopt,nav,rs+i*6,dts+i*2,var+i,
                    svh+i)) {
            continue;
        }
        /* if no precise clock available, use broadcast clock instead */
        if (dts[i*2]==0.0) {
            if (!ephclk(time[i],teph,obs[i].sat,nav,dts+i*2)) continue;
            dts[1+i*2]=0.0;
            *var=SQR(STD_BRDCCLK);
        }
    }
}
#if 0
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*           FP64 *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern UINT8 PntPos(const obsd_t *obs, UINT8 n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol,  FP64 *azel)
{
    prcopt_t opt_=*opt;
    FP64 rs[MAXSAT*6]={0},dts[MAXSAT*2]={0},var[MAXSAT]={0},azel_[MAXSAT*2]={0},resp[MAXSAT]={0};
    SINT32 i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    	
    sol->stat=SOLQ_NONE;
    
    if (n<=0) 	return 0;
    sol->time=obs[0].time;    

    /* satellite positons, velocities and clocks */
    SatPoss(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
	//for (i = 0; i < n; i++)
	//{
	//	printf("\n%d,%.2f %.2f  %.2f",obs[i].sat,rs[0 + i * 6], rs[1 + i * 6], rs[2 + i * 6]);
	//}
    /* estimate receiver position with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
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
