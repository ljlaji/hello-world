/*------------------------------------------------------------------------------
* rtcm3.c : ephemeris functions
*
*          Copyright (C) 2021 by linjun,HuiTian Xpeng, All rights reserved.
* references:rtklib
* history : 2021/05/19 1.0  new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants and macros ------------------------------------------------------*/

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

/* ephemeris selections ------------------------------------------------------*/
static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};

/* variance by ura ephemeris -------------------------------------------------*/
static double var_uraeph(int sys, int ura)
{
    const double ura_value[]={   
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
/* variance by ura ssr (ref [10] table 3.3-1 DF389) --------------------------*/
static double var_urassr(int ura)
{
    double std;
    if (ura<= 0) return SQR(DEFURASSR);
    if (ura>=63) return SQR(5.4665);
    std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
    return SQR(std);
}
/* almanac to satellite position and clock bias --------------------------------
* compute satellite position and clock bias with almanac (gps, galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          alm_t *alm       I   almanac
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
* return : none
* notes  : see ref [1],[7],[8]
*-----------------------------------------------------------------------------*/
extern void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,x,y,sinO,cosO,cosi,mu;
    int n;   
    
    tk=timediff(time,alm->toa);
    
    if (alm->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=0.0;
        return;
    }
    mu=satsys(alm->sat,NULL)==SYS_GAL?MU_GAL:MU_GPS;
    
    M=alm->M0+sqrt(mu/(alm->A*alm->A*alm->A))*tk;
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-alm->e*sin(E)-M)/(1.0-alm->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        return;
    }
    sinE=sin(E); cosE=cos(E);
    u=atan2(sqrt(1.0-alm->e*alm->e)*sinE,cosE-alm->e)+alm->omg;
    r=alm->A*(1.0-alm->e*cosE);
    i=alm->i0;
    O=alm->OMG0+(alm->OMGd-OMGE)*tk-OMGE*alm->toas;
    x=r*cos(u); y=r*sin(u); sinO=sin(O); cosO=cos(O); cosi=cos(i);
    rs[0]=x*cosO-y*cosi*sinO;
    rs[1]=x*sinO+y*cosi*cosO;
    rs[2]=y*sin(i);
    *dts=alm->f0+alm->f1*tk;
}
/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time, const eph_t *eph)
{
    double t,ts;
    int i;
    
    t=ts=timediff(time,eph->toc);
    
    for (i=0;i<2;i++) {
        t=ts-(eph->f0+eph->f1*t+eph->f2*t*t);
    }
    return eph->f0+eph->f1*t+eph->f2*t*t;
}
/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
                    double *var)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;
    
    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=timediff(time,eph->toe);
    
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
    tk=timediff(time,eph->toc);
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;
    
    /* relativity correction */
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);
    
    /* position and clock error variance */
    *var=var_uraeph(sys,eph->sva);
}
/* glonass orbit differential equations --------------------------------------*/
static void deq(const double *x, double *xdot, const double *acc)
{
    double a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);
    
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
static void glorbit(double t, double *x, const double *acc)
{
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;
    
    deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, const geph_t *geph)
{
    double t,ts;
    int i;    
    
    t=ts=timediff(time,geph->toe);
    
    for (i=0;i<2;i++) {
        t=ts-(-geph->taun+geph->gamn*t);
    }
    return -geph->taun+geph->gamn*t;
}
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var)
{
    double t,tt,x[6];
    int i;
    
    t=timediff(time,geph->toe);
    
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
/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return : satellite clock bias (s)
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern double seph2clk(gtime_t time, const seph_t *seph)
{
    double t;
    int i;    
    
    t=timediff(time,seph->t0);
    
    for (i=0;i<2;i++) {
        t-=seph->af0+seph->af1*t;
    }
    return seph->af0+seph->af1*t;
}
/* sbas ephemeris to satellite position and clock bias -------------------------
* compute satellite position and clock bias with sbas ephemeris
* args   : gtime_t time     I   time (gpst)
*          seph_t  *seph    I   sbas ephemeris
*          double  *rs      O   satellite position {x,y,z} (ecef) (m)
*          double  *dts     O   satellite clock bias (s)
*          double  *var     O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var)
{
    double t;
    int i;   
    
    t=timediff(time,seph->t0);
    
    for (i=0;i<3;i++) {
        rs[i]=seph->pos[i]+seph->vel[i]*t+seph->acc[i]*t*t/2.0;
    }
    *dts=seph->af0+seph->af1*t;
    
    *var=var_uraeph(SYS_SBS,seph->sva);
}
/* select ephememeris --------------------------------------------------------*/
static eph_t *seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax,tmin;
    int i,j=-1,sys,sel;    
    
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
/* select glonass ephememeris ------------------------------------------------*/
static geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1;

    
    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        return NULL;
    }
    return nav->geph+j;
}
/* select sbas ephememeris ---------------------------------------------------*/
static seph_t *selseph(gtime_t time, int sat, const nav_t *nav)
{
    double t,tmax=MAXDTOE_SBS,tmin=tmax+1.0;
    int i,j=-1;

    
    for (i=0;i<nav->ns;i++) {
        if (nav->seph[i].sat!=sat) continue;
        if ((t=fabs(timediff(nav->seph[i].t0,time)))>tmax) continue;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (j<0) {
        return NULL;
    }
    return nav->seph+j;
}
/* satellite clock with broadcast ephemeris ----------------------------------*/
static int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  double *dts)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    int sys;
    
    sys=satsys(sat,NULL);
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=seleph(teph,sat,-1,nav))) return 0;
        *dts=eph2clk(time,eph);
    }
    else if (sys==SYS_GLO) {
        if (!(geph=selgeph(teph,sat,-1,nav))) return 0;
        *dts=geph2clk(time,geph);
    }
    else if (sys==SYS_SBS) {
        if (!(seph=selseph(teph,sat,nav))) return 0;
        *dts=seph2clk(time,seph);
    }
    else return 0;
    
    return 1;
}
/* satellite position and clock by broadcast ephemeris -----------------------*/
static int ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  int iode, double *rs, double *dts, double *var, int *svh)
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
#if 0
/* satellite position and clock with sbas correction -------------------------*/
static int satpos_sbas(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                        double *rs, double *dts, double *var, int *svh)
{
    const sbssatp_t *sbs;
    int i;
    
    trace(4,"satpos_sbas: time=%s sat=%2d\n",time_str(time,3),sat);
    
    /* search sbas satellite correciton */
    for (i=0;i<nav->sbssat.nsat;i++) {
        sbs=nav->sbssat.sat+i;
        if (sbs->sat==sat) break;
    }
    if (i>=nav->sbssat.nsat) {
        trace(2,"no sbas correction for orbit: %s sat=%2d\n",time_str(time,0),sat);
        ephpos(time,teph,sat,nav,-1,rs,dts,var,svh);
        *svh=-1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time,teph,sat,nav,sbs->lcorr.iode,rs,dts,var,svh)) return 0;
    
    /* sbas satellite correction (long term and fast) */
    if (sbssatcorr(time,sat,nav,rs,dts,var)) return 1;
    *svh=-1;
    return 0;
}
#endif
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
* notes  : iono-free LC frequencies defined as follows:
*            GPS/QZSS : L1-L2
*            GLONASS  : G1-G2
*            Galileo  : E1-E5b
*            BDS      : B1I-B2I
*            NavIC    : L5-S
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
	double *dant)
{
	const pcv_t *pcv = nav->pcvs + sat - 1;
	double ex[3], ey[3], ez[3], es[3], r[3], rsun[3], gmst, erpv[5] = { 0 }, freq[2];
	double C1, C2, dant1, dant2;
	int i, sys;

	dant[0] = dant[1] = dant[2] = 0.0;

	/* sun position in ecef */
	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, &gmst);

	/* unit vectors of satellite fixed coordinates */
	for (i = 0; i<3; i++) r[i] = -rs[i];
	if (!normv3(r, ez)) return;
	for (i = 0; i<3; i++) r[i] = rsun[i] - rs[i];
	if (!normv3(r, es)) return;
	cross3(ez, es, r);
	if (!normv3(r, ey)) return;
	cross3(ey, ez, ex);

	/* iono-free LC coefficients */
	sys = satsys(sat, NULL);
	if (sys == SYS_GPS || sys == SYS_QZS) { /* L1-L2 */
		freq[0] = FREQ1;
		freq[1] = FREQ2;
	}
	else if (sys == SYS_GLO) { /* G1-G2 */
		freq[0] = sat2freq(sat, CODE_L1C, nav);
		freq[1] = sat2freq(sat, CODE_L2C, nav);
	}
	else if (sys == SYS_GAL) { /* E1-E5b */
		freq[0] = FREQ1;
		freq[1] = FREQ7;
	}
	else if (sys == SYS_CMP) { /* B1I-B2I */
		freq[0] = FREQ1_CMP;
		freq[1] = FREQ2_CMP;
	}
	else if (sys == SYS_IRN) { /* B1I-B2I */
		freq[0] = FREQ5;
		freq[1] = FREQ9;
	}
	else return;

	C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
	C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

	/* iono-free LC */
	for (i = 0; i<3; i++) {
		dant1 = pcv->off[0][0] * ex[i] + pcv->off[0][1] * ey[i] + pcv->off[0][2] * ez[i];
		dant2 = pcv->off[1][0] * ex[i] + pcv->off[1][1] * ey[i] + pcv->off[1][2] * ez[i];
		dant[i] = C1*dant1 + C2*dant2;
	}
}
/* satellite position and clock with ssr correction --------------------------*/
static int satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                      int opt, double *rs, double *dts, double *var, int *svh)
{
    const ssr_t *ssr;
    eph_t *eph;
    double t1,t2,t3,er[3],ea[3],ec[3],rc[3],deph[3],dclk,dant[3]={0},tk;
    int i,sys;
    
    
    ssr=nav->ssr+sat-1;
    
    if (!ssr->t0[0].time) {
        return 0;
    }
    if (!ssr->t0[1].time) {
        return 0;
    }
    /* inconsistency between orbit and clock correction */
    if (ssr->iod[0]!=ssr->iod[1]) {
        *svh=-1;
        return 0;
    }
    t1=timediff(time,ssr->t0[0]);
    t2=timediff(time,ssr->t0[1]);
    t3=timediff(time,ssr->t0[2]);
    
    /* ssr orbit and clock correction (ref [4]) */
    if (fabs(t1)>MAXAGESSR||fabs(t2)>MAXAGESSR) {
        *svh=-1;
        return 0;
    }
    if (ssr->udi[0]>=1.0) t1-=ssr->udi[0]/2.0;
    if (ssr->udi[1]>=1.0) t2-=ssr->udi[1]/2.0;
    
    for (i=0;i<3;i++) deph[i]=ssr->deph[i]+ssr->ddeph[i]*t1;
    dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;
    
    /* ssr highrate clock correction (ref [4]) */
    if (ssr->iod[0]==ssr->iod[2]&&ssr->t0[2].time&&fabs(t3)<MAXAGESSR_HRCLK) {
        dclk+=ssr->hrclk;
    }
    if (norm(deph,3)>MAXECORSSR||fabs(dclk)>MAXCCORSSR) {
        *svh=-1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time,teph,sat,nav,ssr->iode,rs,dts,var,svh)) return 0;
    
    /* satellite clock for gps, galileo and qzss */
    sys=satsys(sat,NULL);
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) {
        if (!(eph=seleph(teph,sat,ssr->iode,nav))) return 0;
        
        /* satellite clock by clock parameters */
        tk=timediff(time,eph->toc);
        dts[0]=eph->f0+eph->f1*tk+eph->f2*tk*tk;
        dts[1]=eph->f1+2.0*eph->f2*tk;
        
        /* relativity correction */
        dts[0]-=2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
    }
    /* radial-along-cross directions in ecef */
    if (!normv3(rs+3,ea)) return 0;
    cross3(rs,rs+3,rc);
    if (!normv3(rc,ec)) {
        *svh=-1;
        return 0;
    }
    cross3(ea,ec,er);
    
    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rs,sat,nav,dant);
    }
    for (i=0;i<3;i++) {
        rs[i]+=-(er[i]*deph[0]+ea[i]*deph[1]+ec[i]*deph[2])+dant[i];
    }
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    dts[0]+=dclk/CLIGHT;
    
    /* variance by ssr ura */
    *var=var_urassr(ssr->ura);
    
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
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *var      O   sat position and clock error variances (m^2)
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
extern void satposs(gtime_t teph, const obsd_t *obs, int n, const nav_t *nav,
                    int ephopt, double *rs, double *dts, double *var, int *svh)
{
    gtime_t time[2*MAXOBS]={{0}};
    double dt,pr;
    int i,j;
    
    
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        for (j=0;j<6;j++) rs [j+i*6]=0.0;
        for (j=0;j<2;j++) dts[j+i*2]=0.0;
        var[i]=0.0; svh[i]=0;
        
        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        
        if (j>=NFREQ) {
            continue;
        }
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);
        
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
/* set selected satellite ephemeris --------------------------------------------
* Set selected satellite ephemeris for multiple ones like LNAV - CNAV, I/NAV -
* F/NAV. Call it before calling satpos(),satposs() to use unselected one.
* args   : int    sys       I   satellite system (SYS_???)
*          int    sel       I   selection of ephemeris
*                                 GPS,QZS : 0:LNAV ,1:CNAV  (default: LNAV)
*                                 GAL     : 0:I/NAV,1:F/NAV (default: I/NAV)
*                                 others  : undefined
* return : none
* notes  : default ephemeris selection for galileo is any.
*-----------------------------------------------------------------------------*/
extern void setseleph(int sys, int sel)
{
    switch (sys) {
        case SYS_GPS: eph_sel[0]=sel; break;
        case SYS_GLO: eph_sel[1]=sel; break;
        case SYS_GAL: eph_sel[2]=sel; break;
        case SYS_QZS: eph_sel[3]=sel; break;
        case SYS_CMP: eph_sel[4]=sel; break;
        case SYS_IRN: eph_sel[5]=sel; break;
        case SYS_SBS: eph_sel[6]=sel; break;
    }
}
/* get selected satellite ephemeris -------------------------------------------
* Get the selected satellite ephemeris.
* args   : int    sys       I   satellite system (SYS_???)
* return : selected ephemeris
*            refer setseleph()
*-----------------------------------------------------------------------------*/
extern int getseleph(int sys)
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
