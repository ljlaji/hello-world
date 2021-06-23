/********************************************************************
* 版权所有（C）2021，广州汇天航空航天科技有限公司
*
* 文件名称：RelPos.c
* 内容摘要：RTK算法主函数
* 其他说明：
* 版本  号：
* 作    者：lj
* 完成日期：2021-06-02
********************************************************************/

#include "RTK.h"

/* ephemeris selections ------------------------------------------------------*/
static int eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};
static const double gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
static const double gst0 []={1999,8,22,0,0,0}; /* galileo system time reference */
static const double bdt0 []={2006,1, 1,0,0,0}; /* beidou time reference */
static const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
    10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
    31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
    46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
    61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
    74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
    88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
    101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
    113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
    126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
    138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};
float64_t g_wavelen[2][3] = {{0.190293672798365,0.244210213424568,0.254828048790854},{0.192039486310276,0.248349369584307,0.236332464604421}};//系统各频点波长
float64_t g_coefWL[3] = { 41.031141868512115, 38.903337899695956, 31.184746877054572 };//宽巷误差因子2*waveWL*waveWL/waveL1/waveL1

/*----------------------------矩阵相关函数-------------------------------*/

/* 矩阵内积 a'*b*/
extern double Dot(const float64_t *a, const float64_t *b, sint16_t n)
{
    float64_t c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* 向量模-----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern float64_t Norm(const float64_t *a, uint16_t n)
{
	return sqrt(Dot(a, a, n));
}

/*矩阵相乘 C[m][k]=A[m][n]*B[n][k]*/
extern void MulMat(const float64_t A[],const float64_t B[],uint16_t iM,uint16_t iN,uint16_t iK,float64_t C[])
{
	uint16_t i,j,t;

	for (i=0;i<iM;i++)
	{
		for (j=0;j<iK;j++)
		{
			C[i+j*iM]=0.0;
			for (t=0;t<iN;t++)
			{
				C[i+j*iM]=C[i+j*iM]+A[i+t*iM]*B[t+j*iN];				
			}
		}
	}
}

/*矩阵求逆 A=B(-1)*/
extern sint32_t InvMat(const float64_t A[],uint32_t m,float64_t B[])
{
	float64_t MatrixInvisDim[MAXDIM*MAXDIM],MatrixInvjsDim[MAXDIM*MAXDIM],MatrixInvtmpDim[MAXDIM*MAXDIM];
	sint32_t i, j, k;
	sint32_t l, u, v;
	float64_t d, p;
	/*if (m <= 0) return -1;*/
	for (i=0;i<m*m;i++) {
		MatrixInvtmpDim[i] = A[i];
	}

	for (k = 0; k <= m - 1; k++) {
		d = 0.0;
		for (i = k; i <= m - 1; i++) {
			for (j = k; j <= m - 1; j++) {
				l = i * m + j;
				p = fabs(MatrixInvtmpDim[l]);
				if (p > d) {
					d = p;
					MatrixInvisDim[k] = i;
					MatrixInvjsDim[k] = j;
				}
			}
		}
		if (fabs(d) < EPSILON) {
			return -1;
		}
		if (MatrixInvisDim[k] != k) {
			for (j = 0; j <= m - 1; j++) {
				u = k * m + j;
				v = (MatrixInvisDim[k] *m + j);
				p = MatrixInvtmpDim[u];
				MatrixInvtmpDim[u] = MatrixInvtmpDim[v];
				MatrixInvtmpDim[v] = p;
			}
		}
		if (MatrixInvjsDim[k] != k) {
			for (i = 0; i <= m - 1; i++) {
				u = i * m + k;
				v = (i * m + MatrixInvjsDim[k]);
				p = MatrixInvtmpDim[u];
				MatrixInvtmpDim[u] = MatrixInvtmpDim[v];
				MatrixInvtmpDim[v] = p;
			}
		}
		l = k * m + k;
		if ((fabs(MatrixInvtmpDim[l])<EPSILON)) MatrixInvtmpDim[l]=EPSILON;
		MatrixInvtmpDim[l] = 1.0/MatrixInvtmpDim[l];
		for (j = 0; j <= m - 1; j++) {
			if (j != k) {
				u = k * m + j;
				MatrixInvtmpDim[u] = MatrixInvtmpDim[u] *MatrixInvtmpDim[l];
			}
		}
		for (i = 0; i <= m - 1; i++) {
			if (i != k) {
				for (j = 0; j <= m - 1; j++) {
					if (j != k) {
						u = i * m + j;
						MatrixInvtmpDim[u] = MatrixInvtmpDim[u] - MatrixInvtmpDim[i *m + k] *MatrixInvtmpDim[k *m + j];
					}
				}
			}
		}
		for (i = 0; i <= m - 1; i++) {
			if (i != k) {
				u = i * m + k;
				MatrixInvtmpDim[u] =  - MatrixInvtmpDim[u] *MatrixInvtmpDim[l];
			}
		}
	}
	for (k = m - 1; k >= 0; k--) {
		if (MatrixInvjsDim[k] != k)
			for (j = 0; j <= m - 1; j++) {
				u = k * m + j;
				v = (MatrixInvjsDim[k] *m + j);
				p = MatrixInvtmpDim[u];
				MatrixInvtmpDim[u] = MatrixInvtmpDim[v];
				MatrixInvtmpDim[v] = p;
			}
			if (MatrixInvisDim[k] != k)
				for (i = 0; i <= m - 1; i++) {
					u = i * m + k;
					v = (i * m + MatrixInvisDim[k]);
					p = MatrixInvtmpDim[u];
					MatrixInvtmpDim[u] = MatrixInvtmpDim[v];
					MatrixInvtmpDim[v] = p;
				}
	}
	for (i=0;i<m*m;i++) {
		B[i] = MatrixInvtmpDim[i];
	}

	return 0;
}
/*单位矩阵 eye(iRow)*/
extern void Eye(float64_t *A,uint32_t iRow)
{
	uint32_t i;
	
	memset(A,0,iRow*iRow*sizeof(float64_t));
	for (i=0;i<iRow;i++)
	{
		A[i+i*iRow]=1.0;
	}
}

/*对角矩阵 B = diag(A)*/
extern void Diag(const float64_t A[],uint32_t iRow,float64_t B[])
{
	uint32_t i;
	memset(B,0,iRow*iRow*sizeof(float64_t));
	for (i=0;i<iRow;i++)
	{
		B[i+i*iRow]=A[i];
	}
}


/*矩阵转置 B=A'*/
extern void TransMat(const float64_t A[],uint32_t iRow,uint32_t iCol,float64_t B[])
{
	uint32_t i,j;

	for (i=0;i<iCol;i++)
	{
		for (j=0;j<iRow;j++)
		{
			B[j*iCol+i]=A[i*iRow+j];
		}
	}
}
/*o矩阵*/
extern void Zeros(float64_t *A,uint16_t n,uint16_t m)
{
	memset(A,0,n*m*sizeof(float64_t));
}
/*矩阵相加 C[iRow][iCol]=A+B*/
extern void AddMat(const float64_t *A,const float64_t *B,uint16_t iRow,uint16_t iCol,float64_t *C)
{
	uint16_t i,j;

	for (i=0;i<iRow;i++)
	{
		for (j=0;j<iCol;j++)
		{
			C[j*iRow+i]=A[j*iRow+i]+B[j*iRow+i];
		}
	}
}
/*矩阵相减 C[iRow][iCol]=A-B*/
extern void SubMat(const float64_t *A,const float64_t *B,uint16_t iRow,uint16_t iCol,float64_t *C)
{
	uint16_t i,j;

	for (i=0;i<iRow;i++)
	{
		for (j=0;j<iCol;j++)
		{
			C[j*iRow+i]=A[j*iRow+i]-B[j*iRow+i];
		}
	}
}
/*打印矩阵A[n][m]，调试用*/
extern void PrintMat(double *A, int n, int m)
{
	int i, j;
	//fprintf(ofp[0],"\r\n");
	printf("\r\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			//fprintf(ofp[0], "%f ", A[j*n + i]);
			printf("%f ", A[j*n + i]);
		}
		//fprintf(ofp[0], "\r\n");
		printf("\r\n");
	}

}
/* 矩阵复制 A[n][m] = B[n][m]-----------------------------------------------------------------*/
extern void CpyMat(double *A, const float64_t *B, uint32_t n, uint32_t m)
{
    memcpy(A,B,sizeof(float64_t)*n*m);
}

/* 矩阵合并 -----------------------------------------------------------------
* 参数   :  int flag         O	合并方式,0:C=[A,B], 1:C=[A
														 B]
*		   double *A        O   矩阵 A (n x m)
*		   double *B        O   矩阵 B (p x q)
*          int    n,m,p,q   O   矩阵A、B的行列数量(0:n=p,1:m=q)
*		   double *C        O   目标矩阵 (0:n x (m+q),1:(n+p)xm)

* 返回值 : 0:error, 1:normal
*-----------------------------------------------------------------------------*/

extern void CombineMat(uint8_t flag,const float64_t *A,const uint16_t n,const uint16_t m,const float64_t *B,const uint16_t p,const uint16_t q,float64_t *C)
{
	uint16_t i,j;
	
	if(flag==0)
	{
		if (n!=p) return 0;
		memset(C,0,n*(m+q)*sizeof(float64_t));
		memcpy(C,A,n*m*sizeof(float64_t));
		memcpy(C+n*m,B,p*q*sizeof(float64_t));		
	}
	if(flag==1)
	{
		if(m!=q) return 0;
		memset(C,0,(n+p)*m*sizeof(float64_t));
		for(i=0;i<(n+p);i++)
		{
			for(j=0;j<m;j++)
			{
				if(i<n)
					C[i+j*(n+p)] = A[i+j*n];
				else
					C[i+j*(n+p)] = B[i-n+j*p];
			}
		}
	}
	
	return 1;
}
/* 矩阵块对角，C[n+p][m+q]=[A[n][m]	O[n][q]
							O[p][m]	B[p][q]]
*-----------------------------------------------------------------------------*/

extern void BlkDiag(const float64_t *A,uint16_t n,uint16_t m,const float64_t *B,uint16_t p, uint16_t q ,float64_t *C)
{
	uint16_t i,j;
	memset(C, 0, (n + p)*(m + q)*sizeof(float64_t));
	for(i=0;i<n+p;i++)
	{
		if(i<n)
			for(j=0;j<m;j++)
				C[i+j*(n+p)] = A[i+j*n];
		else
			for(j=m;j<(m+q);j++)
				C[i+j*(n+p)] = B[i-n+(j-m)*p];
		
	}
}
/*
从A[n][m]中获取B[p][q] ,起始行列为row0,col0; 
*/
extern void SelectMat(const float64_t *A,uint16_t n,uint16_t m, uint16_t row0,uint16_t col0,uint16_t p,uint16_t q,float64_t *B)
{
	uint16_t i,j;
	for(j=0;j<q;j++)
	{
		memcpy(B+j*p,A+row0+(j+col0)*n,sizeof(float64_t)*p);
	}
}
/* least square estimation y[n]=A[n][m]*x[m],Q[n][n]-----------------------------------------------------
* least square estimation by solving normal equation (x=(A'*Q^-1*A)^-1*A'*Q^-1*y,Qx=(A'*Q^-1*A)^-1)
* args   : FP64 *A        I   design matrix (n x m)
*          FP64 *y        I   measurements (n x 1)
		   FP64 *Q		  I	  variance matrix for measurements (nxn)
*          int    n,m       I   number of measurements and parameters   (n>=m)
*          double *x        O   estmated parameters (m x 1)
*          double *Qx        O   esimated parameters covariance matrix (m x m)
* return : status (1:ok,0:error)
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern uint8_t Lsq(const float64_t *A, const float64_t *y, const float64_t *Q,uint16_t n, uint16_t m, float64_t *x,
               float64_t *Qx)
{
    int info;
	float64_t Qinv[MAXDIM*MAXDIM] = { 0 }, A_Qinv[MAXDIM*MAXDIM] = { 0 }, Qxinv[MAXDIM*MAXDIM] = { 0 }, Tmp[MAXDIM*MAXDIM], At[MAXDIM*MAXDIM];

	if (n<m) return 0;
	if(InvMat(Q,n,Qinv)<0) return 0;
	//MatMul("TN",m,n,n,1.0,A,Qinv,0.0,A_Qinv); /* A'*Q^-1 */
 //   MatMul("NN",m,m,n,1.0,A_Qinv,A,0.0,Qxinv);/* A'*Q^-1*A */
	TransMat(A, n, m, At);
	MulMat(At, Qinv, m,n,n,A_Qinv); /* A'*Q^-1 */
	MulMat(A_Qinv, A,m,n,m,Qxinv);/* A'*Q^-1*A */

	if(InvMat(Qxinv,m,Qx)<0) return 0;
	MulMat(Qx, A_Qinv,m,m,n,Tmp);/* (A'*Q^-1*A)^-1*A'*Q^-1 */
	MulMat(Tmp, y,m,n,1, x);/* x=(A'*Q^-1*A)^-1*A'*Q^-1*y */
    return 1;
}	


/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system
* args   : int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern sint16_t SatSys(sint16_t sat, sint16_t *prn)
{
    sint16_t sys=SYS_NONE;
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
    else sat=0;
    if (prn) *prn=sat;
    return sys;
}

/*----------------星历相关函数-------------------------------------------*/
/* variance by ura ephemeris */
static float64_t Var_Uraeph(sint16_t sys, sint16_t ura)
{
    const float64_t ura_value[]={   
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
/* get selected satellite ephemeris 
* Get the selected satellite ephemeris.
* args   : int    sys       I   satellite system (SYS_???)
* return : selected ephemeris
*            refer setseleph()
*-----------------------------------------------------------------------------*/
static sint16_t GetSelEph(sint16_t sys)
{
    switch (sys) {
        case SYS_GPS: return eph_sel[0];
        case SYS_GLO: return eph_sel[1];
        case SYS_GAL: return eph_sel[2];
        case SYS_QZS: return eph_sel[3];
        case SYS_CMP: return eph_sel[4];
        case SYS_IRN: return eph_sel[5];
    }
    return 0;
}

/* select ephememeris --------------------------------------------------------*/
static T_EPH *SelEph(float64_t time, sint16_t sat, sint16_t iode, T_NAV *nav)
{
    float64_t t,tmax,tmin;
    sint16_t i,j=-1,sys,sel;    
    
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
static void Eph2Pos(const float64_t time, const T_EPH *eph,  float64_t *rs,  float64_t *dts,
                     float64_t *var)
{
     float64_t tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
     float64_t xg,yg,zg,sino,coso;
	 
     sint16_t n,sys,prn;
    
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
static void Deq(const float64_t *x, float64_t *xdot, const float64_t *acc)
{
    float64_t a,b,c,r2=Dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);
    
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
static void GlOrbit(float64_t t, float64_t *x, const float64_t *acc)
{
    float64_t k1[6],k2[6],k3[6],k4[6],w[6];
    uint8_t i;
    
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
static void Geph2Pos(float64_t time, const T_GEPH *geph, float64_t *rs, float64_t *dts,
                     float64_t *var)
{
    float64_t t,tt,x[6];
    uint8_t  i;
    
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
static T_GEPH *SelGeph(float64_t time, sint16_t sat, sint16_t iode, const T_NAV *nav)
{
    float64_t t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    sint16_t i,j=-1;

    
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
static sint8_t EphPos(float64_t time, float64_t teph, int sat, const T_NAV *nav,
                  sint16_t iode, float64_t *rs, float64_t *dts, float64_t *var, sint8_t *svh)
{
    T_EPH  *eph;
    T_GEPH *geph;
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
* args   : time     I   time (gpst)
*          teph     I   time to select ephemeris (gpst)
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
extern uint8_t SatPos(float64_t time, float64_t teph, sint16_t sat, uint8_t ephopt,
                  const T_NAV *nav, float64_t *rs, float64_t *dts, float64_t *var,
                  float64_t *svh)
{
    
    *svh=0;
    
    switch (ephopt) {
        case EPHOPT_BRDC  : return EphPos     (time,teph,sat,nav,-1,rs,dts,var,svh);
    }
    *svh=-1;
    return 0;
}
/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : satellite clock does not include relativity correction and tdg  
*-----------------------------------------------------------------------------*/
extern float64_t Eph2Clk(float64_t time, T_EPH *eph)
{
    float64_t t,ts;
    uint8_t i;
    
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
extern float64_t Geph2Clk(float64_t time, T_GEPH *geph)
{
    float64_t t,ts;
    uint8_t i;    
    
    t=ts=time-geph->toe;
    
    for (i=0;i<2;i++) {
        t=ts-(-geph->taun+geph->gamn*t);
    }
    return -geph->taun+geph->gamn*t;
}

/* satellite clock with broadcast ephemeris ----------------------------------*/
static uint8_t EphClk(float64_t time, float64_t teph, sint16_t sat, T_NAV *nav,
                  float64_t *dts)
{
    T_EPH  *eph;
    T_GEPH *geph;
    sint16_t sys;
    
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
/* get group delay parameter (m) ---------------------------------------------*/
static float64_t Gettgd(uint16_t sat, const T_NAV *nav, uint8_t type)
{
    uint16_t i,sys=SatSys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
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
extern void SatPoss(const float64_t teph, const T_OBS *obs, const sint16_t n, const T_NAV *nav, float64_t *rs,  float64_t *dts,  float64_t *var, sint8_t *svh)
{
    float64_t time[MAXOBS]={0};
    float64_t dt,pr;
    sint32_t i,j;
    
    
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
extern void Ecef2Pos(const float64_t *r, float64_t *pos)
{
    float64_t e2=FE_WGS84*(2.0-FE_WGS84),r2=Dot(r,r,2),z,zk,v=RE_WGS84,sinp,x;
	static int i = 0;
	i++;
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
	if (i == 3)
		i = 3;
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
	x = r2 > 1E-12 ? atan2(r[1], r[0]) : 0.0;
	pos[1] = x;
    pos[2]=sqrt(r2+z*z)-v;
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern float64_t GeoDist(const float64_t *rs, const float64_t *rr, float64_t *e)
{
    float64_t r;
    uint8_t i;
    if (Norm(rs,3)<RE_WGS84) return -1.0;
    for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
    r=Norm(e,3);
    for (i=0;i<3;i++) e[i]/=r;
    return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}
/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void Xyz2Enu(const float64_t *pos, float64_t *E)
{
    float64_t sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    
    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}

/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
extern void Ecef2Enu(const float64_t *pos, const float64_t *r, float64_t *e)
{
	float64_t E[9] = {0};
    
    Xyz2Enu(pos,E);
   // MatMul("NN",3,1,3,1.0,E,r,0.0,e);
    MulMat(E,r,3,3,1,e);
}


/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
extern float64_t SatAzel(const float64_t *pos, const float64_t *e, float64_t *azel)
{
    float64_t az=0.0,el=PI/2.0,enu[3];
    
    if (pos[2]>-RE_WGS84) {
        Ecef2Enu(pos,e,enu);
        az=Dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
        if (az<0.0) az+=2*PI;
        el=asin(enu[2]);
    }
    if (azel) {azel[0]=az; azel[1]=el;}
    return el;
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern void Time2epoch(float64_t t, float64_t *ep)
{
    const uint8_t mday[]={ /* # of days in a month */
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    uint64_t days,sec,mon,day;
    
    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(t/86400);
    sec=(int)(t-days*86400);
    for (day=days%1461,mon=0;mon<48;mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t-(uint64_t)t;
}

/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : FP64 standard time (second from 1970)
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern float64_t Epoch2time(const double *ep)
{
    const uint16_t doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    float64_t time=0.0,time1,time2;
    uint16_t days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];
    
    if (year<1970||2099<year||mon<1||12<mon) return time;
    
    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec=(int)floor(ep[5]);
    time1=days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
    time2=ep[5]-sec;
	time = time1+time2;
    return time;
}

/* time to gps time ------------------------------------------------------------
* convert gtime_t struct to week and tow in gps time
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gps time (NULL: no output)
* return : time of week in gps time (s)
*-----------------------------------------------------------------------------*/
extern float64_t Time2gpst(float64_t t, uint16_t *week)
{
    float64_t t0=Epoch2time(gpst0);
    float64_t sec=t-t0;
    uint16_t w=(int)(sec/(86400*7));
    
    if (week) *week=w;
    return (float64_t)(sec-(float64_t)w*86400*7);
}
static uint8_t *obscodes[]={       /* observation code strings */
    
    ""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
    "1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
    "2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
    "6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
    "2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
    "5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
    "6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
};

/* obs code to obs code string -------------------------------------------------
* convert obs code to obs code string
* args   : uint8_t code     I   obs code (CODE_???)
* return : obs code string ("1C","1P","1P",...)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
extern uint8_t *Code2obs(uint8_t code)
{
    if (code<=CODE_NONE||MAXCODE<code) return "";
    return obscodes[code];
}

/* GPS obs code to frequency -------------------------------------------------*/
static sint8_t Code2freq_GPS(uint8_t code, float64_t *freq)
{
    uint8_t *obs=Code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '2': *freq=FREQ2; return 1; /* L2 */
        case '5': *freq=FREQ5; return 2; /* L5 */
    }
    return -1;
}
/* GLONASS obs code to frequency ---------------------------------------------*/
static sint8_t Code2freq_GLO(uint8_t code, uint8_t fcn, float64_t *freq)
{
    uint8_t *obs=Code2obs(code);
    
    if (fcn<-7||fcn>6) return -1;
    
    switch (obs[0]) {
        case '1': *freq=FREQ1_GLO+DFRQ1_GLO*fcn; return 0; /* G1 */
        case '2': *freq=FREQ2_GLO+DFRQ2_GLO*fcn; return 1; /* G2 */
        case '3': *freq=FREQ3_GLO;               return 2; /* G3 */
        case '4': *freq=FREQ1a_GLO;              return 0; /* G1a */
        case '6': *freq=FREQ2a_GLO;              return 1; /* G2a */
    }
    return -1;
}
/* Galileo obs code to frequency ---------------------------------------------*/
static sint8_t Code2freq_GAL(uint8_t code, float64_t *freq)
{
    uint8_t *obs=Code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* E1 */
        case '7': *freq=FREQ7; return 1; /* E5b */
        case '5': *freq=FREQ5; return 2; /* E5a */
        case '6': *freq=FREQ6; return 3; /* E6 */
        case '8': *freq=FREQ8; return 4; /* E5ab */
    }
    return -1;
}
/* QZSS obs code to frequency ------------------------------------------------*/
static sint8_t Code2freq_QZS(uint8_t code, float64_t *freq)
{
    uint8_t *obs=Code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '2': *freq=FREQ2; return 1; /* L2 */
        case '5': *freq=FREQ5; return 2; /* L5 */
        case '6': *freq=FREQ6; return 3; /* L6 */
    }
    return -1;
}
/* BDS obs code to frequency -------------------------------------------------*/
static sint8_t Code2freq_BDS(uint8_t code, float64_t *freq)
{
    uint8_t *obs=Code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1;     return 0; /* B1C */
        case '2': *freq=FREQ1_CMP; return 0; /* B1I */
        case '7': *freq=FREQ2_CMP; return 1; /* B2I/B2b */
        case '5': *freq=FREQ5;     return 2; /* B2a */
        case '6': *freq=FREQ3_CMP; return 3; /* B3 */
        case '8': *freq=FREQ8;     return 4; /* B2ab */
    }
    return -1;
}
/* NavIC obs code to frequency -----------------------------------------------*/
static sint8_t Code2freq_IRN(uint8_t code, float64_t *freq)
{
    uint8_t *obs=Code2obs(code);
    
    switch (obs[0]) {
        case '5': *freq=FREQ5; return 0; /* L5 */
        case '9': *freq=FREQ9; return 1; /* S */
    }
    return -1;
}

/* system and obs code to frequency --------------------------------------------
* convert system and obs code to carrier frequency
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
*          int    fcn       I   frequency channel number for GLONASS
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
extern float64_t Code2freq(uint16_t sys, uint8_t code, uint8_t fcn)
{
    float64_t freq=0.0;
    
    switch (sys) {
        case SYS_GPS: (void)Code2freq_GPS(code,&freq); break;
        case SYS_GLO: (void)Code2freq_GLO(code,fcn,&freq); break;
        case SYS_GAL: (void)Code2freq_GAL(code,&freq); break;
        case SYS_QZS: (void)Code2freq_QZS(code,&freq); break;
        case SYS_CMP: (void)Code2freq_BDS(code,&freq); break;
        case SYS_IRN: (void)Code2freq_IRN(code,&freq); break;
    }
    return freq;
}
/* system and obs code to frequency index --------------------------------------
* convert system and obs code to frequency index
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
* return : frequency index (-1: error)
*                       0     1     2     3     4 
*           --------------------------------------
*            GPS       L1    L2    L5     -     - 
*            GLONASS   G1    G2    G3     -     -  (G1=G1,G1a,G2=G2,G2a)
*            Galileo   E1    E5b   E5a   E6   E5ab
*            QZSS      L1    L2    L5    L6     - 
*            SBAS      L1     -    L5     -     -
*            BDS       B1    B2    B2a   B3   B2ab (B1=B1I,B1C,B2=B2I,B2b)
*            NavIC     L5     S     -     -     - 
*-----------------------------------------------------------------------------*/
extern sint8_t Code2idx(uint16_t sys, uint8_t code)
{
    double freq;
    
    switch (sys) {
        case SYS_GPS: return Code2freq_GPS(code,&freq);
        case SYS_GLO: return Code2freq_GLO(code,0,&freq);
        case SYS_GAL: return Code2freq_GAL(code,&freq);
        case SYS_QZS: return Code2freq_QZS(code,&freq);
        case SYS_CMP: return Code2freq_BDS(code,&freq);
        case SYS_IRN: return Code2freq_IRN(code,&freq);
    }
    return -1;
}

/* satellite and obs code to frequency -----------------------------------------
* convert satellite and obs code to carrier frequency
* args   : int    sat       I   satellite number
*          uint8_t code     I   obs code (CODE_???)
*          nav_t  *nav_t    I   navigation data for GLONASS (NULL: not used)
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
extern float64_t Sat2freq(uint16_t sat, uint8_t code, const T_NAV *nav)
{
    uint16_t i,fcn=0,sys,prn;
    
    sys=SatSys(sat,&prn);
    
    if (sys==SYS_GLO) {
        if (!nav) return 0.0;
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        if (i<nav->ng) {
            fcn=nav->geph[i].frq;
        }
        else if (nav->glo_fcn[prn-1]>0) {
            fcn=nav->glo_fcn[prn-1]-8;
        }
        else return 0.0;
    }
    return Code2freq(sys,code,fcn);
}

/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : FP64 t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *dion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)

* return : none
*-----------------------------------------------------------------------------*/
extern void IonModel(float64_t t, const float64_t *ion, const float64_t *pos,
                       const float64_t *azel,float64_t *dion, float64_t *var)
{
    const float64_t ion_default[]={ /* 2004/1/1 */
        0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    float64_t tt,f,psi,phi,lam,amp,per,x;
    uint16_t week;
    
    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (Norm(ion,8)<=0.0) ion=ion_default;
    
    /* earth centered angle (semi-circle) */
    psi=0.0137/(azel[1]/PI+0.11)-0.022;
    
    /* subionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/PI+psi*cos(azel[0]);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);
    
    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);
    
    /* local time (s) */
    tt=43200.0*lam+Time2gpst(t,&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */
    
    /* slant factor */
    f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);
    
    /* ionospheric delay */
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;

	*dion = CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
	*var=SQR(*dion*ERR_BRDCI);
	
}
extern uint8_t IonCorr(float64_t t, uint16_t sys,T_NAV *nav, const float64_t *pos,const float64_t *azel,float64_t *dion, float64_t *var)
{
	if (sys ==SYS_GPS)
	{
		IonModel(t, nav->ion_gps, pos,azel,dion, var);
		return 1;
	}		
	else if (sys==SYS_CMP)
	{
		IonModel(t, nav->ion_cmp, pos,azel,dion, var);
		return 1;
	}	
	else
	{
		*dion =0;
		*var = SQR(ERR_ION);
		return 0;
		
	}
			
}
/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : FP64 time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double humi      I   relative humidity
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern float64_t TropModel(float64_t time, const float64_t *pos, const float64_t *azel,
                        float64_t humi)
{
    const float64_t temp0=15.0; /* temparature at sea level */
    float64_t hgt,pres,temp,e,z,trph,trpw;
    
    if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;
    
    /* standard atmosphere */
    hgt=pos[2]<0.0?0.0:pos[2];
    
    pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
    temp=temp0-6.5E-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));
    
    /* saastamoninen model */
    z=PI/2.0-azel[1];
    trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
    trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
    return trph+trpw;
}

/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : FP64 time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern uint8_t TropCorr(float64_t time, const float64_t *pos,const float64_t *azel, float64_t *trp, float64_t *var)
{    
      *trp=TropModel(time,pos,azel,REL_HUMI);
      *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
       return 1;

}
					
/* psendorange with code bias correction -------------------------------------*/
static float64_t Prange(const T_OBS *obs, const T_NAV *nav, const T_PRCOPT *opt,float64_t *var)
{
    float64_t P1,P2,gamma,b1,b2;
    uint16_t sat,sys;
    
    sat=obs->sat;
    sys=SatSys(sat,NULL);
    P1=obs->P[0];
    P2=obs->P[1];
    *var=0.0;
    
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;    

    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        *var=SQR(ERR_CBIAS);
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2,G1-G2 */
            gamma=SQR(FREQ1/FREQ2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b */
            gamma=SQR(FREQ1/FREQ7);
            if (GetSelEph(SYS_GAL)) { /* F/NAV */
                P2-=Gettgd(sat,nav,0)-Gettgd(sat,nav,1); /* BGD_E5aE5b */
            }
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* B1-B2 */
            gamma=SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQ1)/FREQ2_CMP);
            if      (obs->code[0]==CODE_L2I) b1=Gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=Gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=Gettgd(sat,nav,2)+Gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            b2=Gettgd(sat,nav,1); /* TGD_B2I/B2bI (m) */
            return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
        }
        else if (sys==SYS_IRN) { /* L5-S */
            gamma=SQR(FREQ5/FREQ9);
            return (P2-gamma*P1)/(1.0-gamma);
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=Gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys==SYS_GLO) { /* G1 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=Gettgd(sat,nav,0); /* -dtaun (m) */
            return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (GetSelEph(SYS_GAL)) b1=Gettgd(sat,nav,0); /* BGD_E1E5a */
            else                    b1=Gettgd(sat,nav,1); /* BGD_E1E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=Gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=Gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=Gettgd(sat,nav,2)+Gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            gamma=SQR(FREQ9/FREQ5);
            b1=Gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
    return P1;
}
/* pseudorange measurement error variance ------------------------------------*/
static float64_t Varerr(const T_PRCOPT *opt, float64_t el, uint16_t sys)
{
    float64_t fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:EFACT_GPS;
    if (el<MIN_EL) el=MIN_EL;
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}					 

/* pseudorange residuals ----------------------
H[NX][nv]
-------------------------------*/
static uint8_t ResCode(uint8_t iter, const T_OBS *obs, uint8_t n, const float64_t *rs,
                   const float64_t *dts, const float64_t *vare, const T_NAV *nav, const float64_t *x, const T_PRCOPT *opt,
                   float64_t *v, float64_t *H, float64_t *var, float64_t *azel, uint8_t *vsat,
                   float64_t *resp, uint8_t *ns)
{
    float64_t time=0;
	float64_t r = 0, freq = 0, dion = 0, dtrp = 0, vmeas = 0, vion = 0, vtrp = 0, rr[3] = { 0 }, pos[3] = {0}, dtr = { 0 }, e[3] = { 0 }, P = 0;
    sint16_t i=0,j=0,nv=0,sat=0,sys=0,mask[NX-3]={0};
    
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    Ecef2Pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=SatSys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            i++;
            continue;
        }
        
        /* geometric distance */
        if ((r=GeoDist(rs+i*6,rr,e))<=0.0) continue;
        
        if (iter>0) {
            /* test elevation mask */
            if (SatAzel(pos,e,azel+i*2)<opt->elmin) continue;            
            /* test SNR mask */
            if(obs[i].SNR[0]*SNR_UNIT< opt->cnrmask) continue;
			if(opt->ionoopt==IONOOPT_IFLC && obs[i].SNR[1]*SNR_UNIT<opt->cnrmask) continue;			
            /* ionospheric correction */
			if(opt->ionoopt!=IONOOPT_IFLC) 
				IonCorr(time, sys,nav, pos,azel+i*2,&dion, &vion);
            if ((freq=Sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
			if(sys==SYS_CMP)
			{
            	dion*=SQR(FREQ1_CMP/freq);
            	vion*=SQR(FREQ1_CMP/freq);			
			}
			else
			{
            	dion*=SQR(FREQ1/freq);
            	vion*=SQR(FREQ1/freq);

			}

			 /* tropospheric correction */
			TropCorr(time, pos,azel+i*2, &dtrp,&vtrp);
        }	
        /* psendorange with code bias correction */
        if ((P=Prange(obs+i,nav,opt,&vmeas))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        /* design matrix H[NX][nv]*/
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }
        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}

        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* variance of pseudorange error */
        var[nv++]=Varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
        
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

extern void Dops(uint8_t ns, const float64_t *azel, float64_t elmin, float64_t *dop)
{
	float64_t H[4 * MAXSAT] = { 0 },Ht[4 * MAXSAT] = { 0 },Q[16] = { 0 }, Qinv[16] = { 0 }, cosel = { 0 }, sinel = {0};
    uint8_t i,n;
    
    for (i=0;i<4;i++) dop[i]=0.0;
    for (i=n=0;i<ns&&i<MAXSAT;i++) {
        if (azel[1+i*2]<elmin||azel[1+i*2]<=0.0) continue;
        cosel=cos(azel[1+i*2]);
        sinel=sin(azel[1+i*2]);
        H[  4*n]=cosel*sin(azel[i*2]);
        H[1+4*n]=cosel*cos(azel[i*2]);
        H[2+4*n]=sinel;
        H[3+4*n++]=1.0;
    }
    if (n<4) return;
    
 //   MatMul("NT",4,4,n,1.0,H,H,0.0,Qinv);
 	TransMat(H,4,n,Ht);
	MulMat(H,Ht,4,n,4,Qinv);
    if (!InvMat(Qinv,4,Q)) {
        dop[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]); /* GDOP */
        dop[1]=SQRT(Q[0]+Q[5]+Q[10]);       /* PDOP */
        dop[2]=SQRT(Q[0]+Q[5]);             /* HDOP */
        dop[3]=SQRT(Q[10]);                 /* VDOP */
    }
}

/* validate solution ---------------------------------------------------------*/
static uint8_t ValSol(const float64_t *azel, const uint8_t *vsat, uint8_t n,
                  const T_PRCOPT *opt, const float64_t *v, uint8_t nv, uint8_t nx)
{
    float64_t azels[MAXOBS*2]={0},dop[4]={0},vv={0};
    uint8_t i=0,ns=0;
        
    /* Chi-square validation of residuals */
    vv=Dot(v,v,nv);
   // if (nv>nx&&vv>chisqr[nv-nx-1])         return 0;
   if (nv>nx&&vv>1000)         return 0;
    /* large GDOP check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    Dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        return 0;
    }
    return 1;
}
 
/* 计算接收机位置----------------------------------------------------
* 参数   :   
*          T_OBS *obs      		I   观测数据
*          UINT8    n         	I   观测数据数量
*          float64_t  *rs      	I   卫星位置、速度，rs[6][n]
*          float64_t *dts     	I  	卫星钟差，rs[2][n]
*          float64_t *vare      I   卫星位置方差，vare[n]
*          sint8_t *svh         I   卫星健康标志，0正常
*          T_NAV  *nav      	I   星历数据
*          T_PRCOPT  *opt      	I  	定位配置
*          T_SOL  *sol      	IO  解算结果
*          float64_t *azel     	IO  卫星方位角,高度角,azel[2][n]
*          uint8_t  *vsat      	IO  卫星应用标志，vsat[n],1使用，0，不使用
*          float64_t *resp     	IO  卫星伪距余量(观测量-估计量-误差),resp[n]
* 返回值 : status(1:正常,0:异常)
*-----------------------------------------------------------------------------*/	

static uint8_t EstPos(const T_OBS *obs, const uint8_t n, const float64_t *rs, const float64_t *dts,
                  const float64_t *vare, const sint8_t *svh, const T_NAV *nav,const T_PRCOPT *opt,
                  T_SOL *sol, float64_t *azel, uint8_t *vsat,float64_t *resp)
{
    float64_t x[NX]={0},dx[NX]={0},Q[(MAXOBS+4)*(MAXOBS+4)]={0},Qx[NX*NX]={0};
	float64_t v[MAXOBS+4]={0},Ht[NX*(MAXOBS+4)]={0},H[NX*(MAXOBS+4)]={0},var[MAXOBS+4]={0},sig=0;
    uint16_t i,j,k,info,stat,nv,ns;
    
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];    
    for (i=0;i<MAXITR;i++) { 

		if (i == 3)
			i = 3;
        /* pseudorange residuals (m) */
		nv = ResCode(i, obs, n, rs,dts, vare,nav,x, opt,v, Ht, var, azel, vsat,resp, &ns);
        if (nv<NX) break;

        /* weighted by Std */
       	for (j=0;j<nv;j++) 
		{
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) Ht[k+j*NX]/=sig;
        }
	  	TransMat(Ht,NX,nv,H);
	  	Eye(Q,nv);	 
		if(!Lsq(H,v,Q,nv,NX,dx,Qx)) break;

        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }
        if (Norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=obs[0].time-x[3]/CLIGHT;
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Qx[j+j*NX];
            sol->qr[3]=(float)Qx[1];    /* cov xy */
            sol->qr[4]=(float)Qx[2+NX]; /* cov yz */
            sol->qr[5]=(float)Qx[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            
            /* validate solution */
            if ((stat=ValSol(azel,vsat,n,opt,v,nv,NX))) 
			{
				sol->stat=SOLQ_SINGLE;
            }
            return stat;
        }
    }
    return 0;
}
/* RAIM FDE (failure detection and exclution) -------------------------------*/
static uint8_t RaimFde(const T_OBS *obs, uint8_t n, const float64_t *rs,
                    const float64_t *dts, const float64_t *vare, const sint8_t *svh,
                    const T_NAV *nav, const T_PRCOPT *opt, T_SOL *sol,
                    float64_t *azel, uint8_t *vsat, float64_t *resp)
{
    T_OBS obs_e[MAXOBS]={0};
    T_SOL sol_e={{0}};
   // char tstr[32],name[16],msg_e[128];
    float64_t rs_e[6*MAXOBS],dts_e[2*MAXOBS],vare_e[MAXOBS],azel_e[2*MAXOBS],resp_e[MAXOBS];
   	float64_t rms_e,rms=100.0;
    sint16_t i,j,k,nvsat,stat=0,svh_e[MAXOBS],vsat_e[MAXOBS],sat=0;
       
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            CpyMat(rs_e +6*k,rs +6*j,6,1);
            CpyMat(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!EstPos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                    vsat_e,resp_e)) {
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);        
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            CpyMat(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
    }
    return stat;
}
/* range rate residuals ------------------------------------------------------*/
static uint8_t Resdop(const T_OBS *obs, uint8_t n, const float64_t *rs, const float64_t *dts,
                  const T_NAV *nav, const float64_t *rr, const float64_t *x,
                  const float64_t *azel, const uint8_t *vsat, float64_t err, float64_t *v,
                  float64_t *H)
{
    float64_t freq,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,sig;
    uint8_t i,j,nv=0;
    
    
    Ecef2Pos(rr,pos); Xyz2Enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        freq=Sat2freq(obs[i].sat,obs[i].code[0],nav);
        
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||Norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        MatMul("TN",3,1,3,1.0,E,a,0.0,e);
        
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=Dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* Std of range rate error (m/s) */
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        
        /* range rate residual (m/s) */
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        
        /* design matrix */
        for (j=0;j<4;j++) {
            H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
        }
        nv++;
    }
    return nv;
}

/* estimate receiver velocity ------------------------------------------------*/
static void EstVel(const T_OBS *obs, uint8_t n, const float64_t *rs, const float64_t *dts,
                   const T_NAV *nav, const T_PRCOPT *opt, T_SOL *sol,
                   const float64_t *azel, const uint8_t *vsat)
{
    float64_t x[4]={0},dx[4]={0},Q[MAXOBS*MAXOBS]={0},Qx[16]={0},v[MAXOBS]={0},H[4*MAXOBS]={0},Ht[4*MAXOBS]={0};
    float64_t err=opt->err[4]; /* Doppler error (Hz) */
    uint8_t i,j,nv;   
    
    for (i=0;i<MAXITR;i++) {
        
        /* range rate residuals (m/s) */
        if ((nv=Resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,Ht))<4) {
            break;
        }
        /* least square estimation */
		TransMat(Ht,4,nv,H);
	  	EYE(Q,nv);	   
		if(!Lsq(H,v,Q,nv,4,dx,Qx)) break;
		       
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            CpyMat(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Qx[0];  /* xx */
            sol->qv[1]=(float)Qx[5];  /* yy */
            sol->qv[2]=(float)Qx[10]; /* zz */
            sol->qv[3]=(float)Qx[1];  /* xy */
            sol->qv[4]=(float)Qx[6];  /* yz */
            sol->qv[5]=(float)Qx[2];  /* zx */
            break;
        }
    }
}

/* 单点定位 ----------------------------------------------------
* 计算接收机位置、速度、钟差
* 参数   :   T_PRCOPT *opt      I   定位选项配置
*          T_OBS *obs      I   观测数据
*          UINT8    n         I   观测数据数量
*          T_NAV  *nav      I   星历数据
*          T_SOL  *sol      IO  解算结果
*          T_SSAT *ssat     IO  解算卫星信息      (NULL: 不输出)
* 返回值 : status(1:正常,0:异常)
*-----------------------------------------------------------------------------*/
extern uint8_t PntPos(const T_PRCOPT *opt,const T_OBS *obs, const uint8_t n, const T_NAV *nav,T_SOL *sol,T_SSAT  *ssat)
{
    float64_t rs[MAXOBS*6]={0},dts[MAXOBS*2]={0},var[MAXOBS]={0},resp[MAXOBS]={0};
	uint8_t i,j,stat;
	sint8_t svh[MAXOBS];
	uint8_t vsat[MAXOBS]={0};
	float64_t azel[MAXOBS*2]={0};
   	
    sol->stat=SOLQ_NONE;
    
    if (n<=0) 	return 0;
    sol->time=obs[0].time;    

    /* 卫星位置、速度、钟差估计 */
    SatPoss(sol->time,obs,n,nav,rs,dts,var,svh);	
    /* 接收机位置估计 */
    stat=EstPos(obs,n,rs,dts,var,svh,nav,opt,sol,azel,vsat,resp);

#if 0    
    /* RAIM 检测 */
    if (!stat&&n>=6) {
        stat=RaimFde(obs,n,rs,dts,var,svh,nav,opt,sol,azel,vsat,resp);
    }
    /* 速度估计 */
    if (stat) {
        EstVel(obs,n,rs,dts,nav,opt,sol,azel,vsat);
    }
#endif
	/*卫星信息存储*/
	if(ssat) 
	{
		float64_t los[COORDIM]={0.0};
		float64_t dis=0.0;
		uint8_t sat  = 0;
								
		for(i=0;i<n;i++)
		{
			sat = obs[i].sat;
			ssat->sat[sat] = sat;
			memcpy(&ssat->pos[sat],rs+i*6,COORDIM*sizeof(float64_t));
			memcpy(&ssat->vel[sat],rs+i*6+3,COORDIM*sizeof(float64_t));
			ssat->ele[sat] = azel[i*2+1];
			ssat->azi[sat] = azel[i*2];
			dis = GeoDist(rs + i * 6, sol->rr, los);
			
			memcpy(&ssat->los[sat],los,COORDIM*sizeof(float64_t));
			ssat->dis[sat] = dis;			
		}
		ssat->n = n;
	}
    return stat;
}
/* 删除给定卫星的观测量 ----------------------------------------------------
* 参数   :   T_OBSS *obss      IO   观测数据
*          uint8_t sat      I   要删除的卫星编号
* 返回值 : 无
*-----------------------------------------------------------------------------*/

extern void DeleteObs(T_OBSS *obss,uint8_t sat)
{
	uint8_t i;
 
	for(i=0;i<obss->n;i++)
	{
		if(obss->data[i].sat==sat)
		{
			if(obss->data[i].rcv==1)  
				obss->nr--;
			else
				obss->nb--;
			memcpy(&obss->data[i], &obss->data[i + 1], sizeof(T_OBS)*(obss->n - i - 1));
			memset(&obss->data[obss->n - 1], 0, sizeof(T_OBS));
			obss->n = obss->n-1;			
		}
	}
}
/* 调整观测值，筛选基准站、移动站共有卫星、删除高度角、载噪比低的卫星
* 参数   :   T_PRCOPT *opt      I   解算配置
*          T_OBSS *obss      IO   观测数据
*          T_SSAT ssat      I   卫星信息
* 返回值 : 无
*-----------------------------------------------------------------------------*/

extern void AdjObss(const T_PRCOPT *opt,T_OBSS *obss,T_SSAT *ssat)
{
	uint8_t i,j;	
	/*search commom satellite*/	
	for(i=0;i<obss->nr;i++) //start from rover station
	{
		for(j = obss->nr;j<obss->n;j++)
		{
			if(obss->data[i].sat == obss->data[j].sat)		break;		
		}
		if (j == obss->n)
		{
			DeleteObs(obss, obss->data[i].sat);
			i--;
		}
	}

	for(i = obss->nr;i<obss->n;i++)  //start from base station
	{
		for(j=0;j<obss->nr;j++)
		{
			if(obss->data[i].sat == obss->data[j].sat)		break;	
		}						
		if (j == obss->nr)
		{
			DeleteObs(obss, obss->data[i].sat);
			i--;
		}
	}
	/*delete satellite low elevation*/
	for(i=0;i<obss->nr;i++)
	{
		if (ssat->ele[obss->data[i].sat] < opt->elmin)
		{
			DeleteObs(obss, obss->data[i].sat);
			i--;
		}
			
	}
	/*delete satellite low cnr*/
	for(i=0;i<obss->nr;i++)
	{
		for(j=0;j<NFREQ;j++)
		{
			if (obss->data[i].SNR[j]*SNR_UNIT < opt->cnrmask && obss->data[i].P[j]!=0)
			{
				//DeleteObs(obss, obss->data[i].sat);
				obss->data[i].P[j] = 0;
				obss->data[i].L[j] = 0;
				i--;
			}
		}			
	}

}
/* 根据双差索引获取卫星的系统和频率，T_DDOBS内部函数
* 参数   :   uint8_t idd      I   双差索引（idd<ddobs.num）
*          uint8_t ddnum[NUM_SYS][NUM_FRE]      I   系统、频点的双差量个数
*          uint8_t isf[3]      0   isf[0]:系统，isf[1]:频率，isf[2]:isf[0]*NUM_FRE+isf[1]
* 返回值 : 无
*-----------------------------------------------------------------------------*/

static uint8_t idd2SysFre(uint8_t idd,uint8_t ddnum[NUM_SYS][NUM_FRE],uint8_t isf[3])
{
	uint8_t i,j,n=0;
	for(i=0;i<NUM_SYS;i++)
		for(j=0;j<NUM_FRE;j++)
			{
				n += ddnum[i][j];
				if (idd<n)
				{
					isf[0]=i;
					isf[1]=j;
					isf[2] = i*NUM_FRE+j;
					return 1;
				}
			}
	return 0;
}

/* 计算双差观测量
* 参数   :   T_PRCOPT opt      I   定位配置
*          T_OBSS obss      I   观测量
*          T_SSAT ssatR      I   移动站卫星信息
*          T_SSAT ssatB      I   基准站卫星信息
*          T_DDOBS ddobs      O   双差观测量
* 返回值 : 无
*-----------------------------------------------------------------------------*/
extern void DDObs(const T_PRCOPT *opt,const T_OBSS *obss,const T_SSAT  *ssatR,const T_SSAT  *ssatB,T_DDOBS *ddobs)
{
	uint8_t i=0,j=0,k=0,ii=0,n = 0;
	uint8_t isf[3]={0};
	uint16_t sat,sys;
	float64_t minele[NUM_SYS][NUM_FRE] = { 0 };
	float64_t dddis, refSDP[NUM_SYS][NUM_FRE] = { 0 }, refSDL[NUM_SYS][NUM_FRE] = {0};
	float64_t lamb,ddmat[COORDIM];
	ddobs->time = obss->data[0].time;
	AdjObss(opt, obss, ssatB);//观测值调整

	/*计算主星*/
	for(i=0;i<NUM_SYS;i++)		
	{
		for(j=0;j<NUM_FRE;j++)
		{
			
			for(k=0;k<obss->nr;k++)
			{
				sat = obss->data[k].sat;
				sys = SatSys(sat, NULL);
				if(!((sys==SYS_GPS && i==0) || (sys == SYS_CMP && i==1)))	continue; //i=0为GPS系统，i=1为北斗系统
				if(fabs(obss->data[k].P[j])<EPSILON || fabs(obss->data[k+obss->nr].P[j])<EPSILON || fabs(obss->data[k].L[j])<EPSILON || fabs(obss->data[k+obss->nr].L[j])<EPSILON)
					continue;
				if (ssatB->ele[sat]>minele[i][j])
				{
					minele[i][j] = ssatB->ele[sat];
					ddobs->refsat[i][j] = sat;
					refSDP[i][j] = obss->data[k].P[j] - obss->data[k+obss->nr].P[j];//主星对应的站间单差值
					refSDL[i][j] = obss->data[k].L[j] - obss->data[k + obss->nr].L[j];
					
				}
			}
		}
	}
	/*双差观测量计算*/
	for(i=0;i<NUM_SYS;i++)		
	{
		for(j=0;j<NUM_FRE;j++)
		{	
			for(k=0;k<obss->nr;k++)
			{
				sat = obss->data[k].sat;
				sys = SatSys(sat, NULL);
				if (!((sys == SYS_GPS && i == 0) || (sys == SYS_CMP && i == 1)))	continue;
				if (sat == ddobs->refsat[i][j]) continue;
				if (!ddobs->refsat[i][j]) continue;
				if(fabs(obss->data[k].P[j])<EPSILON || fabs(obss->data[k+obss->nr].P[j])<EPSILON || fabs(obss->data[k].L[j])<EPSILON || fabs(obss->data[k+obss->nr].L[j])<EPSILON)
					continue;
				ddobs->ddsat[n] = sat;
				ddobs->ddnum[i][j]++;
				ddobs->wavelen[n] = g_wavelen[i][j];
				ddobs->triddL[n] = g_wavelen[i][j]*(obss->data[k].L[j] - obss->data[k + obss->nr].L[j] - refSDL[i][j]);
				ddobs->triddP[n] = obss->data[k].P[j] - obss->data[k + obss->nr].P[j] - refSDP[i][j];
				dddis = ssatR->dis[sat] - ssatB->dis[sat] - (ssatR->dis[ddobs->refsat[i][j]] - ssatB->dis[ddobs->refsat[i][j]]);
				ddobs->ddL[n] =(ddobs->triddL[n] - dddis);
				ddobs->ddP[n] = ddobs->triddP[n] - dddis;
				for (ii = 0; ii < COORDIM; ii++) ddmat[ii] = ssatB->los[ddobs->refsat[i][j]][ii] - ssatB->los[sat][ii];
				memcpy(&ddobs->ddmat[n*COORDIM],&ddmat,COORDIM*sizeof(float64_t));

				n++;
				

			}
		}
	}
	ddobs->num = n;
}

/* 计算宽巷双差观测量
* 参数   :   T_PRCOPT opt      I   定位配置
*          T_OBSS obss      I   观测量
*          T_SSAT ssatR      I   移动站卫星信息
*          T_SSAT ssatB      I   基准站卫星信息
*          T_DDOBS ddobs      O   宽巷双差观测量
* 返回值 : 无
*-----------------------------------------------------------------------------*/
extern void DDObsWL(const T_PRCOPT *opt,const T_OBSS *obss,const T_SSAT  *ssatR,const T_SSAT  *ssatB,T_DDOBS *ddobs)
{
	uint8_t i=0,k=0,ii=0,n = 0;
	uint8_t isf[3]={0};
	uint16_t sat,sys;
	float64_t minele[NUM_SYS] = { 0 };
	float64_t dddis, refSDP[NUM_SYS] = { 0 }, refSDL[NUM_SYS] = {0};
	float64_t lamb,ddmat[COORDIM];
	ddobs->time = obss->data[0].time;
	AdjObss(opt, obss, ssatB);

	/*get reference satellite*/
	for(i=0;i<NUM_SYS;i++)		
	{			
		for(k=0;k<obss->nr;k++)
		{
			sat = obss->data[k].sat;
			sys = SatSys(sat, NULL);
			if(!((sys==SYS_GPS && i==0) || (sys == SYS_CMP && i==1)))	continue;
			if(fabs(obss->data[k].P[0])<EPSILON || fabs(obss->data[k+obss->nr].P[0])<EPSILON || fabs(obss->data[k].L[0])<EPSILON || fabs(obss->data[k+obss->nr].L[0])<EPSILON)
				continue;
			if(fabs(obss->data[k].P[1])<EPSILON || fabs(obss->data[k+obss->nr].P[1])<EPSILON || fabs(obss->data[k].L[1])<EPSILON || fabs(obss->data[k+obss->nr].L[1])<EPSILON)
				continue;

			if (ssatB->ele[sat]>minele[i])
			{
				minele[i] = ssatB->ele[sat];
				ddobs->refsat[i][0] = sat;
				if(i==0)
					refSDP[i] = obss->data[k].P[0] - obss->data[k+obss->nr].P[0];//GPS用L1伪距
				else
					refSDP[i] = obss->data[k].P[1] - obss->data[k+obss->nr].P[1];//BDS用B2伪距
				refSDL[i] = obss->data[k].L[0] - obss->data[k + obss->nr].L[0]-(obss->data[k].L[1] - obss->data[k + obss->nr].L[1]);
				
			}
		}
	}
	/*double-differenced data calculate*/
	for(i=0;i<NUM_SYS;i++)		
	{	
		for(k=0;k<obss->nr;k++)
		{
			sat = obss->data[k].sat;
			sys = SatSys(sat, NULL);
			if (!((sys == SYS_GPS && i == 0) || (sys == SYS_CMP && i == 1)))	continue;
			if (sat == ddobs->refsat[i][0]) continue;
			if (!ddobs->refsat[i][0]) continue;
			
			if(fabs(obss->data[k].P[0])<EPSILON || fabs(obss->data[k+obss->nr].P[0])<EPSILON || fabs(obss->data[k].L[0])<EPSILON || fabs(obss->data[k+obss->nr].L[0])<EPSILON)
				continue;
			if(fabs(obss->data[k].P[1])<EPSILON || fabs(obss->data[k+obss->nr].P[1])<EPSILON || fabs(obss->data[k].L[1])<EPSILON || fabs(obss->data[k+obss->nr].L[1])<EPSILON)
				continue;
			ddobs->ddsat[n] = sat;
			ddobs->ddnum[i][0]++;
			ddobs->wavelen[n] = 1/(1/g_wavelen[i][0]-1/g_wavelen[i][1]);
			ddobs->triddL[n] =ddobs->wavelen[n]*(obss->data[k].L[0] - obss->data[k + obss->nr].L[0] -(obss->data[k].L[1] - obss->data[k + obss->nr].L[1])- refSDL[i]);
			if (i==0)
				ddobs->triddP[n] = obss->data[k].P[0] - obss->data[k + obss->nr].P[0] - refSDP[i];//GPS用L1伪距
			else
				ddobs->triddP[n] = obss->data[k].P[1] - obss->data[k + obss->nr].P[1] - refSDP[i];//BDS用B2伪距

			dddis = ssatR->dis[sat] - ssatB->dis[sat] - (ssatR->dis[ddobs->refsat[i][0]] - ssatB->dis[ddobs->refsat[i][0]]);
			ddobs->ddL[n] =(ddobs->triddL[n] - dddis);
			ddobs->ddP[n] = ddobs->triddP[n] - dddis;
			for (ii = 0; ii < COORDIM; ii++) ddmat[ii] = ssatB->los[ddobs->refsat[i][0]][ii] - ssatB->los[sat][ii];
			memcpy(&ddobs->ddmat[n*COORDIM],&ddmat,COORDIM*sizeof(float64_t));
			
			n++;
			

		}
	}
	ddobs->num = n;
}

/* 高度角随机模型
* 参数   :   float64_t ele      I   卫星高度角
*          float64_t obss      I   初始误差
* 返回值 : 方差值
*-----------------------------------------------------------------------------*/

float64_t VarEle(const float64_t ele, const float64_t sig0)
{
	float64_t var;
	if(ele>30*D2R)
		var = sig0*sig0/sin(ele);
	else
		var = sig0*sig0/sin(ele)/sin(ele);
	return var;
}

/* 双差观测量方差阵计算---------------------------------------------
* 参数   :   T_OBSS obss       I   双差观测量
*          float64_t *sigP   I   伪距观测噪声
*          float64_t *sigL      I   载波相位观测噪声
*          float64_t *eleR      I   移动站高度角
*          float64_t *eleB      I   基准站高度角
*          float64_t *ddQP      O   伪距方差阵
*          float64_t *ddQL      O   载波相位方差阵
* 返回值 : 无
*-----------------------------------------------------------------------------*/

extern void WeightObs(T_DDOBS *ddobs,const float64_t sigP[NUM_SYS][NUM_FRE], const float64_t sigL[NUM_SYS][NUM_FRE],float64_t *eleR, float64_t *eleB, float64_t *ddQP, float64_t *ddQL)

{
	uint8_t i,j,k,n=0,n0=0;
	uint8_t sat,isf[3]={0};
	float64_t ddT[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 }, ddTt[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 };
	float64_t sdvarL[MAXOBS*NUM_FRE]={0},sdvarL_ref[NUM_SYS][NUM_FRE]={0},sdvarP[MAXOBS*NUM_FRE]={0},sdvarP_ref[NUM_SYS][NUM_FRE]={0};
	float64_t sdvarLmat[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0},sdvarPmat[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0};
	float64_t Tmp[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0};
	float64_t varRL,varBL,varRP,varBP;
/*主星观测量单差计算*/
	for(i=0;i<NUM_SYS;i++)
		for(j=0;j<NUM_FRE;j++)
		{
			sat = ddobs->refsat[i][j];
			if(!sat) 
			{
				varRL = 0;varBL = 0;varRP = 0;varBP = 0;
			}
			else
			{
				varRL = VarEle(eleR[sat], sigL[i][j] * g_wavelen[i][j]);
				varBL = VarEle(eleB[sat], sigL[i][j] * g_wavelen[i][j]);
				varRP = VarEle(eleR[sat], sigP[i][j]);
				varBP = VarEle(eleB[sat],sigP[i][j]);
			}			
			sdvarL_ref[i][j] = varRL+varBL;
			sdvarP_ref[i][j] = varRP+varBP;
			sdvarL[n] = sdvarL_ref[i][j];
			sdvarP[n] = sdvarP_ref[i][j];
			n++;			
		}
	n0 = n; //n0为主星数量
/*从星观测量单差计算*/
	for (i = 0; i<ddobs->num; i++)
	{
		idd2SysFre(i, ddobs->ddnum, isf);
		sat = ddobs->ddsat[i];
		varRL = VarEle(eleR[sat], sigL[isf[0]][isf[1]] * g_wavelen[isf[0]][isf[1]]);
		varBL = VarEle(eleB[sat], sigL[isf[0]][isf[1]] * g_wavelen[isf[0]][isf[1]]);
		varRP = VarEle(eleR[sat], sigP[isf[0]][isf[1]]);
		varBP = VarEle(eleB[sat], sigP[isf[0]][isf[1]]);
		sdvarL[n] = varRL+varBL;
		sdvarP[n] = varRP+varBP;

		
		n++;	
	}
	Diag(sdvarL, n, sdvarLmat);
	Diag(sdvarP, n, sdvarPmat);
/*单差转双差矩阵ddT*/

	for (i = 0; i<ddobs->num; i++)
	{
		idd2SysFre(i, ddobs->ddnum, isf);
		ddT[i + (ddobs->num)*isf[2]] = -1.0;
		ddT[i + (ddobs->num)*(i + n0)] = 1.0;
	}
/*双差方差阵计算 ，ddQL = ddT*sdvarLmat*ddTt, ddQP = ddt*sdvarPmat*ddTt ---*/
	TransMat(ddT, ddobs->num, n, ddTt);
	MulMat(ddT, sdvarLmat, ddobs->num, n, n, Tmp);
	MulMat(Tmp, ddTt, ddobs->num, n, ddobs->num, ddQL);
	
	MulMat(ddT, sdvarPmat, ddobs->num, n, n, Tmp);
	MulMat(Tmp, ddTt, ddobs->num, n, ddobs->num, ddQP);
				
}
/* 宽巷双差观测量方差阵计算---------------------------------------------
* 参数   :   T_OBSS obss       I   双差观测量
*          float64_t *sigP   I   伪距观测噪声
*          float64_t *sigL      I   载波相位观测噪声
*          float64_t *eleR      I   移动站高度角
*          float64_t *eleB      I   基准站高度角
*          float64_t *ddQP      O   伪距方差阵
*          float64_t *ddQL      O   载波相位方差阵
* 返回值 : 无
*-----------------------------------------------------------------------------*/

extern void WeightObsWL(T_DDOBS *ddobs,const float64_t sigP[NUM_SYS][NUM_FRE], const float64_t sigL[NUM_SYS][NUM_FRE],float64_t *eleR, float64_t *eleB, float64_t *ddQP, float64_t *ddQL)
{
	uint8_t i, j, k, n = 0, n0 = 0;
	uint8_t sat, isf[3] = { 0 };
	float64_t ddT[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 }, ddTt[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 };
	float64_t sdvarL[MAXOBS*NUM_FRE] = { 0 }, sdvarL_ref[NUM_SYS] = { 0 }, sdvarP[MAXOBS*NUM_FRE] = { 0 }, sdvarP_ref[NUM_SYS] = { 0 };
	float64_t sdvarLmat[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 }, sdvarPmat[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 };
	float64_t Tmp[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE] = { 0 };
	float64_t varRL, varBL, varRP, varBP;
	/*主星观测量单差计算*/
	for (i = 0; i<NUM_SYS; i++)
	{
		sat = ddobs->refsat[i][0];
		if (!sat)
		{
			varRL = 0; varBL = 0; varRP = 0; varBP = 0;
		}
		else
		{
			varRL = VarEle(eleR[sat], sigL[i][0] * g_wavelen[i][0] );
			varBL = VarEle(eleB[sat], sigL[i][0] * g_wavelen[i][0] );
			varRP = VarEle(eleR[sat], sigP[i][0]);
			varBP = VarEle(eleB[sat], sigP[i][0]);
		}
		sdvarL_ref[i] = varRL + varBL;
		sdvarP_ref[i] = varRP + varBP;
		sdvarL[n] = sdvarL_ref[i] * g_coefWL[i];
		sdvarP[n] = sdvarP_ref[i];
		n++;
	}
	n0 = n;
/*从星观测量单差计算*/	
	for (i = 0; i<ddobs->num; i++)
	{
		idd2SysFre(i, ddobs->ddnum, isf);
		sat = ddobs->ddsat[i];
		varRL = VarEle(eleR[sat], sigL[isf[0]][isf[1]] * g_wavelen[isf[0]][isf[1]]);
		varBL = VarEle(eleB[sat], sigL[isf[0]][isf[1]] * g_wavelen[isf[0]][isf[1]]);
		varRP = VarEle(eleR[sat], sigP[isf[0]][isf[1]]);
		varBP = VarEle(eleB[sat], sigP[isf[0]][isf[1]]);
		sdvarL[n] = (varRL + varBL)*g_coefWL[isf[0]];//宽巷载波方差需乘上转换因子
		sdvarP[n] = varRP + varBP;


		n++;
	}
	Diag(sdvarL, n, sdvarLmat);
	Diag(sdvarP, n, sdvarPmat);
	/*单差转双差矩阵ddT*/
	for (i = 0; i<ddobs->num; i++)
	{
		idd2SysFre(i, ddobs->ddnum, isf);
		ddT[i + (ddobs->num)*isf[0]] = -1.0;
		ddT[i + (ddobs->num)*(i + n0)] = 1.0;
	}
	/*双差方差阵计算 ，ddQL = ddT*sdvarLmat*ddTt, ddQP = ddt*sdvarPmat*ddTt ---*/

	TransMat(ddT, ddobs->num, n, ddTt);
	MulMat(ddT, sdvarLmat, ddobs->num, n, n, Tmp);
	MulMat(Tmp, ddTt, ddobs->num, n, ddobs->num, ddQL);

	MulMat(ddT, sdvarPmat, ddobs->num, n, n, Tmp);
	MulMat(Tmp, ddTt, ddobs->num, n, ddobs->num, ddQP);

}
/* 浮点解计算 ----------------------------------------------------
* 最小二乘法计算基线和模糊度浮点解
* 参数  : T_DDOBS *ddobs      I   双差观测数据
*		float64_t *ddQP     I   双差伪距方差阵
*		float64_t *ddQL     I   双差载波相位方差阵
*		float64_t *bfloat	O	基线浮点解		
*       float64_t *afloat	O	模糊度浮点解
*		float64_t *Qbfloat	O	基线浮点解方差阵
*		float64_t *Qafloat	O	模糊度浮点解方差阵
*		float64_t *Qbafloat	O	基线-模糊度浮点解协方差阵
* 返回值 : 0，解算出错，1：正常
*-----------------------------------------------------------------------------*/

uint8_t CalFloatSol(const T_DDOBS *ddobs,const float64_t *ddQP,const float64_t *ddQL,float64_t *bfloat,float64_t *afloat,float64_t *Qbfloat,float64_t *Qafloat,float64_t *Qbafloat)
{
	float64_t y[MAXOBS*NUM_FRE * 2] = { 0 }, x[COORDIM + MAXOBS*NUM_FRE] = {0}, A[MAXOBS*NUM_FRE * 2 * (COORDIM + MAXOBS*NUM_FRE)] = { 0 };
	float64_t O[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE], D[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE], B[MAXOBS*NUM_FRE*COORDIM], tmp1[MAXOBS*NUM_FRE * 2 * COORDIM], tmp2[MAXOBS*NUM_FRE * 2 * MAXOBS*NUM_FRE];
	float64_t Q[MAXOBS*NUM_FRE*2*MAXOBS*NUM_FRE*2],Qx[(COORDIM+MAXOBS*NUM_FRE)*(COORDIM+MAXOBS*NUM_FRE)];

	
	CombineMat(1,ddobs->ddP,ddobs->num,1,ddobs->ddL,ddobs->num,1,y);//y=[ddP;ddL]	
	
	Zeros(O, ddobs->num, ddobs->num);
	Diag(ddobs->wavelen,ddobs->num,D);
	TransMat(ddobs->ddmat, COORDIM, ddobs->num, B);
	CombineMat(1, B, ddobs->num, COORDIM, B, ddobs->num, COORDIM, tmp1);
	CombineMat(1,O,ddobs->num,ddobs->num,D,ddobs->num,ddobs->num,tmp2);
	CombineMat(0, tmp1, ddobs->num * 2, COORDIM, tmp2, ddobs->num * 2, ddobs->num, A);//A = [H,O;H,D]
	
	BlkDiag(ddQP,ddobs->num,ddobs->num,ddQL,ddobs->num,ddobs->num,Q);//Q=[ddQP,O;O,ddQL]
	
	if(!Lsq(A,y,Q, ddobs->num*2,ddobs->num+COORDIM, x,Qx))//y=Ax
		return 0;
	
	SelectMat(x, ddobs->num + COORDIM, 1, 0, 0, COORDIM, 1, bfloat);//bfloat = x[1:3]
	SelectMat(x, ddobs->num + COORDIM, 1,COORDIM, 0, ddobs->num, 1, afloat);//afloat = x[3:end]
	SelectMat(Qx, ddobs->num + COORDIM, ddobs->num + COORDIM, 0, 0, COORDIM, COORDIM, Qbfloat);
	SelectMat(Qx, ddobs->num + COORDIM, ddobs->num + COORDIM, COORDIM, COORDIM, ddobs->num, ddobs->num, Qafloat);
	SelectMat(Qx, ddobs->num + COORDIM, ddobs->num + COORDIM, 0, COORDIM, COORDIM, ddobs->num, Qbafloat);
	return 1;
}

/*条件最小二乘求固定解
公式：bfix = bfloat-Qba*Qaf^-1*(afloat-afix)
* 参数  : uint8_t n      I   模糊度数量
*		float64_t *bfloat     I   基线浮点解
*		float64_t *afloat     I   模糊度浮点解
*		float64_t *afix		I	模糊度固定解		
*		float64_t *Qafloat	I	模糊度浮点解方差阵
*		float64_t *Qbafloat	I	基线-模糊度浮点解协方差阵
*		float64_t *bfix		O	基线固定解

* 返回值 : 无

*/
void FixBaseline(uint8_t n,const float64_t * bfloat,const float64_t*afloat, const float64_t *afix, const float64_t * Qafloat, const float64_t * Qbafloat,float64_t *bfix)
{
	float64_t da[MAXOBS*NUM_FRE]={0},iQaf[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0},Tmp[COORDIM*MAXOBS*NUM_FRE]={0},db[COORDIM]={0};
	uint8_t i;
	InvMat(Qafloat, n, iQaf);
	MulMat(Qbafloat, iQaf, COORDIM, n, n, Tmp);
	SubMat(afloat, afix, n, 1, da);
	MulMat(Tmp, da, COORDIM, n, 1, db);
	SubMat(bfloat,db,COORDIM,1,bfix);			
}
/*RTK定位主函数
* 参数  : T_PRCOPT opt      I   定位配置
*		T_OBSS *obst     I   单历元的观测数据
*		T_NAV *navt     I   导航数据
*		T_RTK *rtk		O	定位结果
* 返回值 : 0，异常 ，1正常
*/
extern uint8_t RtkPos(const T_PRCOPT *opt, const T_OBSS *obst, const T_NAV *navt,T_RTK *rtk)
{

	uint8_t i=0;
	float64_t ddQL[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0}, ddQP[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0};
	float64_t bfloat[COORDIM]={0},afloat[MAXOBS*NUM_FRE]={0},bfix[COORDIM]={0},afix[MAXOBS*NUM_FRE]={0};
	float64_t Qbfloat[COORDIM*COORDIM]={0},Qafloat[MAXOBS*NUM_FRE*MAXOBS*NUM_FRE]={0},Qbafloat[COORDIM*MAXOBS*NUM_FRE]={0};
	float64_t mu[3] = {0};
	float64_t Sqnorm[AMBSOLNUM]={0};
	T_SSAT SsatR = { 0 }, SsatB = { 0 };
	T_DDOBS ddobs = { 0 };

	rtk->time = obst->data[0].time;
	if (!PntPos(opt, obst->data, obst->nr, navt, &rtk->solrover, &SsatR)) //移动站单点定位
		return 0;

	if (!PntPos(opt, obst->data + obst->nr, obst->nb, navt, &rtk->solbase, &SsatB)) //基准站单点定位
		return 0;
	if (opt->modear == MOD_WL) //宽巷组合RTK
	{
		DDObsWL(opt, obst, &SsatR,&SsatB, &ddobs); //宽巷双差观测量获取
		WeightObsWL(&ddobs, opt->sigP, opt->sigL, SsatR.ele, SsatB.ele, ddQP, ddQL); //宽巷观测量噪声方差阵	
	}
	else 
	{

		DDObs(opt, obst, &SsatR, &SsatB, &ddobs); //双差观测量获取
		WeightObs(&ddobs, opt->sigP, opt->sigL, SsatR.ele, SsatB.ele, ddQP, ddQL); //观测量噪声方差阵	
	}

	if(!CalFloatSol(&ddobs,ddQP,ddQL,bfloat,afloat,Qbfloat,Qafloat,Qbafloat))//浮点解
		return 0;	
	if(LambdaFFRT(ddobs.num, AMBSOLNUM, afloat, Qafloat, PROB_FAIL,afix, Sqnorm, mu)) //整周模糊度固定
		return 0;
	FixBaseline(ddobs.num, bfloat, afloat, afix, Qafloat, Qbafloat, bfix);//固定解	

	/*RTK定位结果赋值*/
	for (i = 0; i<COORDIM; i++)
	{
		bfix[i] = rtk->solrover.rr[i] - rtk->solbase.rr[i] + bfix[i];
		rtk->baseline[i] = bfix[i];
		rtk->solrover.rr[i] = rtk->solbase.rr[i] + bfix[i];
	}
	rtk->baselinelen = Norm(rtk->baseline, 3);
	
	rtk->ambfix.time = ddobs.time;
	rtk->ambfix.num = ddobs.num;
	memcpy(rtk->ambfix.refsat,ddobs.refsat,NUM_SYS*NUM_FRE*sizeof(uint8_t));
	memcpy(rtk->ambfix.ddnum,ddobs.ddnum,NUM_SYS*NUM_FRE*sizeof(uint8_t));
	memcpy(rtk->ambfix.ddsat,ddobs.ddsat,ddobs.num*sizeof(uint8_t));
	memcpy(rtk->ambfix.wavelen,ddobs.wavelen,ddobs.num*sizeof(float64_t));
	memcpy(rtk->ambfix.ddamb,afix,ddobs.num*sizeof(float64_t));
	rtk->ambfix.ratioThreshold = mu[0];
	rtk->ambfix.ps = mu[1];
	rtk->ambfix.ratio = mu[2];
	
	if ((rtk->ambfix.ratio-rtk->ambfix.ratioThreshold)>opt->lambdaffrtmask) 
	{
		rtk->solrover.stat = SOLQ_FIX;
		rtk->ambfix.stat = SOLQ_FIX;
	}
		
	else
	{
		rtk->solrover.stat = SOLQ_FLOAT;
		rtk->ambfix.stat = SOLQ_FLOAT;
	}
		
	fprintf(ofp[0], "\r\n%f %d %d %f %f %f %f %f %f %f", ddobs.time, rtk->solrover.stat, ddobs.num, rtk->ambfix.ratioThreshold, rtk->ambfix.ratio, rtk->ambfix.ps,bfix[0], bfix[1], bfix[2], rtk->baselinelen);
					
}