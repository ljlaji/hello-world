/*------------------------------------------------------------------------------
* rtcm.c : rtcm functions
*
*          Copyright (C) 2021 by linjun,HuiTian Xpeng, All rights reserved.
* references:rtklib
* history : 2021/05/19 1.0  new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* function prototypes -------------------------------------------------------*/
extern int decode_rtcm3(rtcm_t *rtcm);

/* constants -----------------------------------------------------------------*/
#define RTCM3PREAMB 0xD3        /* rtcm ver.3 frame preamble */

/* initialize rtcm control -----------------------------------------------------
* initialize rtcm control struct and reallocate memory for observation and
* ephemeris buffer in rtcm control struct
* args   : rtcm_t *raw      IO  rtcm control struct
* return : status (1:ok,0:memory allocation error)
*-----------------------------------------------------------------------------*/
extern int init_rtcm(rtcm_t *rtcm)
{
    gtime_t time0={0};
    obsd_t data0={{0}};
    eph_t  eph0 ={0,-1,-1};
    geph_t geph0={0,-1};
    ssr_t ssr0={{{0}}};
    int i,j;
    
    
    rtcm->staid=rtcm->stah=rtcm->seqno=0;
	rtcm->outtype = 1;
    rtcm->time=rtcm->time_s=time0;
    rtcm->sta.name[0]=rtcm->sta.marker[0]='\0';
    rtcm->sta.antdes[0]=rtcm->sta.antsno[0]='\0';
    rtcm->sta.rectype[0]=rtcm->sta.recver[0]=rtcm->sta.recsno[0]='\0';
    rtcm->sta.antsetup=rtcm->sta.itrf=rtcm->sta.deltype=0;
    for (i=0;i<3;i++) {
        rtcm->sta.pos[i]=rtcm->sta.del[i]=0.0;
    }
    rtcm->sta.hgt=0.0;
    rtcm->dgps=NULL;
    for (i=0;i<MAXSAT;i++) {
        rtcm->ssr[i]=ssr0;
    }
    rtcm->msg[0]=rtcm->msgtype[0]=rtcm->opt[0]='\0';
    for (i=0;i<6;i++) rtcm->msmtype[i][0]='\0';
    rtcm->obsflag=rtcm->ephsat=0;
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ+NEXOBS;j++) {
        rtcm->cp[i][j]=0.0;
        rtcm->lock[i][j]=rtcm->loss[i][j]=0;
        rtcm->lltime[i][j]=time0;
    }
    rtcm->nbyte=rtcm->nbit=rtcm->len=0;
    rtcm->word=0;
    for (i=0;i<100;i++) rtcm->nmsg2[i]=0;
    for (i=0;i<400;i++) rtcm->nmsg3[i]=0;
    
    rtcm->obs.data=NULL;
    rtcm->nav.eph =NULL;
    rtcm->nav.geph=NULL;
    
    /* reallocate memory for observation and ephemeris buffer */
    if (!(rtcm->obs.data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))||
        !(rtcm->nav.eph =(eph_t  *)malloc(sizeof(eph_t )*MAXSAT*2))||
        !(rtcm->nav.geph=(geph_t *)malloc(sizeof(geph_t)*MAXPRNGLO))) {
        free_rtcm(rtcm);
        return 0;
    }
    rtcm->obs.n=0;
    rtcm->nav.n=MAXSAT*2;
	//rtcm->nav.n = 0;
	rtcm->nav.nmax = MAXSAT;
    rtcm->nav.ng=MAXPRNGLO;
    for (i=0;i<MAXOBS   ;i++) rtcm->obs.data[i]=data0;
    for (i=0;i<MAXSAT*2 ;i++) rtcm->nav.eph [i]=eph0;
    for (i=0;i<MAXPRNGLO;i++) rtcm->nav.geph[i]=geph0;
    return 1;
}
/* free rtcm control ----------------------------------------------------------
* free observation and ephemeris buffer in rtcm control struct
* args   : rtcm_t *raw      IO  rtcm control struct
* return : none
*-----------------------------------------------------------------------------*/
extern void free_rtcm(rtcm_t *rtcm)
{
    
    /* free memory for observation and ephemeris buffer */
    free(rtcm->obs.data); rtcm->obs.data=NULL; rtcm->obs.n=0;
    free(rtcm->nav.eph ); rtcm->nav.eph =NULL; rtcm->nav.n=0;
    free(rtcm->nav.geph); rtcm->nav.geph=NULL; rtcm->nav.ng=0;
}
/* input RTCM 3 message from stream --------------------------------------------
* fetch next RTCM 3 message and input a message from byte stream
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          uint8_t data     I   stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 5: input station pos/ant parameters,
*                  10: input ssr messages)
* notes  : before firstly calling the function, time in rtcm control struct has
*          to be set to the approximate time within 1/2 week in order to resolve
*          ambiguity of time in rtcm messages.
*          
*          to specify input options, set rtcm->opt to the following option
*          strings separated by spaces.
*
*          -EPHALL  : input all ephemerides (default: only new)
*          -STA=nnn : input only message with STAID=nnn (default: all)
*          -GLss    : select signal ss for GPS MSM (ss=1C,1P,...)
*          -RLss    : select signal ss for GLO MSM (ss=1C,1P,...)
*          -ELss    : select signal ss for GAL MSM (ss=1C,1B,...)
*          -JLss    : select signal ss for QZS MSM (ss=1C,2C,...)
*          -CLss    : select signal ss for BDS MSM (ss=2I,7I,...)
*          -ILss    : select signal ss for IRN MSM (ss=5A,9A,...)
*          -GALINAV : select I/NAV for Galileo ephemeris (default: all)
*          -GALFNAV : select F/NAV for Galileo ephemeris (default: all)
*
*          supported RTCM 3 messages (ref [7][10][15][16][17][18])
*
*            TYPE       :  GPS   GLONASS Galileo  QZSS     BDS    SBAS    NavIC
*         ----------------------------------------------------------------------
*          OBS COMP L1  : 1001~   1009~     -       -       -       -       -
*              FULL L1  : 1002    1010      -       -       -       -       -
*              COMP L1L2: 1003~   1011~     -       -       -       -       -
*              FULL L1L2: 1004    1012      -       -       -       -       -
*
*          NAV          : 1019    1020    1045**  1044    1042      -     1041
*                           -       -     1046**    -       63*     -       -
*
*          MSM 1        : 1071~   1081~   1091~   1111~   1121~   1101~   1131~
*              2        : 1072~   1082~   1092~   1112~   1122~   1102~   1132~
*              3        : 1073~   1083~   1093~   1113~   1123~   1103~   1133~
*              4        : 1074    1084    1094    1114    1124    1104    1134
*              5        : 1075    1085    1095    1115    1125    1105    1135 
*              6        : 1076    1086    1096    1116    1126    1106    1136 
*              7        : 1077    1087    1097    1117    1127    1107    1137 
*
*          SSR ORBIT    : 1057    1063    1240*   1246*   1258*     -       -
*              CLOCK    : 1058    1064    1241*   1247*   1259*     -       -
*              CODE BIAS: 1059    1065    1242*   1248*   1260*     -       -
*              OBT/CLK  : 1060    1066    1243*   1249*   1261*     -       -
*              URA      : 1061    1067    1244*   1250*   1262*     -       -
*              HR-CLOCK : 1062    1068    1245*   1251*   1263*     -       -
*              PHAS BIAS:   11*     -       12*     13*     14*     -       -
*
*          ANT/RCV INFO : 1007    1008    1033
*          STA POSITION : 1005    1006
*
*          PROPRIETARY  : 4076 (IGS)
*         ----------------------------------------------------------------------
*                            (* draft, ** 1045:F/NAV,1046:I/NAV, ~ only encode)
*
*          for MSM observation data with multiple signals for a frequency,
*          a signal is selected according to internal priority. to select
*          a specified signal, use the input options.
*
*          RTCM 3 message format:
*            +----------+--------+-----------+--------------------+----------+
*            | preamble | 000000 |  length   |    data message    |  parity  |
*            +----------+--------+-----------+--------------------+----------+
*            |<-- 8 --->|<- 6 -->|<-- 10 --->|<--- length x 8 --->|<-- 24 -->|
*            
*-----------------------------------------------------------------------------*/
extern int input_rtcm3(rtcm_t *rtcm, uint8_t data)
{
    
    /* synchronize frame */
    if (rtcm->nbyte==0) {
        if (data!=RTCM3PREAMB) return 0;
        rtcm->buff[rtcm->nbyte++]=data;
        return 0;
    }
    rtcm->buff[rtcm->nbyte++]=data;
    
    if (rtcm->nbyte==3) {
        rtcm->len=getbitu(rtcm->buff,14,10)+3; /* length without parity */
    }
    if (rtcm->nbyte<3||rtcm->nbyte<rtcm->len+3) return 0;
    rtcm->nbyte=0;
    
    /* check parity */
    if (rtk_crc24q(rtcm->buff,rtcm->len)!=getbitu(rtcm->buff,rtcm->len*8,24)) {
        return 0;
    }
    /* decode rtcm3 message */
    return decode_rtcm3(rtcm);
}

/* input RTCM 3 message from file ----------------------------------------------
* fetch next RTCM 3 message and input a messsage from file
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          FILE  *fp        I   file pointer
* return : status (-2: end of file, -1...10: same as above)
* notes  : same as above
*-----------------------------------------------------------------------------*/
extern int input_rtcm3f(rtcm_t *rtcm, FILE *fp)
{
    int i,data=0,ret;    
    
    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) 
			return -2;
        if ((ret=input_rtcm3(rtcm,(uint8_t)data)))
			return ret;
    }
    return 0; /* return at every 4k bytes */
}
extern int rtcm2obs(rtcm_t *rtcm)
{
	if(!rtcm->obsflag){
		return 0;
	}
	
		
	
}