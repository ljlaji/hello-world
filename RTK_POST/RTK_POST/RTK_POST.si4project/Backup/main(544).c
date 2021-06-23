#include "rtklib.h"
#include <process.h>
#include <windows.h>
#include "RelPos.h"

#define RINEX 1
#define HTNEX 0

extern int ObsRead(FILE *fp, obs_t *obss); //读取汇天自定义观测文件
extern int RtcmRead(FILE *fp, obs_t *obss, double *ep);
extern int NavUpdate(FILE *fp, nav_t *nav);
extern void Time2epoch(FP64 t, FP64 *ep);
int iobs = 0;
extern T_PRCOPT opt_default;
int GetEpochObs(obs_t obss, obsd_t * epochobs)
{
	int n = 1, iobs0;
	iobs0 = iobs;

	while (fabs(obss.data[iobs].time.time - obss.data[iobs+1].time.time) < DTTOL)
	{
		iobs++;
		n++;
	}
	memcpy(epochobs, &(obss.data[iobs0]), n*sizeof(obsd_t));
	iobs++;
	return n;
}
void fillNavData(const nav_t nav, T_NAV *navt)
{
	UINT32 i=0;
	navt->n = nav.n;
	navt->nmax = nav.nmax;
	navt->ng = nav.ng;
	navt->ngmax = nav.ngmax;
	for(i=0;i<nav.n;i++)
	{
		navt->eph[i].sat = nav.eph[i].sat;
		navt->eph[i].iode = nav.eph[i].iode;
		navt->eph[i].iodc = nav.eph[i].iodc;
		navt->eph[i].sva = nav.eph[i].sva;
		navt->eph[i].svh = nav.eph[i].svh;
		navt->eph[i].week = nav.eph[i].week;
		navt->eph[i].code = nav.eph[i].code;
		navt->eph[i].toe = (double)nav.eph[i].toe.time+nav.eph[i].toe.sec;
		navt->eph[i].toc = (double)nav.eph[i].toc.time+nav.eph[i].toc.sec;
		navt->eph[i].ttr = (double)nav.eph[i].ttr.time+nav.eph[i].ttr.sec;
		navt->eph[i].A = nav.eph[i].A;
		navt->eph[i].e = nav.eph[i].e;
		navt->eph[i].i0 = nav.eph[i].i0;
		navt->eph[i].OMG0 = nav.eph[i].OMG0;
		navt->eph[i].omg = nav.eph[i].omg;
		navt->eph[i].M0 = nav.eph[i].M0;
		navt->eph[i].deln = nav.eph[i].deln;
		navt->eph[i].OMGd = nav.eph[i].OMGd;
		navt->eph[i].idot = nav.eph[i].idot;
		navt->eph[i].crc = nav.eph[i].crc;
		navt->eph[i].crs = nav.eph[i].crs;
		navt->eph[i].cuc = nav.eph[i].cuc;
		navt->eph[i].cus = nav.eph[i].cus;
		navt->eph[i].cic = nav.eph[i].cic;
		navt->eph[i].cis = nav.eph[i].cis;
		navt->eph[i].toes = nav.eph[i].toes;
		navt->eph[i].fit = nav.eph[i].fit;
		navt->eph[i].f0 = nav.eph[i].f0;
		navt->eph[i].f1 = nav.eph[i].f1;
		navt->eph[i].f2 = nav.eph[i].f2;
		memcpy(navt->eph[i].tgd,nav.eph[i].tgd,sizeof(FP64)*6);

		
	}
	
	
}

void fillObsData(const obsd_t *obs, UINT8 n, T_OBSS *obss)
{
	UINT8 i;
	UINT8 nu=0, nb=0;
	obss->n = n;
	for(i=0;i<n;i++)
	{
		obss->data[i].time = (FP64)obs[i].time.time+ obs[i].time.sec;
		obss->data[i].sat = obs[i].sat;
		if (obs[i].rcv == 1)
			nu++;
		else
			nb++;
			
		obss->data[i].rcv = obs[i].rcv;			
		memcpy(obss->data[i].SNR,obs[i].SNR,NFREQ*sizeof(FP64));
		memcpy(obss->data[i].LLI,obs[i].LLI,NFREQ*sizeof(FP64));
		memcpy(obss->data[i].code,obs[i].code,NFREQ*sizeof(FP64));
		memcpy(obss->data[i].L,obs[i].L,NFREQ*sizeof(FP64));		
		memcpy(obss->data[i].P,obs[i].P,NFREQ*sizeof(FP64));
		memcpy(obss->data[i].D,obs[i].D,NFREQ*sizeof(FP64));		
		
	}
	obss->nu = nu;
	obss->nb = nb;

		
}
int main()
{
	FILE *ifp[3], *ofp[2];
	int i,nsat, nepoch;
	double ep[6] = { 2021, 5, 25, 0, 0, 0 }, rr[3] = { 0 }, pos[3] = { 0 }, enu[3] = { 0 };
	obs_t obss = {0};
	obsd_t epochobs[MAXSAT] = { 0 };
	nav_t nav = { 0 };
	rtk_t rtk = { 0 };
	prcopt_t prcopt = prcopt_default;
	int ifix = 0, isol = 0;

	T_PRCOPT opt = {0};
	opt.cnrmask = 32;
	opt.elmin = 15.0*D2R;
	opt.maxgdop = 30;
	T_OBSS obst = { 0 };
	T_NAV navt = { 0 };
	T_SOL sol = { 0 };


#if RINEX
	char *infiles[3] = { "E:\\Project\\data\\20210525\\gnss_data.obs",
						"E:\\Project\\data\\20210525\\B1145M.21O",
						"E:\\Project\\data\\20210525\\gnss_data.nav" };
	char *outfile = { "E:\\Project\\data\\20210525\\Rinex_result.pos" };

	//char *infiles[3] = { "E:\\Project\\data\\20180319B\\20180319a.obs",
	//	"E:\\Project\\data\\20180319B\\20180319b.obs",
	//	"E:\\Project\\data\\20180319B\\20180319a.nav" };
	//char *outfile = { "E:\\Project\\data\\20180319B\\result.pos" };

	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0;
	ReadRinex(infiles, 3, &obss, &nav, &nepoch);

#elif HTNEX
	char *infiles[3] = { "E:\\Project\\data\\20210525\\ReceivedTofile-TCPCLIENT-2021_5_25_20-38-52.DAT",
		"E:\\Project\\data\\20210525\\range_data(3).txt",//"E:\\Project\\data\\20200525\\range_data(3).txt",
		"E:\\Project\\data\\20210525\\ephem_data(2).txt" };
	char *outfile = { "E:\\Project\\data\\20210525\\Htnex_result.pos" };

	if ((ifp[0] = fopen(infiles[0], "rb")) == NULL) return 0;
	if ((ifp[1] = fopen(infiles[1], "rb")) == NULL) return 0;
	if ((ifp[2] = fopen(infiles[2], "rb")) == NULL) return 0;
	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0;

	if (!NavUpdate(ifp[2], &nav))return 0;
	if (!RtcmRead(ifp[0], &obss, ep)) return 0;
	if (!ObsRead(ifp[1], &obss)) return 0;
	for (i = 0; i < obss.n; i++)
	{
		time2epoch(obss.data[i].time, ep);
		fprintf(ofp[0], "\n RTK %4d%02d%02d %02d %02d %05.2f ", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
		fprintf(ofp[0], "%02d %d %.3f %.3f %.3f P %.3f %.3f %.3f", obss.data[i].sat, obss.data[i].rcv, obss.data[i].L[0], obss.data[i].L[1], obss.data[i].L[2],obss.data[i].P[0], obss.data[i].P[1], obss.data[i].P[2]);
	}
	
#endif

#if 1
	fillNavData(nav, &navt);
	while (1)
	{
		nsat = GetEpochObs(obss, epochobs);
		fillObsData(epochobs, nsat, &obst);

		if (iobs >= obss.n) break;		
		Time2epoch(obst.data[0].time, ep);
		printf("\rRTK-test solution:%4d%02d%02d %02d %02d %05.2f", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
		if ((double)epochobs[0].time.time < 1621946462) continue;
		isol++;
		//if(!rtkpos(&rtk, epochobs, nsat, &nav)) continue;	
		if(!PntPos(&opt, obst.data, obst.nu, &navt, &sol)) continue;
		ifix++;
		fprintf(ofp[0], "\n RTK %4d%02d%02d %02d %02d %05.2f ", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
		fprintf(ofp[0], "%d,%.2f,%.2f,%.2f", sol.stat, sol.rr[0],sol.rr[1], sol.rr[2]);
		
		
	}

	printf("\ntotal epoch :%d\nfix epoch:%d\nfix ratio:%f", isol, ifix, (double)ifix / (double)isol);
	getchar();
	

	/*fclose(ifp[0]);
	fclose(ifp[1]);
	fclose(ifp[2]);*/
#endif
	fclose(ofp[0]);

	return 0;
}