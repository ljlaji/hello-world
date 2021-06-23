#include "rtklib.h"
#include <process.h>
#include <windows.h>

#define RINEX 1
#define HTNEX 0

extern int ObsRead(FILE *fp, obs_t *obss); //读取汇天自定义观测文件
extern int RtcmRead(FILE *fp, obs_t *obss, double *ep);
extern int NavUpdate(FILE *fp, nav_t *nav);

int iobs = 0;
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
	//rtkinit(&rtk, &prcopt);
#if RINEX
	//char *infiles[3] = { "E:\\Project\\data\\20210525\\gnss_data.obs",
	//					"E:\\Project\\data\\20210525\\B1145M.21O",
	//					"E:\\Project\\data\\20210525\\gnss_data.nav" };
	//char *outfile = { "E:\\Project\\data\\20210525\\Rinex_result.pos" };

	char *infiles[3] = { "E:\\Project\\data\\20180319B\\20180319a.obs",
		"E:\\Project\\data\\20180319B\\20180319b.obs",
		"E:\\Project\\data\\20180319B\\20180319a.nav" };
	char *outfile = { "E:\\Project\\data\\20180319B\\result.pos" };

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

#if 0
	while (1)
	{
		nsat = GetEpochObs(obss, epochobs);
		if (iobs >= obss.n) break;		
		time2epoch(epochobs[0].time, ep);
		isol++;
		printf("\rRTK solution:%4d%02d%02d %02d %02d %05.2f", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
	//	if ((double)epochobs[0].time.time < 1621946462) continue;
		if(!rtkpos(&rtk, epochobs, nsat, &nav)) continue;		
		ifix++;
		for (i = 0; i<3; i++) rr[i] = rtk.sol.rr[i] - rtk.rb[i];
		ecef2pos(rtk.rb, pos);
		ecef2enu(pos, rr, enu);
		fprintf(ofp[0], "\n RTK %4d%02d%02d %02d %02d %05.2f ", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
		fprintf(ofp[0], "%d,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f", rtk.sol.stat, rtk.sol.rr[0], rtk.sol.rr[1], rtk.sol.rr[2],enu[0],enu[1],enu[2]);
		
		
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