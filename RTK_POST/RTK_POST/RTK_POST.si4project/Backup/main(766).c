#include "rtklib.h"
#include <process.h>
#include <windows.h>

#define RINEX 0
#define HTNEX 1

int iobs = 0;
rtcm_t rtcm = {0};
extern int ReadRinex(char **files, int n, obs_t *obs, nav_t *nav, int *nepoch);

void ObsUpdate(FILE *fp,obs_t *obs)
{

}
int RtcmUpdate(FILE *fp, obs_t *obs)
{
	int i;
	int ret;
	//lloc(obss)
	while(1)
	{
		ret = input_rtcm3f(&rtcm, fp);
		if(ret ==-2) return 0; 
		else if(ret == 1)
		{
			for(i=0;i<rtcm.obs.n;i++) rtcm.obs.data[i].rcv = 2;				
			memcpy(&(obs->data[0]), &(rtcm.obs.data), sizeof(obsd_t)*rtcm.obs.n);
			obs->n = obs->n+rtcm.obs.n;	
		}														
	}
	
}

int NavUpdate(FILE *fp, nav_t *nav)
{

}
int GetEpochObs(obs_t obss,obsd_t * epochobs)
{
	int n = 1,iobs0;
	iobs0 = iobs;
	
	while (fabs(timediff(obss.data[iobs].time, obss.data[iobs + 1].time)) < DTTOL)
	{
		iobs++;
		n++;
	}
	memcpy(epochobs, &(obss.data[iobs0]),n*sizeof(obsd_t));
	iobs++;
	return n;
}
int main()
{
#if RINEX
	char *infiles[3] = { "E:\\Project\\data\\20200525\\gnss_data2.obs",
						"E:\\Project\\data\\20200525\\B1145M.21O",
						"E:\\Project\\data\\20200525\\gnss_data.nav" };
	char *outfile = { "E:\\Project\\data\\20200525\\result.pos" };

	FILE *ifp[3], *ofp[2];
	int i,nsat, nepoch;
	double ep[6] = {0},rr[3] = { 0 }, pos[3] = { 0 }, enu[3] = { 0 };
	obs_t obss = {0};
	obsd_t epochobs[MAXSAT] = { 0 };
	nav_t nav = { 0 };
	rtk_t rtk = { 0 };	
	prcopt_t prcopt=prcopt_default;

	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0; 
	rtkinit(&rtk, &prcopt);
	ReadRinex(infiles, 3, &obss, &nav, &nepoch);

	while (1)
	{
		nsat = GetEpochObs(obss, epochobs);
		if (iobs >= obss.n) break;
		time2epoch(epochobs[0].time, ep);		
		if(!rtkpos(&rtk, epochobs, nsat, &nav)) continue;		
		for (i = 0; i<3; i++) rr[i] = rtk.sol.rr[i] - rtk.rb[i];
		ecef2pos(rtk.rb, pos);
		ecef2enu(pos, rr, enu);
		fprintf(ofp[0], "\n RTK %4d%02d%02d %02d %02d %05.2f ", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
		fprintf(ofp[0], "%d,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f", rtk.sol.stat, rtk.sol.rr[0], rtk.sol.rr[1], rtk.sol.rr[2],enu[0],enu[1],enu[2]);

	}
	fclose(ofp[0]);
#endif

#if HTNEX
	char *infiles[3] = { "E:\\Project\\data\\20200525\\ReceivedTofile-TCPCLIENT-2021_5_25_20-38-52.DAT",
		"E:\\Project\\data\\20200525\\range_data1",
		"E:\\Project\\data\\20200525\\ephem_data" };
	char *outfile = { "E:\\Project\\data\\20200525\\result.pos" };
	double ep[6] = {2020,5,25,0,0};
	FILE *ifp[3], *ofp[2];
	obs_t obss = { 0 };

	if ((ifp[0] = fopen(infiles[0], "rb")) == NULL) return 0;
	if ((ifp[1] = fopen(infiles[1], "rb")) == NULL) return 0;
	if ((ifp[2] = fopen(infiles[2], "rb")) == NULL) return 0;
	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0;

	init_rtcm(&rtcm);
	rtcm.time = epoch2time(ep);

	RtcmUpdate(ifp[0], &obss);



#endif
	return 0;
}