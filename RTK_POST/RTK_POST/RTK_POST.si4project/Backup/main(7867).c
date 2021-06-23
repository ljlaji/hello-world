#include "rtklib.h"
#include <process.h>
#include <windows.h>

rtcm_t rtcm;

void ObsUpdate(FILE *fp,obs_t *obs)
{

}
int RtcmUpdate(FILE *fp, obs_t *obs)
{
	int i;
	while(1)
	{
		if(input_rtcm3f(&rtcm, fp)==-2) return 0; 
		if(timediff(rtcm.time,obs->data[0].time)>DTTOL)	continue;
		for(i=0;i<rtcm.obs.n;i++)
		{
			rtcm.obs.data[i].rcv = 2;
			obs->data[obs->n+i] = rtcm.obs.data[i];
		}
		obs->n = obs->n+rtcm.obs.n;				
	}
	
}

int NavUpdate(FILE *fp, nav_t *nav)
{

}

int main()
{
	int i, flag;
	char *rtcmfile = { "Data_20201123.dat" };
	char *navfile = {""};
	char *obsfile = { "" };
	char *outfile = { "result.txt" };
	FILE *ifp[3], *ofp[2];

	rtcm_t rtcm = { 0 };



	double timeref[6] = { 2020, 11, 23, 0, 0, 0 };
	obs_t epochobs = {0};
	nav_t nav = {0};
	prcopt_t prcopt=prcopt_default;
	sol_t sol = {0};
	ssat_t ssat[MAXSAT] = {0};
	double azel[MAXSAT] = {0};
	char msg[100];
	double ep[6]={0};
	int feph = 0;

	if ((ifp[0] = fopen(rtcmfile, "rb")) == NULL) return 0;
	if ((ifp[1] = fopen(navfile, "rb")) == NULL) return 0;
	if ((ifp[2] = fopen(obsfile, "rb")) == NULL) return 0;
	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0; 
	
	init_rtcm(&rtcm);
	rtcm.time = epoch2time(timeref);

	
	while (1)
	{
		if(!NavUpdate(ifp[1], &nav)) continue;
		ObsUpdate(ifp[2], &epochobs);
		pntpos(epochobs.data, epochobs.n, &nav, &prcopt, &sol, NULL, NULL, msg);
		RtcmUpdate(ifp[0], &epochobs);

	
		
		
		
		
		if(input_rtcm3f(&rtcm, ifp[0])==-2) break;

		fprintf(ofp, "\n%s", rtcm.msgtype);

		if (rtcm.nmsg3[19] >= 10) Feph = 1;//星历数>10才开启定位

		if(!Feph) continue;

		EpochObs = rtcm.obs;
		Nav = rtcm.nav;

	//	for (i = 0; i < EpochObs.n; i++) printf("\nOBS:%d %.2f,%.2f", EpochObs.data[i].sat, EpochObs.data[i].P[0], EpochObs.data[i].SNR[0] / 1000.0);
	//	for (i = 0; i < EpochNav.nmax; i++) printf("\nNAV:%d %d", i + 1, EpochNav.eph[i].sat);
		/*if (rtcm.msgtype[7] == '7' && rtcm.msgtype[8] == '4')
		{
			time2epoch(rtcm.time, ep);
			
			printf("\rtime:%4d %2d %2d %2d %2d %.2f", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
			pntpos(EpochObs.data, EpochObs.n, &EpochNav, &prcopt, &sol, NULL, NULL, msg);			
			fprintf(ofp, "\nSta:%4d%2d%2d,%2d %2d %.2f,%.2f,%.2f,%.2f", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5], rtcm.sta.pos[0], rtcm.sta.pos[1], rtcm.sta.pos[2]);
			fprintf(ofp, "\nSPP:%.2f,%d,%.2f,%.2f,%.2f\n", (sol.time.time + sol.time.sec), sol.ns, sol.rr[0], sol.rr[1], sol.rr[2]);
		}*/

		
                  
	}

	fclose(ifp);
	fclose(ofp);
	return 0;
}