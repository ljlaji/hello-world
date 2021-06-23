#include "rtklib.h"
#include <process.h>
#include <windows.h>



void ObsUpdate(FILE *fp,obs_t *obs)
{

}
int RtcmUpdate(FILE *fp, obs_t *obs)
{

}

int NavUpdate(FILE *fp, nav_t *nav)
{

}

int main()
{
	int i, flag;
	char *RtcmFile = { "Data_20201123.dat" };
	char *NavFile = {""};
	char *ObsFile = { "" };
	char *OutFile = { "result.txt" };
	FILE *ifp[3], *ofp[2];

	rtcm_t rtcm = { 0 };



	double TimeRef[6] = { 2020, 11, 23, 0, 0, 0 };
	obs_t EpochObs = {0};
	nav_t Nav = {0};
	prcopt_t prcopt=prcopt_default;
	sol_t sol = {0};
	ssat_t ssat[MAXSAT] = {0};
	double Azel[MAXSAT] = {0};
	char msg[100];
	double ep[6]={0};
	int Feph = 0;

	if ((ifp[0] = fopen(RtcmFile, "rb")) == NULL) return 0;
	if ((ifp[1] = fopen(NavFile, "rb")) == NULL) return 0;
	if ((ifp[2] = fopen(ObsFile, "rb")) == NULL) return 0;
	if ((ofp[0] = fopen(OutFile, "w")) == NULL) return 0; 
	
	init_rtcm(&rtcm);
	rtcm.time = epoch2time(TimeRef);

	
	while (1)
	{
		if(!NavUpdate(ifp[1], &Nav)) continue;
		ObsUpdate(ifp[2], &EpochObs);
		RtcmUpdate(ifp[0], &EpochObs);

	
		
		
		
		
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