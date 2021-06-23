#include "rtklib.h"
#include <process.h>
#include <windows.h>

rtcm_t rtcm = {0};
FILE *ifp, *ofp;
int Feph=0;

void RtcmCollect()
{
	while (1)
	{
		if(input_rtcm3f(&rtcm, ifp)==-2) 
		{
			_endthread();
		}
		if (rtcm.nav.n>5) Feph = 1;	
	}
}



int main()
{
	int i, flag;
	char *RtcmFile = { "raw1_2011230232.dat" };
	char *OutFile = { "result.txt" };
	double TimeRef[6] = { 2020, 11, 23, 0, 0, 0 };
	obs_t EpochObs = {0};
	nav_t EpochNav = {0};
	prcopt_t prcopt=prcopt_default;
	sol_t sol = {0};
	ssat_t ssat[MAXSAT] = {0};
	double Azel[MAXSAT] = {0};
	char msg[100];

	if ((ifp = fopen(RtcmFile, "rb")) == NULL) return 0;
	if ((ofp = fopen(OutFile, "w")) == NULL) return 0; 
	
	init_rtcm(&rtcm);
	rtcm.time = epoch2time(TimeRef);

	//_beginthread(RtcmCollect, 0, NULL);
	
	while (1)
	{
		if(input_rtcm3f(&rtcm, ifp)==-2) break;

	//	fprintf(ofp, "\n%s", rtcm.msgtype);

		if (rtcm.nmsg3[19] >= 10) Feph = 1;//星历数>10才开启定位

		if(!Feph) continue;

		EpochObs = rtcm.obs;
		EpochNav = rtcm.nav;

	//	for (i = 0; i < EpochObs.n; i++) printf("\nOBS:%d %.2f,%.2f", EpochObs.data[i].sat, EpochObs.data[i].P[0], EpochObs.data[i].SNR[0] / 1000.0);
	//	for (i = 0; i < EpochNav.nmax; i++) printf("\nNAV:%d %d", i + 1, EpochNav.eph[i].sat);
		if (rtcm.msgtype[7] == '7' && rtcm.msgtype[8] == '4')
		{
			printf(ofp, "\n%s", rtcm.msgtype);
			pntpos(EpochObs.data, EpochObs.n, &EpochNav, &prcopt, &sol, NULL, NULL, msg);
			fprintf(ofp, "\nSta:%.2f,%.2f,%.2f,%.2f", (rtcm.time.time + rtcm.time.sec), rtcm.sta.pos[0], rtcm.sta.pos[1], rtcm.sta.pos[2]);
			fprintf(ofp, "\nSPP:%.2f,%d,%.2f,%.2f,%.2f\n", (sol.time.time + sol.time.sec), sol.ns, sol.rr[0], sol.rr[1], sol.rr[2]);
		}

		
                  
	}

	fclose(ifp);
	fclose(ofp);
	return 0;
}