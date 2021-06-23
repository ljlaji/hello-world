#include "rtklib.h"
#include <process.h>
#include <windows.h>


rtcm_t *rtcm;
FILE *ifp, *ofp;
int Feph=0;//星历存储完成标志

void RtcmCollect()
{
	while (1)
	{
		if(input_rtcm3f(rtcm, ifp)==-2) 
		{
			_endthread();
		}
		if (rtcm->nav.n>5) Feph = 1;	
	}
}



int main()
{
	int i, flag;
	char *RtcmFile = { "raw1_2011230232.dat" };//输入的rtcm文件
	char *OutFile = { "result.txt" };//输出结果文件
	double TimeRef[6] = { 2020, 11, 23, 0, 0, 0 };//参考时间，不需要很准，一周之内即�?
	obs_t EpochObs = {{0}};
	nav_t EpochNav = {{0}};
	prcopt_t prcopt=prcopt_default;
	sol_t sol={{0}};
	ssat_t ssat={{0}};
	double *Azel[MAXSAT];
	char msg[100];

	if ((ifp = fopen(RtcmFile, "rb")) == NULL) return 0;
	if ((ofp = fopen(OutFile, "w")) == NULL) return 0; 

	if (!(rtcm = (rtcm_t *)calloc(1,sizeof(rtcm_t)))) return 0;
	
	init_rtcm(rtcm);
	rtcm->time = epoch2time(TimeRef);

	//_beginthread(RtcmCollect, 0, NULL);
	
	while (1)
	{
		if(input_rtcm3f(rtcm, ifp)==-2) break;
		if (rtcm->nmsg3[19]>100) Feph = 1;	
		if(!Feph) continue;
		EpochObs = rtcm->obs;
		EpochNav = rtcm->nav;

	//	for (i = 0; i < EpochObs.n; i++) printf("\nOBS:%d %.2f,%.2f", EpochObs.data[i].sat, EpochObs.data[i].P[0], EpochObs.data[i].SNR[0] / 1000.0);
	//	for (i = 0; i < EpochNav.nmax; i++) printf("\nNAV:%d %d", i + 1, EpochNav.eph[i].sat);
		pntpos(EpochObs.data, EpochObs.n, &EpochNav, &prcopt, &sol, Azel, &ssat, msg);
		fprintf(ofp,"rtcm");
                  
	}

	fclose(ifp);
	fclose(ofp);
	return 0;
}