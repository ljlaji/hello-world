#include "rtklib.h"
#include <process.h>
#include <windows.h>

extern int ReadRinex(char **files, int n, obs_t *obs, nav_t *nav, int *nepoch);

void ObsUpdate(FILE *fp,obs_t *obs)
{

}
//int RtcmUpdate(FILE *fp, obs_t *obs)
//{
//	int i;
//	while(1)
//	{
//		if(input_rtcm3f(&rtcm, fp)==-2) return 0; 
//		if(timediff(rtcm.time,obs->data[0].time)>DTTOL)	continue;
//		for(i=0;i<rtcm.obs.n;i++)
//		{
//			rtcm.obs.data[i].rcv = 2;
//			obs->data[obs->n+i] = rtcm.obs.data[i];
//		}
//		obs->n = obs->n+rtcm.obs.n;				
//	}
//	
//}

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
	char *infiles[3] = { "E:\\Project\\Rtklib_Post\\data\\20180319B\\20180319a.obs",
		"E:\\Project\\Rtklib_Post\\data\\20180319B\\20180319b.obs",
		"E:\\Project\\Rtklib_Post\\data\\20180319B\\20180319a.nav" };


	rtcm_t rtcm = {0};




	double timeref[6] = { 2018, 03, 19, 0, 0, 0 };
	obs_t obss = {0};

	int nepoch;

	nav_t navs = {0};
	prcopt_t prcopt=prcopt_default;
	sol_t sol = {0};
	ssat_t ssat[MAXSAT] = {0};
	double azel[MAXSAT] = {0};
	char msg[100];
	double ep[6]={0};
	int feph = 0;
	rtk_t rtk={0};



	//if ((ifp[0] = fopen(rtcmfile, "rb")) == NULL) return 0;
	//if ((ifp[1] = fopen(navfile, "rb")) == NULL) return 0;
	//if ((ifp[2] = fopen(obsfile, "rb")) == NULL) return 0;
	//if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0; 
	
	/*init_rtcm(&rtcm);
	rtcm.time = epoch2time(timeref);*/

	ReadRinex(infiles, 3, &obss, &navs, &nepoch);
	rtk.opt = 
	while (1)
	{
		rtkpos(&rtk, epochobs.data, epochobs.n, &nav);	
                  
	}

	fclose(ifp[0]);
	fclose(ifp[1]);
	fclose(ifp[2]);
	fclose(ofp);
	return 0;
}