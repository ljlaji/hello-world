#include "rtklib.h"
#include <process.h>
#include <windows.h>


rtcm_t *rtcm;
FILE *ifp, *ofp;

void RtcmCollect()
{
	while (1)
	{
		input_rtcm3f(rtcm, ifp);
		if (rtcm->)
	}
}



int main()
{

	double ep[6];
	int i, flag;
	char *RtcmFile = { "raw1_2011230232.dat" };
	char *OutFile = { "result.txt" };
	double calendar_ref[6] = { 2020, 11, 23, 0, 0, 0 };//²Î¿¼Ê±¼ä

	if ((ifp = fopen(RtcmFile, "rb")) == NULL)
		return 0;
	if ((ofp = fopen(OutFile, "w")) == NULL)
		return 0; 

	if (!(rtcm = (rtcm_t *)calloc(1,sizeof(rtcm_t))))
		return 0;
	init_rtcm(rtcm);
	rtcm->time = epoch2time(calendar_ref);

	_beginthread(RtcmCollect, 0, NULL);
	
	while (1)
	{
		printf("\n22222222222222222222222222");
		Sleep(2);
#if 0		
		flag = input_rtcm3f(rtcm, ifp);
		if (flag == -2)
			break;
		fprintf(ofp, "\n%s",rtcm->msgtype);
		if(rtcm->nmsg3[19])
		{
			printf("\nrtcm->nmsg3[19]:%d", rtcm->nmsg3[19]);
			fprintf(ofp,"\nrtcm->nmsg3[19]:%d", rtcm->nmsg3[19]);
		}
		if (rtcm->obsflag)
		{
			for (i = 0; i < rtcm->obs.n; i++)
			{
				fprintf(ofp, "\n%2d, %.2f %.2f %.2f, %.2f %.2f %.2f", rtcm->obs.data[i].sat, rtcm->obs.data[i].P[0], rtcm->obs.data[i].P[1], rtcm->obs.data[i].P[2],
					rtcm->obs.data[i].L[0], rtcm->obs.data[i].L[1], rtcm->obs.data[i].L[2]);
			}
		}
		time2epoch(rtcm->time, ep);
	//	printf("\ntime: %.2f %.2f %.2f %.2f %.2f %.2f", ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
#endif
	}

	fclose(ifp);
	fclose(ofp);
	return 0;
}