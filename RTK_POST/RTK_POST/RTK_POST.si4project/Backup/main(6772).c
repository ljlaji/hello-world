#include "rtklib.h"
#include <process.h>
#include <windows.h>

#define RINEX 1
#define HTNEX 0
#define MAXSTRNUM 20000
int iobs = 0;
FILE *ifp[3], *ofp[2];

extern int ReadRinex(char **files, int n, obs_t *obs, nav_t *nav, int *nepoch);
extern int addobsdata(obs_t *obs, const obsd_t *data);

typedef struct split {
	char info[1024];
	int num;
	struct split *next;
}SplitInfo;

void GetSplitString(char* str,char* delim, SplitInfo* splitInfo) //string split
{
	int n = 0;
	char* next = NULL;
	SplitInfo* top = splitInfo;

	char arr[1024];
	//strcpy(arr, str);

	next = strtok(str, delim);
	if (next == NULL) {
		splitInfo->num = 0;
		splitInfo->next = NULL;

		return;
	}

	SplitInfo* in = (SplitInfo*)malloc(sizeof(SplitInfo));
	strcpy(in->info, next);
	in->next = NULL;
	top->next = in;
	top = in;
	n++;

	while (next = strtok(NULL, delim))
	{
		SplitInfo* te = (SplitInfo*)malloc(sizeof(SplitInfo));
		strcpy(te->info, next);
		te->next = NULL;
		top->next = te;
		top = te;
		n++;
	}

	splitInfo->num = n;
}

void SplitFree(SplitInfo* splitInfo) {
	SplitInfo* top = splitInfo->next;
	while (top) {
		SplitInfo* next = top->next;
		free(top);
		top = next;
	}

	splitInfo->next = NULL;
	splitInfo->num = 0;
}
    
    //""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
    //"1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
    //"2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
    //"6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
    //"2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
    //"5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
    //"6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
    
int GetCode(int codeid,int sys)
{
	int code=-1;
	switch (sys)
	{
		case SYS_GPS:
			if(codeid == 0) code = 1;
			else if(codeid==5) code = 9;
			else if (codeid==9) code = 9;
			else if (codeid==14) code = 25;
			else if (codeid==16) code = 1;
			else if (codeid==17) code = 14;
			break;
		case SYS_GLO:	
			if(codeid == 0) code = 1;
			else if(codeid==1) code = 14;
			else if (codeid==5) code = 19;
			break;			
		case SYS_SBS:	
			if(codeid == 0) code = 1;
			else if(codeid==6) code = 14;
			else if (codeid==5) code = 24;
			break;		
		case SYS_GAL:	
			if(codeid == 2) code = 1;
			else if(codeid==6) code = 31;
			else if (codeid==7) code = 32;
			else if(codeid==12) code = 49;
			else if (codeid==17) code = 50;
			else if(codeid==20) code = 50;
			break;
		case SYS_CMP:	
			if(codeid == 0) code = 47;
			else if(codeid==1) code = 40;
			else if (codeid==2) code = 44;
			else if(codeid==4) code = 47;
			else if (codeid==5) code = 40;
			else if(codeid==6) code = 44;
			else if(codeid==7) code = 1;
			else if (codeid==9) code = 40;
			else if(codeid==10) code = 44;
			break;	
		case SYS_QZS:	
			if(codeid == 0) code = 1;
			else if(codeid==14) code = 25;
			else if (codeid==16) code = 1;
			else if(codeid==17) code = 14;
		case SYS_IRN:	
			if(codeid == 0) code = 24;
			break;
	}
	return code;
}
int GetSys(int sysid)
{
	switch(sysid)
	{
	case 0:
		return SYS_GPS;
	case 1:
		return SYS_GLO;
	case 2:
		return SYS_SBS;
	case 3:
		return SYS_GAL;
	case 4:
		return SYS_CMP;
	case 5:
		return SYS_QZS;
	case 6:
		return SYS_IRN;
	default:
		return -1;
	}
}
int ObsRead(FILE *fp,obs_t *obss) //read obs file
{
	obsd_t obs = { 0 };
	char strline[MAXSTRNUM] = "";
	char delim[10] = "\t";
	SplitInfo frame = {0},*split;
	
	int i,j;
	double frame_data[110][27] = {0};
	int nsat = 0;
	char data[1024];

	obsd_t *obs_data;

	int row, col;
	int sys,prn,fre,code,idx;

	if (!(obs_data=(obsd_t *)calloc(1,sizeof(obsd_t)))) return 0;
	while (!feof(fp))
	{
		fgets(strline, MAXSTRNUM, fp);

		GetSplitString(strline, delim, &frame); //get line string
		if (frame.num < 3) continue;
		split = &frame;
		for (i = 0; i < frame.num; i++) //load data by lines
		{
			split = split->next;
			memcpy(data, split->info, sizeof(split->info));
			if (i < 2) continue;
			row = i / 27;
			col = i % 27-2;
			frame_data[row][col] = atof(data);
		}
		nsat = row;
		for (i = 0; i < nsat; i++) //load date by satellites
		{
			obs_data->time = gpst2time(frame_data[i][0], frame_data[i][1]/1000.0);//time			
			if(sys = GetSys((int)frame_data[i][8])<0) return 0;
			prn = (int)frame_data[i][23];
			if (sys==SYS_GLO) 
				prn  = sys-37;
			else if (sys == SYS_QZS) 
				prn  = sys-192;
			obs_data->sat = satno(sys,prn);		
			obs_data->rcv = 1;
			code = GetCode((int)frame_data[i][11],sys);
			idx = code2idx(sys,code);
			obs_data->code[idx]=code;
			obs_data->SNR[idx] = frame_data[i][25];
			//obs_data->LLI[idx] = 
			




			
			for (j = 0; j < 27; j++)
				fprintf(ofp[0], "%.2f ", frame_data[i][j]);
			fprintf(ofp[0], "\n");

		}
		fclose(ofp[0]);
		obss->n = obss + 1;

	}

}
int RtcmRead(FILE *fp, obs_t *obss,double *ep)
{
	int i;
	int ret;
//	double ep[6] = { 2021, 5, 25, 0, 0 };
	rtcm_t rtcm = { 0 };
	init_rtcm(&rtcm);
	rtcm.time = epoch2time(ep);
	while(1)
	{
		ret = input_rtcm3f(&rtcm, fp);
		if(ret ==-2) return 0; 
		else if(ret == 1)
		{
			for(i=0;i<rtcm.obs.n;i++) 
			{
				rtcm.obs.data[i].rcv = 2;	
				addobsdata(obss, rtcm.obs.data+i);
			}			
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
	
#endif

#if HTNEX
	char *infiles[3] = { "E:\\Project\\data\\20200525\\ReceivedTofile-TCPCLIENT-2021_5_25_20-38-52.DAT",
		"E:\\Project\\data\\20200525\\range_data(3).txt",
		"E:\\Project\\data\\20200525\\ephem_data" };
	char *outfile = { "E:\\Project\\data\\20200525\\result.pos" };
	double ep[6] = { 2021,5,25,0,0,0 };


	obs_t obss = { 0 };
	int i;
	gtime_t gtime;

	if ((ifp[0] = fopen(infiles[0], "rb")) == NULL) return 0;
	if ((ifp[1] = fopen(infiles[1], "rb")) == NULL) return 0;
	if ((ifp[2] = fopen(infiles[2], "rb")) == NULL) return 0;
	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0;

	RtcmRead(ifp[0], &obss,ep);
	//for (i = 0; i < obss.n; i++)
	//{
	//	time2epoch(obss.data[i].time, ep);
	//	fprintf(ofp[0], "\n RTK %4d%02d%02d %02d %02d %05.2f ", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
	//	fprintf(ofp[0], "%d %d %.2f %.2f %.2f", obss.data[i].sat, obss.data[i].rcv, obss.data[i].P[0], obss.data[i].P[1], obss.data[i].P[2]);
	//}
	
	ObsRead(ifp[1], &obss);
	fclose(infiles[0]);
	fclose(infiles[1]);
	fclose(infiles[2]);


#endif
	//fclose(ofp[0]);
	return 0;
}