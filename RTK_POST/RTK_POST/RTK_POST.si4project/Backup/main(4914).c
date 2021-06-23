#include "rtklib.h"
#include <process.h>
#include <windows.h>

#define RINEX 1
#define HTNEX 0
#define MAXSTRNUM 20000
#define MAXSPLITNUM 3000
#define MAXSTR 20
int iobs = 0;
FILE *ifp[3], *ofp[2];

extern int ReadRinex(char **files, int n, obs_t *obs, nav_t *nav, int *nepoch);
extern int addobsdata(obs_t *obs, const obsd_t *data);
extern int GetSys(int sysid);
extern int satno(int sys, int prn);

typedef struct split {
	char info[MAXSPLITNUM][MAXSTR];
	int num;
}SplitInfo; //字符串分割结构体

void GetSplitString(char* str,char* delim, SplitInfo* splitInfo) //字符串分割函数
{
	int n = 0;
	char *next = NULL;
	char lines[MAXSTRNUM];

	//strcpy(lines, str);
	next = strtok(str, delim);
	if (next == NULL) {
		splitInfo->num = 0;
		memset(splitInfo->info, 0, sizeof(char)*MAXSPLITNUM * MAXSTR);
		return;
	}	
	
	while (next!=NULL)
	{
		strcpy(splitInfo->info[n], next);	
		n++;
		next = strtok(NULL, delim);
		
	}
	splitInfo->num = n;
	
}
        
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
	int sys = 0;
	switch(sysid)
	{
	case 0:
		sys =  SYS_GPS;
		break;
	case 1:
		sys = SYS_GLO;
		break;
	case 2:
		sys = SYS_SBS;
		break;
	case 3:
		sys = SYS_GAL;
		break;
	case 4:
		sys = SYS_CMP;
		break;
	case 5:
		sys = SYS_QZS;
		break;
	case 6:
		sys = SYS_IRN;
		break;
	default:
		sys = -1;
	}
	return sys;
}
/*frame_data[sat][sat_data] 1:27:
周数	周内秒	跟踪状态	通道号	相位锁定	校验已知	伪码锁定	相关器类型	卫星系统	预留	分组标志	信号类型	预留	L1为首要通道
0		1		2			3		4			5			6			7			8			9		10			11			12		13
载波相位测量值+半周	滤波器指示	PRN锁定标志	通道分配	多普勒	伪距	ADR	伪距标准差	ADR标准差	卫星号	锁定时间	载噪比		GLONASS频带数
14					15			16			17			18		19		20	 21			22			23		24			25			26
参考：《UG016_数据通信接口协议_北云科技》
*/
int ObsRead(FILE *fp,obs_t *obss) //读取汇天自定义观测文件
{
	obsd_t obs = { 0 };
	char strline[MAXSTRNUM] = "";
	char delim[10] = "\t ";
	SplitInfo frame = { 0 };
	int i,j,k=0;
	double frame_data[110][27] = {0};
	int nsat = 0;
	char data[1024];
	char split_data[4000][20] = {0};

	obsd_t *obs_data;

	int row, col;
	int sys,prn,code,idx;
	double fre,lamb;
	double adr,adr_roll,max_value=8388608.0;

	if (!(obs_data=(obsd_t *)calloc(1,sizeof(obsd_t)))) return 0;
	while (!feof(fp))
	{
		fgets(strline, MAXSTRNUM, fp); //获得观测数据每行数据（每秒）
		//printf("\n%d", ++k);
		GetSplitString(strline, delim, &frame); //行数据分割为字符串

		if (frame.num < 3) continue; 
		for (i = 0; i < frame.num; i++) //存储行数据到变量frame_data
		{
			strcpy(data, frame.info[i]);
			if (i < 2) continue;
			row = i / 27;
			col = i % 27-2;
			frame_data[row][col] = atof(data);
		}
		nsat = row;
		for (i = 0; i < nsat; i++) //存储每秒数据到obss
		{
			obs_data->time = gpst2time(frame_data[i][0], frame_data[i][1]/1000.0);//时间
			sys = GetSys(frame_data[i][8]);//系统
			if (sys < 0) continue;
			prn = (int)frame_data[i][23];//prn
			if (sys==SYS_GLO) 
				prn  = sys-37;
			else if (sys == SYS_QZS) 
				prn  = sys-192;
			obs_data->sat = satno(sys,prn);	//卫星号
			obs_data->rcv = 1;//接收机号
			for (j = 0; j < obss->n; j++) //判断是否为同一时刻、同一卫星不同频率数据
			{
				if (obss->data[j].sat == obs_data->sat && obss->data[j].rcv == obs_data->rcv &&fabs(timediff(obss->data[j].time, obs_data->time)) < DTTOL) //if the same satellite,different freq
				{
					code = GetCode((int)frame_data[i][11], sys);//频率代码
					idx = code2idx(sys, code);//频率id
					obss->data[j].code[idx] = code;//频率代码
					obss->data[j].SNR[idx] = frame_data[i][25];//SNR
					obss->data[j].LLI[idx] = 0;
					obss->data[j].P[idx] = frame_data[i][19];//伪距
					obss->data[j].D[idx] = frame_data[i][18];//多普勒
					fre = sat2freq(obs_data->sat, obs_data->code[idx], NULL);
					lamb = CLIGHT / fre;
					adr = frame_data[i][20];
					adr_roll = (obs_data->P[idx] / lamb + adr) / max_value;
					adr_roll = ROUND(adr_roll);
					obss->data[j].L[idx] = adr - max_value*adr_roll;//载波相位（求解参考《UG016_数据通信接口协议_北云科技》）
					break;
				}
			}
			if (j == obss->n) //同一时刻、新的卫星数据
			{
				code = GetCode((int)frame_data[i][11], sys);
				idx = code2idx(sys, code);
				obs_data->code[idx] = code;
				obs_data->SNR[idx] = frame_data[i][25];
				obs_data->LLI[idx] = 0;
				obs_data->P[idx] = frame_data[i][19];
				obs_data->D[idx] = frame_data[i][18];
				fre = sat2freq(obs_data->sat, obs_data->code[idx], NULL);
				lamb = CLIGHT / fre;
				adr = frame_data[i][20];
				adr_roll = (obs_data->P[idx] / lamb + adr) / max_value;
				adr_roll = ROUND(adr_roll);
				obs_data->L[idx] = adr - max_value*adr_roll;
				addobsdata(obss, obs_data);
			}	
		}
	}
	sortobs(obss);
	free(obs_data);
	return 1;
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
		if(ret ==-2) return 1; 
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
	char strline[1000] = "";
	char delim[10] = "\t ";
	SplitInfo frame = { 0 };

	while (!feof(fp))
	{
		fgets(strline, 1000, fp); //获得观测数据每行数据（每秒）
		GetSplitString(strline, delim, &frame); //行数据分割为字符串

	}
	return 1;
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
		"E:\\Project\\data\\20200525\\range_data(3).txt",//"E:\\Project\\data\\20200525\\range_data(3).txt",
		"E:\\Project\\data\\20200525\\ephem_data(2).txt" };
	char *outfile = { "E:\\Project\\data\\20200525\\result.pos" };
	double ep[6] = { 2021,5,25,0,0,0 };


	obs_t obss = { 0 };
	nav_t nav = { 0 };
	int i;
	gtime_t gtime;

	if ((ifp[0] = fopen(infiles[0], "rb")) == NULL) return 0;
	if ((ifp[1] = fopen(infiles[1], "rb")) == NULL) return 0;
	if ((ifp[2] = fopen(infiles[2], "rb")) == NULL) return 0;
	if ((ofp[0] = fopen(outfile, "w")) == NULL) return 0;

	if (!NavUpdate(ifp[2], &nav))return 0;
	if (!RtcmRead(ifp[0], &obss, ep)) return 0;
	if (!ObsRead(ifp[1], &obss)) return 0;

	//for (i = 0; i < obss.n; i++)
	//{
	//	time2epoch(obss.data[i].time, ep);
	//	fprintf(ofp[0], "\n RTK %4d%02d%02d %02d %02d %05.2f ", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], ep[5]);
	//	fprintf(ofp[0], "%d %d %.2f %.2f %.2f", obss.data[i].sat, obss.data[i].rcv, obss.data[i].P[0], obss.data[i].P[1], obss.data[i].P[2]);
	//}
	

	fclose(infiles[0]);
	fclose(infiles[1]);
	fclose(infiles[2]);


#endif
	//fclose(ofp[0]);
	return 0;
}