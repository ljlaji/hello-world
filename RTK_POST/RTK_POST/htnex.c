#include "rtklib.h"

#define MAXSTRNUM 20000
#define MAXSPLITNUM 3000
#define MAXSTR 20

extern int ReadRinex(char **files, int n, obs_t *obs, nav_t *nav, int *nepoch);
extern int add_eph(nav_t *nav, const eph_t *eph);
extern int addobsdata(obs_t *obs, const obsd_t *data);
extern int GetSys(int sysid);
extern int satno(int sys, int prn);
extern int uraindex(double value);

typedef struct split {
	char info[MAXSPLITNUM][MAXSTR];
	int num;
}SplitInfo; //字符串分割结构体

void GetSplitString(char* str, char* delim, SplitInfo* splitInfo) //字符串分割函数
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

	while (next != NULL)
	{
		strcpy(splitInfo->info[n], next);
		n++;
		next = strtok(NULL, delim);

	}
	splitInfo->num = n;

}
/*
""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E",   0- 9 
"1A", "1B", "1X", "1Z", "2C", "2D", "2S", "2L", "2X", "2P",  10-19 
"2W", "2Y", "2M", "2N", "5I", "5Q", "5X", "7I", "7Q", "7X",  20-29 
"6A", "6B", "6C", "6X", "6Z", "6S", "6L", "8L", "8Q", "8X",  30-39 
"2I", "2Q", "6I", "6Q", "3I", "3Q", "3X", "1I", "1Q", "5A",  40-49 
"5B", "5C", "9A", "9B", "9C", "9X", "1D", "5D", "5P", "5Z",	 50-59
"6E", "7D", "7P", "7Z", "8D", "8P", "4A", "4B", "4X", ""     60-69
 */
int GetCode(int codeid, int sys)
{
	int code = 0;
	switch (sys)
	{
	case SYS_GPS:
		if (codeid == 0) code = 1; //L1C/A
		//else if (codeid == 5) code = 19;//L2P
		//else if (codeid == 9) code = 19;//L2P加密
		//else if (codeid == 14) code = 25;//L5Q
		//else if (codeid == 16) code = 1;//L1C
		else if (codeid == 17) code = 14;//L2C
		break;
	case SYS_GLO:
		if (codeid == 0) code = 1;
		else if (codeid == 1) code = 14;
		else if (codeid == 5) code = 19;
		break;
	case SYS_SBS:
		if (codeid == 0) code = 1;
		else if (codeid == 6) code = 14;
		else if (codeid == 5) code = 24;
		break;
	case SYS_GAL:
		if (codeid == 2) code = 1;
		else if (codeid == 6) code = 31;
		else if (codeid == 7) code = 32;
		else if (codeid == 12) code = 49;
		else if (codeid == 17) code = 50;
		else if (codeid == 20) code = 50;
		break;
	case SYS_CMP:
		if (codeid == 0) code = 47; //B1D1
		else if (codeid == 1) code = 27;//B2D1(C7I)
	//	else if (codeid == 2) code = 44;//B3D1
		else if (codeid == 4) code = 40;//B1D2(C2I)
		else if (codeid == 5) code = 27;//B2D2(C7I)
	//	else if (codeid == 6) code = 44;//B3D2
	//	else if (codeid == 7) code = 1;//B1C
	//	else if (codeid == 9) code = 40;//B2a
	//	else if (codeid == 10) code = 44;//B2b
		break;
	case SYS_QZS:
		if (codeid == 0) code = 1;
		else if (codeid == 14) code = 25;
		else if (codeid == 16) code = 1;
		else if (codeid == 17) code = 14;
	case SYS_IRN:
		if (codeid == 0) code = 24;
		break;
	}
	return code;
}
int GetSys(int sysid)
{
	int sys = 0;
	switch (sysid)
	{
	case 0:
		sys = SYS_GPS;
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
int ObsRead(FILE *fp, obs_t *obss) //读取汇天自定义观测文件
{
	obsd_t obs = { 0 };
	char strline[MAXSTRNUM] = "";
	char delim[10] = "\t ";
	SplitInfo frame = { 0 };
	int i, j, k = 0;
	double frame_data[110][27] = { 0 };
	int nsat = 0;
	char data[1024];
	char split_data[4000][20] = { 0 };

	obsd_t *obs_data;

	int row, col;
	int sys, prn, code, idx;
	double fre, lamb;
	double adr, adr_roll, max_value = 8388608.0;

	if (!(obs_data = (obsd_t *)calloc(1, sizeof(obsd_t)))) return 0;
	printf("\n obs data reading...\n");
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
			col = i % 27 - 2;
			frame_data[row][col] = atof(data);
		}
		nsat = row;
		for (i = 0; i < nsat; i++) //存储每秒数据到obs
		{
			memset(obs_data, 0, sizeof(obsd_t));
			obs_data->time = gpst2time(frame_data[i][0], frame_data[i][1] / 1000.0);//时间


			sys = GetSys(frame_data[i][8]);//系统
			if (!(sys == SYS_GPS || sys == SYS_CMP)) continue;

			prn = (int)frame_data[i][23];//prn
			if (sys == SYS_GLO)
				prn = sys - 37;
			else if (sys == SYS_QZS)
				prn = sys - 192;
			if (!(obs_data->sat = satno(sys, prn))) continue;	//卫星号
			
			if (obs_data->sat == 34)
				obs_data->sat = obs_data->sat;
			obs_data->rcv = 1;//接收机号


			for (j = 0; j < obss->n; j++) //判断是否为同一时刻、同一卫星不同频率数据
			{
				if (obss->data[j].sat == obs_data->sat && obss->data[j].rcv == obs_data->rcv &&fabs(timediff(obss->data[j].time, obs_data->time)) < DTTOL) //if the same satellite,different freq
				{
					if (!(code = GetCode((int)frame_data[i][11], sys))) break;//频率代码
					idx = code2idx(sys, code);//频率id
					obss->data[j].code[idx] = code;//频率代码
					obss->data[j].SNR[idx] = frame_data[i][25]/SNR_UNIT;//SNR
					obss->data[j].LLI[idx] = 0;
					obss->data[j].P[idx] = frame_data[i][19];//伪距
					obss->data[j].D[idx] = frame_data[i][18];//多普勒
					fre = sat2freq(obss->data[j].sat, obss->data[j].code[idx], NULL);
					lamb = CLIGHT / fre;
					adr = frame_data[i][20];
					adr_roll = (obss->data[j].P[idx] / lamb + adr) / max_value;
					adr_roll = ROUND(adr_roll);
					obss->data[j].L[idx] = -(adr - max_value*adr_roll);//载波相位（求解参考《UG016_数据通信接口协议_北云科技》）
					break;
				}
			}
			if (j == obss->n) //同一时刻、新的卫星数据
			{
				//code = GetCode((int)frame_data[i][11], sys);
				if (!(code = GetCode((int)frame_data[i][11], sys))) continue;//频率代码
				idx = code2idx(sys, code);
				obs_data->code[idx] = code;
				obs_data->SNR[idx] = frame_data[i][25] / SNR_UNIT;
				obs_data->LLI[idx] = 0;
				obs_data->P[idx] = frame_data[i][19];
				obs_data->D[idx] = frame_data[i][18];
				fre = sat2freq(obs_data->sat, obs_data->code[idx], NULL);
				lamb = CLIGHT / fre;
				adr = frame_data[i][20];
				adr_roll = (obs_data->P[idx] / lamb + adr) / max_value;
				adr_roll = ROUND(adr_roll);
				obs_data->L[idx] = -(adr - max_value*adr_roll);
				addobsdata(obss, obs_data);
			}
			if (fabs((double)obs_data->time.time - 1621946462) < 0.001 && obs_data->sat == 2)
				obs_data->time.time = obs_data->time.time;
		}
	}
	sortobs(obss);
	free(obs_data);
	return 1;
}
int RtcmRead(FILE *fp, obs_t *obss, double *ep)
{
	int i;
	int ret;
	//	double ep[6] = { 2021, 5, 25, 0, 0 };
	rtcm_t rtcm = { 0 };
	init_rtcm(&rtcm);
	rtcm.time = epoch2time(ep);
	printf("\n rtcm data reading...");
	while (1)
	{
		ret = input_rtcm3f(&rtcm, fp);
		if (ret == -2) return 1;
		else if (ret == 1)
		{
			for (i = 0; i<rtcm.obs.n; i++)
			{
				rtcm.obs.data[i].rcv = 2;
				addobsdata(obss, rtcm.obs.data + i);
			}
		}
	}

}
/*
GPS
0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32
GPSEPHEM PRN	TOW	Health	IODE1	IODE2	WN	Z WN	Toe	A	dn	M0	ecc	w	cuc	cus	crc	crs	cic	cis	i0	IDOT	Omega0	Omega	IODC	Toc	Tgd	a0	a1	a2	AS	N	URA
BDS

0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28
GPSEPHEM PRN	WN	URA	Health	tgd1	tgd2	AODC	toc	a0	a1	a2	AODE	toe	RootA	ecc	w	dn	M0	Omega0	Omega	i0	IDOT	cuc	cus	crc	crs	cic	cis

*/
int NavUpdate(FILE *fp, nav_t *nav)
{
	char strline[1000] = "";
	char delim[10] = "\t ";
	SplitInfo frame = { 0 };
	eph_t eph;
	int n = 0;
	int sys, prn;
	double toe, toc, ttr;
	printf("\n nav data reading...");
	while (!feof(fp))
	{
		fgets(strline, 1000, fp); //获得观测数据每行数据（每秒）
		GetSplitString(strline, delim, &frame); //行数据分割为字符串
		memset(&eph, 0, sizeof(eph_t));
		if (!strcmp(frame.info[0], "GPSEPHEM"))
		{
			
			sys = SYS_GPS;
			prn = atof(frame.info[1]);
			eph.sat = satno(sys, prn);
			eph.iode = atof(frame.info[4]);
			eph.iodc = atof(frame.info[24]);
			eph.sva = uraindex(atof(frame.info[32]));
			eph.svh = atof(frame.info[3]);
			eph.week = atof(frame.info[6]);
			eph.code = 0;
			eph.flag = 0;
			toe = atof(frame.info[8]);
			toc = atof(frame.info[25]);
			ttr = atof(frame.info[2]);
			eph.toe = gpst2time(eph.week, toe);
			eph.toc = gpst2time(eph.week, toc);
			eph.ttr = gpst2time(eph.week, ttr);
			eph.A = atof(frame.info[9]);
			eph.e = atof(frame.info[12]);
			eph.i0 = atof(frame.info[20]);
			eph.OMG0 = atof(frame.info[22]);
			eph.omg = atof(frame.info[13]);
			eph.M0 = atof(frame.info[11]);
			eph.deln = atof(frame.info[10]);
			eph.OMGd = atof(frame.info[23]);
			eph.idot = atof(frame.info[21]);
			eph.crc = atof(frame.info[16]);
			eph.crs = atof(frame.info[17]);
			eph.cuc = atof(frame.info[14]);
			eph.cus = atof(frame.info[15]);
			eph.cic = atof(frame.info[18]);
			eph.cis = atof(frame.info[19]);
			eph.toes = toe;
			eph.fit = 0;
			eph.f0 = atof(frame.info[27]);
			eph.f1 = atof(frame.info[28]);
			eph.f2 = atof(frame.info[29]);
			eph.tgd[0] = atof(frame.info[26]);
			eph.Adot = 0;
			eph.ndot = 0;
			add_eph(nav, &eph);
		}
		else if (!strcmp(frame.info[0], "BDSEPHEM"))
		{
			sys = SYS_CMP;
			prn = atof(frame.info[1]);
			eph.sat = satno(sys, prn);
			eph.iode = atof(frame.info[12]);
			eph.iodc = atof(frame.info[7]);
			eph.sva = uraindex(atof(frame.info[3]));
			eph.svh = atof(frame.info[4]);
			eph.week = atof(frame.info[2]);
			eph.code = 0;
			eph.flag = 0;
			toe = atof(frame.info[13]);
			toc = atof(frame.info[8]);
			ttr = atof(frame.info[8]);
			eph.toe = bdt2gpst(bdt2time(eph.week, toe));
			eph.toc = bdt2gpst(bdt2time(eph.week, toc));
			eph.ttr = bdt2gpst(bdt2time(eph.week, ttr));
			eph.A = atof(frame.info[14]);
			eph.A = eph.A*eph.A;
			eph.e = atof(frame.info[15]);
			eph.i0 = atof(frame.info[21]);
			eph.OMG0 = atof(frame.info[19]);
			eph.omg = atof(frame.info[16]);
			eph.M0 = atof(frame.info[18]);
			eph.deln = atof(frame.info[17]);
			eph.OMGd = atof(frame.info[20]);
			eph.idot = atof(frame.info[22]);
			eph.crc = atof(frame.info[25]);
			eph.crs = atof(frame.info[26]);
			eph.cuc = atof(frame.info[23]);
			eph.cus = atof(frame.info[24]);
			eph.cic = atof(frame.info[27]);
			eph.cis = atof(frame.info[28]);
			eph.toes = toe;
			eph.fit = 0;
			eph.f0 = atof(frame.info[9]);
			eph.f1 = atof(frame.info[10]);
			eph.f2 = atof(frame.info[11]);
			eph.tgd[0] = atof(frame.info[5]);
			eph.tgd[1] = atof(frame.info[6]);
			eph.Adot = 0;
			eph.ndot = 0;
			add_eph(nav, &eph);
		}



	}
	uniqnav(nav);
	return 1;
}
