#ifndef GLOBALS_H
#define GLOBALS_H
#include<QVector>
#include "QtSql/QSqlDatabase"
#include "QSqlQuery"
#include "QDebug"
#include "QSqlQueryModel"

#include "qcustomplot.h"

#define NETWORK_LOC "https://paloz.marum.de/QAnalySeriesWASM/laskar.db3"

extern QVector<double> x_d, y_d, x_d1, x_t, y_t, x_tm, y_tm, y_dn, y_tn; //dn and tn - for normalized graphs;
extern QVector<int> xind_tm, yind_tm, xind_tm_orig, yind_tm_orig;
extern int n_d, n_t, n_tm, ind0, ind1, ind2, screenw, screenh;
extern double psd_min, psd_max;
extern QString fileName,fType;
extern QString fileContent;
extern QSqlDatabase db;
//extern QSqlQuery query;
extern QTemporaryFile tmpFile;
extern QString laskarFileName;
extern QString exampleProjectName;
extern QString astroFileURL;

struct DataStr
{
    QVector<double> x;
    QVector<double> y;

} ;
extern struct DataStr data_orig;
extern struct DataStr target_orig;
extern struct DataStr tm_orig;

struct Orbit
{
    double time;
    double ecc;
    double eps;
    double varpi;
};

typedef		int			CommandT;
typedef		int			PaneIDT;

const CommandT	cmd_None					    = 9999;
const CommandT	cmd_OpenData					= 10000;
const CommandT	cmd_InsertGroup					= 10001;
const CommandT	cmd_SaveData					= 10002;
const CommandT	cmd_Info						= 10003;
const CommandT	cmd_Draw						= 10004;

const CommandT	cmd_DaylyInso					= 10010;
const CommandT	cmd_MeanInso					= 10011;
const CommandT	cmd_Eccentricity				= 10012;
const CommandT	cmd_Obliquity					= 10013;
const CommandT	cmd_Precession					= 10014;
const CommandT	cmd_MoreInso					= 10015;
const CommandT	cmd_PrecessionAngle				= 10016;
const CommandT	cmd_Noise						= 10017;

const CommandT	cmd_Correlations				= 10020;
const CommandT	cmd_Sampling					= 10021;
const CommandT	cmd_Smoothing					= 10022;
//consCommandTnt	cmd_Run							= 10023;
const CommandT	cmd_SSA							= 10023;
const CommandT	cmd_Fitting						= 10024;
const CommandT	cmd_Filtering					= 10025;
const CommandT	cmd_Princ_Compon				= 10026;
const CommandT	cmd_Periodogram					= 10027;
const CommandT	cmd_BTukey						= 10028;
const CommandT	cmd_MaxEntropy					= 10029;
const CommandT	cmd_MTM							= 10030;
const CommandT	cmd_Stats						= 10031;

const CommandT	cmd_Calder						= 10035;
const CommandT	cmd_Imbrie						= 10036;
const CommandT	cmd_Paillard_1998				= 10037;
const CommandT	cmd_Paillard_2004				= 10038;

const CommandT	cmd_Linage						= 10040;
const CommandT	cmd_LinageNext					= 10041;
const CommandT	cmd_LinageRevert				= 10042;
const CommandT	cmd_LinageSwitch				= 10043;
const CommandT	cmd_FindAscendant				= 10044;
const CommandT	cmd_Splinage					= 10045;
const CommandT	cmd_AgeScale					= 10046;
const CommandT	cmd_Combine						= 10047;
const CommandT	cmd_AddFiltered					= 10048;

const CommandT	cmd_LinageSwitch_first			= 10100;
const CommandT	cmd_LinageSwitch_last			= 10199;		//	at most 100 compPairs !
const CommandT	cmd_work_Windows_first			= 10200;
const CommandT	cmd_work_Windows_last			= 10299;		//	at most 100 windows !
const CommandT	cmd_draw_Windows_first			= 10300;
const CommandT	cmd_draw_Windows_last			= 10399;		//	at most 100 windows !

const PaneIDT	k_squareWindow		= 1;
const PaneIDT	k_BartlettWindow	= 2;
const PaneIDT	k_TukeyWindow		= 3;
const PaneIDT	k_ParzenWindow		= 4;
const PaneIDT	k_WelchWindow		= 5;

#endif // GLOBALS_H
