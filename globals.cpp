#include "globals.h"
#include "qcustomplot.h"
QVector<double> x_d, y_d, x_d1, x_t, y_t, x_tm, y_tm, y_dn, y_tn; //dn and tn - for normalized graphs;
QVector<int> xind_tm, yind_tm, xind_tm_orig, yind_tm_orig;
int n_d, n_t, n_tm, ind0, ind1, ind2, screenw, screenh;
double psd_min, psd_max;
QString fileName,fType;
QString fileContent;

struct DataStr data_orig;
struct DataStr target_orig;
struct DataStr tm_orig;

QSqlDatabase db;
QTemporaryFile tmpFile(qApp);

QString laskarFileName = "laskar.db3";
QString exampleProjectName = "exampleProject.qaproj";
QString astroFileURL = NETWORK_LOC;
    //"https://paloz.marum.de/QAnalySeriesWASM/laskar.db3";
