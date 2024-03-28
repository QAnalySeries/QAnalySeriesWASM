#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"

#include "correlation.h"
#include "sedrate.h"
#include "tuning.h"

#include <QtCore/qtimer.h>

#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QSqlDriver>
#include <QTemporaryFile>
#include <QStandardPaths>


#include <QtNetwork/QNetworkRequest>
#include <QtNetwork/QNetworkAccessManager>
#include <QtNetwork/QNetworkReply>
#include <QtNetwork/QtNetwork>
//#include "qhtml5file.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void recalculateData();
    void replotAll();

private:
    Ui::MainWindow *ui;
    QCustomPlot *correlation;
    QCustomPlot *sedrate;
    QCustomPlot *tuning;

private slots:
    void mousePress();
    void mouseWheel();
    void onSedrateDrop();

    void reverseDataYButton_clicked();
    void reverseDataXButton_clicked();
    void reverseTargetYButton_clicked();
    void reverseTargetXButton_clicked();

    void on_actionOpen_Example_Dataset_triggered();
    void actionOpen_Laskar_triggered(); //what does this do again?
    void actionETP_Target_triggered();
    void actionInsolation_La04_triggered();
    void on_actionSSA_triggered();
    void actionData_Processor_triggered();
    void actionSpectral_Analysis_triggered();
    void actionFetchAstroData_triggered();

    
    void actionAbout_triggered();

    void actionSave_Project_as_triggered();
    void actionOpen_Project_triggered();

    void actionOpen_Data_triggered();
    //void actionOpen_Data2_triggered();

    void actionSave_new_Data_as_triggered();
    void actionSave_Data_with_ages_triggered();

    void actionOpen_Target_triggered();
    void actionSave_new_Target_as_triggered();

    void actionOpen_TimeModel_triggered();
    void actionSave_TimeModel_as_triggered();

    void actionBackToOriginal_clicked();




    void setCurrentIndex();
    void toggleMarkersCheckbox_stateChanged(int);
    void toggleCorrelationCheckbox_stateChanged(int);
    void showPointToolTip(QMouseEvent *event);


    void initLaskarDB();

public slots:
    void dataLoadCallback();
    void targetLoadCallback();
    void tmLoadCallback();

    void etpCallback();
    void laskarCallback();
    void insolationLa04Callback();
    void openDataCallback(const QString &fileName, const QByteArray &fileContentArray, const MainWindow *window);
    void dataProcessingCallback();
    void ssaProcessingCallback();
    void spectralProcessingCallback();
    /*
//    void on_pushButton_clicked();
//    void on_pushButton_4_clicked();
//    void on_pushButton_3_clicked();


    void on_actionExit_triggered();
    void on_pushButton_5_clicked();
    void on_pushButton_6_clicked();
    void on_pushButton_7_clicked();
    void on_pushButton_8_clicked();
    void on_checkBox_2_stateChanged(int arg1);*/
};

void readFile(QWidget *, MainWindow *);
#endif // MAINWINDOW_H
