#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include "readcsv.h"
#include "globals.h"
#include "qcustomplot.h"

#include "spectraldialog.h"

#include "mynr.h"
#include "laskar.h"
#include "dataproc.h"
#include "ssa.h"
#include "etp.h"
#include "insolation.h"

#include <QDir>
#include <QDebug>
#include <QVectorIterator>


static QVector<QCPItemLine*> arrows; // initiates an array of pointers on arrows between data and target tie points
static QVector<QCPItemText*> depths; // initiates an array of pointers on depth text labels of tuned data
static QString fName = "qa_untitled.txt";
static QString prName = "untitled.qaprj";
static int command_line;


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{

    ui->setupUi(this);
    n_tm=0; n_d = 0; n_t = 0;
    x_tm.clear(); y_tm.clear();
    tm_orig.x.clear();tm_orig.y.clear();

    // connect reverse buttons

    connect(ui->actionAbout, SIGNAL(triggered()), this, SLOT(actionAbout_triggered()));

    connect(ui->ReverseTargetXButton, SIGNAL(clicked()), this, SLOT(reverseTargetXButton_clicked()));
    connect(ui->ReverseTargetYButton, SIGNAL(clicked()), this, SLOT(reverseTargetYButton_clicked()));
    connect(ui->ReverseDataXButton, SIGNAL(clicked()), this, SLOT(reverseDataXButton_clicked()));
    connect(ui->ReverseDataYButton, SIGNAL(clicked()), this, SLOT(reverseDataYButton_clicked()));
    connect(ui->BackToOriginal, SIGNAL(clicked()), this, SLOT(actionBackToOriginal_clicked()));
    connect(ui->actionOpen_Project, SIGNAL(triggered()), this, SLOT(actionOpen_Project_triggered()));
    connect(ui->actionSave_Project_as, SIGNAL(triggered()), this, SLOT(actionSave_Project_as_triggered()));

    connect(ui->actionOpen_Data, SIGNAL(triggered()), this, SLOT(actionOpen_Data_triggered()));
    connect(ui->actionSave_new_Data_as, SIGNAL(triggered()), this, SLOT(actionSave_new_Data_as_triggered()));
    connect(ui->actionSave_Data_with_ages, SIGNAL(triggered()), this, SLOT(actionSave_Data_with_ages_triggered()));

    connect(ui->actionOpen_TimeModel, SIGNAL(triggered()), this, SLOT(actionOpen_TimeModel_triggered()));
    connect(ui->actionSave_TimeModel_as, SIGNAL(triggered()), this, SLOT(actionSave_TimeModel_as_triggered()));


    connect(ui->actionOpen_Target, SIGNAL(triggered()), this, SLOT(actionOpen_Target_triggered()));
    connect(ui->actionSave_new_Target_as, SIGNAL(triggered()), this, SLOT(actionSave_new_Target_as_triggered()));

    connect(ui->actionFetch_Astro_Data, SIGNAL(triggered()), this, SLOT(actionFetchAstroData_triggered()));

    connect(ui->actionETP_Target, SIGNAL(triggered()), this, SLOT(actionETP_Target_triggered()));
    connect(ui->actionOpen_Laskar, SIGNAL(triggered()), this, SLOT(actionOpen_Laskar_triggered()));
    connect(ui->actionInsolation_La04, SIGNAL(triggered()), this, SLOT(actionInsolation_La04_triggered()));

    connect(ui->MarkersCheckbox, SIGNAL(stateChanged(int)), this, SLOT(toggleMarkersCheckbox_stateChanged(int)));
    connect(ui->CorrelationCheckbox, SIGNAL(stateChanged(int)), this, SLOT(toggleCorrelationCheckbox_stateChanged(int)));

    connect(ui->actionSpectral_Analysis, SIGNAL(triggered()), this, SLOT(actionSpectral_Analysis_triggered()));
    connect(ui->actionData_Processor, SIGNAL(triggered()), this, SLOT(actionData_Processor_triggered()));


    // setup correlation
    QList<QCPAxis *> axlist;
    axlist << ui->correlation->xAxis <<  ui->correlation->xAxis2 ;

    ui->correlation->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->correlation->axisRect()->setRangeDragAxes(axlist);
    ui->correlation->axisRect()->setRangeZoomAxes(axlist);
    ui->correlation->addGraph(ui->correlation->xAxis2, ui->correlation->yAxis2);
    ui->correlation->xAxis2->setLabel("Data, x");
    ui->correlation->yAxis2->setLabel("Data, y");
    ui->correlation->graph(0)->rescaleAxes();
    ui->correlation->graph(0)->setName("Data");
    ui->correlation->graph(0)->setSelectable(QCP::stSingleData);
    ui->correlation->xAxis2->setVisible(true);
    ui->correlation->yAxis2->setVisible(true);
    ui->correlation->legend->setVisible(true);
    ui->correlation->addGraph(ui->correlation->xAxis, ui->correlation->yAxis);
    ui->correlation->graph(1)->setPen(QPen(Qt::red));
    ui->correlation->xAxis->setLabel("Target, x");
    ui->correlation->yAxis->setLabel("Target, y");
    ui->correlation->graph(1)->rescaleAxes();
    ui->correlation->graph(1)->setName("Target");
    ui->correlation->graph(1)->setSelectable(QCP::stSingleData);

    if (ui->MarkersCheckbox->isChecked()) {
        ui->correlation->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::white, 5));
        ui->correlation->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::red, Qt::white, 5));
    }
    else {
        ui->correlation->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone));
        ui->correlation->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone));
    }

    ui->correlation->show();


    axlist.clear();
    axlist << ui->sedrate->xAxis;
    ui->sedrate->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->sedrate->axisRect()->setRangeDragAxes(axlist);
    ui->sedrate->axisRect()->setRangeZoomAxes(axlist);
    ui->sedrate->addGraph(ui->sedrate->xAxis, ui->sedrate->yAxis);
    ui->sedrate->xAxis->setLabel("Target, x");
    ui->sedrate->yAxis->setLabel("Sed. rate");
    ui->sedrate->yAxis2->setLabel(" ");

    ui->sedrate->graph(0)->rescaleAxes();
    ui->sedrate->graph(0)->setName("Sedimentation rate");
    ui->sedrate->xAxis->setVisible(true);
    ui->sedrate->show();
    axlist.clear();

    axlist << ui->tuning->xAxis <<  ui->tuning->xAxis2 ;
    ui->tuning->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iSelectItems);
    ui->tuning->axisRect()->setRangeDragAxes(axlist);
    ui->tuning->axisRect()->setRangeZoomAxes(axlist);
    ui->tuning->addGraph(ui->tuning->xAxis2, ui->tuning->yAxis2);
    ui->tuning->xAxis2->setLabel("Data, x");


    ui->tuning->yAxis2->setRange(-1.2, 1.2);
    ui->tuning->graph(0)->setName("Data");
    ui->tuning->graph(0)->setPen(QPen(Qt::blue));
    ui->tuning->graph(0)->setSelectable(QCP::stSingleData);
    ui->tuning->xAxis2->setVisible(true);
    ui->tuning->yAxis2->setVisible(true);
    ui->tuning->yAxis->setVisible(true);
    ui->tuning->legend->setVisible(true);
    ui->tuning->addGraph(ui->tuning->xAxis, ui->tuning->yAxis);
    ui->tuning->graph(1)->setPen(QPen(Qt::red));
    ui->tuning->xAxis->setLabel("Target, x");
    ui->tuning->yAxis->setLabel("Target, y");
    ui->tuning->yAxis2->setLabel("Data, y");

    ui->tuning->yAxis->setRange(-1.2, 1.2);
    ui->tuning->graph(1)->setName("Target");
    ui->tuning->graph(1)->setPen(QPen(Qt::red));
    ui->tuning->graph(1)->setSelectable(QCP::stSingleData);
    connect(ui->tuning, &QCustomPlot::mouseMove, this, &MainWindow::showPointToolTip);

    // Selection markers
    ui->tuning->addGraph(ui->tuning->xAxis2, ui->tuning->yAxis2);
    ui->tuning->graph(2)->setPen(QPen(Qt::black));
    ui->tuning->graph(2)->setLineStyle(QCPGraph::lsNone);
    ui->tuning->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCrossCircle, 10));
    ui->tuning->addGraph(ui->tuning->xAxis, ui->tuning->yAxis);
    ui->tuning->graph(3)->setPen(QPen(Qt::black));
    ui->tuning->graph(3)->setLineStyle(QCPGraph::lsNone);
    ui->tuning->graph(3)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCrossCircle, 10));
    ui->tuning->legend->setVisible(false);

    if (ui->MarkersCheckbox->isChecked()) {
        ui->tuning->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
        ui->tuning->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
    }
    else {
        ui->tuning->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone));
        ui->tuning->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone));
    }


    ui->tuning->show();

    connect(ui->tuning, &QCustomPlot::selectionChangedByUser, this, &MainWindow::mousePress);
    connect(ui->tuning, &QCustomPlot::mouseWheel, this, &MainWindow::mouseWheel);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void readFile(QWidget *parent, MainWindow *mainWindow) // fType: d (data), t (target), or tm (time model)
{
    Readcsv *readcsv = new Readcsv(mainWindow);
    QObject::connect(readcsv, SIGNAL(accepted()), readcsv, SLOT(buttonBox_accepted()));
    readcsv->setModal(true);
    readcsv->open();
}

void MainWindow::mousePress()
{
    QCPDataSelection selection0 = ui->tuning->graph(0)->selection();
    if (!selection0.isEmpty()) {
        ind0 = selection0.dataRange().begin();
        ui->tuning->xAxis2->setLabel(QString("Data, x: %1, y:  %2").arg(x_d[ind0]).arg(y_d[ind0]));
        ui->tuning->graph(2)->data()->clear();
        ui->tuning->graph(2)->addData(x_d1[ind0], y_dn[ind0]);

    }

    QCPDataSelection selection1 = ui->tuning->graph(1)->selection();
    if (!selection1.isEmpty()){
        ind1 = selection1.dataRange().begin();
        ui->tuning->xAxis->setLabel(QString("Target, x: %1, y:  %2").arg(x_t[ind1]).arg(y_t[ind1]));
        ui->tuning->graph(3)->data()->clear();
        ui->tuning->graph(3)->addData(x_t[ind1], y_tn[ind1]);

    }

    if(QApplication::keyboardModifiers() & Qt::ShiftModifier)  // recalculate data
    {
        int i = 0;
        while (i<n_tm && x_tm[i] < x_d[ind0]) ++i;
        if (i==n_tm || i==0  ||  ((x_t[ind1]-y_tm[i-1])*(x_t[ind1]-y_tm[i])<0.0 && (x_d[ind0]-x_tm[i-1])*(x_d[ind0]-x_tm[i])<0.0 )) // check for Time Model consistantcy
        {
            x_tm.insert(i, x_d[ind0]);
            y_tm.insert(i, x_t[ind1]);
            xind_tm.insert(i, ind0);
            yind_tm.insert(i, ind1);
            n_tm++;
            recalculateData();
            ui->tuning->graph(2)->data()->clear();
            ui->tuning->graph(3)->data()->clear();}

    }

    if(QApplication::keyboardModifiers() & Qt::AltModifier)
    {
        bool selection2=false;
        for (int i = 0; i < n_tm; i++) {
            if (arrows[i]->selected()) {
                ind2=i; selection2=true;
            }
        }

        if (selection2){
            x_tm.remove(ind2);
            y_tm.remove(ind2);
            xind_tm.remove(ind2);
            yind_tm.remove(ind2);
            n_tm--;
            if (n_tm==0) {
                arrows.clear();
                depths.clear();
                for (int i = ui->tuning->itemCount()-1; i>=0 ; i--) {
                    ui->tuning->removeItem(i); // clear arrows and depths
                }
            }
            if (n_tm>0) {
                recalculateData();
            }
        }
    }
}

void MainWindow::mouseWheel()
{
    QList<QCPAxis *> axlist;
    if(QApplication::keyboardModifiers() & Qt::ShiftModifier) {
        axlist << ui->tuning->xAxis <<  ui->tuning->xAxis2 << ui->tuning->yAxis <<  ui->tuning->yAxis2 ;}
    else {
        axlist.clear();
        axlist << ui->tuning->xAxis <<  ui->tuning->xAxis2 ;
    }

    ui->tuning->axisRect()->setRangeDragAxes(axlist);
    ui->tuning->axisRect()->setRangeZoomAxes(axlist);

}

void MainWindow::recalculateData()
{
    double r;

    if (n_d<2 || n_t<2) {QMessageBox msgBox;
        msgBox.setText("Download Data and Target first!");
        msgBox.show();
        return;}

    if (n_tm == 0) {x_d1=x_d; return;}

    if (x_tm.first()<x_d.first() || x_tm.last()>x_d.last() || y_tm.first()<x_t.first() || y_tm.last()>x_t.last() ||
            x_tm.first()>x_d.last() || x_tm.last()<x_d.first() || y_tm.first()>x_t.last() || y_tm.last()<x_t.first()) {
        QMessageBox msgBox;
        msgBox.setText("Time Modes is inconsistent with Data or Target");
        msgBox.setModal(true);
        msgBox.show();
        return;
    }


    arrows.clear();
    depths.clear();

    for (int i = ui->tuning->itemCount()-1; i>=0 ; i--) {
        ui->tuning->removeItem(i); // clear arrows and depths
    }

    if (n_tm==1) {
        arrows<< new QCPItemLine(ui->tuning);
        arrows[0]->start->setAxes(ui->tuning->xAxis2,ui->tuning->yAxis2);
        arrows[0]->start->setCoords(x_d1[xind_tm[0]], y_dn[xind_tm[0]]);
        arrows[0]->end->setCoords(x_t[yind_tm[0]], y_tn[yind_tm[0]]);
        arrows[0]->setPen(Qt::DashDotLine);
        x_d1=x_d;
        return;
    }

    x_d1.resize(x_d.length());

    mynr::linterpv(x_tm, y_tm, x_d, x_d1); // faster but x_d must be in increasing order...
    if (n_tm==2) {
        ui->correlation->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
        ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
        ui->tuning->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
        ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    }

    ui->tuning->graph(0)->data()->clear();
    ui->tuning->graph(0)->setData(x_d1, y_dn);
    ui->tuning->graph(1)->data()->clear();
    ui->tuning->graph(1)->setData(x_t, y_tn);

    ui->correlation->graph(0)->data()->clear();
    ui->correlation->graph(0)->setData(x_d1, y_d);
    ui->correlation->replot();

    ui->sedrate->graph(0)->data()->clear();
    for (int i=1; i<n_tm; i++) {
        ui->sedrate->graph(0)->addData(y_tm[i-1],(x_tm[i]-x_tm[i-1])/(y_tm[i]-y_tm[i-1]));
        ui->sedrate->graph(0)->addData(y_tm[i],(x_tm[i]-x_tm[i-1])/(y_tm[i]-y_tm[i-1]));
    }
    ui->sedrate->rescaleAxes();
    ui->sedrate->xAxis->setRange(ui->tuning->xAxis->range().lower, ui->tuning->xAxis->range().upper);
    ui->sedrate->replot();

    // plot arrows and depths labels
    for (int i=0; i < n_tm; i++) {//ui->tuning->graph(2)->addData(y_tm[i], 0.0);
        arrows<< new QCPItemLine(ui->tuning);
        arrows[i]->start->setCoords(x_d1[xind_tm[i]], y_dn[xind_tm[i]]);
        arrows[i]->end->setCoords(x_t[yind_tm[i]], y_tn[yind_tm[i]]);

        //  Depth
        depths<< new QCPItemText(ui->tuning);
        depths[i]->position->setCoords(x_t[yind_tm[i]], 1.15); // place position at center/top of axis rect
        depths[i]->setText(QString::number(x_d[xind_tm[i]],'g',3));

    }

    // Correlation
    QVector<double> y_t1;
    mynr::linterpv(x_t, y_t, x_d1, y_t1);
    r = mynr::corr(y_d, y_t1);
    if (ui->CorrelationCheckbox->isChecked()) {
        ui->correlationLabel->setText(QString::number(r,'g',3));
    } else {ui->correlationLabel->setText("   ");}
    ui->tuning->replot();
}

void MainWindow::replotAll()
{
    double minv, maxv;
    if (n_d > 0) {
        minv = *std::min_element(y_d.constBegin(), y_d.constEnd());
        maxv = *std::max_element(y_d.constBegin(), y_d.constEnd());
        y_dn=y_d;// normalized data to [0 1];
        for (int i=0; i<n_d; i++) {y_dn[i]=(y_d[i]-minv)/(maxv-minv)+0.1;};
        //normalized = (x-min(x))/(max(x)-min(x));
        ui->correlation->graph(0)->data()->clear();
        ui->correlation->graph(0)->setData(x_d1,y_d);
        ui->correlation->graph(0)->rescaleAxes();
        ui->correlation->setVisible(true);
        ui->correlation->replot();

        ui->tuning->graph(0)->data()->clear();
        ui->tuning->graph(0)->setData(x_d1,y_dn);

        ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
        ui->tuning->yAxis2->setRange(-1.2 , 1.2);
        ui->tuning->setVisible(true);
        ui->tuning->replot();

        ui->sedrate->rescaleAxes();
        ui->sedrate->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
        ui->sedrate->replot();
    }

    if (n_t > 0) {
        minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
        maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());
        y_tn=y_t;// normalized target to [-1 0];
        for (int i=0; i<n_t; i++) {y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;};
        ui->correlation->graph(1)->data()->clear();
        ui->correlation->graph(1)->setData(x_t,y_t);
        ui->correlation->graph(1)->rescaleAxes();
        ui->correlation->setVisible(true);
        ui->correlation->replot();

        ui->tuning->graph(1)->data()->clear();
        ui->tuning->graph(1)->setData(x_t,y_tn);
        ui->tuning->xAxis->setRange(x_t[0], x_t[n_t-1]);
        ui->tuning->yAxis->setRange(-1.2 , 1.2);
        ui->tuning->setVisible(true);
        ui->tuning->replot();
    }

    ui->sedrate->graph(0)->data()->clear();
    for (int i=1; i<n_tm; i++) {

        ui->sedrate->graph(0)->addData(y_tm[i-1],(x_tm[i]-x_tm[i-1])/(y_tm[i]-y_tm[i-1]));
        ui->sedrate->graph(0)->addData(y_tm[i],(x_tm[i]-x_tm[i-1])/(y_tm[i]-y_tm[i-1]));
    }

}

void MainWindow::onSedrateDrop()
{
    recalculateData();
    ui->tuning->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->replot();
    ui->correlation->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
    ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->correlation->replot();
}

void MainWindow::reverseTargetXButton_clicked()
{
    if (ui->reverseTargetXLabel->text() == "Normal") {
      ui->reverseTargetXLabel->setText("Reversed");
    } else {
      ui->reverseTargetXLabel->setText("Normal");
    }
    ui->reverseTargetXLabel->repaint();

    for (int i=0; i<n_t; ++i) {
        x_t[i]=x_t[i]*(-1.0);
    }
    for (int i=0; i<floor(n_t/2.0);++i) {
        std::swap(x_t[i],x_t[n_t-1-i]);
        std::swap(y_t[i],y_t[n_t-1-i]);
    }
    n_tm = 0;
    x_tm.clear(); y_tm.clear();
    depths.clear(); arrows.clear();

    replotAll();
}

void MainWindow::reverseTargetYButton_clicked()
{
    if (ui->reverseTargetYLabel->text() == "Normal") {
        ui->reverseTargetYLabel->setText("Reversed");
    } else {
         ui->reverseTargetYLabel->setText("Normal");
    }

    ui->reverseTargetYLabel->repaint();

    for (int i=0; i<n_t; ++i) y_t[i]=y_t[i]*(-1.0);
    replotAll();
    if (n_tm > 1) {
        recalculateData();
        ui->correlation->xAxis->setRange(x_t[0], x_t[n_t-1]);
        ui->correlation->xAxis2->setRange(x_t[0], x_t[n_t-1]);
        ui->correlation->replot();
        ui->tuning->xAxis->setRange(x_t[0], x_t[n_t-1]);
        ui->tuning->xAxis2->setRange(x_t[0], x_t[n_t-1]);
        ui->tuning->replot();
    }
}

void MainWindow::reverseDataXButton_clicked()
{
    if (ui->reverseDataXLabel->text() == "Normal") {
        ui->reverseDataXLabel->setText("Reversed");
    } else {
         ui->reverseDataXLabel->setText("Normal");
    }

    ui->reverseDataXLabel->repaint();

    for (int i=0; i<n_d; ++i) {
         x_d[i]=x_d[i]*(-1.0); x_d1[i]=x_d1[i]*(-1.0);
    }
    for (int i=0; i<floor(n_d/2.0);++i) {
        std::swap(x_d[i],x_d[n_d-1-i]);
        std::swap(x_d1[i],x_d1[n_d-1-i]);
        std::swap(y_d[i],y_d[n_d-1-i]);
    }
    n_tm = 0;
    x_tm.clear(); y_tm.clear();
    depths.clear(); arrows.clear();

    replotAll();
}

void MainWindow::reverseDataYButton_clicked()
{
    if (ui->reverseDataYLabel->text() == "Normal") {
        ui->reverseDataYLabel->setText("Reversed");
    } else {
         ui->reverseDataYLabel->setText("Normal");
    }

    ui->reverseDataYLabel->repaint();

    for (int i=0; i<n_d; ++i) {
         y_d[i] = y_d[i] * (-1.0);
    }

    replotAll();
    if (n_tm > 1) {recalculateData();
        ui->correlation->xAxis->setRange(x_t[0], x_t[n_t-1]);
        ui->correlation->xAxis2->setRange(x_t[0], x_t[n_t-1]);
        ui->correlation->replot();
        ui->tuning->xAxis->setRange(x_t[0], x_t[n_t-1]);
        ui->tuning->xAxis2->setRange(x_t[0], x_t[n_t-1]);
        ui->tuning->replot();
    }
}

void MainWindow::actionBackToOriginal_clicked()  // back to original
{
    ui->reverseTargetXLabel->setText("Normal"); ui->reverseTargetXLabel->repaint();
    ui->reverseTargetYLabel->setText("Normal"); ui->reverseTargetXLabel->repaint();
    ui->reverseDataXLabel->setText("Normal"); ui->reverseTargetXLabel->repaint();
    ui->reverseDataXLabel->setText("Normal"); ui->reverseTargetXLabel->repaint();

    x_d = data_orig.x;
    y_d = data_orig.y;
    x_t = target_orig.x;
    y_t = target_orig.y;
    x_tm = tm_orig.x;
    y_tm = tm_orig.y;
    xind_tm = xind_tm_orig;
    yind_tm = yind_tm_orig;
    n_d = x_d.length();
    n_t = x_t.length();
    n_tm = x_tm.length();
    recalculateData();
    replotAll();
    ui->tuning->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->replot();
    ui->correlation->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
    ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->correlation->replot();
}


void MainWindow::on_actionOpen_Example_Dataset_triggered()
{
    QFile exampleProjectFile = QFile(":/Resources/"+exampleProjectName, this);


    if (exampleProjectFile.open(QIODevice::ReadOnly)) {
        qDebug() << "exampleProjectFileOK = true, fileName = " << exampleProjectFile.fileName();

        QByteArray content = exampleProjectFile.readAll();
        int bytesRead = content.length();
        //qDebug() << "exampleProject Bytes read = " << bytesRead ;

        QString fileContent = QString(content);

        fType = "d";
        double x,y;

        arrows.clear();
        depths.clear();

        x_d.clear();y_d.clear();x_d1.clear();y_dn.clear();x_t.clear();y_t.clear();y_tn.clear();
        x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();n_tm=0;

        QTextStream infile(&fileContent);

        infile >> n_d >> n_t >> n_tm;
        for (int i=0; i<n_d; ++i)
        {   infile >> x >> y;
            x_d.append(x); y_d.append(y);
        }
        for (int i=0; i<n_t; ++i)
        {   infile >> x >> y;
            x_t.append(x); y_t.append(y);
        }
        for (int i=0; i<n_tm; ++i)
        {   infile >> x >> y;
            x_tm.append(x); y_tm.append(y);
        }
        for (int i=0; i<n_tm; ++i)
        {   infile >> x >> y;
            xind_tm.append(x); yind_tm.append(y);
        }

        double minv, maxv;
        minv = *std::min_element(y_d.constBegin(), y_d.constEnd());
        maxv = *std::max_element(y_d.constBegin(), y_d.constEnd());
        y_dn=y_d;// normalized data to [0 1];
        for (int i=0; i<n_d; i++) {
            y_dn[i]=(y_d[i]-minv)/(maxv-minv)+0.1;
        }
        minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
        maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());
        y_tn=y_t;// normalized target to [-1 0];
        for (int i=0; i<n_t; i++) {
            y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;
        }

        recalculateData();
        replotAll();

        if (n_tm>=2) {
            ui->tuning->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
            ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
            ui->tuning->replot();
            ui->correlation->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
            ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
            ui->correlation->replot();
        }

        tm_orig.x = x_tm;
        tm_orig.y = y_tm;
        xind_tm_orig = xind_tm;
        yind_tm_orig = yind_tm;
        data_orig.x=x_d;
        data_orig.y=y_d;
        target_orig.x=x_t;
        target_orig.y=y_t;

        exampleProjectFile.close();
        return;
    } else { // file open error
        qDebug() << "error opening file" << exampleProjectName << ": " << exampleProjectFile.errorString();
    }


    x_d.clear();y_d.clear();x_d1.clear();y_dn.clear();
    x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();n_tm=0;
    x_t.clear();y_t.clear();x_d1.clear();y_tn.clear();

    arrows.clear();
    depths.clear();
    for (int i = ui->tuning->itemCount()-1; i>=0 ; i--) {
        ui->tuning->removeItem(i); // clear arrows
    }

    // Create sinusiod data
    double t,x;
    const double PI  = 3.141592653589793238463;
    t = 0.0;

    for (int i=0; i<1000; ++i)
    {
        x = sin(2*PI*t*10+7.0*sin(2*PI*t))*sin(3*PI*t);
        x_d.append(t); y_d.append(x);
        t = t+0.001;
    }

    n_d = x_d.length();
    x_d1=x_d;

    // Create target
    t = 0.0;
    for (int i=0; i<1000; ++i)
    {
        x = sin(2*PI*t*10);
        x_t.append(t+10); y_t.append(x);
        t = t+0.001;
    }

    n_t=y_t.length();
    data_orig.x=x_d;
    data_orig.y=y_d;
    target_orig.x=x_t;
    target_orig.y=y_t;

    replotAll();
}

void MainWindow::showPointToolTip(QMouseEvent *event)
{

    if (n_d==0 || n_t==0) {
        return;
    }
    int it;

    double y = ui->tuning->yAxis->pixelToCoord(event->pos().y());

    if (y<0.0) {
        double x = ui->tuning->xAxis->pixelToCoord(event->pos().x());
        it = ui->tuning->graph(1)->findBegin(x,true);
        ui->tuning->xAxis->setLabel(QString("Target, x: %1, y:  %2").arg(x_t[it]).arg(y_t[it]));
        ui->tuning->replot();
    } else {
        double x = ui->tuning->xAxis2->pixelToCoord(event->pos().x());
        it = ui->tuning->graph(0)->findBegin(x,true);
        ui->tuning->xAxis2->setLabel(QString("Data, x: %1, y:  %2").arg(x_d[it]).arg(y_d[it]));
        ui->tuning->replot();
    }
}

void MainWindow::actionAbout_triggered()
{
    QMessageBox *msgBox = new QMessageBox(this);
    msgBox->setText("QAnalySeries\n H.Pälike, S. Kotov, 2018-2023 (C)\n https://www.marum.de");
    msgBox->open();
}

void MainWindow::actionOpen_Project_triggered()
{
    QFileDialog::getOpenFileContent(tr("QAS Project Files (*.*)"), [this](const QString &fName, const QByteArray &fData){
        qDebug() << fName<<" opened";
        qDebug() << fData.size() <<" bytes";
        if(fName.isEmpty()){
            qDebug() << "File does not exist: "<<fName;
            qDebug() << fData.size() <<" bytes";
            return;
        } else {
            qDebug() << fName<<" opened";
            qDebug() << fData.size() <<" bytes";

            fileName = fName;
            fileContent = QString(fData);

            fType = "d";
            qDebug() << "about to call openProject from dataFileContentReady";

            //readFile(ui->centralwidget, this);
            qDebug() << "returned from openProject in dataFileContentReady";

            double x,y;

            arrows.clear();
            depths.clear();

            x_d.clear();y_d.clear();x_d1.clear();y_dn.clear();x_t.clear();y_t.clear();y_tn.clear();
            x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();n_tm=0;

            QTextStream infile(&fileContent);

            infile >> n_d >> n_t >> n_tm;
            for (int i=0; i<n_d; ++i)
            {   infile >> x >> y;
                x_d.append(x); y_d.append(y);
            }
            for (int i=0; i<n_t; ++i)
            {   infile >> x >> y;
                x_t.append(x); y_t.append(y);
            }
            for (int i=0; i<n_tm; ++i)
            {   infile >> x >> y;
                x_tm.append(x); y_tm.append(y);
            }
            for (int i=0; i<n_tm; ++i)
            {   infile >> x >> y;
                xind_tm.append(x); yind_tm.append(y);
            }

            double minv, maxv;
            minv = *std::min_element(y_d.constBegin(), y_d.constEnd());
            maxv = *std::max_element(y_d.constBegin(), y_d.constEnd());
            y_dn=y_d;// normalized data to [0 1];
            for (int i=0; i<n_d; i++) {y_dn[i]=(y_d[i]-minv)/(maxv-minv)+0.1;};
            minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
            maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());
            y_tn=y_t;// normalized target to [-1 0];
            for (int i=0; i<n_t; i++) {y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;};

            recalculateData();
            replotAll();
            if (n_tm>=2) {
                ui->tuning->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
                ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
                ui->tuning->replot();
                ui->correlation->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
                ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
                ui->correlation->replot();
            }

            tm_orig.x = x_tm;
            tm_orig.y = y_tm;
            xind_tm_orig = xind_tm;
            yind_tm_orig = yind_tm;
            data_orig.x=x_d;
            data_orig.y=y_d;
            target_orig.x=x_t;
            target_orig.y=y_t;
        }
    }, NULL);
}

void MainWindow::actionSave_Project_as_triggered()
{
    if (command_line == 0) {
        QString outString;
        QTextStream stream(&outString);
        stream << n_d << Qt::endl <<  n_t << Qt::endl << n_tm << Qt::endl ;
        for (int i=0; i<n_d; i++){
            stream << x_d[i] << "\t" << y_d[i] << Qt::endl;
        }
        for (int i=0; i<n_t; i++){
            stream << x_t[i] << "\t" << y_t[i] << Qt::endl;
        }
        for (int i=0; i<n_tm; i++){
            stream << x_tm[i] << "\t" << y_tm[i] << Qt::endl;
        }
        for (int i=0; i<n_tm; i++){
            stream << xind_tm[i] << "\t" << yind_tm[i] << Qt::endl;
        }

        //qDebug() << "actionSave_Project_as_triggered output data = " << outString;
        QFileDialog::saveFileContent(outString.toUtf8(), "myQAProject.qaproj");
    }
}


void MainWindow::actionOpen_Data_triggered()
{
    QFileDialog::getOpenFileContent(tr("Data TextFiles (*.*)"), [this](const QString &fNameD, const QByteArray &fDataD){
        if(fNameD.isEmpty()){
            qDebug() << "File does not exist: "<<fNameD;
            return;
        } else {
            qDebug() << fNameD <<" opened";
            qDebug() << fDataD.size() <<" bytes";

            fileName = fNameD;
            fileContent = QString(fDataD);

            fType = "d";
            qDebug() << "about to call readFile from dataFileContentReady";

            readFile(ui->centralwidget, this);
            qDebug() << "returned from readFile in dataFileContentReady";

        }

    }, nullptr);
}

void MainWindow::openDataCallback(const QString &fileName, const QByteArray &fileContentArray, const MainWindow *window) {

            printf("filename is %s.", fileName.toUtf8().toStdString().c_str());
//            qDebug() << fileName << " opened";
            fileContent = QString(fileContentArray);
//            qDebug() << fileContent.size() <<" bytes";

            fType = "d";
//            qDebug() << "about to call readFile from dataFileContentReady";
//            qDebug() << "returned from readFile in dataFileContentReady";

            arrows.clear();
            depths.clear();
            //replotAll();
            fileContent = "";
            //fileName = "";
            n_d=y_d.length();
            data_orig.x=x_d;
            data_orig.y=y_d;
}

void MainWindow::actionSave_new_Data_as_triggered()
{
    QString outString;
    QTextStream stream(&outString);
    stream << "dataAbscissa" << "\t" << "dataOrdinate" << Qt::endl;
    for (int i=0; i<n_d; i++) {
        stream << x_d[i] << "\t" << y_d[i] << Qt::endl;
    }
//    qDebug() << "actionSave_new_Data_as_triggered output data = " << outString;
    QFileDialog::saveFileContent(outString.toUtf8(), "newData.txt");
}

void MainWindow::actionSave_Data_with_ages_triggered()
{
    QString outString;
    QTextStream stream(&outString);
    if (n_tm<2) {
        QMessageBox *msgBox = new QMessageBox(this);
        msgBox->setText("Create at least 2 tie points");
        msgBox->show();
        return;
    }

    QVector<double> x_d2;
    x_d2.resize(data_orig.x.length());
    mynr::linterpv(x_tm, y_tm, data_orig.x, x_d2);

    stream << "dataAbscissa/Depth" << "\t" << "Age" << "\t" << "dataOrdinate/Value" << Qt::endl;

    for (int i=0; i<data_orig.x.length(); i++) {
        stream << data_orig.x[i] << "\t" << x_d2[i] << "\t" << data_orig.y[i] << Qt::endl;
    }
    qDebug() << "actionSave_Data_with_ages_triggered output data = " << outString;
    QFileDialog::saveFileContent(outString.toUtf8(), "newDataAged.txt");
}




void MainWindow::actionOpen_Target_triggered()
{
    QFileDialog::getOpenFileContent(tr("Target TextFiles (*.*)"), [this](const QString &fNameT, const QByteArray &fDataT){
        if(fNameT.isEmpty()){
            qDebug() << "File does not exist: "<< fNameT;
            return;
        } else {
            qDebug() << fNameT <<" opened";
            qDebug() << fDataT.size() <<" bytes";

            fileName = fNameT;
            fileContent = QString(fDataT);

            fType = "t";
            qDebug() << "about to call readFile from targetFileContentReady";

            readFile(ui->centralwidget, this);
            qDebug() << "returned from readFile in targetFileContentReady";
        }
    }, nullptr);

}

void MainWindow::actionSave_new_Target_as_triggered()
{
    if (command_line == 0) {
        QString outString;
        QTextStream stream(&outString);
        stream << "targetAbscissa" << "\t" << "targetOrdinate" << Qt::endl;

        for (int i=0; i<n_t; i++) {
            stream << x_t[i] << "\t" << y_t[i] << Qt::endl;
        }
        qDebug() << "actionSave_new_Target_as_triggered output data = " << outString;
        QFileDialog::saveFileContent(outString.toUtf8(), "newTarget.txt");
    }
}

void MainWindow::actionOpen_TimeModel_triggered(){
    if (n_d<2 || n_t<2) {
        QMessageBox * msgBox = new QMessageBox(this);
        msgBox->setText("Import Data and Target first!");
        msgBox->setModal(true);
        msgBox->show();
        return;
    }

    QFileDialog::getOpenFileContent(tr("TimeModel TextFiles (*.*)"), [this](const QString &fNameTM, const QByteArray &fDataTM){
        if(fNameTM.isEmpty()){
            qDebug() << "File does not exist: "<< fNameTM;
            return;
        } else {
//            qDebug() << fNameTM <<" opened";
//            qDebug() << fDataTM.size() <<" bytes";

            fileName = fNameTM;
            fileContent = QString(fDataTM);

            fType = "tm";
 //           qDebug() << "about to call readFile from timemodelFileContentReady";

            readFile(ui->centralwidget, this);
 //           qDebug() << "returned from readFile in timemodelFileContentReady";
            // now in tmLoadCallback

        }
    }, nullptr);
}

void MainWindow::actionSave_TimeModel_as_triggered(){
    QString outString;
    QTextStream stream(&outString);
    stream << "From_Data" << "\t" << "To_Target" << Qt::endl;
    for (int i=0; i<n_tm; i++){
        stream << x_tm[i] << "\t" << y_tm[i] << Qt::endl;
    }
//    qDebug() << "actionSave_TimeModel_as_triggered output data = " << outString;
    QFileDialog::saveFileContent(outString.toUtf8(), "newTimeModel.txt");
}



void MainWindow::actionFetchAstroData_triggered() {
    QUrl url(astroFileURL);
    QNetworkRequest request(url);
    QNetworkAccessManager* manager = new QNetworkAccessManager();
    QObject::connect(manager, &QNetworkAccessManager::finished, [this](QNetworkReply *reply) {
        qDebug() << "download finished";
        if (reply->error()) {
            qDebug() << "error" << reply->errorString();
        } else {
            QByteArray payload = reply->readAll();
//            qDebug() << "size" << payload.size();

            QByteArray hash = QCryptographicHash::hash(payload, QCryptographicHash::Sha256);
//            qDebug() << "sha256" << hash.toHex();
            QFile laskarFile = QFile(laskarFileName, this);

            if (laskarFile.open(QIODevice::ReadWrite)) {
//                qDebug() << "laskarFileOK = true, fileName = " << laskarFile.fileName();

                int bytesWritten = laskarFile.write(payload);
                laskarFile.flush();
                laskarFile.close();
//                qDebug() << "laskarFile Bytes written = " << bytesWritten ;
            } else { // file open error
//                qDebug() << "error opening file" << laskarFileName << ": " << laskarFile.errorString();
                return;
            }

            db = QSqlDatabase::addDatabase("QSQLITE");
            db.setDatabaseName(laskarFileName);
            bool success = db.open();
//            qDebug() << "db.open == success " << success;
            ui->actionOpen_Laskar->setEnabled(true);
            ui->actionInsolation_La04->setEnabled(true);
            ui->actionETP_Target->setEnabled(true);
//            qDebug() << "enabled Astro Menus ";
        }

    });

//    qDebug() << "download begin";
    manager->get(request);
}

void loadDatabaseFromNetwork(QSqlDatabase& db)
{
    // Download the SQLite database file from the network
    QNetworkAccessManager manager;
    QNetworkReply* reply = manager.get(QNetworkRequest(QUrl(NETWORK_LOC)));

    // Wait for the download to finish
    QEventLoop loop;
    QObject::connect(reply, &QNetworkReply::finished, &loop, &QEventLoop::quit);
    loop.exec();

    // Check if the download was successful
    if (reply->error() != QNetworkReply::NoError) {
        qDebug() << "Error downloading the database:" << reply->errorString();
        return;
    }

    // Save the downloaded data to a temporary file
    QString tempFileName = QStandardPaths::writableLocation(QStandardPaths::TempLocation) + "/temp_database.sqlite";
    QFile tempFile(tempFileName);
    if (!tempFile.open(QIODevice::WriteOnly)) {
        qDebug() << "Error creating temporary file:" << tempFile.errorString();
        return;
    }
    tempFile.write(reply->readAll());
    tempFile.close();

    // Load the SQLite database from the temporary file
    db = QSqlDatabase::addDatabase("QSQLITE", "in_memory_db");
    db.setDatabaseName(":memory:");
    if (!db.open()) {
        qDebug() << "Error opening in-memory database:" << db.lastError().text();
        return;
    }

    // Read the contents of the temporary file
    QFile databaseFile(tempFileName);
    if (!databaseFile.open(QIODevice::ReadOnly)) {
        qDebug() << "Error opening database file:" << databaseFile.errorString();
        return;
    }
    QByteArray databaseData = databaseFile.readAll();
    databaseFile.close();

    QSqlQuery query(db);
    query.prepare("ATTACH DATABASE ':memory:' AS temp_db");
    if (!query.exec()) {
        qDebug() << "Error attaching in-memory database:" << query.lastError().text();
        return;
    }
    query.prepare("SELECT sql FROM temp_db.sqlite_master WHERE type='table'");
    if (!query.exec()) {
        qDebug() << "Error retrieving table SQL statements:" << query.lastError().text();
        return;
    }
    while (query.next()) {
        QString createTableQuery = query.value(0).toString();
        createTableQuery.replace("CREATE TABLE", "CREATE TABLE IF NOT EXISTS");
        if (!query.exec(createTableQuery)) {
            qDebug() << "Error creating table:" << query.lastError().text();
            return;
        }
    }
    query.prepare("DETACH DATABASE temp_db");
    if (!query.exec()) {
        qDebug() << "Error detaching in-memory database:" << query.lastError().text();
        return;
    }

    // Clean up the temporary file
    tempFile.remove();
}

void MainWindow::initLaskarDB() {

    db = QSqlDatabase::addDatabase("QSQLITE");
    tmpFile.setFileTemplate("XXXXXX.db3");

    if (tmpFile.open()) {
        QString tmp_filename=tmpFile.fileName();
        qDebug() << "temporary" << tmp_filename;

        QFile file(":Resources/laskar.db3");
        if (file.open(QIODevice::ReadOnly)) {
            tmpFile.write(file.readAll());
        } else {
            qDebug() << "cannot find Resource.\n";
        }

        tmpFile.close();
    }

    db.setDatabaseName(tmpFile.fileName());
    qDebug() << "db = " << db;

    int success = db.open();
    qDebug() << "db.open == success " << success;
}

void MainWindow::actionOpen_Laskar_triggered()
{
//    qDebug() << "in actionOpen_Laskar_triggered";
    arrows.clear();
    depths.clear();

    for (int i = tuning->itemCount()-1; i>=0 ; i--) tuning->removeItem(i); // clear arrows
    Laskar *laskar = new Laskar(this);

    QObject::connect(laskar, SIGNAL(accepted()), this, SLOT(laskarCallback()));

    laskar->setModal(true);
//    qDebug() << "in actionOpen_Laskar_triggered: calling open()";
    laskar->open();
}

void MainWindow::actionETP_Target_triggered()
{
    arrows.clear();
    depths.clear();
    for (int i = ui->tuning->itemCount()-1; i>=0 ; i--) {
        tuning->removeItem(i); // clear arrows
    }

    ETP *etp = new ETP(this);

    QObject::connect(etp, SIGNAL(accepted()), this, SLOT(etpCallback()));

    etp->setModal(true);
    etp->open();
}

void MainWindow::actionInsolation_La04_triggered(){
    arrows.clear();
    depths.clear();
    for (int i = ui->tuning->itemCount()-1; i>=0 ; i--) {
        ui->tuning->removeItem(i); // clear arrows
    }

    Insolation *insolation = new Insolation(this);

    QObject::connect(insolation, SIGNAL(accepted()), this, SLOT(insolationLa04Callback()));

    insolation->setModal(true);
    insolation->open();
}


void MainWindow::dataLoadCallback()
{
    arrows.clear();
    depths.clear();
    for (int i = this->ui->tuning->itemCount()-1; i>=0 ; i--) {
        this->ui->tuning->removeItem(i); // clear arrows
    }

    replotAll();
    fileContent = "";
    fileName = "";

    n_d=y_d.length();
    data_orig.x=x_d;
    data_orig.y=y_d;
}

void MainWindow::targetLoadCallback() {
    target_orig.x=x_t;
    target_orig.y=y_t;
    n_t=y_t.length();

    arrows.clear();
    depths.clear();
    for (int i = this->ui->tuning->itemCount()-1; i>=0 ; i--) {
        this->ui->tuning->removeItem(i); // clear arrows
    }
    replotAll();
    fileContent = "";
    fileName = "";
}

void MainWindow::tmLoadCallback()
{
//    qDebug() << "in tmLoadCallback" << Qt::endl;
    arrows.clear();
    depths.clear();
    for (int i = this->ui->tuning->itemCount()-1; i>=0 ; i--) {
        this->ui->tuning->removeItem(i); // clear arrows
    }

    double x,y, interp;
    QVector<double> xtmp(2), ytmp(2);

    x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();
    int i=0;

    for (int j = 0; j<tm_orig.x.length();j++)
    {
        x = tm_orig.x[j]; y = tm_orig.y[j];
        x_tm.append(x); y_tm.append(y);

        i=0;
        while (i<x_d.length() && x>x_d[i]) {
            i++; xind_tm.append(i); //
        }

        if (i>(x_d.length()-1)) {
            xtmp[0]=x_d[i-2]; xtmp[1]=x_d[i-1];
            ytmp[0]=y_d[i-2]; ytmp[1]=y_d[i-1];
        }
        else if (i==0  && x!=x_d[i]) {
            xtmp[0]=x_d[0]; xtmp[1]=x_d[1];
            ytmp[0]=y_d[0]; ytmp[1]=y_d[1];
        }
        else if (i>0 && i<=(x_d.length()-1) && x!=x_d[i]) {
            xtmp[0]=x_d[i-1]; xtmp[1]=x_d[i];
            ytmp[0]=y_d[i-1]; ytmp[1]=y_d[i];}

        if (i>(x_d.length()-1) || x!=x_d[i]) {
            x_d.insert(i,x);
            mynr::linterp(xtmp,ytmp,x_d[i],interp);
            y_d.insert(i,interp);
        }

        i=0;
        while (i<x_t.length() && y>x_t[i]) {
            i++; yind_tm.append(i);
        }

        if (i>(x_t.length()-1)) {
            xtmp[0]=x_t[i-2]; xtmp[1]=x_t[i-1];
            ytmp[0]=y_t[i-2]; ytmp[1]=y_t[i-1];}
        else if (i==0 && y!=x_t[i] ) {
            xtmp[0]=x_t[0]; xtmp[1]=x_t[1];
            ytmp[0]=y_t[0]; ytmp[1]=y_t[1];
        }
        else if (i>0 && i<=(x_t.length()-1) && y!=x_t[i]) {
            xtmp[0]=x_t[i-1]; xtmp[1]=x_t[i];
            ytmp[0]=y_t[i-1]; ytmp[1]=y_t[i];
        }

        if (i>(x_t.length()-1) || y!=x_t[i]) {
            x_t.insert(i,y);
            mynr::linterp(xtmp,ytmp,x_t[i],interp);
            y_t.insert(i,interp);}

    }


    n_tm=x_tm.length();
    n_d=x_d.length();
    n_t=x_t.length();

    y_dn=y_d;// normalized target to [-1 0];
    double minv, maxv;

    minv = *std::min_element(y_d.constBegin(), y_d.constEnd());
    maxv = *std::max_element(y_d.constBegin(), y_d.constEnd());

    for (int i=0; i<n_d; i++) {
        y_dn[i]=(y_d[i]-minv)/(maxv-minv)+0.1;
    }

    y_tn=y_t;// normalized target to [-1 0];

    minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
    maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());

    for (int i=0; i<n_t; i++) {
        y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;
    }

    recalculateData();
    replotAll();

    ui->tuning->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->replot();
    ui->correlation->xAxis->setRange(x_d1[0], x_d1[n_d-1]);
    ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->correlation->replot();

    tm_orig.x = x_tm;
    tm_orig.y = y_tm;
    xind_tm_orig = xind_tm;
    yind_tm_orig = yind_tm;
    data_orig.x=x_d;
    data_orig.y=y_d;
    target_orig.x=x_t;
    target_orig.y=y_t;
}

void MainWindow::setCurrentIndex()
{
    printf("in setCurrentIndex\n");
}

void MainWindow::etpCallback()
{
    qDebug() << "in etpCallback\n";
    if (n_t == 0) {
        qDebug() << "in etpCallback: n_t == 0";
        return;
    }
    n_t=y_t.length();
//    qDebug() << "in etpCallback: n_t == " << n_t;
    ui->correlation->graph(1)->data()->clear();
    ui->correlation->graph(1)->setData(x_t,y_t);
    ui->correlation->graph(1)->rescaleAxes();
    ui->correlation->setVisible(true);
    ui->correlation->replot();

    y_tn=y_t;// normalized target to [-1 0];
    double minv, maxv;

    minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
    maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());

    for (int i=0; i<n_t; i++) {
        y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;
    };
    ui->tuning->graph(1)->data()->clear();
    ui->tuning->graph(1)->setData(x_t,y_tn);
    ui->tuning->xAxis->setRange(x_t[0], x_t[n_t-1]);
    ui->tuning->setVisible(true);
    ui->tuning->replot();
    target_orig.x=x_t;
    target_orig.y=y_t;
}

void MainWindow::laskarCallback()
{
//    qDebug() << "in laskarsCallback\n";
    n_t=y_t.length();
    if (n_t == 0) {
        qDebug() << "in laskarCallback: n_t == 0";
        return;
    }
//    qDebug() << "in laskarCallback: n_t == " << n_t;

    ui->correlation->graph(1)->data()->clear();
    ui->correlation->graph(1)->setData(x_t,y_t);
    ui->correlation->graph(1)->rescaleAxes();
    ui->correlation->setVisible(true);
    ui->correlation->replot();

    y_tn=y_t;// normalized target to [-1 0];
    double minv, maxv;

    minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
    maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());

    for (int i=0; i<n_t; i++) {
        y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;
    }

    ui->tuning->graph(1)->data()->clear();
    ui->tuning->graph(1)->setData(x_t,y_tn);
    ui->tuning->xAxis->setRange(x_t[0], x_t[n_t-1]);
    ui->tuning->setVisible(true);
    ui->tuning->replot();
    target_orig.x=x_t;
    target_orig.y=y_t;
}

void MainWindow::insolationLa04Callback(){
    n_t=y_t.length();
    if (n_t == 0) {
        return;
    }

//    qDebug() << "in insolationLa04Callback\n";
    n_t=y_t.length();

    if (n_t == 0) {
        qDebug() << "in insolationLa04Callback: n_t == 0";
        return;
    }
//    qDebug() << "in insolationLa04Callback: n_t == " << n_t;

    ui->correlation->graph(1)->data()->clear();
    ui->correlation->graph(1)->setData(x_t,y_t);
    ui->correlation->graph(1)->rescaleAxes();
    ui->correlation->setVisible(true);
    ui->correlation->replot();

    y_tn=y_t;// normalized target to [-1 0];
    double minv, maxv;

    minv = *std::min_element(y_t.constBegin(), y_t.constEnd());
    maxv = *std::max_element(y_t.constBegin(), y_t.constEnd());

    for (int i=0; i<n_t; i++) {
        y_tn[i]=(y_t[i]-minv)/(maxv-minv)-1.1;
    };
    ui->tuning->graph(1)->data()->clear();
    ui->tuning->graph(1)->setData(x_t,y_tn);
    ui->tuning->xAxis->setRange(x_t[0], x_t[n_t-1]);
    ui->tuning->setVisible(true);
    ui->tuning->replot();
    target_orig.x=x_t;
    target_orig.y=y_t;
}

void MainWindow::toggleMarkersCheckbox_stateChanged(int state)
{
    if (ui->MarkersCheckbox->isChecked()) {
        ui->tuning->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
        ui->tuning->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));
        ui->correlation->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::white, 5));
        ui->correlation->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::red, Qt::white, 5));

    }
    else {
        ui->tuning->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot));
        ui->tuning->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot));
        ui->correlation->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot));
        ui->correlation->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDot));

    }
    ui->tuning->replot();
    ui->correlation->replot();

}

void MainWindow::toggleCorrelationCheckbox_stateChanged(int state)
{
   recalculateData();
}


void MainWindow::actionData_Processor_triggered()
{
    DataProc *dataproc;
    dataproc = new DataProc(this);
    dataproc->setModal(true);

    QObject::connect(dataproc, SIGNAL(accepted()), this, SLOT(dataProcessingCallback()));
    auto dialogClosed = [=](int code) {
        Q_UNUSED(code);
        delete dataproc;
    };
    connect(dataproc, &DataProc::finished, dialogClosed);
    dataproc->show();

}

void MainWindow::on_actionSSA_triggered()
{
    SSA *ssa;
    ssa = new SSA(this);
    ssa->setModal(true);
    QObject::connect(ssa, SIGNAL(accepted()), this, SLOT(ssaProcessingCallback()));

    ssa->show();

}


void MainWindow::ssaProcessingCallback(){
    n_d = x_d.length();
    n_t = x_t.length();
    n_tm = 0;
    x_tm.clear(); y_tm.clear();
    depths.clear(); arrows.clear();
    recalculateData();
    replotAll();

    ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
}

void MainWindow::dataProcessingCallback()
{
    n_d = x_d.length();
    x_d1 = x_d; //should I ?
    n_tm = 0;
    x_tm.clear(); y_tm.clear();
    depths.clear(); arrows.clear();
    recalculateData();
    replotAll();

    ui->correlation->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
    ui->tuning->xAxis2->setRange(x_d1[0], x_d1[n_d-1]);
}


void MainWindow::spectralProcessingCallback() {
//    qDebug() << "in spectralProcessingCallback" << Qt::endl;
}

void MainWindow::actionSpectral_Analysis_triggered() {
    if ((x_d.count() == 0) && (x_t.count() == 0)) {
        QMessageBox *msgBox = new QMessageBox(this);
        msgBox->setText("Please load Data and or Target first");
        msgBox->setModal(true);
        msgBox->show();
        return;
    }

    SpectralDialog * spectralDialog;
    spectralDialog = new SpectralDialog(this);

    auto dialogClosed = [=](int code) {
        Q_UNUSED(code);
        delete spectralDialog;
    };

    spectralDialog->show();
};





