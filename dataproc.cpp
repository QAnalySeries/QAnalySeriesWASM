#include "dataproc.h"
#include "ui_dataproc.h"
#include "globals.h"
#include "mynr.h"

QVector<double> xx, yy;

DataProc::DataProc(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DataProc)
{
    ui->setupUi(this);

    connect(ui->trimData, &QPushButton::pressed, this, &DataProc::on_trimData_clicked);
    connect(ui->comboBox_2, &QComboBox::currentTextChanged, this, &DataProc::on_comboBox_2_currentIndexChanged);

    ui->customPlot->addGraph(0);
    ui->customPlot->xAxis->setLabel("x");
    ui->customPlot->yAxis->setLabel("y");
    ui->customPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, Qt::white, 5));
    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->comboBox_7->addItem("1"); ui->comboBox_7->addItem("2"); ui->comboBox_7->addItem("3");
    ui->comboBox_7->setCurrentIndex(2);
    ui->comboBox->addItem("Taner Filter"); ui->comboBox->addItem("Envelope");
    ui->comboBox_2->addItem("Data"); ui->comboBox_2->addItem("Target");


    ui->spectrumPlot->addGraph(0);
    ui->spectrumPlot->xAxis->setLabel("Frequency (cycles/unit)");
    ui->spectrumPlot->yAxis->setLabel("PSD");
    ui->spectrumPlot->addGraph();
    ui->spectrumPlot->graph(1)->setPen(QPen(Qt::red));// taner filter
    ui->spectrumPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    connect(ui->customPlot, &QCustomPlot::mouseWheel, this, &DataProc::mouseWheel);

    if (ui->comboBox_2->currentText() == "Data") {
        xx = x_d;
        yy = y_d;
    } else {
        xx = x_t;
        yy = y_t;
    }

    setupPlot();
    ui->customPlot->replot();
}

DataProc::~DataProc()
{
    delete ui;
}

void DataProc::on_saveResultsAndCloseWindow_clicked()
{
//    qDebug() << "in on_saveResultsAndCloseWindow_clicked" << Qt::endl;
    if (ui->comboBox_2->currentText() == "Data") {
        x_d = xx;
        y_d = yy;
        n_d = x_d.length();
    } else {
        x_t = xx;
        y_t = yy;
        n_t = x_t.length();
    }
    DataProc::close();
}


void DataProc::setupPlot() {
    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
}

void DataProc::on_pushButton_5_clicked() { // detrend
    int n_d;
    n_d = xx.length();

    if (n_d < 2) {QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;}

    double a, b;
    mynr::lr(xx, yy, a, b); // calculate linear regression y = a + b *x

    double shift = a + b*xx[0];

    for (int i=0; i<n_d; i++) {
        yy[i] = yy[i] - (a + b*xx[i]) + shift;
        }

        ui->customPlot->graph(0)->data()->clear();
        ui->customPlot->graph(0)->setData(xx, yy);
        ui->customPlot->graph(0)->rescaleAxes();
        ui->customPlot->replot();

}

void DataProc::on_pushButton_6_clicked() { // back
    xx.clear(); yy.clear();
    if (ui->comboBox_2->currentText() == "Data") {
        xx = x_d;
        yy = y_d;
    } else {
        xx = x_t;
        yy = y_t;
    }

    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}

void DataProc::on_pushButton_7_clicked() { // remove outliers
    int n_d;
    int n_sigma = 1;
    double mean, sigma, var;

    n_d=xx.length();

    if (n_d < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    switch ( ui->comboBox_7->currentIndex() )
    {
    case 0:
        n_sigma = 1;
        break;
    case 1:
        n_sigma = 2;
        break;
    case 2:
        n_sigma = 3;
        break;
    }


    mynr::meanandvar(yy, mean, var);
    sigma = n_sigma * sqrt(var);

    for (int i=n_d-1; i>=0; i--) {if (yy[i]>(mean+sigma) || yy[i]<(mean-sigma)) {xx.removeAt(i); yy.removeAt(i);}}
    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}

void DataProc::on_trimData_clicked() { // select subset from data
    QVector<double> x_tmp, y_tmp;
    double x1, x2, y1, y2;

    x_tmp=xx; y_tmp=yy; n_d=xx.length();

    if (n_d < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    bool ok;

    if (ui->lineEdit->text().isEmpty())   {
        x1 = xx.first();
    } else {
        x1 = ui->lineEdit->text().toDouble(&ok);
        if (!ok) {
            QMessageBox msgBox;
            msgBox.setText("Invalid number format");
            msgBox.show();
            return;}
    }
    if (ui->lineEdit_2->text().isEmpty()) {
        x2 = xx.last();} else {x2 = ui->lineEdit_2->text().toDouble(&ok);
        if (!ok) {
            QMessageBox msgBox;
            msgBox.setText("Invalid number format");
            msgBox.show();
            return;
        }
    }

    if (ui->lineEdit_4->text().isEmpty()) {
        y1 = *std::min_element(yy.constBegin(), yy.constEnd());
    }
    else {
        y1 = ui->lineEdit_4->text().toDouble(&ok);
        if (!ok) {
            QMessageBox msgBox;
            msgBox.setText("Invalid number format");
            msgBox.show();
            return;}
    }
    if (ui->lineEdit_5->text().isEmpty()) {
        y2 = *std::max_element(yy.constBegin(), yy.constEnd());
    }
    else {
        y2 = ui->lineEdit_5->text().toDouble(&ok);
        if (!ok) {
            QMessageBox msgBox;
            msgBox.setText("Invalid number format");
            msgBox.show();
            return;
        }
    }

    xx.clear(); yy.clear();

    for (int i=0; i<n_d; i++) {
        if (x_tmp[i] >= x1 && x_tmp[i] <= x2 && y_tmp[i] >= y1 && y_tmp[i] <= y2) {
            xx.append(x_tmp[i]);
            yy.append(y_tmp[i]);
        }
    }

    n_d = xx.length();
    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}

void DataProc::on_pushButton_8_clicked() // resample evenly
{
    int n_d;
    QVector<double> x_tmp, y_tmp;
    double step, x, y, xmax;
    bool ok;
    step=ui->lineEdit_3->text().toDouble(&ok);

    if (!ok) {
        QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }

    n_d=xx.length(); xmax=xx.last(); x=xx.first();

    if (n_d < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    x_tmp.append(xx[0]); y_tmp.append(yy[0]);
    for (int i=1; i<n_d; i++) {
        if (xx[i]>x_tmp.last())  { x_tmp.append(xx[i]); y_tmp.append(yy[i]);}
    }

    xmax=x_tmp.last();
    xx.clear(); yy.clear();

    while (x<=xmax) {
      xx.append(x);
      mynr::linterp(x_tmp, y_tmp, x, y);
      yy.append(y);
      x+=step;
    }

    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}

void DataProc::on_pushButton_9_clicked() { // remove jumps (n*sigma(sec. derivatives))

    int n_d;
    int n_sigma = 1;
    double mean, var, sigma;
    QVector<double> d2;

    n_d=xx.length();

    if (n_d < 5) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    switch ( ui->comboBox_7->currentIndex() )
    {
    case 0:
        n_sigma = 1;
        break;
    case 1:
        n_sigma = 2;
        break;
    case 2:
        n_sigma = 3;
        break;

    }

    // calculate array of second derivatives
    for (int i=1; i<(n_d-1); i++) {
        if (xx[i]==xx[i-1] || xx[i+1]==xx[i-1] || xx[i+1]==xx[i]) {
            d2.append(0.0);
        } else {
            d2.append(  2.0*yy[i-1]/(xx[i]-xx[i-1])/(xx[i+1]-xx[i-1])
                      - 2.0*yy[i]/(xx[i+1]-xx[i])/(xx[i]-xx[i-1])
                      + 2.0*yy[i+1]/(xx[i+1]-xx[i])/(xx[i+1]-xx[i-1])
                     );
        }
    }


    mynr::meanandvar(d2, mean, var);
    sigma = n_sigma*sqrt(var);

    for (int i=n_d-3; i>0; i--) {
        if (d2[i]>(mean+sigma) || d2[i]<(mean-sigma)) {
            xx.removeAt(i+1); yy.removeAt(i+1);
        }
    }

    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}


void DataProc::on_tabWidget_currentChanged(int index) {// Filter tab activation
    double dt;
    if  (ui->tabWidget->currentIndex()!=1) {
        return;
    }

    if (yy.length() < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    dt = std::abs(xx.last() - xx.first())/(xx.length()-1.0);
    QVector <double> x_psd, psd, filtv;
    psd = mynr::periodogram(yy);
    psd_min = *std::min_element(psd.constBegin()+1, psd.constEnd());
    psd_max = *std::max_element(psd.constBegin()+1, psd.constEnd());


    for (int i=0;i<psd.length();++i) {
        x_psd.append((float) i/yy.length()/dt);
    }
    int halfpsd = floor((float) x_psd.length()/2.0);
    ui->spectrumPlot->graph(0)->data()->clear();
    ui->spectrumPlot->graph(0)->setData(x_psd, psd);
    ui->spectrumPlot->yAxis->setRange(psd_min, psd_max);
    ui->spectrumPlot->xAxis->setRange(0, x_psd[psd.length()-1]);

    std::vector<double> filt = mynr::taner(yy, dt, 0, x_psd[halfpsd], 10000000000);
    filtv = QVector<double>(filt.begin(), filt.end());
    filtv.resize(x_psd.length());
    ui->spectrumPlot->graph(1)->data()->clear();
    ui->spectrumPlot->graph(1)->setData(x_psd, filtv);
    ui->spectrumPlot->replot();

    ui->lineEdit_6->setText("0.0");
    ui->lineEdit_7->setText(QString::number(x_psd[halfpsd]));

}

void DataProc::on_pushButton_11_clicked()  { // normalize data
    if (yy.length() < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }
    int n=yy.length();
    double mean, var;
    mynr::meanandvar(yy, mean, var);
    var = sqrt(var);  //standard deviation

    for (int i=0;i<n;++i) {
        yy[i] = (yy[i]-mean)/var;
    }
    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}

void DataProc::on_pushButton_12_clicked()  // filter
{
    if (yy.length() < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    bool ok;
    double fl = ui->lineEdit_6->text().toDouble(&ok);
    if (!ok) {
        QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }

    double fh = ui->lineEdit_7->text().toDouble(&ok);

    if (!ok) {
        QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }

    if (fl >= fh) {
        QMessageBox msgBox;
        msgBox.setText("fl <= fh");
        msgBox.show();
        return;
    }

    double dt = std::abs(xx.last() - xx.first())/(xx.length()-1.0);
    QVector<double> yv;
    QVector <double> x_psd, psd;
    std::vector<double> filt;

    if (ui->comboBox->currentText() == "Envelope") {
        yy = mynr::envelope(yy);
    } else {
        filt = mynr::taner(yy, dt, fl, fh, 10000000000);
        yy = mynr::filter(yy, filt);
    }

    psd = mynr::periodogram(yy);
    psd_min = *std::min_element(psd.constBegin()+1, psd.constEnd());
    psd_max = *std::max_element(psd.constBegin()+1, psd.constEnd());

    for (int i=0;i<psd.length();++i) {
        x_psd.append((float) i/yy.length()/dt);
    }
    ui->spectrumPlot->graph(0)->data()->clear();
    ui->spectrumPlot->graph(0)->setData(x_psd, psd);
    ui->spectrumPlot->graph(0)->rescaleAxes();

    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->spectrumPlot->yAxis->setRange(psd_min, psd_max);


    QVector<double> filtv = QVector<double>(filt.begin(), filt.end());
    filtv.resize(x_psd.length());
    ui->spectrumPlot->graph(1)->data()->clear();
    ui->spectrumPlot->graph(1)->setData(x_psd, filtv);
    ui->spectrumPlot->replot();
    ui->customPlot->replot();
}

void DataProc::on_pushButton_13_clicked() // design filter
{
    if (yy.length() < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }

    bool ok;
    double fl = ui->lineEdit_6->text().toDouble(&ok);
    if (!ok) {
        QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }

    double fh = ui->lineEdit_7->text().toDouble(&ok);

    if (!ok) {
        QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }
    if (fl >= fh) {
        QMessageBox msgBox;
        msgBox.setText("fl <= fh");
        msgBox.show();
        return;
    }

    double  dt = std::abs(xx.last() - xx.first())/(xx.length()-1.0);
 //   QVector<double> yv;
    QVector <double> x_psd, psd;

    std::vector<double> filt = mynr::taner(yy, dt, fl, fh, 10000000000);
    psd = mynr::periodogram(yy);
    for (int i=0;i<psd.length();++i) {
        x_psd.append((float) i/yy.length()/dt);
    }
    QVector<double> filtv = QVector<double>(filt.begin(), filt.end());
    filtv.resize(x_psd.length());
    ui->spectrumPlot->graph(1)->data()->clear();
    ui->spectrumPlot->graph(1)->setData(x_psd, filtv);
    ui->spectrumPlot->yAxis->setRange(psd_min, psd_max);
    ui->spectrumPlot->replot();
    ui->customPlot->replot();
}

void DataProc::on_checkBox_clicked()
{
    if (ui->checkBox->isChecked()) {
        ui->spectrumPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    } else {
        ui->spectrumPlot->yAxis->setScaleType(QCPAxis::stLinear);
    }
    ui->spectrumPlot->yAxis->setRange(psd_min, psd_max);
    ui->spectrumPlot->replot();
}

void DataProc::on_pushButton_2_clicked() //cancel
{
    DataProc::close();
}

void DataProc::on_comboBox_2_currentIndexChanged(const QString &arg1)
{
//    qDebug() << "in on_comboBox_2_currentIndexChanged" << Qt::endl;
    if (ui->comboBox_2->currentText() == "Data") {
        xx = data_orig.x;
        yy = data_orig.y;
    } else {
        xx = target_orig.x;
        yy = target_orig.y;
    }
    setupPlot();
    ui->customPlot->replot();
}

void DataProc::mouseWheel()
{
    QList<QCPAxis *> axlist;
    if(QApplication::keyboardModifiers() & Qt::ShiftModifier) {
        axlist << ui->customPlot->xAxis <<  ui->customPlot->xAxis << ui->customPlot->yAxis;
    } else {
        axlist.clear(); axlist << ui->customPlot->xAxis;
    }

    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iSelectItems);
    ui->customPlot->axisRect()->setRangeDragAxes(axlist);
    ui->customPlot->axisRect()->setRangeZoomAxes(axlist);

}

void DataProc::on_pushButton_3_clicked()
{
    if (yy.length() < 2) {
        QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;
    }
    int n=yy.length();

    double miny = *std::min_element(yy.constBegin(), yy.constEnd());

    for (int i=0;i<n;++i) {
        yy[i] = qLn(yy[i]+1.0-miny);
    }
    ui->customPlot->graph(0)->data()->clear();
    ui->customPlot->graph(0)->setData(xx, yy);
    ui->customPlot->graph(0)->rescaleAxes();
    ui->customPlot->replot();
}
