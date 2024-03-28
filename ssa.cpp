#include "ssa.h"
#include "ui_ssa.h"
#include "globals.h"
#include "mynr.h"
#include "math.h"

int delay, nn;
QVector<double>  yy1_rec, eigval, eof, xx1, yy1;
QVector<QVector<double>> eigvec, rcs;

SSA::SSA(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SSA)
{
    ui->setupUi(this);

    delay = ui->lineEdit->text().toInt(); // delay window size

    ui->xyPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->xyPlot->addGraph(0);
    ui->xyPlot->addGraph();
    ui->xyPlot->xAxis->setLabel("Depth/Time");

    ui->evalPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    ui->evalPlot->addGraph(0);
    ui->evalPlot->xAxis->setLabel("N");
    ui->evalPlot->yAxis->setLabel("Eigen Values");
    ui->evalPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));

    ui->eofPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->eofPlot->legend->setVisible(true);
    ui->eofPlot->xAxis->setLabel("N");
    ui->eofPlot->yAxis->setLabel("EOFs");

    ui->spectrumPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->spectrumPlot->xAxis->setLabel("Norm. freq.");
    ui->spectrumPlot->yAxis->setLabel("SD");

    ui->comboBox->addItem("Data"); ui->comboBox->addItem("Target");

    connect(ui->xyPlot, &QCustomPlot::mouseWheel, this, &SSA::mouseWheel);

    if (ui->comboBox->currentText() == "Data") {
        xx1 = x_d;
        yy1 = y_d;
        nn = n_d;
    } else {
        xx1 = x_t;
        yy1 = y_t;
        nn = n_t;
    }



    ui->xyPlot->graph(0)->setData(xx1,yy1);
    ui->xyPlot->graph(0)->rescaleAxes();
    ui->xyPlot->replot();
}

SSA::~SSA()
{
    delete ui;
}

void SSA::on_pushButton_3_clicked()  // save results and close window
{
    if (ui->comboBox->currentText() == "Data") {
        x_d = xx1;
        y_d = yy1_rec;
        n_d = x_d.length();
    } else {
        x_t = xx1;
        y_t = yy1_rec;
        n_t = x_t.length();
    }

    SSA::close();
    accept();
}


void SSA::on_pushButton_clicked() // Run SSA
{
    if (yy1.length() < 2) {QMessageBox msgBox;
        msgBox.setText("Get data first!");
        msgBox.show();
        return;}
    delay = ui->lineEdit_2->text().toInt(); // delay window size

    QVector<double>   x, y;
    mynr::ssa(yy1, delay, eigval, eigvec, rcs); // SSA
    for (int i=0; i < eigval.length(); ++i) {
        x<<i;
        y<<eigval[i];
    }

    ui->evalPlot->graph(0)->setData(x,y);
    ui->evalPlot->graph(0)->rescaleAxes();
    ui->evalPlot->replot();
}

void SSA::on_pushButton_2_clicked() // Reconstruct
{
    if (rcs.length() < 2) {
        QMessageBox msgBox;
        msgBox.setText("Run SSA first!");
        msgBox.show();
        return;
    }

    ui->eofPlot->clearGraphs();
    ui->spectrumPlot->clearGraphs();
    QStringList selList;
    selList = ui->lineEdit->text().split(QRegularExpression("\\W+"), Qt::SkipEmptyParts);
    QVector<double>   x, y, freq, psd;
    yy1_rec.clear(); yy1_rec.resize(nn);


    QStringList colors = {"red", "green","blue","black", "cyan"};
    x.clear(); y.clear(); freq.clear(); psd.clear();

    int col = 0;
    for (int i=0; i<delay; ++i) {
        x<<i;
    }
    foreach (QString sel, selList) {
        y = eigvec.at(sel.toInt());
        ui->eofPlot->addGraph();
        ui->eofPlot->graph()->setData(x,y);
        ui->eofPlot->graph()->setName(sel);
        ui->eofPlot->graph()->setPen(QPen(colors.at(col)));


        psd = mynr::periodogram(rcs[sel.toInt()]);
        freq.clear();
        for (int i=0; i<psd.length(); ++i) freq<< (float) i/psd.length()*0.5;

        ui->spectrumPlot->addGraph();
        ui->spectrumPlot->graph()->setData(freq,psd);
        ui->spectrumPlot->graph()->setName(sel);
        ui->spectrumPlot->graph()->setPen(QPen(colors.at(col)));


        for (int i=0; i<nn;++i) yy1_rec[i] = yy1_rec[i] + rcs[sel.toInt()][i];
        ++col;
        if (col==colors.size()) {
            col=0;
        }
    }
    ui->eofPlot->graph()->rescaleAxes();
    ui->eofPlot->replot();
    ui->spectrumPlot->graph()->rescaleAxes();
    ui->spectrumPlot->replot();


    ui->xyPlot->graph(1)->setData(xx1,yy1_rec);
    ui->xyPlot->graph(1)->rescaleAxes();
    ui->xyPlot->graph(1)->setPen(QPen(Qt::red));
    ui->xyPlot->replot();
}

void SSA::on_checkBox_stateChanged(int arg1)
{
    if (ui->checkBox->isChecked()) {
        ui->spectrumPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else { ui->spectrumPlot->yAxis->setScaleType(QCPAxis::stLinear);
    }
    ui->spectrumPlot->replot();
}

void SSA::on_pushButton_4_clicked() // back
{
    xx1.clear(); yy1.clear();
    if (ui->comboBox->currentText() == "Data") {
        xx1 = x_d;
        yy1 = y_d;
    } else {
        xx1 = x_t;
        yy1 = y_t;
    }

    ui->xyPlot->graph(0)->data()->clear();
    ui->xyPlot->graph(0)->setData(xx1, yy1);
    ui->xyPlot->graph(0)->rescaleAxes();
    ui->xyPlot->graph(1)->data()->clear();
    ui->xyPlot->replot();
}

void SSA::on_pushButton_5_clicked()  // cancel
{
    SSA::close();
}

void SSA::on_comboBox_currentIndexChanged(const int arg1)
{
    if (ui->comboBox->currentText() == "Data") {
        xx1 = x_d;
        yy1 = y_d;
        nn = n_d;
    } else {
        xx1 = x_t;
        yy1 = y_t;
        nn = n_t;
    }
    setupPlot();
    ui->xyPlot->replot();
}

void SSA::setupPlot()
{
    ui->xyPlot->graph(0)->setData(xx1,yy1);
    ui->xyPlot->graph(0)->rescaleAxes();
    ui->xyPlot->replot();

}

void SSA::mouseWheel()
{
    QList<QCPAxis *> axlist;
    if(QApplication::keyboardModifiers() & Qt::ShiftModifier) {
        axlist << ui->xyPlot->xAxis <<  ui->xyPlot->yAxis ;}
    else {axlist.clear(); axlist << ui->xyPlot->xAxis ;}

    ui->xyPlot->axisRect()->setRangeDragAxes(axlist);
    ui->xyPlot->axisRect()->setRangeZoomAxes(axlist);

}
