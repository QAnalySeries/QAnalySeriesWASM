#include <QApplication>
#include <QStatusBar>
#include <QToolTip>
#include <QDebug>

#include "globals.h"

#include "tuning.h"
#include "ui_tuning.h"
#include "qcustomplot.h"
#include "mainwindow.h"



Tuning::Tuning(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Tuning)
{
    ui->setupUi(this);

    Tuning::setGeometry(screenw*0.15,screenh*0.65,screenw*0.85,screenh*0.35);

    tuneGraphic = new QCustomPlot();
    ui->gridLayout->addWidget(tuneGraphic,0,0,1,1);

    setAcceptDrops(true);
    connect(tuneGraphic, &QCustomPlot::mouseMove, this, &Tuning::showPointToolTip);
}

Tuning::~Tuning()
{
    delete ui;
}

void Tuning::dragEnterEvent(QDragEnterEvent *e)
{
    if (e->mimeData()->hasUrls()) {
        e->acceptProposedAction();
    }
}

void Tuning::dropEvent(QDropEvent *e)
{
    foreach (const QUrl &url, e->mimeData()->urls()) {
        QString fileName = url.toLocalFile();
        QPoint p = e->position().toPoint();
        QSize s = ui->widget->size();
        QFile file(fileName);
        if(!file.exists()){
            qDebug() << "File does not exist: "<<fileName;return;
        }else{
            qDebug() << fileName<<" opened";
        }

        QString fType;

        if (p.y() < s.height()/2.0) {fType = "d";} else {fType = "t";}
        readFile(ui->widget, nullptr);

        emit myDropSignal();
    }
}

void Tuning::showPointToolTip(QMouseEvent *event)
{

    if (n_d==0 || n_t==0) return;
    int it;
    double y = tuneGraphic->yAxis->pixelToCoord(event->pos().y());

    if (y<0.0) {double x = tuneGraphic->xAxis->pixelToCoord(event->pos().x());
        it = tuneGraphic->graph(1)->findBegin(x,true);
        tuneGraphic->xAxis->setLabel(QString("Target, x: %1, y:  %2").arg(x_t[it]).arg(y_t[it]));
        tuneGraphic->replot();}
    else {double x = tuneGraphic->xAxis2->pixelToCoord(event->pos().x());
        it = tuneGraphic->graph(0)->findBegin(x,true);
        tuneGraphic->xAxis2->setLabel(QString("Data, x: %1, y:  %2").arg(x_d[it]).arg(y_d[it]));
        tuneGraphic->replot();
    }
    //setToolTip(QString("%1 , %2").arg(x).arg(y));
}
