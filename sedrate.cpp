#include <QApplication>

#include "sedrate.h"
#include "ui_sedrate.h"
#include "globals.h"
#include "qcustomplot.h"
#include "mynr.h"
#include "globals.h"
#include "mainwindow.h"


Sedrate::Sedrate(QWidget *parent) : QDialog(parent), ui(new Ui::Sedrate)
{
    ui->setupUi(this);

    Sedrate::setGeometry(screenw*0.15,screenh*0.35,screenw*0.85,screenh*0.30);
    sedGraphic = new QCustomPlot();
    ui->gridLayout->addWidget(sedGraphic,0,0,1,1);

    setAcceptDrops(true);
}

Sedrate::~Sedrate()
{
    delete ui;
}

void Sedrate::dragEnterEvent(QDragEnterEvent *e)
{
    if (e->mimeData()->hasUrls()) {
        e->acceptProposedAction();
    }
}

void Sedrate::dropEvent(QDropEvent *e)
{
    if (n_d<2 || n_t<2) {
        QMessageBox msgBox;
        msgBox.setText("Download Data and Target first!");
        msgBox.show();
        return;
    }

    foreach (const QUrl &url, e->mimeData()->urls()) {
        fileName = url.toLocalFile();
        QFile file(fileName);
        if(!file.exists()){
            qDebug() << "File does not exist: "<<fileName;return;
        } else {
            qDebug() << fileName<<" opened";
        }

        double x,y, interp;
        QVector<double> xtmp(2), ytmp(2);

        x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();
        int i=0;

        fType = "tm";
        readFile(this, nullptr);



        for (int j = 0; j<tm_orig.x.length();j++) {
            x = tm_orig.x[j]; y = tm_orig.y[j];
            x_tm.append(x); y_tm.append(y);

            i=0;
            while (i<x_d.length() && x>x_d[i]) {
                i++; xind_tm.append(i); //
            }

            if (i>(x_d.length()-1)) {
                xtmp[0]=x_d[i-2];
                xtmp[1]=x_d[i-1];
                ytmp[0]=y_d[i-2];
                ytmp[1]=y_d[i-1];
            } else if (i==0  && x!=x_d[i])
            {xtmp[0]=x_d[0];
                xtmp[1]=x_d[1];
                ytmp[0]=y_d[0];
                ytmp[1]=y_d[1];
            } else if (i>0 && i<=(x_d.length()-1) && x!=x_d[i]) {
                xtmp[0]=x_d[i-1];
                xtmp[1]=x_d[i];
                ytmp[0]=y_d[i-1];
                ytmp[1]=y_d[i];
            }

            if (i>(x_d.length()-1) || x!=x_d[i]) {
                x_d.insert(i,x);
                mynr::linterp(xtmp,ytmp,x_d[i],interp);
                y_d.insert(i,interp);
            }

            i=0;
            while (i<x_t.length() && y>x_t[i]) {
                i++;
                yind_tm.append(i);
            }

            if (i>(x_t.length()-1)) {
                xtmp[0]=x_t[i-2];
                xtmp[1]=x_t[i-1];
                ytmp[0]=y_t[i-2];
                ytmp[1]=y_t[i-1];
            } else if (i==0 && y!=x_t[i] ) {
                xtmp[0]=x_t[0];
                xtmp[1]=x_t[1];
                ytmp[0]=y_t[0];
                ytmp[1]=y_t[1];
            } else if (i>0 && i<=(x_t.length()-1) && y!=x_t[i]) {
                xtmp[0]=x_t[i-1];
                xtmp[1]=x_t[i];
                ytmp[0]=y_t[i-1];
                ytmp[1]=y_t[i];
            }

            if (i>(x_t.length()-1) || y!=x_t[i]) {
                x_t.insert(i,y);
                mynr::linterp(xtmp,ytmp,x_t[i],interp);
                y_t.insert(i,interp);
            }
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

        tm_orig.x = x_tm;
        tm_orig.y = y_tm;
        xind_tm_orig = xind_tm;
        yind_tm_orig = yind_tm;
        data_orig.x=x_d;
        data_orig.y=y_d;
        target_orig.x=x_t;
        target_orig.y=y_t;

        emit myDropSignal();
    }
}
