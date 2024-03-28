#include "correlation.h"
#include "ui_correlation.h"
#include "globals.h"
#include "qcustomplot.h"
#include <QApplication>

Correlation::Correlation(QWidget *parent) : QDialog(parent), ui(new Ui::Correlation) {
    ui->setupUi(this);
    Correlation::setGeometry(screenw*0.15,screenh*0.05,screenw*0.85,screenh*0.30);
    corrGraphic = new QCustomPlot();
    ui->gridLayout->addWidget(corrGraphic,0,0,1,1);
}

Correlation::~Correlation() {
    delete ui;
}
