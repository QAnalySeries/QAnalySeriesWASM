#include "readcsv.h"
#include <ui_readcsv.h>
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <QMessageBox>
#include <QDebug>
#include "globals.h"
#include "mainwindow.h"

QRegularExpression sep_vars("\\s*,\\s*|\\s*;\\s*|\\s+");
QRegularExpression sep_lines("\u2029|\\r\\n|\\r|\\n");
QString oneline;
QStringList lines;


Readcsv::Readcsv(MainWindow *mainWindow) :
    QDialog(mainWindow),
    ui(new Ui::Readcsv)
{
    ui->setupUi(this);
    _mainWindow = mainWindow;
    qDebug() << "Readcsv started";

    if (!fileContent.isEmpty()) {
 //       qDebug() << "Readcsv in fileContent not empty branch";
        oneline = fileContent;
    } else if (!fileName.isEmpty()) {
        qDebug() << "Readcsv in fileContent empty branch";

        QFile inputFile(fileName);

        if (inputFile.open(QIODevice::ReadOnly)) {
          QTextStream in(&inputFile);
          oneline = in.readAll();
          inputFile.close();
        }
    }
    oneline = oneline.trimmed();
//    qDebug() << "read oneline: " << oneline << "\n";
    lines = oneline.split(sep_lines);

    QString line = lines[0];
    line = line.trimmed();
    QStringList list = line.split(sep_vars);

    bool ok = true;
            //str_ok = true; // first line are numbers

    double fl;
    foreach(QString str, list){
        str.toDouble(&ok);
        if (!ok) {
          str_ok = false;
        }
    }

    QStringList headerLabels;

    for (int i=0; i<list.length(); ++i){
        if (str_ok) {
          ui->comboBox->addItem(QString::number(i+1));
          ui->comboBox_2->addItem(QString::number(i+1));
          headerLabels << QString::number(i+1);
        } else {
          ui->comboBox->addItem(list[i]);
          ui->comboBox_2->addItem(list[i]);
          headerLabels << list[i];
        }
    }

    ui->comboBox->setCurrentIndex(0);
    ui->comboBox_2->setCurrentIndex(1);

    ui->tableWidget->setColumnCount(list.length());
    ui->tableWidget->setHorizontalHeaderLabels(headerLabels);


    if(!str_ok) lines.removeFirst();

    for(int i = 0; i<lines.length();i++) {
        line = lines[i];
        line = line.trimmed();
        list = line.split(sep_vars);
        ui->tableWidget->insertRow(i);
        for(int j=0; j<list.length();j++) {
            ui->tableWidget->setItem(i,j,new QTableWidgetItem(list[j]));
        };
    }
//    qDebug() << "Readcsv ended" << Qt::endl;;
}

Readcsv::~Readcsv()
{
    delete ui;
}


void Readcsv::buttonBox_accepted()
{
    QVector<double> x, y;
    bool ok;
//    qDebug() << "Readcsv buttonBox_accepted started";
//    qDebug() << "Line length" << lines.length();
       for(int i=0;i<lines.length();i++) {
              QString line = lines[i];
 //             qDebug() << "i= " << i << "  " << line;
              QStringList list = line.trimmed().split(sep_vars);
 //             qDebug() << "list.length " << list.length();
              if (list.length() <= 1 || list[1] == "") {
                qDebug() << "Skipping line " << line;
                continue;
              }
              //qDebug() << "currentIndex X" << ui->comboBox->currentIndex();
              //qDebug() << "currentIndex Y" << ui->comboBox_2->currentIndex();
              x.append(list.at(ui->comboBox->currentIndex()).toDouble(&ok));
              if (!ok) {
                QMessageBox msgBox;
                  msgBox.setText("Invalid number format, line: "+QString::number(i+1));
                  msgBox.show();
                  return;
              }
              y.append(list.at(ui->comboBox_2->currentIndex()).toDouble(&ok));
              if (!ok) {
                  QMessageBox msgBox;
                  msgBox.setText("Invalid number format, line: "+QString::number(i+1));
                  msgBox.show();
                  return;
              }

              int xLength = x.length();
              if (xLength>2) {
                if(x[xLength-1] < x[xLength-2]) {
                  QMessageBox msgBox;
                  msgBox.setText("x must be in increasing order!");
                  msgBox.show();
                  return;
                }
              }
           }




    if (fType == "d") {
        x_d.clear();y_d.clear();x_d1.clear();y_dn.clear();
        x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();n_tm=0;

        x_d = x;
        y_d = y;
        n_d=y_d.length();
        x_d1=x_d;
        data_orig.x=x_d;
        data_orig.y=y_d;
    }

    if (fType == "t") {
        x_t.clear();y_t.clear();y_tn.clear();
        x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();n_tm=0;
        x_t = x;
        y_t = y;
        n_t=y_t.length();
        target_orig.x=x_t;
        target_orig.y=y_t;
    }

    if (fType == "tm") {
        tm_orig.x=x;
        tm_orig.y=y;
    }

    this->setModal(false);

    this->close();
    //QCoreApplication::processEvents();
    qDebug() << "Readcsv buttonBox_accepted ended";

    if (fType == "tm") {
      this->_mainWindow->tmLoadCallback();
    } else if (fType == "t") {
      this->_mainWindow->targetLoadCallback();
    } else if (fType == "d") {
      this->_mainWindow->dataLoadCallback();
    }
    this->_mainWindow = nullptr;
}
