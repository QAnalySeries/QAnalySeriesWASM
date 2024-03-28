#include "etp.h"
#include "ui_etp.h"
#include "QtSql/QSqlDatabase"
#include "QSqlQuery"
#include "QDebug"
#include "QSqlQueryModel"
#include "QSqlError"
#include "QVector"
#include "QFile"
#include "globals.h"
#include <QMessageBox>
#include "QDir"
#include <QtMath>
#include "mynr.h"

ETP::ETP(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ETP)
{
  
    ui->setupUi(this);

    ui->comboBox->addItem("e*sin(w)");
    ui->comboBox->addItem("sin(w)");
    ui->comboBox->setCurrentIndex(0);

    QSqlDatabase db = QSqlDatabase::database();
    bool ok = db.open();
 //   qDebug() << "etp db open = " << ok;

    QSqlQuery query = QSqlQuery(db);

    query.exec("SELECT min(t) FROM La2004ecc");
    query.first();
    ui->lineEdit->setText(query.value(0).toString());
    query.clear();
    query.exec("SELECT max(t) FROM La2004ecc");
    query.first();
    ui->lineEdit_2->setText(query.value(0).toString());
}

ETP::~ETP()
{
    delete ui;
}

void ETP::on_cancelButton_clicked() // Cancel button
{
    ETP::close();
}

void ETP::on_okButton_clicked()  // Ok button
{
    QVector <double> ecc, obl, pr;
    double mean, var, eWt, tWt, pWt;

    bool ok;
    ui->lineEdit->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;}
    ui->lineEdit_2->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;}
    ui->lineEdit_3->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;}
    ui->lineEdit_4->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;}
    ui->lineEdit_5->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;}

    eWt = ui->lineEdit_3->text().toDouble();
    tWt = ui->lineEdit_4->text().toDouble();
    pWt = ui->lineEdit_5->text().toDouble();


    x_tm.clear();
    y_tm.clear();
    xind_tm.clear();
    yind_tm.clear();
    x_t.clear();
    y_t.clear();
    n_t=0;
    n_tm=0;

    QSqlDatabase db = QSqlDatabase::database();
    if (db.open()) {
//      //qDebug() << "etp db open = " << dbok;

      QSqlQuery query = QSqlQuery(db);
      // Eccentricity
      query.exec("SELECT t, s FROM La2004ecc where t between " + ui->lineEdit->text()+" and " + ui->lineEdit_2->text());
      while (query.next())
      {
        x_t.append(query.value("t").toDouble());
        ecc.append(query.value("s").toDouble());
      }
      n_t = x_t.length();
      query.clear();

      // Obliquity
      query.exec("SELECT s FROM La2004obl where t between "+ui->lineEdit->text()+" and "+ui->lineEdit_2->text());
      while (query.next())   obl.append(query.value("s").toDouble());
      query.clear();

      // Precession
      query.exec("SELECT s FROM La2004pibar where t between "+ui->lineEdit->text()+" and "+ui->lineEdit_2->text());
      while (query.next())  pr.append(qSin(query.value("s").toDouble()));
      query.clear();
      if (ui->comboBox->currentIndex() == 0) {
        for (int i=0;i<n_t;++i) {
            pr[i] = pr[i]*ecc[i];
        }
      }

      // Normalisation
      if (ui->checkBox->isChecked()) {
       mynr::meanandvar(ecc, mean, var);
       var = sqrt(var);  //standard deviation
       for (int i=0;i<n_t;++i) {
            ecc[i] = (ecc[i]-mean)/var;
       }

       mynr::meanandvar(obl, mean, var);
       var = sqrt(var);  //standard deviation
       for (int i=0;i<n_t;++i) {
            obl[i] = (obl[i]-mean)/var;
       }

       mynr::meanandvar(pr, mean, var);
       var = sqrt(var);  //standard deviation
       for (int i=0;i<n_t;++i) {
            pr[i] = (pr[i]-mean)/var;
       }
      }

      y_t.resize(n_t);
      for (int i=0;i<n_t;++i) {
        y_t[i] = ecc[i]*eWt+obl[i]*tWt+pr[i]*pWt;
      }
    } else {
      qDebug() << "opening db NOT OK: " << db.lastError();
    }

    ETP::close();
    accept();
}
