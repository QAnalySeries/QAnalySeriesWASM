#include "laskar.h"
#include "ui_laskar.h"
#include "QDebug"
#include "QtSql/QSqlDatabase"
#include "QSqlQuery"
#include "QDebug"
#include "QSqlQueryModel"
#include "QSqlError"
#include "QVector"
#include "QFile"
#include "globals.h"
#include <QMessageBox>
#include <QtMath>
#include "QDir"

Laskar::Laskar(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Laskar)
{
    ui->setupUi(this);

    qDebug() << "in Laskar Init";
    QSqlDatabase db = QSqlDatabase::database();
    if (db.open()) {
      qDebug() << "laskar db open = ok";
      QSqlQuery query = QSqlQuery(db);
      query.clear();
      query.exec("SELECT * FROM sqlite_master WHERE type='table'");

      ui->comboBox->clear();
      while (query.next())
      {
          QString item = query.value("tbl_name").toString();
          qDebug() << "adding :" << item;
         ui->comboBox->addItem(item);
      }

      ui->comboBox->setCurrentIndex(0);

      query.clear();
      query.exec("SELECT min(t) FROM " + ui->comboBox->itemText(0));
      query.first();
      ui->lineEdit->setText(query.value(0).toString());
      query.clear();
      query.exec("SELECT max(t) FROM " + ui->comboBox->itemText(0));
      query.first();
      ui->lineEdit_2->setText(query.value(0).toString());
    } else {
      qDebug() << "laskar db open = NOT OK: " << db.lastError();
    }
}

Laskar::~Laskar()
{
    delete ui;
}


void Laskar::on_cancelButton_clicked()
{
    qDebug() << "in Laskar::on_cancelButton_clicked";
    Laskar::close();
}

void Laskar::on_okButton_clicked()
{
    qDebug() << "in Laskar::on_OKButton_clicked";

    x_tm.clear(); y_tm.clear(); xind_tm.clear(); yind_tm.clear();
    bool ok;
    ui->lineEdit->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }
    ui->lineEdit_2->text().toDouble(&ok);
    if (!ok) {QMessageBox msgBox;
        msgBox.setText("Invalid number format");
        msgBox.show();
        return;
    }

    QSqlDatabase db = QSqlDatabase::database();
    if (db.open()) {
      qDebug() << "laskar db open = OK";
      QSqlQuery query = QSqlQuery(db);
      qDebug() << "laskar db selected  = " << ui->comboBox->currentText();
      query.exec("SELECT t, s FROM "+ ui->comboBox->currentText()+" where t between "+ui->lineEdit->text()+" and "+ui->lineEdit_2->text());

      x_t.clear();y_t.clear(); n_t=0; n_tm=0;
      while (query.next()) {
        x_t.append(query.value("t").toDouble());
        y_t.append(query.value("s").toDouble());
      }
    } else {
      qDebug() << "laskar db open = NOT OK: " << db.lastError();
    }
    Laskar::close();
    accept();
}

void Laskar::on_comboBox_currentIndexChanged(const int arg1)
{
    return;
    QSqlDatabase db = QSqlDatabase::database();
    bool ok2 = db.open();
    qDebug() << "laskar db2 open = " << ok2;
    QSqlQuery query = QSqlQuery(db);

    query.clear();
    query.exec("SELECT min(t) FROM " + ui->comboBox->currentText());
    query.first();
    ui->lineEdit->setText(query.value(0).toString());
    query.clear();
    query.exec("SELECT max(t) FROM " + ui->comboBox->currentText());
    query.first();
    ui->lineEdit_2->setText(query.value(0).toString());
}
