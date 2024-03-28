#include "insolation.h"
#include "ui_insolation.h"
#include "globals.h"
#include "math.h"
#include "QDir"
#include "QtSql/QSqlDatabase"
#include "QSqlQuery"
#include "QDebug"
#include "QSqlQueryModel"
#include "QSqlError"
#include <QMessageBox>

// Algorithm is based on PalInsol R package after Michel Crucifix, 2016
// (https://bitbucket.org/mcrucifix/insol)

Insolation::Insolation(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Insolation)
{
    ui->setupUi(this);

    QSqlDatabase db = QSqlDatabase::database();
    if (db.open()) {
    //qDebug() << "insolation db open = " << ok;

      QSqlQuery query = QSqlQuery(db);

      query.exec("SELECT min(t) FROM La2004ecc");
      query.first();
      ui->lineEdit->setText(query.value(0).toString());
      query.clear();
      query.exec("SELECT max(t) FROM La2004ecc");
      query.first();
      ui->lineEdit_2->setText(query.value(0).toString());
      ui->comboBox->addItem("Daily Ins.");
      ui->comboBox->addItem("Integrated Ins.");
      query.clear();
    } else {
      qDebug() << "insolation db open = NOT OK: " << db.lastError();
    }
}

Insolation::~Insolation()
{
    delete ui;
}


double Insolation::Insol(Orbit orb, double lon, double lat, double S0)
{
    double  varpi = orb.varpi;
    double  eps = orb.eps;
    double  ecc = orb.ecc;

    double  nu = lon - varpi;
    double  rho = (1.0-ecc*ecc)/(1.0+ecc*cos(nu));
    double  sindelta = sin(eps)*sin(lon);
    double  cosdelta = sqrt(1-sindelta*sindelta);
    double  sinlatsindelta = sin(lat)*sindelta;
    double  coslatcosdelta = cos(lat)*cosdelta;
    double  cosH0 = std::fmin(std::fmax(-1.0,-sinlatsindelta/coslatcosdelta),1.0);
    double  sinH0 = sqrt(1.0-cosH0*cosH0);
    double  H0 = acos(cosH0);
    return  S0/(M_PI*rho*rho)*(H0*sinlatsindelta+coslatcosdelta*sinH0);
 }

void Insolation::on_cancelButton_clicked()
{
    qDebug() << "Insolation Cancel Button clicked";
    Insolation::close();
}

void Insolation::on_okButton_clicked()
{
    qDebug() << "Insolation OK Button clicked";

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

    Orbit orb;
    QVector<double> la_time, la_ecc, la_eps, la_varpi;
    int t1, t2;

    t1 =  ui->lineEdit->text().toDouble();
    t2 =  ui->lineEdit_2->text().toDouble();
    x_t.clear(); y_t.clear();

    QSqlDatabase db = QSqlDatabase::database();
    if (db.open()) {
      qDebug() << "insolLa04 db open = OK";
      QSqlQuery query = QSqlQuery(db);
      query.exec("SELECT t, s FROM  La2004pibar  where t between " + QString::number(t1) +" and " + QString::number(t2));
      query.first();

      while (query.next())
      {
          la_time.append(query.value("t").toDouble());
          la_varpi.append(query.value("s").toDouble());
      }

      query.clear();
      query.exec("SELECT s FROM  La2004ecc  where t between "+ QString::number(t1) +" and " + QString::number(t2));
      query.first();

      while (query.next())
      {
          la_ecc.append(query.value("s").toDouble());
      }

      query.clear();
      query.exec("SELECT s FROM  La2004obl  where t between " + QString::number(t1) + " and " + QString::number(t2));
      query.first();

      while (query.next())
      {
          la_eps.append(query.value("s").toDouble());
      }

      // Calculate insolation

      double ins, lon1, lon2, lat, S0;

      lon1 = ui->lineEdit_4 ->text().toDouble()*M_PI/180.0;
      lon2 = ui->lineEdit_6 ->text().toDouble()*M_PI/180.0;
      lat = ui->lineEdit_5->text().toDouble()*M_PI/180.0;
      S0 = ui->lineEdit_3->text().toDouble();

      for (int i=0; i<la_time.length(); i++) {
          orb.time = la_time[i];
          orb.ecc = la_ecc[i];
          orb.eps = la_eps[i];
          orb.varpi = fmod(qAbs((la_varpi[i]+M_PI)), (2.0*M_PI));
          x_t.append(la_time[i]);
          if (ui->comboBox->currentIndex()==0) {ins = Insol(orb, lon1, lat, S0);} else {
              ins = Insol_l1l2(orb, lon1, lon2, lat, S0);
          }


          y_t.append(ins);
      }
    } else {
      qDebug() << "Insolation OK Button db open NOT OK: " << db.lastError();
    }

    Insolation::close();
    accept();
}

// time increment corresponding a tsl increment
double Insolation::dtdnu(Orbit orb, double lon)
{
    double  varpi = orb.varpi;
    double  ecc = orb.ecc;
    double nu = lon - varpi;
    double xec = ecc*ecc;
    double rho = (1.0-xec)/(1.0+ecc*cos(nu));
    return rho*rho/sqrt(1.0-xec);
}

double Insolation::Insol_l1l2(Orbit orb, double l1, double l2, double lat, double S0)
{
    if (l2<l1) {
      l2=l2+2*M_PI;
    }
    double  Dl0 = fmod((l2-l1), (2.0*M_PI));

    if (Dl0 == 0.0) {
      Dl0=2.0*M_PI;
    }
    int N =  ceil(Dl0*180.0/M_PI);
    double dl = Dl0/(double)N;
    QVector<double> L;
    for (int i=0; i< N+1; i++) {
      L.append(l1+i*dl);
    }

    QVector<double> is; is.resize(L.length());
    for (int i=0; i<is.length();i++) {
      is[i] = Insol(orb,L[i], lat, S0)*dtdnu(orb, L[i])*180.0/M_PI;
    }
    double XCORR = 86.4 *  365.24219876 / 360.0;
    double sum=0.0;
    for (int i=1; i<N;++i) {
      sum=sum+is[i];
    }
    sum = (sum + 0.5*is[0]+0.5*is[N])*dl*XCORR;
    return sum;
}

void Insolation::on_comboBox_currentIndexChanged(const int arg1)
{
    if ( ui->comboBox->currentIndex()==0) {
        ui->lineEdit_6->setEnabled(false);
        ui->label_6->setEnabled(false);
        ui->lineEdit_4->setText("90");
    } else {
        ui->lineEdit_6->setEnabled(true);
        ui->label_6->setEnabled(true);
        ui->lineEdit_4->setText("0");
        ui->lineEdit_6->setText("360");
    }
}
