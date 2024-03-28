#ifndef INSOLATION_H
#define INSOLATION_H

#include <QDialog>
#include "globals.h"

namespace Ui {
class Insolation;
}

class Insolation : public QDialog
{
    Q_OBJECT

public:
    explicit Insolation(QWidget *parent = nullptr);
    double Insol(Orbit orb, double lon, double lat, double S0);
    double Insol_l1l2(Orbit orb, double l1, double l2, double lat, double S0);
    double dtdnu(Orbit orb, double lon);
    ~Insolation();

private slots:
    void on_cancelButton_clicked();

    void on_okButton_clicked();

    void on_comboBox_currentIndexChanged(const int arg1);

private:
    Ui::Insolation *ui;
};

#endif // INSOLATION_H
