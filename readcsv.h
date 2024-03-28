#ifndef READCSV_H
#define READCSV_H

#include <QDialog>
#include "mainwindow.h"

namespace Ui {
class Readcsv;
}

class Readcsv : public QDialog
{
    Q_OBJECT

public:
    explicit Readcsv(MainWindow *mainWindow = nullptr);

    ~Readcsv();

public slots:

    void buttonBox_accepted();

private:
    Ui::Readcsv *ui;
    bool str_ok = true;
    MainWindow* _mainWindow;

};

#endif // READCSV_H
