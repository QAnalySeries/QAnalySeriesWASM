#ifndef SEDRATE_H
#define SEDRATE_H
#include <qcustomplot.h>
#include <QDialog>

namespace Ui {
class Sedrate;
}

class Sedrate : public QDialog
{
    Q_OBJECT

public:
    explicit Sedrate(QWidget *parent = 0);
    ~Sedrate();
    QCustomPlot *sedGraphic;

    void dragEnterEvent(QDragEnterEvent *e);
    void dropEvent(QDropEvent *e);

signals:
        void myDropSignal();

private:
    Ui::Sedrate *ui;

};

#endif // SEDRATE_H
