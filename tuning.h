#ifndef TUNING_H
#define TUNING_H
#include <qcustomplot.h>
#include <QDialog>


namespace Ui {
class Tuning;
}

class Tuning : public QDialog
{
    Q_OBJECT

public:
    explicit Tuning(QWidget *parent = 0);
    ~Tuning();
    QCustomPlot *tuneGraphic;
    QLabel *data_selected;
    QLabel *target_selected;
    void dragEnterEvent(QDragEnterEvent *e);
    void dropEvent(QDropEvent *e);
    void showPointToolTip(QMouseEvent *event);
//    bool event( QEvent *event );

signals:
        void myDropSignal();


private:
    Ui::Tuning *ui;
};

#endif // TUNING_H
