#ifndef SSA_H
#define SSA_H

#include <QDialog>

namespace Ui {
class SSA;
}

class SSA : public QDialog
{
    Q_OBJECT

public:
    explicit SSA(QWidget *parent = 0);
    ~SSA();
    void setupPlot();

private slots:
    void on_pushButton_3_clicked();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_checkBox_stateChanged(int arg1);


    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_comboBox_currentIndexChanged(const int arg1);

    void mouseWheel();

private:
    Ui::SSA *ui;
};

#endif // SSA_H
