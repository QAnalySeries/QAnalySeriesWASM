#ifndef DATAPROC_H
#define DATAPROC_H

#include <QDialog>

namespace Ui {
class DataProc;
}

class DataProc : public QDialog
{
    Q_OBJECT

public:
    explicit DataProc(QWidget *parent = 0);
    ~DataProc();

    void setupPlot();
private slots:
    void on_saveResultsAndCloseWindow_clicked();

    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_trimData_clicked();

    void on_pushButton_8_clicked();

    void on_pushButton_9_clicked();

   void on_tabWidget_currentChanged(int index);

    void on_pushButton_11_clicked();

    void on_pushButton_12_clicked();

    void on_pushButton_13_clicked();

    void on_checkBox_clicked();

    void on_pushButton_2_clicked();

    void on_comboBox_2_currentIndexChanged(const QString &arg1);

    void mouseWheel();

    void on_pushButton_3_clicked();

private:
    Ui::DataProc *ui;
};

#endif // DATAPROC_H
