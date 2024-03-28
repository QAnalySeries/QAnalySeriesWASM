#ifndef LASKAR_H
#define LASKAR_H

#include <QDialog>

namespace Ui {
class Laskar;
}

class Laskar : public QDialog
{
    Q_OBJECT

public:
    explicit Laskar(QWidget *parent = 0);
    ~Laskar();



private slots:
    void on_cancelButton_clicked();

    void on_okButton_clicked();

    void on_comboBox_currentIndexChanged(const int arg1);

private:
    Ui::Laskar *ui;
};

#endif // LASKAR_H
