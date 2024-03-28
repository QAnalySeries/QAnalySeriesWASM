#ifndef ETP_H
#define ETP_H

#include <QDialog>

namespace Ui {
class ETP;
}

class ETP : public QDialog
{
    Q_OBJECT

public:
    explicit ETP(QWidget *parent = nullptr);
    ~ETP();

signals:
    void etpOKButtonPressed();

private slots:
    void on_cancelButton_clicked();

    void on_okButton_clicked();

private:
    Ui::ETP *ui;
};

#endif // ETP_H
